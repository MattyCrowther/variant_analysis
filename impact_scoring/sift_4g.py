import subprocess
from pathlib import Path
import gzip
import shutil
import tempfile
import pysam
import pandas as pd

def run(vcf_gz_path, database, final_output_vcf):
    vcf_gz_path = Path(vcf_gz_path).resolve()
    database = Path(database).resolve()
    final_output_vcf = Path(final_output_vcf).resolve()
    output_dir = final_output_vcf.parent.resolve()
    jar_path = Path("impact_scoring/SIFT4G_Annotator/SIFT4G_Annotator.jar").resolve()

    # Step 1: Uncompress to temp VCF
    with tempfile.NamedTemporaryFile("wt", suffix=".vcf", delete=False) as tmp_vcf:
        with gzip.open(vcf_gz_path, "rt") as f_in:
            shutil.copyfileobj(f_in, tmp_vcf)
        tmp_vcf_path = Path(tmp_vcf.name)

    # Step 2: Run SIFT4G with proper -r output directory
    subprocess.run([
        "java", "-jar", str(jar_path),
        "-c",
        "-i", str(tmp_vcf_path),
        "-d", str(database),
        "-r", str(output_dir)
    ], check=True)

    # Step 3: Locate most recent *_SIFTpredictions.vcf file
    vcf_candidates = sorted(output_dir.glob("*_SIFTpredictions.vcf"), key=lambda p: p.stat().st_mtime, reverse=True)
    if not vcf_candidates:
        raise FileNotFoundError("SIFT4G did not generate a predictions VCF in the output directory.")
    
    sift_output = vcf_candidates[0]

    # Step 4: Move/rename to final output
    shutil.move(str(sift_output), str(final_output_vcf))

    # Optional: clean up temp file
    tmp_vcf_path.unlink()


def parse_sift_scores(vcf_path):
    vcf_path = Path(vcf_path)
    results = []

    # Decide how to open
    open_func = gzip.open if vcf_path.suffix == ".gz" else open

    with open_func(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
            info = fields[7]

            sift_score = None
            sift_pred = None
            for entry in info.split(";"):
                if entry.startswith("SIFT4G="):
                    sift_score = float(entry.split("=")[1])
                elif entry.startswith("SIFT4G_pred="):
                    sift_pred = entry.split("=")[1]

            results.append({
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
                "sift_score": sift_score,
                "sift_prediction": sift_pred
            })

    return results


def filter_by_score():
    parsed = parse_sift_scores("sift4g_output/sift_scored.vcf.gz")

    # Keep only predicted damaging variants
    return [v for v in parsed if 
            v["sift_score"] is not None 
            and v["sift_score"] < 0.05]


def write_scores_to_tsv(score_list, output_path):
    """
    Save parsed SIFT score entries to a TSV file.
    """
    df = pd.DataFrame(score_list)
    df.to_csv(output_path, sep="\t", index=False)
    return output_path


def write_filtered_vcf(input_vcf_gz, filtered_variants, output_vcf_gz):
    """
    Write a new VCF containing only the filtered variants (by pos/ref/alt match).
    """
    vcf_in = pysam.VariantFile(input_vcf_gz)
    vcf_out = pysam.VariantFile(output_vcf_gz, "w", header=vcf_in.header)

    key_set = {(v["chrom"], v["pos"], v["ref"], v["alt"]) for v in filtered_variants}

    for rec in vcf_in.fetch():
        if (rec.contig, rec.pos, rec.ref, rec.alts[0]) in key_set:
            vcf_out.write(rec)

    vcf_out.close()
    return output_vcf_gz


def fix_vcf_header(vcf_path):
    vcf_path = Path(vcf_path)
    fixed_lines = []

    with vcf_path.open("r") as f:
        for line in f:
            if line.startswith("##SIFT_Threshold:"):
                continue  # skip invalid header
            fixed_lines.append(line)

    fixed_path = vcf_path.with_name(vcf_path.stem + "_fixed.vcf")
    with fixed_path.open("w") as out:
        out.writelines(fixed_lines)

    return fixed_path