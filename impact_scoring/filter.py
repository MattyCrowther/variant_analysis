import gzip
import subprocess
from pathlib import Path
import tempfile

def is_snv(ref, alt):
    return len(ref) == 1 and len(alt) == 1 and ref != alt

def is_biallelic(alt):
    return "," not in alt

def filter_SNV_biallelic(vcf, output):
    """
    Extracts biallelic SNVs from a gzipped VCF.
    Writes uncompressed VCF to disk, then bgzips + indexes.
    """
    # Step 1: Write uncompressed .vcf to temp file
    with tempfile.NamedTemporaryFile("wt", delete=False, suffix=".vcf") as tmp:
        tmp_path = Path(tmp.name)
        with gzip.open(vcf, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    tmp.write(line)
                    continue
                fields = line.strip().split("\t")
                ref, alt = fields[3], fields[4]
                if is_snv(ref, alt) and is_biallelic(alt):
                    tmp.write(line)

    # Step 2: bgzip + tabix
    subprocess.run(["bgzip", "-c", str(tmp_path)], stdout=open(output, "wb"), check=True)
    subprocess.run(["tabix", "-p", "vcf", str(output)], check=True)

    # Clean up
    tmp_path.unlink()

def filter_for_missense(vcf, output):
    """
    Extracts missense variants from a gzipped VCF.
    Assumes `ANN=` field from SnpEff is present.
    Writes uncompressed VCF to temp, then bgzips + indexes.
    """
    with tempfile.NamedTemporaryFile("wt", delete=False, suffix=".vcf") as tmp:
        tmp_path = Path(tmp.name)
        with gzip.open(vcf, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    tmp.write(line)
                    continue
                info = line.split("\t")[7]
                ann_field = next((x for x in info.split(";") if x.startswith("ANN=")), None)
                if ann_field and "missense_variant" in ann_field:
                    tmp.write(line)

    subprocess.run(["bgzip", "-c", str(tmp_path)], stdout=open(output, "wb"), check=True)
    subprocess.run(["tabix", "-p", "vcf", str(output)], check=True)

    tmp_path.unlink()
