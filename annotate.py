import subprocess
from pathlib import Path
import os
import shutil
import gzip

def remove_vcf_fields(
    vcf_in: Path, fields_to_remove: str, vcf_out: Path = None
) -> Path:
    """
    Remove unwanted INFO/FORMAT fields from a VCF using bcftools annotate.

    Parameters:
        vcf_in (Path): Input VCF file (.vcf.gz)
        fields_to_remove (str): Comma-separated list of fields to remove (e.g. "INFO/OLD_VARIANT,FORMAT/DP4")
        vcf_out (Path, optional): Output path. If None, auto-generates based on input name.

    Returns:
        Path: Path to the tidied VCF (.vcf.gz)
    """
    if vcf_out is None:
        base = vcf_in.name.replace(".vcf.gz", "").replace(".vcf", "")
        vcf_out = vcf_in.parent / f"{base}_tidy.vcf.gz"

    subprocess.run(
        [
            "bcftools",
            "annotate",
            "--remove",
            fields_to_remove,
            "-Oz",
            "-o",
            str(vcf_out),
            str(vcf_in),
        ],
        check=True,
    )

    return vcf_out




def decompress_if_needed(src_path, dest_path):
    """
    Decompress a .gz file if needed, else copy as-is.
    Ensures the final file is plain text, as required by SnpEff.
    """
    src_path = Path(src_path)
    dest_path = Path(dest_path)

    if src_path.suffix == ".gz":
        with gzip.open(src_path, "rt") as fin, open(dest_path, "wt") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        shutil.copy(src_path, dest_path)


def build_snpeff_db(
    db_dir,
    key,
    gtf_file,
    reference_file,
    config_path="snpeff.config",
    genome_label=None,
):
    """
    Create and build a custom SnpEff genome database.
    Accepts compressed or uncompressed GTF files.

    Args:
        db_dir (str): Path to parent directory that will hold the genome subfolder
        key (str): Genome key (e.g. 'SCEREVISIAE_YEAST') — must match subfolder name
        gtf_file (str): Path to GTF annotation file (.gtf or .gtf.gz)
        reference_file (str): Path to FASTA reference genome
        config_path (str): Path to SnpEff config file to write/update
        genome_label (str): Optional label for the genome (default = same as key)
    """
    genome_label = genome_label or key
    genome_path = os.path.join(db_dir, key)
    os.makedirs(genome_path, exist_ok=True)
    
    # Copy reference genome into genome folder for build
    shutil.copy(reference_file, os.path.join(genome_path, "sequences.fa"))

    # Copy GTF file with correct output name
    gtf_dest = "genes.gtf.gz" if str(gtf_file).endswith(".gz") else "genes.gtf"
    shutil.copy(gtf_file, os.path.join(genome_path, gtf_dest))

    # Write genome config
    with open(config_path, "w") as config_file:
        config_file.write(f"{key}.genome : {genome_label}\n")

    # Run snpEff build
    subprocess.run([
        "snpEff", "build",
        "-v",
        "-dataDir", db_dir,
        "-c", config_path,
        "-noCheckCds",
        "-noCheckProtein",
        key
    ], check=True)



def annotate_vcf_with_snpeff(
    input_vcf,
    output_vcf,
    genome_key="SCEREVISIAE_YEAST",
    config_path="snpeff.config",
    data_dir="snpeff",
    verbose=True,
):
    """
    Annotate a VCF using a custom SnpEff genome database.

    Args:
        input_vcf (str): Path to input VCF file (.vcf or .vcf.gz)
        output_vcf (str): Path to write annotated output VCF
        genome_key (str): Genome key defined in SnpEff config
        config_path (str): Path to snpEff.config file
        data_dir (str): Directory containing your custom genome
        verbose (bool): If True, run SnpEff with '-v' for progress info
    """
    cmd = [
        "snpEff", "ann",
        "-noDownload", 
        "-dataDir", data_dir,
        "-c", config_path,
    ]
    if verbose:
        cmd.append("-v")
    cmd.append(genome_key)
    cmd.append(input_vcf)

    with open(output_vcf, "w") as out_f:
        subprocess.run(cmd, check=True, stdout=out_f)
