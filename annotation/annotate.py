import subprocess
import os
import shutil


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
