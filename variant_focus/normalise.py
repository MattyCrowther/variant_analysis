import subprocess
from pathlib import Path

def validate_ref_alleles(vcf_in: Path, fasta_ref: Path, vcf_out: Path = None) -> Path:
    """
    Validate that REF alleles in a VCF match the reference genome.

    Parameters:
        vcf_in (Path): Input VCF file (.vcf or .vcf.gz)
        fasta_ref (Path): Reference genome FASTA file
        vcf_out (Path, optional): Output VCF path. If None, auto-generates based on input name.

    Returns:
        Path: Path to the validated VCF (.vcf.gz)
    """
    if vcf_out is None:
        base = vcf_in.name.replace(".vcf.gz", "").replace(".vcf", "")
        vcf_out = vcf_in.parent / f"{base}_validated.vcf.gz"

    subprocess.run([
        "bcftools", "norm",
        "-f", str(fasta_ref),
        "-c", "s",
        "-Oz",
        "-o", str(vcf_out),
        str(vcf_in)
    ], check=True)

    return vcf_out
