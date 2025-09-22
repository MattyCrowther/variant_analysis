import subprocess
from pathlib import Path

def sort_vcf(vcf_in: Path, vcf_out: Path = None) -> Path:
    """
    Sorts a VCF file by chromosome and position using bcftools sort.

    Parameters:
        vcf_in (Path): Input VCF file (.vcf.gz)
        vcf_out (Path, optional): Output sorted VCF. If None, auto-generates.

    Returns:
        Path: Path to the sorted VCF file (.vcf.gz)
    """
    if vcf_out is None:
        base = vcf_in.name.replace(".vcf.gz", "").replace(".vcf", "")
        vcf_out = vcf_in.parent / f"{base}_sorted.vcf.gz"

    subprocess.run([
        "bcftools", "sort",
        "-Oz",
        "-o", str(vcf_out),
        str(vcf_in)
    ], check=True)

    return vcf_out
