import subprocess
from pathlib import Path

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