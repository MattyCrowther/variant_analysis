import subprocess
from pathlib import Path

def index_vcf(vcf_file: Path) -> Path:
    """
    Indexes a compressed VCF file using tabix.

    Parameters:
        vcf_file (Path): Compressed VCF file (.vcf.gz)

    Returns:
        Path: Path to the created index file (.tbi)
    """
    subprocess.run(["tabix", "-p", "vcf", str(vcf_file)], check=True)
    return vcf_file.with_suffix(vcf_file.suffix + ".tbi")
