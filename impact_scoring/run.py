import os
import subprocess
from pathlib import Path
from impact_scoring.filter import filter_SNV_biallelic
from impact_scoring.filter import filter_for_missense
from impact_scoring.sift_4g import run as sift_run
from impact_scoring.sift_4g import parse_sift_scores
from impact_scoring.sift_4g import write_scores_to_tsv
from impact_scoring.sift_4g import write_filtered_vcf
from impact_scoring.sift_4g import fix_vcf_header
def run(vcf_filename, output_vcf_path, database_dir, write_tsv=None, write_vcf=None):
    """
    Run the full SIFT4G scoring pipeline:
    - Prepare VCF (compress/index)
    - Filter SNVs + biallelic
    - Extract missense variants
    - Run SIFT4G
    - Parse scores
    - (Optional) Write TSV
    - (Optional) Write filtered VCF of damaging SNVs

    Args:
        vcf_filename: Path to annotated input VCF
        output_vcf_path: Path to final output VCF scored by SIFT
        database_dir: Path to SIFT4G DB dir
        write_tsv: Optional path to write parsed scores TSV
        write_filtered_vcf: Optional path to output VCF of deleterious-only SNVs
    Returns:
        List of parsed variant dicts
    """
    # Step 1: Prepare VCF
    vcf_filename = _ensure_bgzipped_and_indexed(vcf_filename)

    # Step 2: Filter biallelic SNVs
    SNV_biallelic_out = Path("SNV_biallelic.vcf.gz")
    filter_SNV_biallelic(vcf_filename, SNV_biallelic_out)

    # Step 3: Filter missense only
    missense_out = Path("biallelic_missense.vcf.gz")
    filter_for_missense(SNV_biallelic_out, missense_out)
    os.remove(SNV_biallelic_out)

    # Step 4: Run SIFT4G
    sift_run(missense_out, database_dir, output_vcf_path)
    os.remove(missense_out)

    output_vcf_path = fix_vcf_header(output_vcf_path)

    # Step 5: Parse predictions
    parsed = parse_sift_scores(output_vcf_path)

    # Step 6 (Optional): Write TSV of parsed scores
    if write_tsv:
        write_scores_to_tsv(parsed, write_tsv)

    # Step 7 (Optional): Write filtered VCF of damaging SNVs
    if write_vcf:
        damaging = [v for v in parsed if v["sift_score"] is not None and v["sift_score"] < 0.05]
        write_filtered_vcf(output_vcf_path, damaging, write_vcf)

    return parsed


def _ensure_bgzipped_and_indexed(vcf_path):
    """
    Ensure a VCF is bgzipped and indexed with tabix.
    If not already bgzipped, compresses and indexes the VCF.
    Returns path to the .vcf.gz file.
    """
    vcf_path = Path(vcf_path)

    if vcf_path.suffix != ".gz":
        # Compress to .vcf.gz
        gz_path = vcf_path.with_suffix(".vcf.gz")
        subprocess.run(["bgzip", "-c", str(vcf_path)], stdout=open(gz_path, "wb"), check=True)
        subprocess.run(["tabix", "-p", "vcf", str(gz_path)], check=True)
        vcf_path.unlink()  # Remove original uncompressed
        return gz_path
    else:
        # Ensure index file exists
        tbi_path = vcf_path.with_name(vcf_path.name + ".tbi")
        if not tbi_path.exists():
            subprocess.run(["tabix", "-p", "vcf", str(vcf_path)], check=True)
        return vcf_path


