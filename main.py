import os
from pathlib import Path

from downloader import download
from utility import decompress_gzip
from converter import gtf_to_bed, filter_bed_to_vcf
from normalise import validate_ref_alleles
from annotate import remove_vcf_fields
from annotate import build_snpeff_db
from annotate import annotate_vcf_with_snpeff
from sort import sort_vcf
from index import index_vcf
from interpreter import parse_annotated_vcf

from variant_qc import (
    load_variants,
    filter_reliable_snvs,
    export_vcf,
    plot_depth_and_ab,
    stratify_by_impact,
    plot_impact_qc
)

DATA_DIR = Path("vcf")
REF_DIR = Path("ref")
STORAGE_DIR = Path("storage")

vcf_input = DATA_DIR / "variants_raw.vcf.gz"
fasta_ref = REF_DIR / "Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"
gtf_url = "https://ftp.ensembl.org/pub/release-109/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.109.gtf.gz"
gtf_gz_path = STORAGE_DIR / "Saccharomyces_cerevisiae.R64-1-1.109.gtf.gz"
bed_output_path = STORAGE_DIR / "Saccharomyces_cerevisiae.CDS.bed"
final_vcf_path = STORAGE_DIR / "yeast_final.vcf.gz"

def download_gtp(gtf_url, gtf_gz_path):
    if not os.path.isfile(gtf_gz_path):
        gtf_gz_path = download(gtf_url, gtf_gz_path)
    return decompress_gzip(gtf_gz_path)

def focus_vcf(gtf_file):
    # Step 2: Convert GTF to BED (CDS only)
    bed_file = gtf_to_bed(gtf_file, bed_output=bed_output_path)

    # Step 3: Filter VCF by CDS regions
    filtered_vcf = filter_bed_to_vcf(vcf_input, bed_file)

    # Step 4: Validate reference alleles
    validated_vcf = validate_ref_alleles(filtered_vcf, fasta_ref)
    filtered_vcf.unlink()

    # Step 5: Remove unnecessary INFO/FORMAT fields
    cleaned_vcf = remove_vcf_fields(validated_vcf, "INFO/OLD_VARIANT,FORMAT/DP4")
    validated_vcf.unlink()

    # Step 6: Sort and compress
    final_vcf = sort_vcf(cleaned_vcf, final_vcf_path)

    # Step 7: Index
    index_vcf(final_vcf)

    return final_vcf


def annotate_vcf(cleaned_vcf,output_vcf):
    build_snpeff_db(
        db_dir="snpeff",
        key="SCEREVISIAE_YEAST",
        gtf_file=gtf_gz_path,
        reference_file=fasta_ref,
        config_path="snpeff.config",
        genome_label="Saccharomyces_cerevisiae_R64-1-1"
    )

    annotate_vcf_with_snpeff(cleaned_vcf,output_vcf)

# === Paths ===
vcf_path = "vcf/variants.ann.vcf"
output_vcf = "vcf/variants_reliable_only.vcf"

raw_qc_dir = "qc/raw"
reliable_qc_dir = "qc/reliable"
os.makedirs(raw_qc_dir, exist_ok=True)
os.makedirs(reliable_qc_dir, exist_ok=True)

# === Step 1: Load annotated VCF ===
vcf_raw = list(load_variants(vcf_path))  # Cache in memory

# === Step 2: Filter reliable variants ===
reliable_snvs = filter_reliable_snvs(vcf_raw, min_dp=10, min_ab=0.2)

# === Step 3: Export to new VCF ===
export_vcf(reliable_snvs, output_vcf, load_variants(vcf_path))  # Reload here to pass original header

# === Step 4: Plot QC metrics for raw SNVs ===
plot_depth_and_ab(vcf_raw, raw_qc_dir)
impact_ab, impact_dp = stratify_by_impact(vcf_raw)
plot_impact_qc(impact_ab, impact_dp, raw_qc_dir)

# === Step 5: Plot QC metrics for reliable SNVs ===
plot_depth_and_ab(reliable_snvs, reliable_qc_dir)
impact_ab_rel, impact_dp_rel = stratify_by_impact(reliable_snvs)
plot_impact_qc(impact_ab_rel, impact_dp_rel, reliable_qc_dir)
