import os
from technical_reliability.variant_qc import (
    load_variants,
    filter_reliable_snvs,
    export_vcf,
    plot_depth_and_ab,
    stratify_by_impact,
    plot_impact_qc
)

def run(annotated_variants_file,output):
    # === Paths ===
    raw_qc_dir = "annotation_prio_qc/raw"
    reliable_qc_dir = "annotation_prio_qc/reliable"
    os.makedirs(raw_qc_dir, exist_ok=True)
    os.makedirs(reliable_qc_dir, exist_ok=True)

    # === Step 1: Load annotated VCF ===
    vcf_raw = list(load_variants(annotated_variants_file))  # Cache in memory

    # === Step 2: Filter reliable variants ===
    reliable_snvs = filter_reliable_snvs(vcf_raw, min_dp=10, min_ab=0.2)

    # === Step 3: Export to new VCF ===
    export_vcf(reliable_snvs, output, load_variants(annotated_variants_file))

    # === Step 4: Plot QC metrics for raw SNVs ===
    plot_depth_and_ab(vcf_raw, raw_qc_dir)
    impact_ab, impact_dp = stratify_by_impact(vcf_raw)
    plot_impact_qc(impact_ab, impact_dp, raw_qc_dir)

    # === Step 5: Plot QC metrics for reliable SNVs ===
    plot_depth_and_ab(reliable_snvs, reliable_qc_dir)
    impact_ab_rel, impact_dp_rel = stratify_by_impact(reliable_snvs)
    plot_impact_qc(impact_ab_rel, impact_dp_rel, reliable_qc_dir)

    return output