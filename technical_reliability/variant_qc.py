from cyvcf2 import VCF, Writer
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import os


def load_variants(vcf_path):
    """Load variants from a VCF file using cyvcf2.

    Parameters:
        vcf_path (str): Path to the input VCF file.

    Returns:
        cyvcf2.VCF: A VCF object that can be iterated over to access variants.
    """
    return VCF(vcf_path)


def calculate_ab(ad):
    """Calculate allele balance from the AD (allele depth) field.

    Parameters:
        ad (list[int]): A list of allele depths (e.g., [ref_count, alt_count]).

    Returns:
        float: Allele balance (alt / (ref + alt)). Returns 0.0 if undefined.
    """
    if ad is None or sum(ad) == 0:
        return 0.0
    return ad[1] / sum(ad)


def is_reliable(
    variant, min_dp=10, min_ab=0.2, 
    min_mq=40, allowed_filters={"PASS", ".", None}
):
    """Check whether a variant passes technical reliability filters.

    Parameters:
        variant (cyvcf2.Variant): The variant record.
        min_dp (int): Minimum required read depth.
        min_ab (float): Minimum allele balance (alt / total).
        min_mq (int): Minimum mapping quality.
        allowed_filters (set[str]): Set of FILTER values considered acceptable.

    Returns:
        bool: True if variant passes all reliability checks, False otherwise.
    """
    dp = variant.format("DP")[0]
    ad = variant.format("AD")[0]
    ab = calculate_ab(ad)
    filt = variant.FILTER

    if dp < min_dp or ab < min_ab:
        return False
    if filt not in allowed_filters:
        return False
    if check_strand_bias(variant):
        return False
    if low_mapping_quality(variant, threshold=min_mq):
        return False

    return True


def filter_reliable_snvs(
    vcf, min_dp=10, min_ab=0.2, min_mq=40, allowed_filters={"PASS", ".", None}
):
    """Filter a VCF for SNVs passing reliability thresholds.

    Parameters:
        vcf (cyvcf2.VCF): The input VCF object.
        min_dp (int): Minimum read depth threshold.
        min_ab (float): Minimum allele balance.
        min_mq (int): Minimum mapping quality.
        allowed_filters (set[str]): Acceptable FILTER field values.

    Returns:
        list[cyvcf2.Variant]: List of variants that meet reliability criteria.
    """
    return [
        v
        for v in vcf
        if is_reliable(
            v,
            min_dp=min_dp,
            min_ab=min_ab,
            min_mq=min_mq,
            allowed_filters=allowed_filters,
        )
    ]


def export_vcf(variants, output_path, template_vcf):
    """Export a list of variant records to a new VCF file.

    Parameters:
        variants (list[cyvcf2.Variant]): Variants to write.
        output_path (str): Path to the output VCF file (can be .vcf or .vcf.gz).
        template_vcf (cyvcf2.VCF): Template VCF used to preserve header/meta info.

    Returns:
        None
    """
    writer = Writer(output_path, template_vcf)
    for var in variants:
        writer.write_record(var)
    writer.close()


def check_strand_bias(variant):
    """Detect strand bias using SAF and SAR INFO fields.

    Parameters:
        variant (cyvcf2.Variant): Variant record with SAF and SAR fields.

    Returns:
        bool: True if all ALT reads are on a single strand, False otherwise.
    """
    saf = variant.INFO.get("SAF")
    sar = variant.INFO.get("SAR")
    return saf == 0 or sar == 0 if saf is not None and sar is not None else False


def low_mapping_quality(variant, threshold=40):
    """Check if a variant has low mapping quality.

    Parameters:
        variant (cyvcf2.Variant): Variant record.
        threshold (int): Minimum acceptable MQ value.

    Returns:
        bool: True if MQ is below threshold, False otherwise.
    """
    mq = variant.INFO.get("MQ")
    return mq is not None and mq < threshold


def stratify_by_impact(vcf):
    """Group variants by functional impact for QC visualization.

    Parameters:
        vcf (cyvcf2.VCF): VCF object with functional annotations (ANN field).

    Returns:
        tuple:
            - impact_ab (dict[str, list[float]]): Allele balance by impact.
            - impact_dp (dict[str, list[int]]): Read depth by impact.
    """
    impact_ab = defaultdict(list)
    impact_dp = defaultdict(list)

    for var in vcf:
        ann_field = var.INFO.get("ANN")
        if not ann_field:
            continue
        annotations = ann_field.split(",")
        impact = annotations[0].split("|")[2]  # IMPACT field

        # Safely flatten DP and AD arrays
        dp_raw = var.format("DP")
        ad_raw = var.format("AD")

        if dp_raw is None or ad_raw is None:
            continue

        dp = dp_raw[0] if isinstance(dp_raw[0], (int, float)) else dp_raw[0][0]
        ad = ad_raw[0] if isinstance(ad_raw[0], (int, float)) else ad_raw[0]

        ab = calculate_ab(ad)
        impact_ab[impact].append(ab)
        impact_dp[impact].append(dp)

    return impact_ab, impact_dp


def plot_depth_and_ab(vcf, out_dir):
    """Plot histograms of read depth and allele balance across variants.

    Parameters:
        vcf (cyvcf2.VCF or iterable): Variants to plot.
        out_dir (str): Directory to save histogram images.

    Returns:
        None
    """
    depths = []
    balances = []

    for v in vcf:
        dp_array = v.format("DP")  # Usually shape (1, 1) or (1,)
        ad_array = v.format("AD")  # Usually shape (1, 2) or more

        if dp_array is not None and len(dp_array) > 0:
            dp = dp_array[0]
            if isinstance(dp, (list, tuple)):
                dp = dp[0]
            depths.append(dp)

        if ad_array is not None and len(ad_array[0]) >= 2:
            ad = ad_array[0]
            ab = calculate_ab(ad)
            balances.append(ab)

    # Ensure depths and balances are flat scalars
    depths = [int(d[0]) if isinstance(d, (list, tuple)) else int(d) for d in depths]
    balances = [float(b) for b in balances]

    # Plot Read Depth Histogram
    plt.figure(figsize=(8, 4))
    plt.hist(depths, bins=range(0, 600, 20), color="gray", edgecolor="black")
    plt.title("Read Depth Distribution")
    plt.xlabel("DP")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "read_depth_distribution.png"))
    plt.close()

    # Plot Allele Balance Histogram
    plt.figure(figsize=(8, 4))
    plt.hist(balances, bins=30, color="teal", edgecolor="black")
    plt.title("Allele Balance Distribution")
    plt.xlabel("AB")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "allele_balance_distribution.png"))
    plt.close()


def plot_impact_qc(impact_ab, impact_dp, out_dir):
    """Generate boxplots of read depth and allele balance grouped by impact class.

    Parameters:
        impact_ab (dict[str, list[float]]): AB values stratified by impact.
        impact_dp (dict[str, list[int]]): DP values stratified by impact.
        out_dir (str): Directory to save the plots.

    Returns:
        None
    """
    # Filter non-empty impact groups for depth
    dp_items = sorted((k, v) for k, v in impact_dp.items() if v)
    ab_items = sorted((k, v) for k, v in impact_ab.items() if v)
    # Unzip
    dp_labels, dp_data = zip(*dp_items)
    ab_labels, ab_data = zip(*ab_items)

    # Plot Read Depth
    plt.figure(figsize=(6, 4))
    plt.boxplot(dp_data, labels=dp_labels, patch_artist=True)
    plt.title("Read Depth by Variant Impact")
    plt.ylabel("DP")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "impact_vs_depth.png"))
    plt.close()

    # Plot Allele Balance
    plt.figure(figsize=(6, 4))
    plt.boxplot(ab_data, labels=ab_labels, patch_artist=True)
    plt.title("Allele Balance by Variant Impact")
    plt.ylabel("AB")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "impact_vs_ab.png"))
    plt.close()
