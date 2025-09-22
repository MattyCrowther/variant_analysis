from collections import Counter
from collections import defaultdict


def variant_key(annotation):
    """
    Return a key uniquely identifying a variant (chrom, pos, ref, alt).
    """
    return (
        annotation["Chrom"],
        annotation["Pos"],
        annotation["Ref"],
        annotation["Alt"],
    )

def summarize_variant_consequences(annotations):
    """
    Count the number of times each consequence appears.
    """
    return Counter(a["Annotation"] for a in annotations)

def count_variants_by_impact_per_gene(annotations, impact_level="HIGH"):
    """
    Return a count of variants per gene with the given impact level.
    """
    counter = Counter()
    for a in annotations:
        if a["Impact"] == impact_level:
            counter[a["Gene_Name"]] += 1
    return counter

def write_gene_list_to_file(gene_list, output_path):
    """
    Write a list of genes to a file.
    """
    with open(output_path, "w") as f:
        for gene in sorted(gene_list):
            f.write(gene + "\n")

def filter_annotations(annotations, field, allowed_values):
    """
    Generic filter: keep annotations where `annotation[field]` is in allowed_values.
    """
    return [a for a in annotations if a.get(field) in allowed_values]

def filter_by_impact(annotations, impact_levels={"HIGH", "MODERATE"}):
    """
    Return annotations with an impact level in the given set.
    """
    return filter_annotations(annotations, "Impact", impact_levels)

def filter_by_consequence(annotations, consequences):
    """
    Return annotations with a consequence (effect) in the given set.
    """
    return filter_annotations(annotations, "Annotation", consequences)

def group_annotations_by_gene(annotations):
    """
    Group annotations into a dict keyed by gene name.
    """
    gene_to_annots = defaultdict(list)
    for a in annotations:
        gene_to_annots[a["Gene_Name"]].append(a)
    return gene_to_annots

def top_genes_by_variant_count(annotations, n=10):
    """
    Return the top N genes with the most annotations.
    Groups annotations by gene name and ranks them by count.
    """
    gene_to_annots = defaultdict(list)
    for a in annotations:
        gene_to_annots[a["Gene_Name"]].append(a)

    return sorted(gene_to_annots.items(), key=lambda x: len(x[1]), reverse=True)[:n]

def count_multi_transcript_variants(annotations):
    """
    Count how many unique variant-allele pairs appear in multiple transcripts.
    """
    variant_to_transcripts = defaultdict(set)

    for a in annotations:
        variant_to_transcripts[variant_key(a)].add(a["Feature_ID"])

    return sum(1 for t in variant_to_transcripts.values() if len(t) > 1)

def collapse_to_most_severe_annotation_per_variant(annotations, coding_only=False):
    """
    Collapse transcript-level annotations by selecting the most severe one per variant.

    Args:
        annotations (List[dict]): List of annotations parsed from ANN fields.
        coding_only (bool): If True, restrict to protein-coding transcripts.

    Returns:
        List[dict]: One annotation per variant (most severe).
    """
    impact_rank = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}

    def variant_key(a):
        return (a["Chrom"], a["Pos"], a["Ref"], a["Alt"])

    # Optionally restrict to coding transcripts
    if coding_only:
        annotations = [
            a for a in annotations if a.get("Transcript_BioType") == "protein_coding"
        ]

    # Group annotations by variant
    variant_to_annots = defaultdict(list)
    for ann in annotations:
        key = variant_key(ann)
        variant_to_annots[key].append(ann)

    # Keep the most severe annotation for each variant
    collapsed = []
    for ann_list in variant_to_annots.values():
        best = max(ann_list, key=lambda a: impact_rank.get(a.get("Impact", ""), -1))
        collapsed.append(best)

    return collapsed

def classify_protein_effects(annotations):
    """
    Classify protein-level consequences using HGVS.p field.

    Returns:
        Counter: Counts of 'missense', 'nonsense', 'synonymous', 'unknown', 'other'
    """
    effect_types = []

    for ann in annotations:
        p = ann.get("HGVS.p", "")

        if not p or p == "p.?":
            effect = "unknown"
        elif "Ter" in p or "*" in p:
            effect = "nonsense"
        elif "=" in p:
            effect = "synonymous"
        elif p.startswith("p.") and len(p) > 5:
            effect = "missense"
        else:
            effect = "other"

        effect_types.append(effect)

    return Counter(effect_types)

def extract_amino_acid_positions(annotations):
    """
    Extract numeric amino acid positions from HGVS.p for positional analysis.

    Returns:
        List[int]: Parsed amino acid positions (if available)
    """
    positions = []

    for ann in annotations:
        p = ann.get("HGVS.p", "")

        if p.startswith("p.") and len(p) > 5:
            digits = "".join(c for c in p if c.isdigit())
            if digits:
                positions.append(int(digits))

    return positions

def histogram_amino_acid_regions(positions, bin_size=100):
    """
    Summarize positional distribution of variants across protein length bins.

    Args:
        positions (List[int])
        bin_size (int): Amino acid window size for bucketing

    Returns:
        dict: Region buckets and variant counts
    """
    region_counts = Counter((pos // bin_size) * bin_size for pos in positions)
    return dict(sorted(region_counts.items()))

def filter_protein_altering_variants(annotations):
    """
    Filter annotations to retain only protein-altering variants.

    Returns:
        List[dict]: Subset of annotations with coding effects.
    """
    relevant = {"missense_variant", "stop_gained", "start_lost", "frameshift_variant"}

    return [a for a in annotations if a.get("Annotation") in relevant]
