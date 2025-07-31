from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import os

def plot_consequence_distribution(annotations, output_dir="plots"):
    """
    Save a horizontal bar chart of variant consequence types.
    """
    os.makedirs(output_dir, exist_ok=True)
    consequence_counts = Counter(a["Annotation"] for a in annotations)
    labels, values = zip(*consequence_counts.most_common())

    plt.figure(figsize=(8, 5))
    plt.barh(labels, values)
    plt.xlabel("Variant Count")
    plt.title("Consequence Type Distribution")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "consequence_distribution.png"))
    plt.close()


def plot_impact_distribution(annotations, output_dir="plots"):
    """
    Save a bar chart of SnpEff impact categories.
    """
    os.makedirs(output_dir, exist_ok=True)
    impact_counts = Counter(a["Impact"] for a in annotations)
    labels, values = zip(*impact_counts.items())

    plt.figure(figsize=(6, 4))
    plt.bar(labels, values)
    plt.ylabel("Variant Count")
    plt.title("Impact Level Distribution")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "impact_distribution.png"))
    plt.close()



def plot_top_genes_by_variant_count(annotations, n=10, output_dir="plots"):
    """
    Save a bar chart of the top N genes with most variant annotations.
    """
    os.makedirs(output_dir, exist_ok=True)
    gene_hits = defaultdict(int)
    for a in annotations:
        gene = a["Gene_Name"]
        if gene:
            gene_hits[gene] += 1

    top_genes = sorted(gene_hits.items(), key=lambda x: x[1], reverse=True)[:n]
    labels, values = zip(*top_genes)

    plt.figure(figsize=(8, 5))
    plt.barh(labels, values)
    plt.xlabel("Variant Count")
    plt.title("Top Genes by Variant Hit Count")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "top_genes_by_count.png"))
    plt.close()



def plot_amino_acid_position_distribution(annotations, bins=20, output_dir="plots"):
    """
    Save a histogram of amino acid positions from HGVS.p field.
    """
    os.makedirs(output_dir, exist_ok=True)
    positions = []

    for a in annotations:
        p = a.get("HGVS.p", "")
        if p.startswith("p.") and any(c.isdigit() for c in p):
            digits = "".join(c for c in p if c.isdigit())
            if digits:
                positions.append(int(digits))

    if positions:
        plt.figure(figsize=(8, 4))
        plt.hist(positions, bins=bins)
        plt.xlabel("Amino Acid Position")
        plt.ylabel("Variant Count")
        plt.title("Distribution of Amino Acid Positions")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "amino_acid_position_histogram.png"))
        plt.close()



