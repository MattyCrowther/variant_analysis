import gzip

def load_phenotype_genes(path):
    """Load phenotype-associated genes as a set of STANDARD names."""
    genes = set()
    with open(path) as f:
        for line in f:
            if line.startswith("Feature") or line.strip() == "":
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                gene = parts[0]  # Standard gene name
                if gene:
                    genes.add(gene)
            else:
                print(f"Skipping malformed line: {line.strip()}")
    return genes


def load_gene_name_map(mapping_path):
    """Returns a dict: SYSTEMATIC -> STANDARD"""
    name_map = {}
    with open(mapping_path) as f:
        for line in f:
            if line.startswith("Feature") or line.strip() == "":
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                systematic = parts[0]
                standard = parts[3]
                if standard:
                    name_map[systematic] = standard
    return name_map


def annotate_with_phenotype_tags(vcf_in, vcf_out, known_genes, name_map):
    """Tag variants with PHENO_HIT=YES if any ANN genes match known phenotype-associated genes."""
    with gzip.open(vcf_in, "rt") as fin, gzip.open(vcf_out, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            fields = line.strip().split("\t")
            info = fields[7]
            ann = next((x for x in info.split(";") if x.startswith("ANN=")), "")
            tag = "PHENO_HIT=NO"

            if ann:
                entries = ann.replace("ANN=", "").split(",")
                genes = {
                    name_map.get(entry.split("|")[3], None)
                    for entry in entries
                }
                genes = {g for g in genes if g is not None}
                if known_genes & genes:
                    tag = "PHENO_HIT=YES"

            fields[7] = info + ";" + tag
            fout.write("\t".join(fields) + "\n")
