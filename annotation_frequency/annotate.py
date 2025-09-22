import gzip

def load_known_variants(vcf_path):
    """Build a set of variant keys from reference VCF: chrom:pos:ref:alt"""
    known = set()
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
            key = f"{chrom}:{pos}:{ref}:{alt}"
            known.add(key)
    return known


def annotate_with_frequency(vcf_in_path, vcf_out_path, known_variants):
    with gzip.open(vcf_in_path, "rt") as fin, gzip.open(vcf_out_path, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.strip().split("\t")
            chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
            key = f"{chrom}:{pos}:{ref}:{alt}"

            info = fields[7]
            tag = "FREQ=SEEN" if key in known_variants else "FREQ=NOVEL"
            fields[7] = info + ";" + tag
            fout.write("\t".join(fields) + "\n")

