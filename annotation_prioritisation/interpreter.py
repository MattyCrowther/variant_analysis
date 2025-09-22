import gzip


def open_vcf(filepath):
    """
    Open a VCF file that may be gzipped or plain text.
    """
    if filepath.endswith(".gz"):
        return gzip.open(filepath, "rt")
    return open(filepath, "r")


def extract_snpeff_annotations(vcf_path):
    """
    Parse annotated VCF and extract functional consequence and gene name from SnpEff 'ANN' field.

    Args:
        vcf_path (str): Path to a .vcf or .vcf.gz file annotated with SnpEff

    Returns:
        List of (consequence, gene_name) tuples
    """
    annotations = []

    with open_vcf(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue  # Skip headers

            fields = line.strip().split("\t")
            info_field = fields[7]

            for entry in info_field.split(";"):
                if entry.startswith("ANN="):
                    ann_values = entry[len("ANN=") :].split(",")
                    for ann in ann_values:
                        ann_fields = ann.split("|")
                        if len(ann_fields) >= 4:
                            consequence = ann_fields[1]
                            gene_name = ann_fields[3]
                            annotations.append((consequence, gene_name))
                    break  # One ANN per INFO field

    return annotations


def parse_annotated_vcf(vcf_path):
    """
    Parse a SnpEff-annotated VCF file (optionally gzipped) and extract ANN fields.

    Args:
        vcf_path (str): Path to a .vcf or .vcf.gz file.

    Returns:
        List[dict]: Parsed annotations. Each dict includes keys like:
            'Allele', 'Annotation', 'Impact', 'Gene_Name', 'Gene_ID', etc.,
            plus: 'Chrom', 'Pos', 'Ref', 'Alt'
    """
    import gzip

    ann_fields = [
        "Allele",
        "Annotation",
        "Impact",
        "Gene_Name",
        "Gene_ID",
        "Feature_Type",
        "Feature_ID",
        "Transcript_BioType",
        "Rank",
        "HGVS.c",
        "HGVS.p",
        "cDNA.pos / cDNA.length",
        "CDS.pos / CDS.length",
        "AA.pos / AA.length",
        "Distance",
        "Errors/Warnings/Info",
    ]

    annotations = []
    open_func = gzip.open if vcf_path.endswith(".gz") else open

    with open_func(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            info = fields[7]

            ann_entries = next(
                (entry[4:] for entry in info.split(";") if entry.startswith("ANN=")),
                None,
            )

            if ann_entries:
                for raw_ann in ann_entries.split(","):
                    parts = raw_ann.split("|")
                    ann_record = {
                        key: parts[i] if i < len(parts) else None
                        for i, key in enumerate(ann_fields)
                    }

                    # Add positional + allele data
                    ann_record.update({
                        "Chrom": chrom,
                        "Pos": pos,
                        "Ref": ref,
                        "Alt": alt,
                    })

                    annotations.append(ann_record)

    return annotations

