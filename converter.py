import subprocess
from pathlib import Path


def gtf_to_bed(
    gtf_path: Path, feature_type: str = "CDS", bed_output: Path = None
) -> Path:
    """
    Extract a specified feature type from a GTF file and convert to BED format.

    Parameters:
        gtf_path (Path): Path to the input GTF file.
        feature_type (str): GTF feature to extract (e.g., "CDS", "exon", "gene").
        bed_output (Path, optional): Output BED file path. If not provided,
                                     it will be generated automatically.

    Returns:
        Path: Path to the generated BED file.
    """
    # Generate output path if not provided
    if bed_output is None:
        out_stem = gtf_path.stem.split(".")[0]  # Just the base name
        bed_output = gtf_path.parent / f"{out_stem}.{feature_type}.bed"

    bed_output.parent.mkdir(parents=True, exist_ok=True)

    # AWK command to extract the desired feature
    awk_cmd = f"""awk '$3 == "{feature_type}" {{ print $1 "\\t" $4-1 "\\t" $5 }}'"""

    with bed_output.open("w") as out:
        gffread = subprocess.Popen(
            ["gffread", str(gtf_path), "-T", "-o-"], stdout=subprocess.PIPE
        )
        awk = subprocess.Popen(
            ["bash", "-c", awk_cmd], stdin=gffread.stdout, stdout=subprocess.PIPE
        )
        sort = subprocess.Popen(
            ["sort", "-k1,1", "-k2,2n"], stdin=awk.stdout, stdout=out
        )
        sort.communicate()

    return bed_output


def filter_bed_to_vcf(vcf_input: Path, bed_file: Path, vcf_output: Path = None) -> Path:
    """
    Subset a VCF file to regions defined in a BED file using bcftools.

    Parameters:
        vcf_input (Path): Path to the input VCF (.vcf.gz, must be indexed).
        bed_file (Path): Path to BED file defining regions to retain.
        vcf_output (Path, optional): Output path for filtered VCF. If not provided,
                                     a name is generated based on the inputs.

    Returns:
        Path: Path to the filtered and indexed VCF.
    """
    if vcf_output is None:
        base = vcf_input.stem.replace(".vcf", "")
        region = bed_file.stem
        vcf_output = vcf_input.parent / f"{base}_in_{region}.vcf.gz"

    vcf_input = Path(vcf_input)
    bed_file = Path(bed_file)
    vcf_output.parent.mkdir(parents=True, exist_ok=True)

    # Run bcftools filtering
    subprocess.run(
        [
            "bcftools",
            "view",
            "-R",
            str(bed_file),
            "-Oz",
            "-o",
            str(vcf_output),
            str(vcf_input),
        ],
        check=True,
    )

    # Index the output VCF
    subprocess.run(["tabix", "-p", "vcf", str(vcf_output)], check=True)

    return vcf_output
