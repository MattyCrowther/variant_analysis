import shutil
from annotation.annotate import build_snpeff_db
from annotation.annotate import annotate_vcf_with_snpeff


def run(cleaned_vcf,gtf_gz_path,reference,output_vcf,key,genome_label):
    '''
    shutil.copy(
        reference,
        f"snpeff/{key}/sequences.fa"
    )

    shutil.copy(
        gtf_gz_path,
        f"snpeff/{key}/genes.gtf"
    )
    '''

    with open("snpeff.config", "w") as config_file:
        print("writing")
        config_file.write(f"{key}.genome : {reference}\n")

    build_snpeff_db(
        db_dir="snpeff",
        key=key,
        gtf_file=gtf_gz_path,
        reference_file=reference,
        config_path="snpeff.config",
        genome_label=genome_label
    )
    annotate_vcf_with_snpeff(cleaned_vcf,output_vcf,genome_key=key)
    return output_vcf