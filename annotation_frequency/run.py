from pathlib import Path
import subprocess
import gzip

from annotation_frequency.download import download
from annotation_frequency.annotate import load_known_variants
from annotation_frequency.annotate import annotate_with_frequency

from annotation_frequency.phenotype_tag import load_phenotype_genes
from annotation_frequency.phenotype_tag import annotate_with_phenotype_tags
from annotation_frequency.phenotype_tag import load_gene_name_map

def run(ref_url, reference_vcf, query_vcf, phenotype_url=None):

    download(ref_url, reference_vcf)


    known_variants = load_known_variants(reference_vcf)

    freq_anno_file = Path(str(query_vcf).replace(".vcf", "_freq.vcf.gz"))
    annotate_with_frequency(
        vcf_in_path=query_vcf,
        vcf_out_path=freq_anno_file,
        known_variants=known_variants
    )


    if phenotype_url is None:
        return freq_anno_file


    pheno_out = Path(phenotype_url.split("/")[-1])
    download(phenotype_url, pheno_out)

    known_genes = load_phenotype_genes(pheno_out)
    gene_map = load_gene_name_map("gene_literature.tab")

    pheno_anno_file = Path(str(freq_anno_file).replace("_freq.vcf.gz", "_freq_pheno.vcf"))
    annotate_with_phenotype_tags(freq_anno_file, pheno_anno_file, known_genes,gene_map)

    return pheno_anno_file






    