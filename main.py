import os
from pathlib import Path

from downloader import download
from utility import decompress_gzip
from impact_scoring.run import run as impact_scoring_run
from annotation.snpeff import run as annotate_run
from technical_reliability.run import run as tech_run
from impact_scoring.run import run as impact_run
from annotation_frequency.run import run as frequency_run
VCF_DIR = Path("vcf")
REF_DIR = Path("ref")
STORAGE_DIR = Path("storage")

vcf_input = STORAGE_DIR/VCF_DIR/"variants_raw.vcf.gz"

fasta_ref = STORAGE_DIR/REF_DIR / "Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"
gtf_url = "https://ftp.ensembl.org/pub/release-109/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.109.gtf.gz"
gtf_gz_path = STORAGE_DIR / "Saccharomyces_cerevisiae.R64-1-1.109.gtf.gz"
bed_output_path = STORAGE_DIR / "Saccharomyces_cerevisiae.CDS.bed"
final_vcf_path = STORAGE_DIR / "yeast_final.vcf.gz"

def download_gtp(gtf_url, gtf_gz_path):
    if not os.path.isfile(gtf_gz_path):
        gtf_gz_path = download(gtf_url, gtf_gz_path)
    return decompress_gzip(gtf_gz_path)


'''

key = "SCEREVISIAE_YEAST"
genome_label = "Saccharomyces_cerevisiae_R64-1-1"
annotate_run(vcf_input,gtf_gz_path,fasta_ref,annotated_vcf,key,genome_label)
'''

'''
impact_db = Path("impact_scoring/sift4g_db/R64-1-1.23")
annotated_vcf = STORAGE_DIR/VCF_DIR/"technical_filter_vcf.vcf.gz"
technical_filter_vcf = STORAGE_DIR/VCF_DIR/"impact_prio_vcf.vcf"
write_tsv="sift_scores.tsv"
write_filtered_vcf="storage/vcf/damaging_only.vcf.gz"
impact_run(annotated_vcf,technical_filter_vcf,impact_db,
           write_tsv=write_tsv,write_vcf=write_filtered_vcf)
'''

output_path = Path("storage/vcf/reference/saccharomyces_cerevisiae.vcf.gz")
qry_vcf = Path("storage/vcf/technical_filter_vcf.vcf.gz")
url = "https://ftp.ensembl.org/pub/release-109/variation/vcf/saccharomyces_cerevisiae/saccharomyces_cerevisiae.vcf.gz"
phenotype_url = "http://sgd-archive.yeastgenome.org/curation/literature/phenotype_data.tab"
frequency_run(url,output_path,qry_vcf,phenotype_url)



