from variant_focus.converter import gtf_to_bed, filter_bed_to_vcf
from variant_focus.normalise import validate_ref_alleles
from variant_focus.annotate import remove_vcf_fields
from variant_focus.sort import sort_vcf
from variant_focus.index import index_vcf

def focus_vcf(gtf_file,vcf_input,reference,vcf_output,bed_output):
    # Step 2: Convert GTF to BED (CDS only)
    bed_file = gtf_to_bed(gtf_file, bed_output=bed_output)

    # Step 3: Filter VCF by CDS regions
    filtered_vcf = filter_bed_to_vcf(vcf_input, bed_file)

    # Step 4: Validate reference alleles
    validated_vcf = validate_ref_alleles(filtered_vcf, reference)
    filtered_vcf.unlink()

    # Step 5: Remove unnecessary INFO/FORMAT fields
    cleaned_vcf = remove_vcf_fields(validated_vcf, "INFO/OLD_VARIANT,FORMAT/DP4")
    validated_vcf.unlink()

    # Step 6: Sort and compress
    final_vcf = sort_vcf(cleaned_vcf, vcf_output)

    # Step 7: Index
    index_vcf(final_vcf)

    return final_vcf