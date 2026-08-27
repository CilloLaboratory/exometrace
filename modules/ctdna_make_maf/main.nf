process CTDNA_VCF_TO_MAF {
    label 'cpu_small'
    publishDir "${params.results_dir}/ctdna/cohort", mode: 'copy'
    container { params.qc_container }

    input:
    tuple val(meta), path(annotated_vcf), path(comparison_tsv), path(umi_qc_tsv), val(min_family_support)

    output:
    tuple val(meta), path("${meta.patient_id}.ctdna_mutations_high_sensitivity.maf.tsv"), path("${meta.patient_id}.ctdna_mutations_high_confidence.maf.tsv")

    script:
    """
    python3 ${projectDir}/scripts/ctdna_vcf_to_maf.py \
      --input ${annotated_vcf} \
      --comparison ${comparison_tsv} \
      --umi-qc ${umi_qc_tsv} \
      --min-family-support ${min_family_support} \
      --high-sensitivity-output ${meta.patient_id}.ctdna_mutations_high_sensitivity.maf.tsv \
      --high-confidence-output ${meta.patient_id}.ctdna_mutations_high_confidence.maf.tsv
    """

    stub:
    """
    python3 ${projectDir}/scripts/ctdna_vcf_to_maf.py \
      --input ${annotated_vcf} \
      --comparison ${comparison_tsv} \
      --umi-qc ${umi_qc_tsv} \
      --min-family-support ${min_family_support} \
      --high-sensitivity-output ${meta.patient_id}.ctdna_mutations_high_sensitivity.maf.tsv \
      --high-confidence-output ${meta.patient_id}.ctdna_mutations_high_confidence.maf.tsv
    """
}
