process MAKE_MAF {
    label 'cpu_small'
    container { params.qc_container }
    publishDir "${params.results_dir}/cohort", mode: 'copy'

    input:
    tuple val(meta), path(annotated_vcf), path(comparison_tsv)

    output:
    tuple val(meta), path("${meta.patient_id}.somatic_mutations.maf")

    script:
    """
    python3 ${projectDir}/scripts/vep_to_maf.py \
      --input ${annotated_vcf} \
      --comparison ${comparison_tsv} \
      --output ${meta.patient_id}.somatic_mutations.maf
    """

    stub:
    """
    python3 ${projectDir}/scripts/vep_to_maf.py \
      --input ${annotated_vcf} \
      --comparison ${comparison_tsv} \
      --output ${meta.patient_id}.somatic_mutations.maf
    """
}
