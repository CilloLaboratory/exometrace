process SOMATIC_COMPARE {
    label 'cpu_small'
    container { params.qc_container }
    publishDir "${params.results_dir}/somatic/comparison", mode: 'copy'

    input:
    tuple val(meta), path(deepsomatic_vcf), path(mutect2_vcf)

    output:
    tuple val(meta), path("${meta.patient_id}.variant_concordance.tsv"), path("${meta.patient_id}.variant_counts.tsv")

    script:
    """
    python3 ${projectDir}/scripts/compare_somatic_vcfs.py \
      --deepsomatic ${deepsomatic_vcf} \
      --mutect2 ${mutect2_vcf} \
      --output-table ${meta.patient_id}.variant_concordance.tsv \
      --output-counts ${meta.patient_id}.variant_counts.tsv
    """

    stub:
    """
    python3 ${projectDir}/scripts/compare_somatic_vcfs.py \
      --deepsomatic ${deepsomatic_vcf} \
      --mutect2 ${mutect2_vcf} \
      --output-table ${meta.patient_id}.variant_concordance.tsv \
      --output-counts ${meta.patient_id}.variant_counts.tsv
    """
}
