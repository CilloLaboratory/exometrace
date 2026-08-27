process ARM_LEVEL_CNV {
    label 'cpu_small'
    container { params.qc_container }
    publishDir "${params.results_dir}/cohort", mode: 'copy', pattern: "*.arm_level_cnv_long.tsv"
    publishDir "${params.results_dir}/cohort", mode: 'copy', pattern: "*.arm_level_cnv_matrix.tsv"

    input:
    tuple val(meta), path(call_cns), val(reference_config), val(arms_bed)

    output:
    tuple val(meta), path("${meta.patient_id}.arm_level_cnv_long.tsv"), path("${meta.patient_id}.arm_level_cnv_matrix.tsv")

    script:
    """
    python3 ${projectDir}/scripts/summarize_cnv.py \
      --patient-id ${meta.patient_id} \
      --segments ${call_cns} \
      --arms-bed ${arms_bed} \
      --output-long ${meta.patient_id}.arm_level_cnv_long.tsv \
      --output-matrix ${meta.patient_id}.arm_level_cnv_matrix.tsv \
      --deep-deletion-log2 ${params.copy_number_arm_deep_del ?: -1.1} \
      --loss-log2 ${params.copy_number_arm_loss ?: -0.3} \
      --gain-log2 ${params.copy_number_arm_gain ?: 0.2} \
      --amplification-log2 ${params.copy_number_arm_amp ?: 0.7}
    """

    stub:
    """
    python3 ${projectDir}/scripts/summarize_cnv.py \
      --patient-id ${meta.patient_id} \
      --segments ${call_cns} \
      --arms-bed ${arms_bed} \
      --output-long ${meta.patient_id}.arm_level_cnv_long.tsv \
      --output-matrix ${meta.patient_id}.arm_level_cnv_matrix.tsv
    """
}
