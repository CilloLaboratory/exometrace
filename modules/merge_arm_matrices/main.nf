process MERGE_ARM_MATRICES {
    label 'cpu_small'
    container { params.qc_container }
    publishDir "${params.results_dir}/cohort", mode: 'copy'

    input:
    path matrices

    output:
    path "arm_level_cnv_matrix.tsv"

    script:
    """
    python3 ${projectDir}/scripts/merge_arm_matrices.py \
      --inputs ${matrices.join(' ')} \
      --output arm_level_cnv_matrix.tsv
    """
}
