process GENOMIC_FEATURES {
    label 'cpu_small'
    container { params.qc_container }
    publishDir "${params.results_dir}/cohort", mode: 'copy'

    input:
    path tmb_tables
    path purity_tables
    path signature_tables
    path clonality_tables
    path driver_matrices
    path arm_matrices

    output:
    path "genomic_features.tsv"

    script:
    """
    python3 ${projectDir}/scripts/merge_features.py \
      --tmb ${tmb_tables.join(' ')} \
      --purity ${purity_tables.join(' ')} \
      --signatures ${signature_tables.join(' ')} \
      --clonality ${clonality_tables.join(' ')} \
      --drivers ${driver_matrices.join(' ')} \
      --arm-matrices ${arm_matrices.join(' ')} \
      --output genomic_features.tsv
    """
}
