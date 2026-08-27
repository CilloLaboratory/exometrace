process MERGE_CTDNA_TABLES {
    label 'cpu_small'
    container { params.qc_container }
    publishDir "${params.results_dir}/ctdna/cohort", mode: 'copy'

    input:
    path tables
    val output_name
    val union_by_first_column

    output:
    path "${output_name}"

    script:
    """
    python3 ${projectDir}/scripts/merge_tables.py \
      --inputs ${tables.join(' ')} \
      ${union_by_first_column ? '--union-by-first-column' : ''} \
      --output ${output_name}
    """
}

process MERGE_CTDNA_QC_TABLES {
    label 'cpu_small'
    container { params.qc_container }
    publishDir "${params.results_dir}/ctdna/qc", mode: 'copy'

    input:
    path tables
    val output_name
    val union_by_first_column

    output:
    path "${output_name}"

    script:
    """
    python3 ${projectDir}/scripts/merge_tables.py \
      --inputs ${tables.join(' ')} \
      ${union_by_first_column ? '--union-by-first-column' : ''} \
      --output ${output_name}
    """
}
