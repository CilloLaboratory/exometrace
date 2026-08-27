process MERGE_TABLES {
    label 'cpu_small'
    container { params.qc_container }
    publishDir "${params.results_dir}/cohort", mode: 'copy'

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
