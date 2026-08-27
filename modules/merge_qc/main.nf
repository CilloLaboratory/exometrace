process MERGE_SAMPLE_QC {
    label 'cpu_small'
    container { params.qc_container }
    publishDir "${params.results_dir}/cohort", mode: 'copy'

    input:
    path qc_files

    output:
    path 'sample_qc.tsv'

    script:
    """
    python3 ${projectDir}/scripts/merge_sample_qc.py \
      --inputs ${qc_files.join(' ')} \
      --output sample_qc.tsv
    """

    stub:
    """
    cp ${qc_files[0]} sample_qc.tsv
    """
}
