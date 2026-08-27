process COLLECT_VERSIONS {
    label 'cpu_small'
    publishDir "${params.results_dir}/provenance", mode: 'copy'

    output:
    path 'software_versions.tsv'

    script:
    """
    bash ${projectDir}/scripts/collect_versions.sh software_versions.tsv
    """
}
