process VALIDATE_REFERENCE {
    label 'cpu_small'
    publishDir "${params.results_dir}/provenance", mode: 'copy'

    input:
    tuple val(reference_config), val(pipeline_config)

    output:
    path 'reference_manifest.tsv'

    script:
    """
    python3 ${projectDir}/scripts/validate_reference.py \
      --config ${reference_config} \
      --pipeline-config ${pipeline_config} \
      --emit-manifest reference_manifest.tsv
    """
}
