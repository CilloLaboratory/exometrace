process CONTAINER_MANIFEST {
    label 'cpu_small'
    publishDir "${params.results_dir}/provenance", mode: 'copy'

    input:
    val containers_config

    output:
    path "containers.tsv"

    script:
    """
    python3 ${projectDir}/scripts/generate_containers_manifest.py \
      --config ${containers_config} \
      --output containers.tsv
    """
}

process PIPELINE_PARAMETERS {
    label 'cpu_small'
    publishDir "${params.results_dir}/provenance", mode: 'copy'

    input:
    val samplesheet
    val reference_config

    output:
    path "pipeline_parameters.yaml"

    script:
    """
    python3 ${projectDir}/scripts/dump_pipeline_params.py \
      --samplesheet ${samplesheet} \
      --reference-config ${reference_config} \
      --results-dir ${params.results_dir} \
      --output pipeline_parameters.yaml
    """
}

process MULTIQC_STUB {
    label 'cpu_small'
    publishDir "${params.results_dir}/qc/multiqc", mode: 'copy'

    input:
    path sample_qc

    output:
    path "multiqc_report.html"

    script:
    """
    python3 ${projectDir}/scripts/generate_multiqc_stub.py \
      --sample-qc ${sample_qc} \
      --output multiqc_report.html
    """
}

process EXECUTION_REPORT {
    label 'cpu_small'
    publishDir "${params.results_dir}/provenance", mode: 'copy'

    input:
    path sample_qc
    path maf
    path tmb
    path purity
    path signatures
    path clonality
    path drivers
    path arm_matrix
    path features
    path software_versions
    path containers
    path reference_manifest
    path pipeline_parameters

    output:
    path "execution_report.html"

    script:
    """
    python3 ${projectDir}/scripts/generate_execution_report.py \
      --sample-qc ${sample_qc} \
      --maf ${maf} \
      --tmb ${tmb} \
      --purity ${purity} \
      --signatures ${signatures} \
      --clonality ${clonality} \
      --drivers ${drivers} \
      --arm-matrix ${arm_matrix} \
      --features ${features} \
      --software-versions ${software_versions} \
      --containers ${containers} \
      --reference-manifest ${reference_manifest} \
      --pipeline-parameters ${pipeline_parameters} \
      --output execution_report.html
    """
}
