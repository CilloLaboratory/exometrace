process SIGNATURES {
    label 'cpu_small'
    publishDir "${params.results_dir}/signatures", mode: 'copy'
    container { params.signature_container }

    input:
    tuple val(meta), path(maf), val(reference_config), val(pipeline_config)

    output:
    tuple val(meta), path("${meta.patient_id}.signature_exposures.tsv"), path("${meta.patient_id}.signature_qc.tsv")

    script:
    """
    python3 ${projectDir}/scripts/calculate_signature_exposure.py \
      --patient-id ${meta.patient_id} \
      --maf ${maf} \
      --reference-config ${reference_config} \
      --pipeline-config ${pipeline_config} \
      --cpu ${task.cpus} \
      --output-exposure ${meta.patient_id}.signature_exposures.tsv \
      --output-qc ${meta.patient_id}.signature_qc.tsv
    """

    stub:
    """
    python3 ${projectDir}/scripts/calculate_signature_exposure.py \
      --patient-id ${meta.patient_id} \
      --maf ${maf} \
      --reference-config ${reference_config} \
      --pipeline-config ${pipeline_config} \
      --stub \
      --output-exposure ${meta.patient_id}.signature_exposures.tsv \
      --output-qc ${meta.patient_id}.signature_qc.tsv
    """
}
