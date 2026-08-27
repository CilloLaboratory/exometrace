process CLONALITY {
    label 'cpu_small'
    container { params.qc_container }
    publishDir "${params.results_dir}/clonality", mode: 'copy'

    input:
    tuple val(meta), path(maf), path(purity_tsv), path(allele_specific_tsv), val(pipeline_config)

    output:
    tuple val(meta), path("${meta.patient_id}.ccf.tsv"), path("${meta.patient_id}.clonality_summary.tsv")

    script:
    """
    python3 ${projectDir}/scripts/calculate_ccf.py \
      --patient-id ${meta.patient_id} \
      --maf ${maf} \
      --purity ${purity_tsv} \
      --segments ${allele_specific_tsv} \
      --pipeline-config ${pipeline_config} \
      --output-ccf ${meta.patient_id}.ccf.tsv \
      --output-summary ${meta.patient_id}.clonality_summary.tsv
    """

    stub:
    """
    python3 ${projectDir}/scripts/calculate_ccf.py \
      --patient-id ${meta.patient_id} \
      --maf ${maf} \
      --purity ${purity_tsv} \
      --segments ${allele_specific_tsv} \
      --pipeline-config ${pipeline_config} \
      --output-ccf ${meta.patient_id}.ccf.tsv \
      --output-summary ${meta.patient_id}.clonality_summary.tsv
    """
}
