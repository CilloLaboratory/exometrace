process TMB {
    label 'cpu_small'
    container { params.qc_container }
    publishDir "${params.results_dir}/cohort", mode: 'copy'

    input:
    tuple val(meta), path(maf), path(callable_mb)

    output:
    tuple val(meta), path("${meta.patient_id}.tmb.tsv")

    script:
    """
    python3 ${projectDir}/scripts/calculate_tmb.py \
      --maf ${maf} \
      --callable-mb ${callable_mb} \
      --patient-id ${meta.patient_id} \
      --tumor-sample ${meta.tumor_sample} \
      --output ${meta.patient_id}.tmb.tsv
    """

    stub:
    """
    python3 ${projectDir}/scripts/calculate_tmb.py \
      --maf ${maf} \
      --callable-mb ${callable_mb} \
      --patient-id ${meta.patient_id} \
      --tumor-sample ${meta.tumor_sample} \
      --output ${meta.patient_id}.tmb.tsv
    """
}
