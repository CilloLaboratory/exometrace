process DRIVER_ANNOTATION {
    label 'cpu_small'
    container { params.qc_container }
    publishDir "${params.results_dir}/cohort", mode: 'copy'

    input:
    tuple val(meta), path(maf), val(reference_config), path(census), path(hotspots)

    output:
    tuple val(meta), path("${meta.patient_id}.drivers_long.tsv"), path("${meta.patient_id}.driver_matrix.tsv")

    script:
    """
    python3 ${projectDir}/scripts/annotate_drivers.py \
      --maf ${maf} \
      --census ${census} \
      --hotspots ${hotspots} \
      --output-long ${meta.patient_id}.drivers_long.tsv \
      --output-matrix ${meta.patient_id}.driver_matrix.tsv
    """

    stub:
    """
    python3 ${projectDir}/scripts/annotate_drivers.py \
      --maf ${maf} \
      --census ${census} \
      --hotspots ${hotspots} \
      --output-long ${meta.patient_id}.drivers_long.tsv \
      --output-matrix ${meta.patient_id}.driver_matrix.tsv
    """
}
