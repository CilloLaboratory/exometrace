process VALIDATE_SAMPLESHEET {
    label 'cpu_small'
    publishDir "${params.results_dir}/provenance", mode: 'copy'

    input:
    val samplesheet

    output:
    path 'validated.samplesheet.csv'

    script:
    """
    python3 ${projectDir}/scripts/validate_samplesheet.py \
      --samplesheet ${samplesheet} \
      --check-files \
      --emit-normalized validated.samplesheet.csv
    """
}
