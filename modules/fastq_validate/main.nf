process FASTQ_VALIDATE {
    label 'cpu_small'
    cache 'deep'
    publishDir "${params.results_dir}/qc/fastq_validation", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2), path(bait_bed)

    output:
    tuple val(meta), path("${meta.sample_id}.fastq_validation.tsv"), path(r1), path(r2), path(bait_bed)

    script:
    """
    python3 ${projectDir}/scripts/validate_fastq_pairs.py \
      --sample-id ${meta.sample_id} \
      --r1 ${r1} \
      --r2 ${r2} \
      --output ${meta.sample_id}.fastq_validation.tsv
    """

    stub:
    """
    cat > ${meta.sample_id}.fastq_validation.tsv <<'EOF'
    sample\treads_r1\treads_r2\tstatus
    ${meta.sample_id}\t2\t2\tPASS
    EOF
    """
}
