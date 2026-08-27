process ALIGN_BWA_MEM2 {
    label 'cpu_large'
    publishDir "${params.results_dir}/alignment", mode: 'copy'
    container { params.align_container }

    input:
    tuple val(meta), path(r1), path(r2), path(bait_bed), val(reference_config), val(ref_fasta)

    output:
    tuple val(meta), path("${meta.sample_id}.sorted.bam"), path("${meta.sample_id}.sorted.bam.bai"), path(bait_bed), val(reference_config)

    script:
    """
    SORT_THREADS=${task.cpus > 8 ? 8 : task.cpus}
    ${params.bwa_binary} mem \
      -t ${task.cpus} \
      -R '@RG\\tID:${meta.sample_id}\\tSM:${meta.sample_id}\\tLB:${meta.sample_id}_lib1\\tPL:ILLUMINA\\tPU:${meta.patient_id}' \
      ${ref_fasta} \
      ${r1} \
      ${r2} \
    | samtools sort -@ "\$SORT_THREADS" -o ${meta.sample_id}.sorted.bam
    samtools index ${meta.sample_id}.sorted.bam
    samtools quickcheck ${meta.sample_id}.sorted.bam
    """

    stub:
    """
    touch ${meta.sample_id}.sorted.bam ${meta.sample_id}.sorted.bam.bai
    """
}
