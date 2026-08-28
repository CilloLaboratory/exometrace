process MARK_DUPLICATES {
    label 'cpu_medium'
    publishDir "${params.results_dir}/alignment", mode: 'copy'
    container { params.gatk_container }

    input:
    tuple val(meta), path(bam), path(bai), path(bait_bed), val(reference_config)

    output:
    tuple val(meta), path("${meta.sample_id}.markdup.bam"), path("${meta.sample_id}.markdup.bam.bai"), path("${meta.sample_id}.duplicate_metrics.txt"), path(bait_bed), val(reference_config)

    script:
    def gatk_heap_gb = Math.max(2, task.memory.toGiga().intValue() - 4)
    """
    gatk --java-options "-Xms1g -Xmx${gatk_heap_gb}g" MarkDuplicates \
      -I ${bam} \
      -O ${meta.sample_id}.markdup.bam \
      -M ${meta.sample_id}.duplicate_metrics.txt \
      --CREATE_INDEX true
    if [[ ! -f ${meta.sample_id}.markdup.bam.bai && -f ${meta.sample_id}.markdup.bai ]]; then
      ln -sf ${meta.sample_id}.markdup.bai ${meta.sample_id}.markdup.bam.bai
    fi
    samtools quickcheck ${meta.sample_id}.markdup.bam
    """

    stub:
    """
    touch ${meta.sample_id}.markdup.bam ${meta.sample_id}.markdup.bam.bai
    cat > ${meta.sample_id}.duplicate_metrics.txt <<'EOF'
    ## METRICS CLASS	picard.sam.DuplicationMetrics
    LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	SECONDARY_OR_SUPPLEMENTARY_RDS	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE
    ${meta.sample_id}	0	100	0	0	0	10	0	0.100000	1000
    EOF
    """
}
