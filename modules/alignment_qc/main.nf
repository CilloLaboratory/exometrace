process ALIGNMENT_QC {
    label 'cpu_medium'
    publishDir "${params.results_dir}/qc/alignment", mode: 'copy'
    publishDir "${params.results_dir}/qc/coverage", mode: 'copy', pattern: "*.mosdepth*"
    container { params.qc_container }

    input:
    tuple val(meta), path(bam), path(bai), path(dup_metrics), path(bait_bed), val(reference_config), path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict), path(target_intervals)

    output:
    tuple val(meta), path("${meta.sample_id}.qc.tsv")

    script:
    def gatk_heap_gb = Math.max(2, task.memory.toGiga().intValue() - 4)
    """
    samtools flagstat ${bam} > ${meta.sample_id}.flagstat.txt
    samtools stats ${bam} > ${meta.sample_id}.stats.txt
    export JAVA_TOOL_OPTIONS="-Xms1g -Xmx${gatk_heap_gb}g"
    gatk --java-options "-Xms1g -Xmx${gatk_heap_gb}g" CollectAlignmentSummaryMetrics \
      -R ${ref_fasta} \
      -I ${bam} \
      -O ${meta.sample_id}.alignment_summary_metrics.txt
    gatk --java-options "-Xms1g -Xmx${gatk_heap_gb}g" CollectInsertSizeMetrics \
      -I ${bam} \
      -O ${meta.sample_id}.insert_size_metrics.txt \
      -H ${meta.sample_id}.insert_size_histogram.pdf
    if [[ ! -f ${meta.sample_id}.insert_size_metrics.txt ]]; then
      cat > ${meta.sample_id}.insert_size_metrics.txt <<'EOF'
MEDIAN_INSERT_SIZE
NA
EOF
    fi
    gatk --java-options "-Xms1g -Xmx${gatk_heap_gb}g" CollectHsMetrics \
      -R ${ref_fasta} \
      -I ${bam} \
      -BAIT_INTERVALS ${target_intervals} \
      -TARGET_INTERVALS ${target_intervals} \
      -O ${meta.sample_id}.hs_metrics.txt
    mosdepth --by ${bait_bed} --threads ${task.cpus} ${meta.sample_id}.mosdepth ${bam}
    python3 ${projectDir}/scripts/summarize_alignment_qc.py \
      --patient-id ${meta.patient_id} \
      --sample ${meta.sample_id} \
      --sample-type ${meta.sample_type} \
      --flagstat ${meta.sample_id}.flagstat.txt \
      --duplication-metrics ${dup_metrics} \
      --insert-metrics ${meta.sample_id}.insert_size_metrics.txt \
      --hs-metrics ${meta.sample_id}.hs_metrics.txt \
      --output ${meta.sample_id}.qc.tsv
    """

    stub:
    """
    touch ${meta.sample_id}.mosdepth.regions.bed.gz ${meta.sample_id}.mosdepth.summary.txt
    cat > ${meta.sample_id}.qc.tsv <<'EOF'
    patient_id\tsample\tsample_type\ttotal_reads\tmapped_reads\tmapping_rate\tproper_pair_rate\tduplicate_rate\tmean_target_coverage\tmedian_target_coverage\tpct_target_10x\tpct_target_20x\tpct_target_30x\tpct_target_50x\tpct_target_100x\tmean_insert_size
    ${meta.patient_id}\t${meta.sample_id}\t${meta.sample_type}\t100\t98\t0.980000\t0.950000\t0.100000\t120.000000\t110.000000\t0.990000\t0.970000\t0.950000\t0.900000\t0.800000\t250.000000
    EOF
    """
}
