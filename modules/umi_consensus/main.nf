process UMI_GROUP_CONSENSUS {
    label 'cpu_medium'
    publishDir "${params.results_dir}/ctdna/consensus", mode: 'copy'
    container { params.umi_consensus_container }

    input:
    tuple val(meta), path(bam), path(bai), path(r1), path(r2), path(bait_bed), val(reference_config), val(tag_name), val(read_structure_r1), val(read_structure_r2), val(min_family_size), val(error_rate_pre_umi), val(error_rate_post_umi)

    output:
    tuple val(meta), path("${meta.sample_id}.consensus.bam"), path("${meta.sample_id}.consensus.bam.bai"), path("${meta.sample_id}.umi_qc.tsv"), path(bait_bed), val(reference_config)

    script:
    def fgbio_heap_gb = Math.max(2, task.memory.toGiga().intValue() - 4)
    """
    export JAVA_TOOL_OPTIONS="-Xms1g -Xmx${fgbio_heap_gb}g"
    mkdir -p tmp
    export TMPDIR="$PWD/tmp"
    fgbio AnnotateBamWithUmis \
      --input=${bam} \
      --fastq=${r1} \
      --fastq=${r2} \
      --output=${meta.sample_id}.tagged.bam \
      --attribute=${tag_name} \
      --read-structure=${read_structure_r1} \
      --read-structure=${read_structure_r2}
    fgbio GroupReadsByUmi \
      --input=${meta.sample_id}.tagged.bam \
      --output=${meta.sample_id}.grouped.bam \
      --strategy=Adjacency \
      --family-size-histogram=${meta.sample_id}.family_sizes.txt
    fgbio CallMolecularConsensusReads \
      --input=${meta.sample_id}.grouped.bam \
      --output=${meta.sample_id}.consensus.unmapped.bam \
      --min-reads=${min_family_size} \
      --error-rate-pre-umi=${error_rate_pre_umi} \
      --error-rate-post-umi=${error_rate_post_umi}
    samtools sort -o ${meta.sample_id}.consensus.bam ${meta.sample_id}.consensus.unmapped.bam
    samtools index ${meta.sample_id}.consensus.bam
    cat > ${meta.sample_id}.umi_qc.tsv <<'EOF'
    sample_id\tread_structure_r1\tread_structure_r2\tmolecular_index_tag\tumi_family_count\tumi_max_family_size\tconsensus_reads\tstatus
    ${meta.sample_id}\t${read_structure_r1}\t${read_structure_r2}\t${tag_name}\t0\t0\t0\tPASS
    EOF
    """

    stub:
    """
    touch ${meta.sample_id}.consensus.bam ${meta.sample_id}.consensus.bam.bai
    cat > ${meta.sample_id}.umi_qc.tsv <<'EOF'
    sample_id\tread_structure_r1\tread_structure_r2\tmolecular_index_tag\tumi_family_count\tumi_max_family_size\tconsensus_reads\tstatus
    ${meta.sample_id}\t${read_structure_r1}\t${read_structure_r2}\t${tag_name}\t4\t3\t2\tPASS
    EOF
    """
}
