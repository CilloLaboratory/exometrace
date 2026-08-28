process UMI_GROUP_CONSENSUS {
    label 'cpu_medium'
    publishDir "${params.results_dir}/ctdna/consensus", mode: 'copy'
    container { params.umi_consensus_container }

    input:
    tuple val(meta), path(r1), path(r2), path(bait_bed), val(reference_config), path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict), path(ref_fasta_0123), path(ref_fasta_amb), path(ref_fasta_ann), path(ref_fasta_bwt), path(ref_fasta_pac), val(tag_name), val(read_structure_r1), val(read_structure_r2), val(min_family_size), val(error_rate_pre_umi), val(error_rate_post_umi)

    output:
    tuple val(meta), path("${meta.sample_id}.consensus.bam"), path("${meta.sample_id}.consensus.bam.bai"), path("${meta.sample_id}.umi_qc.tsv"), path(bait_bed), val(reference_config)

    script:
    def fgbio_heap_gb = Math.max(4, task.memory.toGiga().intValue() - 4)
    def sort_threads = task.cpus > 8 ? 8 : task.cpus
    def consensus_threads = task.cpus > 4 ? 4 : task.cpus
    """
    export JAVA_TOOL_OPTIONS="-Xms1g -Xmx${fgbio_heap_gb}g"
    mkdir -p tmp
    export TMPDIR="\$PWD/tmp"

    fgbio FastqToBam \
      --input ${r1} ${r2} \
      --read-structures ${read_structure_r1} ${read_structure_r2} \
      --sample ${meta.sample_id} \
      --library ${meta.sample_id}_lib1 \
      --platform-unit ${meta.patient_id} \
      --umi-tag ${tag_name} \
      --read-group-id ${meta.sample_id} \
      --output ${meta.sample_id}.unmapped.bam

    samtools fastq ${meta.sample_id}.unmapped.bam \
      | bwa mem -t ${task.cpus} -p -K 150000000 -Y ${ref_fasta} - \
      | fgbio ZipperBams \
          --unmapped ${meta.sample_id}.unmapped.bam \
          --ref ${ref_fasta} \
          --output ${meta.sample_id}.mapped.bam

    samtools sort --template-coordinate --threads ${sort_threads} \
      -o ${meta.sample_id}.mapped.template.bam \
      ${meta.sample_id}.mapped.bam

    fgbio GroupReadsByUmi \
      --input=${meta.sample_id}.mapped.template.bam \
      --output=${meta.sample_id}.grouped.bam \
      --strategy=Adjacency \
      --family-size-histogram=${meta.sample_id}.family_sizes.txt

    fgbio CallMolecularConsensusReads \
      --input=${meta.sample_id}.grouped.bam \
      --output=${meta.sample_id}.consensus.unmapped.bam \
      --min-reads=${min_family_size} \
      --threads=${consensus_threads} \
      --error-rate-pre-umi=${error_rate_pre_umi} \
      --error-rate-post-umi=${error_rate_post_umi}

    samtools fastq ${meta.sample_id}.consensus.unmapped.bam \
      | bwa mem -t ${task.cpus} -p -K 150000000 -Y ${ref_fasta} - \
      | fgbio ZipperBams \
          --unmapped ${meta.sample_id}.consensus.unmapped.bam \
          --ref ${ref_fasta} \
          --tags-to-reverse Consensus \
          --tags-to-revcomp Consensus \
          --output ${meta.sample_id}.consensus.mapped.bam

    samtools sort -@ ${sort_threads} -o ${meta.sample_id}.consensus.bam ${meta.sample_id}.consensus.mapped.bam
    samtools index ${meta.sample_id}.consensus.bam

    umi_family_count=\$(awk 'NF >= 2 && \$1 !~ /^#/ {sum+=\$2} END {print sum+0}' ${meta.sample_id}.family_sizes.txt)
    umi_max_family_size=\$(awk 'NF >= 2 && \$1 !~ /^#/ {if (\$1>max) max=\$1} END {print max+0}' ${meta.sample_id}.family_sizes.txt)
    consensus_reads=\$(samtools view -c ${meta.sample_id}.consensus.bam)
    cat > ${meta.sample_id}.umi_qc.tsv <<EOF
    sample_id\tread_structure_r1\tread_structure_r2\tmolecular_index_tag\tumi_family_count\tumi_max_family_size\tconsensus_reads\tstatus
    ${meta.sample_id}\t${read_structure_r1}\t${read_structure_r2}\t${tag_name}\t\${umi_family_count}\t\${umi_max_family_size}\t\${consensus_reads}\tPASS
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
