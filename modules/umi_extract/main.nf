process UMI_TEMPLATE_TRIM {
    label 'cpu_small'
    cache 'deep'
    publishDir "${params.results_dir}/ctdna/qc", mode: 'copy', pattern: "*.template_trim.tsv"

    input:
    tuple val(meta), path(r1), path(r2), path(bait_bed), val(read_structure_r1), val(read_structure_r2)

    output:
    tuple val(meta), path("${meta.sample_id}.trimmed_R1.fastq.gz"), path("${meta.sample_id}.trimmed_R2.fastq.gz"), path("${meta.sample_id}.template_trim.tsv"), path(bait_bed)

    script:
    """
    python3 ${projectDir}/scripts/trim_ctdna_template.py \
      --r1 ${r1} \
      --r2 ${r2} \
      --read-structure-r1 ${read_structure_r1} \
      --read-structure-r2 ${read_structure_r2} \
      --output-r1 ${meta.sample_id}.trimmed_R1.fastq.gz \
      --output-r2 ${meta.sample_id}.trimmed_R2.fastq.gz \
      --qc-output ${meta.sample_id}.template_trim.tsv \
      --sample-id ${meta.sample_id}
    """

    stub:
    """
    cp ${r1} ${meta.sample_id}.trimmed_R1.fastq.gz
    cp ${r2} ${meta.sample_id}.trimmed_R2.fastq.gz
    cat > ${meta.sample_id}.template_trim.tsv <<'EOF'
    sample_id\tread_pairs\tread_structure_r1\tread_structure_r2\ttemplate_trim_bases_r1\ttemplate_trim_bases_r2\tstatus
    ${meta.sample_id}\t2\t${read_structure_r1}\t${read_structure_r2}\t6\t6\tPASS
    EOF
    """
}
