process FASTQC {
    label 'cpu_small'
    publishDir "${params.results_dir}/qc/fastqc", mode: 'copy'
    container { params.qc_container }

    input:
    tuple val(meta), path(r1), path(r2), path(bait_bed)

    output:
    tuple val(meta), path("${meta.sample_id}_R1_fastqc.html"), path("${meta.sample_id}_R1_fastqc.zip"), path("${meta.sample_id}_R2_fastqc.html"), path("${meta.sample_id}_R2_fastqc.zip"), path(r1), path(r2), path(bait_bed)

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${r1} ${r2}
    r1_base=\$(basename ${r1} .fastq.gz)
    r2_base=\$(basename ${r2} .fastq.gz)
    if [[ "\${r1_base}_fastqc.html" != "${meta.sample_id}_R1_fastqc.html" ]]; then
      mv "\${r1_base}_fastqc.html" ${meta.sample_id}_R1_fastqc.html
    fi
    if [[ "\${r1_base}_fastqc.zip" != "${meta.sample_id}_R1_fastqc.zip" ]]; then
      mv "\${r1_base}_fastqc.zip" ${meta.sample_id}_R1_fastqc.zip
    fi
    if [[ "\${r2_base}_fastqc.html" != "${meta.sample_id}_R2_fastqc.html" ]]; then
      mv "\${r2_base}_fastqc.html" ${meta.sample_id}_R2_fastqc.html
    fi
    if [[ "\${r2_base}_fastqc.zip" != "${meta.sample_id}_R2_fastqc.zip" ]]; then
      mv "\${r2_base}_fastqc.zip" ${meta.sample_id}_R2_fastqc.zip
    fi
    """

    stub:
    """
    touch ${meta.sample_id}_R1_fastqc.zip ${meta.sample_id}_R2_fastqc.zip
    cat > ${meta.sample_id}_R1_fastqc.html <<'EOF'
    <html><body>stub fastqc r1</body></html>
    EOF
    cat > ${meta.sample_id}_R2_fastqc.html <<'EOF'
    <html><body>stub fastqc r2</body></html>
    EOF
    """
}
