process CFSNV_STD_PREP {
    label 'cpu_large'
    publishDir "${params.results_dir}/ctdna/cfsnv", mode: 'copy'
    container { params.cfsnv_container }

    input:
    tuple val(meta), path(plasma_r1), path(plasma_r2), path(wbc_r1), path(wbc_r2), path(bait_bed), val(reference_config), val(ref_fasta), val(snp_database)

    output:
    tuple val(meta), path("${meta.patient_id}.plasma.recal.bam"), path("${meta.patient_id}.plasma.recal.bam.bai"), path("${meta.patient_id}.wbc.recal.bam"), path("${meta.patient_id}.wbc.recal.bam.bai"), path("${meta.patient_id}.ctdna_sample_qc.tsv")

    script:
    """
    Rscript ${projectDir}/scripts/cfsnv_wrapper.R STDprep \
      --fastq1 ${plasma_r1} \
      --fastq2 ${plasma_r2} \
      --reference ${ref_fasta} \
      --snp-database ${snp_database} \
      --sample-id ${meta.patient_id}.tumor \
      --output-dir .
    Rscript ${projectDir}/scripts/cfsnv_wrapper.R STDprep \
      --fastq1 ${wbc_r1} \
      --fastq2 ${wbc_r2} \
      --reference ${ref_fasta} \
      --snp-database ${snp_database} \
      --sample-id ${meta.patient_id}.normal \
      --output-dir .
    mv ${meta.patient_id}.tumor.recal.bam ${meta.patient_id}.plasma.recal.bam
    mv ${meta.patient_id}.normal.recal.bam ${meta.patient_id}.wbc.recal.bam
    samtools index ${meta.patient_id}.plasma.recal.bam
    samtools index ${meta.patient_id}.wbc.recal.bam
    cat > ${meta.patient_id}.ctdna_sample_qc.tsv <<'EOF'
    patient_id\tplasma_sample\twbc_sample\tplasma_bam\twbc_bam\tstatus
    ${meta.patient_id}\t${meta.tumor_sample}\t${meta.normal_sample}\t${meta.patient_id}.plasma.recal.bam\t${meta.patient_id}.wbc.recal.bam\tPASS
    EOF
    """

    stub:
    """
    touch ${meta.patient_id}.plasma.recal.bam ${meta.patient_id}.plasma.recal.bam.bai
    touch ${meta.patient_id}.wbc.recal.bam ${meta.patient_id}.wbc.recal.bam.bai
    cat > ${meta.patient_id}.ctdna_sample_qc.tsv <<'EOF'
    patient_id\tplasma_sample\twbc_sample\tplasma_bam\twbc_bam\tstatus
    ${meta.patient_id}\t${meta.tumor_sample}\t${meta.normal_sample}\t${meta.patient_id}.plasma.recal.bam\t${meta.patient_id}.wbc.recal.bam\tPASS
    EOF
    """
}
