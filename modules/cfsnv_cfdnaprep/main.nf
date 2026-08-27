process CFSNV_CFDNA_PREP {
    label 'cpu_large'
    publishDir "${params.results_dir}/ctdna/cfsnv", mode: 'copy'
    container { params.cfsnv_container }

    input:
    tuple val(meta), path(plasma_r1), path(plasma_r2), path(bait_bed), val(reference_config), val(ref_fasta), val(snp_database)

    output:
    tuple val(meta), path("${meta.patient_id}.plasma.extendedFrags.recal.bam"), path("${meta.patient_id}.plasma.extendedFrags.recal.bam.bai"), path("${meta.patient_id}.plasma.notCombined.recal.bam"), path("${meta.patient_id}.plasma.notCombined.recal.bam.bai")

    script:
    """
    Rscript ${projectDir}/scripts/cfsnv_wrapper.R cfDNAprep \
      --fastq1 ${plasma_r1} \
      --fastq2 ${plasma_r2} \
      --reference ${ref_fasta} \
      --snp-database ${snp_database} \
      --sample-id ${meta.patient_id} \
      --output-dir .
    mv ${meta.patient_id}.extendedFrags.recal.bam ${meta.patient_id}.plasma.extendedFrags.recal.bam
    mv ${meta.patient_id}.notCombined.recal.bam ${meta.patient_id}.plasma.notCombined.recal.bam
    samtools index ${meta.patient_id}.plasma.extendedFrags.recal.bam
    samtools index ${meta.patient_id}.plasma.notCombined.recal.bam
    """

    stub:
    """
    touch ${meta.patient_id}.plasma.extendedFrags.recal.bam ${meta.patient_id}.plasma.extendedFrags.recal.bam.bai
    touch ${meta.patient_id}.plasma.notCombined.recal.bam ${meta.patient_id}.plasma.notCombined.recal.bam.bai
    """
}
