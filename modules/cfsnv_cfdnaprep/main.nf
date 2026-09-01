process CFSNV_CFDNA_PREP {
    label 'cpu_large'
    publishDir "${params.results_dir}/ctdna/cfsnv", mode: 'copy'
    container { params.cfsnv_container }

    input:
    tuple val(meta), path(plasma_r1), path(plasma_r2), path(bait_bed), val(reference_config), path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict), path(ref_fasta_amb), path(ref_fasta_ann), path(ref_fasta_bwt), path(ref_fasta_pac), path(ref_fasta_sa), path(snp_database), path(snp_database_tbi)

    output:
    tuple val(meta), path("${meta.patient_id}.plasma.extendedFrags.recal.bam"), path("${meta.patient_id}.plasma.extendedFrags.recal.bam.bai"), path("${meta.patient_id}.plasma.notCombined.recal.bam"), path("${meta.patient_id}.plasma.notCombined.recal.bam.bai")

    script:
    """
    mkdir -p tmp r_libs
    export TMPDIR="\$PWD/tmp"
    export CFSNV_R_LIB_ROOT="\$PWD/r_libs"
    export CFSNV_JAVA="/opt/conda/bin/java"
    export CFSNV_PICARD_JAR="/usr/local/share/cfsnv-tools/picard.jar"
    export CFSNV_GATK_JAR="/usr/local/share/cfsnv-tools/GenomeAnalysisTK.jar"
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
