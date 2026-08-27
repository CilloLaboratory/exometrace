process COMPARE_CTDNA_CALLSETS {
    label 'cpu_small'
    publishDir "${params.results_dir}/ctdna/comparison", mode: 'copy'
    container { params.qc_container }

    input:
    tuple val(meta), path(cfsnv_vcf), path(mutect2_vcf)

    output:
    tuple val(meta), path("${meta.patient_id}.ctdna_call_concordance.tsv")

    script:
    """
    python3 ${projectDir}/scripts/compare_ctdna_callsets.py \
      --cfsnv ${cfsnv_vcf} \
      --mutect2 ${mutect2_vcf} \
      --output ${meta.patient_id}.ctdna_call_concordance.tsv
    """

    stub:
    """
    python3 ${projectDir}/scripts/compare_ctdna_callsets.py \
      --cfsnv ${cfsnv_vcf} \
      --mutect2 ${mutect2_vcf} \
      --output ${meta.patient_id}.ctdna_call_concordance.tsv
    """
}
