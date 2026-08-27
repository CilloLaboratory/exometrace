process FILTER_MUTECT {
    label 'cpu_small'
    publishDir "${params.results_dir}/somatic/mutect2", mode: 'copy'
    container { params.gatk_container }

    input:
    tuple val(meta), path(unfiltered_vcf), path(unfiltered_vcf_tbi), path(stats_file), path(f1r2_tar), path(tumor_pileups), path(normal_pileups), val(reference_config), val(ref_fasta)

    output:
    tuple val(meta), path("${meta.patient_id}.mutect2.filtered.vcf.gz"), path("${meta.patient_id}.mutect2.filtered.vcf.gz.tbi"), path("${meta.patient_id}.contamination.table"), path("${meta.patient_id}.segments.table")

    script:
    """
    gatk CalculateContamination \
      -I ${tumor_pileups} \
      -matched ${normal_pileups} \
      -O ${meta.patient_id}.contamination.table \
      --tumor-segmentation ${meta.patient_id}.segments.table

    gatk LearnReadOrientationModel \
      -I ${f1r2_tar} \
      -O ${meta.patient_id}.artifact-priors.tar.gz

    gatk FilterMutectCalls \
      -V ${unfiltered_vcf} \
      -R ${ref_fasta} \
      --stats ${stats_file} \
      --contamination-table ${meta.patient_id}.contamination.table \
      --tumor-segmentation ${meta.patient_id}.segments.table \
      --orientation-bias-artifact-priors ${meta.patient_id}.artifact-priors.tar.gz \
      -O ${meta.patient_id}.mutect2.filtered.vcf.gz
    bcftools view ${meta.patient_id}.mutect2.filtered.vcf.gz > /dev/null
    """

    stub:
    """
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${meta.tumor_sample}\t${meta.normal_sample}\nchr1\t5\t.\tA\tT\t50\tPASS\t.\tGT:DP\t0/1:31\t0/0:29\nchr1\t8\t.\tC\tG\t45\tPASS\t.\tGT:DP\t0/1:27\t0/0:26\n' > ${meta.patient_id}.mutect2.filtered.vcf
    gzip -c ${meta.patient_id}.mutect2.filtered.vcf > ${meta.patient_id}.mutect2.filtered.vcf.gz
    touch ${meta.patient_id}.mutect2.filtered.vcf.gz.tbi
    printf 'sample\tcontamination\n${meta.tumor_sample}\t0.01\n' > ${meta.patient_id}.contamination.table
    printf 'contig\tstart\tend\tminor_allele_fraction\nchr1\t1\t12\t0.45\n' > ${meta.patient_id}.segments.table
    """
}
