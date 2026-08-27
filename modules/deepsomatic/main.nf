process DEEPSOMATIC {
    label 'gpu'
    publishDir "${params.results_dir}/somatic/deepsomatic", mode: 'copy'
    container { params.deepsomatic_container }

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(bait_bed), val(reference_config), val(ref_fasta)

    output:
    tuple val(meta), path("${meta.patient_id}.somatic.vcf.gz"), path("${meta.patient_id}.somatic.vcf.gz.tbi")

    script:
    """
    mkdir -p .mplconfig .cache
    export MPLCONFIGDIR=\$PWD/.mplconfig
    export XDG_CACHE_HOME=\$PWD/.cache

    run_deepsomatic \
      --model_type=${meta.ffpe == 'true' ? 'FFPE_WES' : 'WES'} \
      --ref=${ref_fasta} \
      --reads_normal=${normal_bam} \
      --reads_tumor=${tumor_bam} \
      --regions=${bait_bed} \
      --output_vcf=${meta.patient_id}.somatic.vcf.gz \
      --output_gvcf=${meta.patient_id}.somatic.g.vcf.gz \
      --sample_name_tumor=${meta.tumor_sample} \
      --sample_name_normal=${meta.normal_sample} \
      --num_shards=${task.cpus} \
      --logging_dir=logs \
      --intermediate_results_dir intermediate_results_dir
    bcftools view ${meta.patient_id}.somatic.vcf.gz > /dev/null
    """

    stub:
    """
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${meta.tumor_sample}\t${meta.normal_sample}\nchr1\t5\t.\tA\tT\t60\tPASS\t.\tGT:DP\t0/1:30\t0/0:28\n' > ${meta.patient_id}.somatic.vcf
    gzip -c ${meta.patient_id}.somatic.vcf > ${meta.patient_id}.somatic.vcf.gz
    touch ${meta.patient_id}.somatic.vcf.gz.tbi
    """
}
