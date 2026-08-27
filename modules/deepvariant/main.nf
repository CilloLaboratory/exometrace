process DEEPVARIANT {
    label 'gpu'
    publishDir "${params.results_dir}/germline/deepvariant", mode: 'copy'
    container { params.deepvariant_container }

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(bait_bed), val(reference_config), val(ref_fasta)

    output:
    tuple val(meta), path("${meta.sample_id}.deepvariant.vcf.gz"), path("${meta.sample_id}.deepvariant.vcf.gz.tbi"), path("${meta.sample_id}.deepvariant.g.vcf.gz")

    script:
    """
    mkdir -p .mplconfig .cache
    export MPLCONFIGDIR=\$PWD/.mplconfig
    export XDG_CACHE_HOME=\$PWD/.cache

    run_deepvariant \
      --model_type=WES \
      --ref=${ref_fasta} \
      --reads=${normal_bam} \
      --regions=${bait_bed} \
      --output_vcf=${meta.sample_id}.deepvariant.vcf.gz \
      --output_gvcf=${meta.sample_id}.deepvariant.g.vcf.gz \
      --num_shards=${task.cpus} \
      --logging_dir=logs
    bcftools view ${meta.sample_id}.deepvariant.vcf.gz > /dev/null
    """

    stub:
    """
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${meta.sample_id}\nchr1\t3\t.\tG\tA\t55\tPASS\t.\tGT:DP\t0/1:25\n' > ${meta.sample_id}.deepvariant.vcf
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${meta.sample_id}\nchr1\t3\t.\tG\t<NON_REF>\t.\tPASS\tEND=12\tGT:DP\t0/0:25\n' > ${meta.sample_id}.deepvariant.g.vcf
    gzip -c ${meta.sample_id}.deepvariant.vcf > ${meta.sample_id}.deepvariant.vcf.gz
    gzip -c ${meta.sample_id}.deepvariant.g.vcf > ${meta.sample_id}.deepvariant.g.vcf.gz
    touch ${meta.sample_id}.deepvariant.vcf.gz.tbi
    """
}
