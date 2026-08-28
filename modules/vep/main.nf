process VEP {
    label 'cpu_medium'
    publishDir "${params.results_dir}/annotation", mode: 'copy'
    container { params.vep_container }

    input:
    tuple val(meta), path(input_vcf), val(reference_config), val(callset_name), path(ref_fasta), path(ref_fasta_fai), path(vep_cache)

    output:
    tuple val(meta), path("${meta.patient_id}.${callset_name}.vep.vcf.gz")

    script:
    """
    vep \
      --input_file ${input_vcf} \
      --format vcf \
      --output_file ${meta.patient_id}.${callset_name}.vep.vcf.gz \
      --vcf \
      --compress_output bgzip \
      --offline \
      --cache \
      --dir_cache ${vep_cache} \
      --assembly GRCh38 \
      --fasta ${ref_fasta} \
      --symbol \
      --canonical \
      --protein \
      --biotype \
      --hgvs \
      --variant_class \
      --af_gnomad \
      --force_overwrite
    """

    stub:
    """
    printf '##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|HGVSc|HGVSp">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${meta.tumor_sample}\t${meta.normal_sample}\nchr1\t5\t.\tA\tT\t60\tPASS\tCSQ=T|missense_variant|MODERATE|TP53|ENSG00000141510|Transcript|ENST00000269305|protein_coding|ENST00000269305.8:c.215C>T|ENSP00000269305.4:p.Pro72Leu\tGT:DP:AD\t0/1:30:20,10\t0/0:28:28,0\n' | gzip -c > ${meta.patient_id}.${callset_name}.vep.vcf.gz
    """
}
