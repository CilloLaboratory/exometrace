process MUTECT2 {
    label 'cpu_medium'
    publishDir "${params.results_dir}/somatic/mutect2", mode: 'copy'
    container { params.gatk_container }

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(bait_bed), val(reference_config), path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict), path(germline_resource), path(germline_resource_tbi), path(common_snps), path(common_snps_tbi), path(pon), path(pon_tbi), path(target_intervals)

    output:
    tuple val(meta), path("${meta.patient_id}.mutect2.unfiltered.vcf.gz"), path("${meta.patient_id}.mutect2.unfiltered.vcf.gz.tbi"), path("${meta.patient_id}.mutect2.unfiltered.vcf.gz.stats"), path("${meta.patient_id}.f1r2.tar.gz"), path("${meta.patient_id}.tumor.pileups.table"), path("${meta.patient_id}.normal.pileups.table"), val(reference_config), path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict)

    script:
    def gatk_heap_gb = Math.max(2, task.memory.toGiga().intValue() - 4)
    """
    export JAVA_TOOL_OPTIONS="-Xms1g -Xmx${gatk_heap_gb}g"
    gatk --java-options "-Xms1g -Xmx${gatk_heap_gb}g" Mutect2 \
      -R ${ref_fasta} \
      -I ${tumor_bam} \
      -I ${normal_bam} \
      -tumor ${meta.tumor_sample} \
      -normal ${meta.normal_sample} \
      --germline-resource ${germline_resource} \
      --panel-of-normals ${pon} \
      -L ${target_intervals} \
      --f1r2-tar-gz ${meta.patient_id}.f1r2.tar.gz \
      -O ${meta.patient_id}.mutect2.unfiltered.vcf.gz

    gatk --java-options "-Xms1g -Xmx${gatk_heap_gb}g" IndexFeatureFile \
      -I ${meta.patient_id}.mutect2.unfiltered.vcf.gz

    gatk --java-options "-Xms1g -Xmx${gatk_heap_gb}g" GetPileupSummaries \
      -I ${tumor_bam} \
      -V ${common_snps} \
      -L ${target_intervals} \
      -O ${meta.patient_id}.tumor.pileups.table

    gatk --java-options "-Xms1g -Xmx${gatk_heap_gb}g" GetPileupSummaries \
      -I ${normal_bam} \
      -V ${common_snps} \
      -L ${target_intervals} \
      -O ${meta.patient_id}.normal.pileups.table
    """

    stub:
    """
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${meta.tumor_sample}\t${meta.normal_sample}\nchr1\t5\t.\tA\tT\t50\tPASS\t.\tGT:DP\t0/1:31\t0/0:29\nchr1\t8\t.\tC\tG\t45\tPASS\t.\tGT:DP\t0/1:27\t0/0:26\n' > ${meta.patient_id}.mutect2.unfiltered.vcf
    gzip -c ${meta.patient_id}.mutect2.unfiltered.vcf > ${meta.patient_id}.mutect2.unfiltered.vcf.gz
    touch ${meta.patient_id}.mutect2.unfiltered.vcf.gz.tbi
    echo "stub stats" > ${meta.patient_id}.mutect2.unfiltered.vcf.gz.stats
    tar -czf ${meta.patient_id}.f1r2.tar.gz ${meta.patient_id}.mutect2.unfiltered.vcf
    printf 'contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\nchr1\t5\t20\t10\t0\t0.05\n' > ${meta.patient_id}.tumor.pileups.table
    printf 'contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\nchr1\t5\t25\t0\t0\t0.05\n' > ${meta.patient_id}.normal.pileups.table
    """
}
