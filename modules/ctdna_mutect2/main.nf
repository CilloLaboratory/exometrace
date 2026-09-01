process CTDNA_MUTECT2 {
    label 'cpu_medium'
    publishDir "${params.results_dir}/ctdna/mutect2", mode: 'copy'
    container { params.gatk_container }

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(bait_bed), val(reference_config), path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict), path(germline_resource), path(germline_resource_tbi), path(common_snps), path(common_snps_tbi), path(pon), path(pon_tbi), path(target_intervals)

    output:
    tuple val(meta), path("${meta.patient_id}.ctdna.mutect2.unfiltered.vcf.gz"), path("${meta.patient_id}.ctdna.mutect2.unfiltered.vcf.gz.tbi"), path("${meta.patient_id}.ctdna.mutect2.unfiltered.vcf.gz.stats"), path("${meta.patient_id}.ctdna.f1r2.tar.gz"), path("${meta.patient_id}.ctdna.tumor.pileups.table"), path("${meta.patient_id}.ctdna.normal.pileups.table"), val(reference_config), path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict)

    script:
    def gatk_heap_gb = Math.max(2, task.memory.toGiga().intValue() - 4)
    def gatk_cmd = params.gatk_local_jar ? "java -Xms1g -Xmx${gatk_heap_gb}g -jar ${params.gatk_local_jar}" : "gatk --java-options \"-Xms1g -Xmx${gatk_heap_gb}g\""
    """
    export JAVA_TOOL_OPTIONS="-Xms1g -Xmx${gatk_heap_gb}g"
    export JAVA_OPTS="-Xms1g -Xmx${gatk_heap_gb}g"
    export GATK_JAVA_OPTIONS="-Xms1g -Xmx${gatk_heap_gb}g"
    ${gatk_cmd} Mutect2 \
      -R ${ref_fasta} \
      -I ${tumor_bam} \
      -I ${normal_bam} \
      -tumor ${meta.tumor_sample} \
      -normal ${meta.normal_sample} \
      --germline-resource ${germline_resource} \
      --panel-of-normals ${pon} \
      -L ${target_intervals} \
      --f1r2-tar-gz ${meta.patient_id}.ctdna.f1r2.tar.gz \
      -O ${meta.patient_id}.ctdna.mutect2.unfiltered.vcf.gz
    ${gatk_cmd} IndexFeatureFile -I ${meta.patient_id}.ctdna.mutect2.unfiltered.vcf.gz
    ${gatk_cmd} GetPileupSummaries -I ${tumor_bam} -V ${common_snps} -L ${target_intervals} -O ${meta.patient_id}.ctdna.tumor.pileups.table
    ${gatk_cmd} GetPileupSummaries -I ${normal_bam} -V ${common_snps} -L ${target_intervals} -O ${meta.patient_id}.ctdna.normal.pileups.table
    """

    stub:
    """
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${meta.tumor_sample}\t${meta.normal_sample}\nchr1\t5\t.\tA\tT\t50\tPASS\t.\tGT:DP:AD\t0/1:31:20,11\t0/0:29:29,0\nchr1\t11\t.\tG\tC\t45\tPASS\t.\tGT:DP:AD\t0/1:18:15,3\t0/0:20:20,0\n' > ${meta.patient_id}.ctdna.mutect2.unfiltered.vcf
    gzip -c ${meta.patient_id}.ctdna.mutect2.unfiltered.vcf > ${meta.patient_id}.ctdna.mutect2.unfiltered.vcf.gz
    touch ${meta.patient_id}.ctdna.mutect2.unfiltered.vcf.gz.tbi
    echo "stub stats" > ${meta.patient_id}.ctdna.mutect2.unfiltered.vcf.gz.stats
    tar -czf ${meta.patient_id}.ctdna.f1r2.tar.gz ${meta.patient_id}.ctdna.mutect2.unfiltered.vcf
    printf 'contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\nchr1\t5\t20\t11\t0\t0.05\n' > ${meta.patient_id}.ctdna.tumor.pileups.table
    printf 'contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\nchr1\t5\t25\t0\t0\t0.05\n' > ${meta.patient_id}.ctdna.normal.pileups.table
    """
}

process CTDNA_FILTER_MUTECT2 {
    label 'cpu_small'
    publishDir "${params.results_dir}/ctdna/mutect2", mode: 'copy'
    container { params.gatk_container }

    input:
    tuple val(meta), path(unfiltered_vcf), path(unfiltered_vcf_tbi), path(stats_file), path(f1r2_tar), path(tumor_pileups), path(normal_pileups), val(reference_config), path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict)

    output:
    tuple val(meta), path("${meta.patient_id}.ctdna.mutect2.filtered.vcf.gz"), path("${meta.patient_id}.ctdna.mutect2.filtered.vcf.gz.tbi"), path("${meta.patient_id}.ctdna.contamination.table"), path("${meta.patient_id}.ctdna.segments.table")

    script:
    def gatk_heap_gb = Math.max(2, task.memory.toGiga().intValue() - 2)
    def gatk_cmd = params.gatk_local_jar ? "java -Xms1g -Xmx${gatk_heap_gb}g -jar ${params.gatk_local_jar}" : "gatk --java-options \"-Xms1g -Xmx${gatk_heap_gb}g\""
    """
    export JAVA_TOOL_OPTIONS="-Xms1g -Xmx${gatk_heap_gb}g"
    export JAVA_OPTS="-Xms1g -Xmx${gatk_heap_gb}g"
    export GATK_JAVA_OPTIONS="-Xms1g -Xmx${gatk_heap_gb}g"
    ${gatk_cmd} CalculateContamination \
      -I ${tumor_pileups} \
      -matched ${normal_pileups} \
      -O ${meta.patient_id}.ctdna.contamination.table \
      --tumor-segmentation ${meta.patient_id}.ctdna.segments.table
    ${gatk_cmd} LearnReadOrientationModel -I ${f1r2_tar} -O ${meta.patient_id}.ctdna.artifact-priors.tar.gz
    ${gatk_cmd} FilterMutectCalls \
      -V ${unfiltered_vcf} \
      -R ${ref_fasta} \
      --stats ${stats_file} \
      --contamination-table ${meta.patient_id}.ctdna.contamination.table \
      --tumor-segmentation ${meta.patient_id}.ctdna.segments.table \
      --orientation-bias-artifact-priors ${meta.patient_id}.ctdna.artifact-priors.tar.gz \
      -O ${meta.patient_id}.ctdna.mutect2.filtered.vcf.gz
    bcftools view ${meta.patient_id}.ctdna.mutect2.filtered.vcf.gz > /dev/null
    """

    stub:
    """
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${meta.tumor_sample}\t${meta.normal_sample}\nchr1\t5\t.\tA\tT\t50\tPASS\t.\tGT:DP:AD\t0/1:31:20,11\t0/0:29:29,0\nchr1\t11\t.\tG\tC\t45\tPASS\t.\tGT:DP:AD\t0/1:18:15,3\t0/0:20:20,0\n' > ${meta.patient_id}.ctdna.mutect2.filtered.vcf
    gzip -c ${meta.patient_id}.ctdna.mutect2.filtered.vcf > ${meta.patient_id}.ctdna.mutect2.filtered.vcf.gz
    touch ${meta.patient_id}.ctdna.mutect2.filtered.vcf.gz.tbi
    printf 'sample\tcontamination\n${meta.tumor_sample}\t0.01\n' > ${meta.patient_id}.ctdna.contamination.table
    printf 'contig\tstart\tend\tminor_allele_fraction\nchr1\t1\t12\t0.45\n' > ${meta.patient_id}.ctdna.segments.table
    """
}

process CTDNA_CFSNV_CALL {
    label 'cpu_medium'
    publishDir "${params.results_dir}/ctdna/cfsnv", mode: 'copy'
    container { params.cfsnv_container }

    input:
    tuple val(meta), path(plasma_std_bam), path(plasma_std_bai), path(wbc_std_bam), path(wbc_std_bai), path(plasma_extended_bam), path(plasma_extended_bai), path(plasma_notcombined_bam), path(plasma_notcombined_bai), path(bait_bed), path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict), path(snp_database), path(snp_database_tbi), path(blocked_positions_vcf), path(blocked_positions_tbi), val(min_hold_support), val(min_pass_support)

    output:
    tuple val(meta), path("${meta.patient_id}.cfsnv.vcf.gz"), path("${meta.patient_id}.cfsnv.vcf.gz.tbi")

    script:
    """
    mkdir -p tmp r_libs
    export TMPDIR="\$PWD/tmp"
    export CFSNV_R_LIB_ROOT="\$PWD/r_libs"
    export CFSNV_JAVA="/opt/conda/bin/java"
    export CFSNV_PICARD_JAR="/usr/local/share/cfsnv-tools/picard.jar"
    Rscript ${projectDir}/scripts/cfsnv_wrapper.R DetectMuts \
      --tumor-bam ${plasma_std_bam} \
      --normal-bam ${wbc_std_bam} \
      --extended-bam ${plasma_extended_bam} \
      --not-combined-bam ${plasma_notcombined_bam} \
      --targets ${bait_bed} \
      --reference ${ref_fasta} \
      --snp-database ${snp_database} \
      --blocked-positions ${blocked_positions_vcf} \
      --sample-id ${meta.patient_id} \
      --tumor-sample ${meta.tumor_sample} \
      --normal-sample ${meta.normal_sample} \
      --min-hold-support ${min_hold_support} \
      --min-pass-support ${min_pass_support} \
      --output ${meta.patient_id}.cfsnv.vcf.gz
    """

    stub:
    """
    printf '##fileformat=VCFv4.2\n##INFO=<ID=UMI_FAMILY_COUNT,Number=1,Type=Integer,Description="Consensus family support">\n##INFO=<ID=UMI_MAX_FAMILY_SIZE,Number=1,Type=Integer,Description="Largest observed family size">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${meta.tumor_sample}\t${meta.normal_sample}\nchr1\t5\t.\tA\tT\t65\tPASS\tUMI_FAMILY_COUNT=4;UMI_MAX_FAMILY_SIZE=3\tGT:DP:AD\t0/1:28:20,8\t0/0:32:32,0\nchr1\t8\t.\tC\tG\t55\tPASS\tUMI_FAMILY_COUNT=2;UMI_MAX_FAMILY_SIZE=2\tGT:DP:AD\t0/1:20:18,2\t0/0:30:30,0\n' > ${meta.patient_id}.cfsnv.vcf
    gzip -c ${meta.patient_id}.cfsnv.vcf > ${meta.patient_id}.cfsnv.vcf.gz
    touch ${meta.patient_id}.cfsnv.vcf.gz.tbi
    """
}
