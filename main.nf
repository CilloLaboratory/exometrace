nextflow.enable.dsl = 2

def readConfigValue(String configPath, String key) {
    def command = ['python3', "${projectDir}/scripts/read_yaml_value.py", '--config', configPath, '--key', key]
    def process = command.execute()
    def stdout = new StringBuffer()
    def stderr = new StringBuffer()
    process.consumeProcessOutput(stdout, stderr)
    process.waitFor()
    if (process.exitValue() != 0) {
        error "Failed to read '${key}' from ${configPath}: ${stderr.toString().trim()}"
    }
    stdout.toString().trim()
}

include { VALIDATE_SAMPLESHEET } from './modules/validate_samplesheet'
include { VALIDATE_REFERENCE } from './modules/validate_reference'
include { COLLECT_VERSIONS } from './modules/collect_versions'
include { FASTQ_VALIDATE } from './modules/fastq_validate'
include { FASTQC } from './modules/fastqc'
include { UMI_TEMPLATE_TRIM } from './modules/umi_extract'
include { CFSNV_STD_PREP } from './modules/cfsnv_stdprep'
include { CFSNV_CFDNA_PREP } from './modules/cfsnv_cfdnaprep'
include { UMI_GROUP_CONSENSUS } from './modules/umi_consensus'
include { CTDNA_MUTECT2; CTDNA_FILTER_MUTECT2; CTDNA_CFSNV_CALL } from './modules/ctdna_mutect2'
include { COMPARE_CTDNA_CALLSETS } from './modules/ctdna_compare'
include { CTDNA_VCF_TO_MAF } from './modules/ctdna_make_maf'
include { MERGE_CTDNA_TABLES as MERGE_CTDNA_CONCORDANCE } from './modules/ctdna_merge_tables'
include { MERGE_CTDNA_TABLES as MERGE_CTDNA_HIGH_SENSITIVITY } from './modules/ctdna_merge_tables'
include { MERGE_CTDNA_TABLES as MERGE_CTDNA_HIGH_CONFIDENCE } from './modules/ctdna_merge_tables'
include { MERGE_CTDNA_QC_TABLES as MERGE_CTDNA_SAMPLE_QC } from './modules/ctdna_merge_tables'
include { MERGE_CTDNA_QC_TABLES as MERGE_CTDNA_UMI_QC } from './modules/ctdna_merge_tables'
include { ALIGN_BWA_MEM2 } from './modules/align'
include { MARK_DUPLICATES } from './modules/markduplicates'
include { ALIGNMENT_QC } from './modules/alignment_qc'
include { MERGE_SAMPLE_QC } from './modules/merge_qc'
include { DEEPSOMATIC } from './modules/deepsomatic'
include { MUTECT2 } from './modules/mutect2'
include { FILTER_MUTECT } from './modules/filter_mutect'
include { DEEPVARIANT } from './modules/deepvariant'
include { SOMATIC_COMPARE } from './modules/somatic_compare'
include { CNVKIT } from './modules/cnvkit'
include { PURITY_PLOIDY } from './modules/purity_ploidy'
include { ARM_LEVEL_CNV } from './modules/arm_level_cnv'
include { MERGE_TABLES as MERGE_PURITY_PLOIDY } from './modules/merge_tables'
include { MERGE_TABLES as MERGE_ALLELE_SPECIFIC } from './modules/merge_tables'
include { MERGE_TABLES as MERGE_ARM_LEVEL_LONG } from './modules/merge_tables'
include { MERGE_TABLES as MERGE_MAF } from './modules/merge_tables'
include { MERGE_TABLES as MERGE_DRIVERS_LONG } from './modules/merge_tables'
include { MERGE_TABLES as MERGE_DRIVER_MATRIX } from './modules/merge_tables'
include { MERGE_TABLES as MERGE_TMB } from './modules/merge_tables'
include { MERGE_TABLES as MERGE_SIGNATURE_EXPOSURES } from './modules/merge_tables'
include { MERGE_TABLES as MERGE_CLONALITY_SUMMARY } from './modules/merge_tables'
include { MERGE_ARM_MATRICES } from './modules/merge_arm_matrices'
include { VEP as VEP_DEEPSOMATIC } from './modules/vep'
include { VEP as VEP_MUTECT2 } from './modules/vep'
include { VEP as VEP_CFSNV } from './modules/vep'
include { VEP as VEP_CTDNA_MUTECT2 } from './modules/vep'
include { MAKE_MAF } from './modules/make_maf'
include { DRIVER_ANNOTATION } from './modules/driver_annotation'
include { CALLABLE_TERRITORY } from './modules/callable_territory'
include { TMB } from './modules/tmb'
include { SIGNATURES } from './modules/signatures'
include { CLONALITY } from './modules/clonality'
include { GENOMIC_FEATURES } from './modules/genomic_features'
include { CONTAINER_MANIFEST } from './modules/report'
include { PIPELINE_PARAMETERS } from './modules/report'
include { MULTIQC_STUB } from './modules/report'
include { EXECUTION_REPORT } from './modules/report'

workflow {
    main:
    if (!params.samplesheet) {
        error 'Missing required parameter: --samplesheet'
    }

    sample_sheet = Channel.value(file(params.samplesheet, checkIfExists: true).toString())
    reference_path = file(params.reference_config, checkIfExists: true)
    reference_config = Channel.value(reference_path.toString())
    pipeline_config_path = file(params.pipeline_config, checkIfExists: true).toString()

    analysis_mode = readConfigValue(pipeline_config_path, 'analysis.mode')
    reference_fasta = readConfigValue(reference_path.toString(), 'reference.fasta')
    reference_fasta_path = file(reference_fasta, checkIfExists: true)
    reference_fasta_index_path = file("${reference_fasta}.fai", checkIfExists: true)
    reference_fasta_dict_path = file(reference_fasta.replaceFirst(/\.fa(sta)?$/, '.dict'), checkIfExists: true)
    reference_bwa_index_0123 = file("${reference_fasta}.0123", checkIfExists: true)
    reference_bwa_index_amb = file("${reference_fasta}.amb", checkIfExists: true)
    reference_bwa_index_ann = file("${reference_fasta}.ann", checkIfExists: true)
    reference_bwa_index_bwt = file("${reference_fasta}.bwt.2bit.64", checkIfExists: true)
    reference_bwa_index_pac = file("${reference_fasta}.pac", checkIfExists: true)
    target_intervals = readConfigValue(reference_path.toString(), 'intervals.interval_list')
    target_intervals_path = file(target_intervals, checkIfExists: true)
    germline_resource = readConfigValue(reference_path.toString(), 'gatk.germline_resource')
    germline_resource_path = file(germline_resource, checkIfExists: true)
    germline_resource_index_path = file("${germline_resource}.tbi", checkIfExists: true)
    common_snps = readConfigValue(reference_path.toString(), 'gatk.common_snps')
    common_snps_path = file(common_snps, checkIfExists: true)
    common_snps_index_path = file("${common_snps}.tbi", checkIfExists: true)
    panel_of_normals = readConfigValue(reference_path.toString(), 'gatk.panel_of_normals')
    panel_of_normals_path = file(panel_of_normals, checkIfExists: true)
    panel_of_normals_index_path = file("${panel_of_normals}.tbi", checkIfExists: true)
    vep_cache = readConfigValue(reference_path.toString(), 'vep.cache_dir')
    vep_cache_path = file(vep_cache, checkIfExists: true)
    reference_build = readConfigValue(reference_path.toString(), 'reference.build')
    chromosome_arms = readConfigValue(reference_path.toString(), 'annotations.chromosome_arms')
    chromosome_arms_path = file(chromosome_arms, checkIfExists: true)
    cancer_gene_census = readConfigValue(reference_path.toString(), 'annotations.cancer_gene_census')
    cancer_gene_census_path = file(cancer_gene_census, checkIfExists: true)
    hotspots = readConfigValue(reference_path.toString(), 'annotations.hotspots')
    hotspots_path = file(hotspots, checkIfExists: true)
    facets_cval = readConfigValue(pipeline_config_path, 'copy_number.purity_ploidy.facets.cval')
    facets_ndepth = readConfigValue(pipeline_config_path, 'copy_number.purity_ploidy.facets.ndepth')
    facets_ndepthmax = readConfigValue(pipeline_config_path, 'copy_number.purity_ploidy.facets.ndepthmax')
    callable_tumor_min_depth = readConfigValue(pipeline_config_path, 'callable.tumor_min_depth')
    callable_normal_min_depth = readConfigValue(pipeline_config_path, 'callable.normal_min_depth')
    callable_min_mapping_quality = readConfigValue(pipeline_config_path, 'callable.min_mapping_quality')
    callable_min_base_quality = readConfigValue(pipeline_config_path, 'callable.min_base_quality')

    validated_samplesheet = VALIDATE_SAMPLESHEET(sample_sheet)
    validated_reference = VALIDATE_REFERENCE(Channel.value(tuple(reference_path.toString(), pipeline_config_path)))
    provenance = COLLECT_VERSIONS()
    container_manifest = CONTAINER_MANIFEST(file("${projectDir}/config/containers.yaml").toString())
    pipeline_parameters = PIPELINE_PARAMETERS(file(params.samplesheet, checkIfExists: true).toString(), reference_path.toString())

    sample_fastqs = validated_samplesheet
        .splitCsv(header: true)
        .flatMap { row ->
            def sex = row.sex ?: 'NA'
            def ffpe = row.FFPE ? row.FFPE.toString().toLowerCase() : 'false'
            [
                tuple(
                    [patient_id: row.patient_id, sample_id: row.tumor_sample, sample_type: 'tumor', sample_role: 'plasma', sex: sex, ffpe: ffpe, tumor_sample: row.tumor_sample, normal_sample: row.normal_sample],
                    file(row.tumor_r1, checkIfExists: true),
                    file(row.tumor_r2, checkIfExists: true),
                    file(row.bait_bed, checkIfExists: true)
                ),
                tuple(
                    [patient_id: row.patient_id, sample_id: row.normal_sample, sample_type: 'normal', sample_role: 'wbc', sex: sex, ffpe: ffpe, tumor_sample: row.tumor_sample, normal_sample: row.normal_sample],
                    file(row.normal_r1, checkIfExists: true),
                    file(row.normal_r2, checkIfExists: true),
                    file(row.bait_bed, checkIfExists: true)
                )
            ]
        }

    if (analysis_mode == 'ctdna_umi') {
        ctdna_read_structure_r1 = readConfigValue(pipeline_config_path, 'ctdna.umi.read_structure_r1')
        ctdna_read_structure_r2 = readConfigValue(pipeline_config_path, 'ctdna.umi.read_structure_r2')
        ctdna_umi_tag_name = readConfigValue(pipeline_config_path, 'ctdna.umi.tag_name')
        ctdna_min_family_size = readConfigValue(pipeline_config_path, 'ctdna.umi.consensus.min_family_size')
        ctdna_error_rate_pre_umi = readConfigValue(pipeline_config_path, 'ctdna.umi.consensus.error_rate_pre_umi')
        ctdna_error_rate_post_umi = readConfigValue(pipeline_config_path, 'ctdna.umi.consensus.error_rate_post_umi')
        ctdna_min_hold = readConfigValue(pipeline_config_path, 'ctdna.cfsnv.min_hold')
        ctdna_min_pass = readConfigValue(pipeline_config_path, 'ctdna.cfsnv.min_pass')
        ctdna_emit_high_sensitivity = readConfigValue(pipeline_config_path, 'ctdna.outputs.emit_high_sensitivity').toLowerCase() == 'true'
        ctdna_emit_high_confidence = readConfigValue(pipeline_config_path, 'ctdna.outputs.emit_high_confidence').toLowerCase() == 'true'
        ctdna_blocked_positions_vcf = readConfigValue(reference_path.toString(), 'cfsnv.blocked_positions_vcf')
        ctdna_blocked_positions_vcf_path = file(ctdna_blocked_positions_vcf, checkIfExists: true)
        ctdna_blocked_positions_index_path = file("${ctdna_blocked_positions_vcf}.tbi", checkIfExists: true)

        fastq_validated = FASTQ_VALIDATE(sample_fastqs)
        fastqc_inputs = fastq_validated.map { meta, validation_tsv, r1, r2, bait_bed ->
            tuple(meta, r1, r2, bait_bed)
        }
        fastqc_results = FASTQC(fastqc_inputs)

        umi_template_trim_inputs = fastqc_results.map { meta, r1_html, r1_zip, r2_html, r2_zip, r1, r2, bait_bed ->
            tuple(meta, r1, r2, bait_bed, ctdna_read_structure_r1, ctdna_read_structure_r2)
        }
        template_trimmed = UMI_TEMPLATE_TRIM(umi_template_trim_inputs)

        ctdna_pairs = template_trimmed
            .map { meta, r1, r2, template_trim_qc, bait_bed ->
                tuple(meta.patient_id, [meta: meta, r1: r1, r2: r2, template_trim_qc: template_trim_qc, bait_bed: bait_bed])
            }
            .groupTuple()
            .map { patient_id, records ->
                def plasma = records.find { it.meta.sample_role == 'plasma' }
                def wbc = records.find { it.meta.sample_role == 'wbc' }
                if (!plasma || !wbc) {
                    error "Missing plasma/WBC pair for patient ${patient_id}"
                }
                def pairMeta = [
                    patient_id: patient_id,
                    tumor_sample: plasma.meta.sample_id,
                    normal_sample: wbc.meta.sample_id,
                    ffpe: plasma.meta.ffpe ?: 'false'
                ]
                tuple(pairMeta, plasma.r1, plasma.r2, wbc.r1, wbc.r2, plasma.bait_bed, plasma.template_trim_qc, wbc.template_trim_qc)
            }

        cfsnv_stdprep_inputs = ctdna_pairs.map { meta, plasma_r1, plasma_r2, wbc_r1, wbc_r2, bait_bed, plasma_template_trim_qc, wbc_template_trim_qc ->
            tuple(meta, plasma_r1, plasma_r2, wbc_r1, wbc_r2, bait_bed, reference_path.toString(), reference_fasta_path, reference_fasta_index_path, reference_fasta_dict_path, reference_bwa_index_0123, reference_bwa_index_amb, reference_bwa_index_ann, reference_bwa_index_bwt, reference_bwa_index_pac, common_snps_path, common_snps_index_path)
        }
        cfsnv_stdprep = CFSNV_STD_PREP(cfsnv_stdprep_inputs)

        cfsnv_cfdnaprep_inputs = ctdna_pairs.map { meta, plasma_r1, plasma_r2, wbc_r1, wbc_r2, bait_bed, plasma_template_trim_qc, wbc_template_trim_qc ->
            tuple(meta, plasma_r1, plasma_r2, bait_bed, reference_path.toString(), reference_fasta_path, reference_fasta_index_path, reference_fasta_dict_path, reference_bwa_index_0123, reference_bwa_index_amb, reference_bwa_index_ann, reference_bwa_index_bwt, reference_bwa_index_pac, common_snps_path, common_snps_index_path)
        }
        cfsnv_cfdnaprep = CFSNV_CFDNA_PREP(cfsnv_cfdnaprep_inputs)

        consensus_inputs = fastqc_results.map { meta, r1_html, r1_zip, r2_html, r2_zip, r1, r2, bait_bed ->
            tuple(meta, r1, r2, bait_bed, reference_path.toString(), reference_fasta_path, reference_fasta_index_path, reference_fasta_dict_path, reference_bwa_index_0123, reference_bwa_index_amb, reference_bwa_index_ann, reference_bwa_index_bwt, reference_bwa_index_pac, ctdna_umi_tag_name, ctdna_read_structure_r1, ctdna_read_structure_r2, ctdna_min_family_size, ctdna_error_rate_pre_umi, ctdna_error_rate_post_umi)
        }
        consensus_outputs = UMI_GROUP_CONSENSUS(consensus_inputs)

        consensus_pairs = consensus_outputs
            .map { meta, bam, bai, umi_qc, bait_bed, ref_cfg ->
                tuple(meta.patient_id, [meta: meta, bam: bam, bai: bai, umi_qc: umi_qc, bait_bed: bait_bed, ref_cfg: ref_cfg])
            }
            .groupTuple()
            .map { patient_id, records ->
                def plasma = records.find { it.meta.sample_role == 'plasma' }
                def wbc = records.find { it.meta.sample_role == 'wbc' }
                def pairMeta = [
                    patient_id: patient_id,
                    tumor_sample: plasma.meta.sample_id,
                    normal_sample: wbc.meta.sample_id,
                    ffpe: plasma.meta.ffpe ?: 'false'
                ]
                tuple(pairMeta, plasma, wbc)
            }

        ctdna_mutect2_inputs = consensus_pairs.map { meta, plasma, wbc ->
            tuple(meta, plasma.bam, plasma.bai, wbc.bam, wbc.bai, plasma.bait_bed, plasma.ref_cfg, reference_fasta_path, reference_fasta_index_path, reference_fasta_dict_path, germline_resource_path, germline_resource_index_path, common_snps_path, common_snps_index_path, panel_of_normals_path, panel_of_normals_index_path, target_intervals_path)
        }
        ctdna_mutect2_raw = CTDNA_MUTECT2(ctdna_mutect2_inputs)
        ctdna_mutect2_filtered = CTDNA_FILTER_MUTECT2(ctdna_mutect2_raw)

        ctdna_cfsnv_inputs = cfsnv_stdprep
            .join(cfsnv_cfdnaprep)
            .join(ctdna_pairs.map { meta, plasma_r1, plasma_r2, wbc_r1, wbc_r2, bait_bed, plasma_template_trim_qc, wbc_template_trim_qc -> tuple(meta, bait_bed, plasma_template_trim_qc) })
            .map { meta, plasma_std_bam, plasma_std_bai, wbc_std_bam, wbc_std_bai, sample_qc_tsv, plasma_extended_bam, plasma_extended_bai, plasma_notcombined_bam, plasma_notcombined_bai, bait_bed, plasma_template_trim_qc ->
                tuple(meta, plasma_std_bam, plasma_std_bai, wbc_std_bam, wbc_std_bai, plasma_extended_bam, plasma_extended_bai, plasma_notcombined_bam, plasma_notcombined_bai, bait_bed, reference_fasta_path, reference_fasta_index_path, reference_fasta_dict_path, common_snps_path, common_snps_index_path, ctdna_blocked_positions_vcf_path, ctdna_blocked_positions_index_path, ctdna_min_hold, ctdna_min_pass)
            }
        ctdna_cfsnv_calls = CTDNA_CFSNV_CALL(ctdna_cfsnv_inputs)

        ctdna_comparison_inputs = ctdna_cfsnv_calls
            .join(ctdna_mutect2_filtered)
            .map { meta, cfsnv_vcf, cfsnv_tbi, mutect_vcf, mutect_tbi, contamination, segments ->
                tuple(meta, cfsnv_vcf, mutect_vcf)
            }
        ctdna_comparison = COMPARE_CTDNA_CALLSETS(ctdna_comparison_inputs)

        cfsnv_vep_inputs = ctdna_cfsnv_calls.map { meta, vcf, tbi ->
            tuple(meta, vcf, reference_path.toString(), 'cfsnv', reference_fasta_path, reference_fasta_index_path, vep_cache_path)
        }
        ctdna_mutect2_vep_inputs = ctdna_mutect2_filtered.map { meta, vcf, tbi, contamination, segments ->
            tuple(meta, vcf, reference_path.toString(), 'ctdna_mutect2', reference_fasta_path, reference_fasta_index_path, vep_cache_path)
        }
        cfsnv_annotated = VEP_CFSNV(cfsnv_vep_inputs)
        ctdna_mutect2_annotated = VEP_CTDNA_MUTECT2(ctdna_mutect2_vep_inputs)

        plasma_umi_qc_by_patient = ctdna_pairs.map { meta, plasma_r1, plasma_r2, wbc_r1, wbc_r2, bait_bed, plasma_template_trim_qc, wbc_template_trim_qc ->
            tuple(meta.patient_id, plasma_template_trim_qc)
        }
        ctdna_maf_inputs = cfsnv_annotated
            .join(ctdna_comparison)
            .join(plasma_umi_qc_by_patient)
            .map { meta, annotated_vcf, comparison_tsv, umi_qc_tsv ->
                tuple(meta, annotated_vcf, comparison_tsv, umi_qc_tsv, ctdna_min_pass)
            }
        ctdna_maf_rows = CTDNA_VCF_TO_MAF(ctdna_maf_inputs)

        ctdna_sample_qc = MERGE_CTDNA_SAMPLE_QC(cfsnv_stdprep.map { meta, plasma_bam, plasma_bai, wbc_bam, wbc_bai, sample_qc -> sample_qc }.collect(), 'ctdna_sample_qc.tsv', false)
        ctdna_umi_qc = MERGE_CTDNA_UMI_QC(consensus_outputs.map { meta, bam, bai, umi_qc, bait_bed, ref_cfg -> umi_qc }.collect(), 'ctdna_umi_qc.tsv', false)
        ctdna_call_concordance = MERGE_CTDNA_CONCORDANCE(ctdna_comparison.map { meta, comparison_tsv -> comparison_tsv }.collect(), 'ctdna_call_concordance.tsv', false)

        if (ctdna_emit_high_sensitivity) {
            MERGE_CTDNA_HIGH_SENSITIVITY(ctdna_maf_rows.map { meta, hs, hc -> hs }.collect(), 'ctdna_mutations_high_sensitivity.maf.tsv', false)
        }
        if (ctdna_emit_high_confidence) {
            MERGE_CTDNA_HIGH_CONFIDENCE(ctdna_maf_rows.map { meta, hs, hc -> hc }.collect(), 'ctdna_mutations_high_confidence.maf.tsv', false)
        }
    } else {
        if (analysis_mode != 'tumor_normal_wes') {
            error "Unsupported analysis.mode '${analysis_mode}'. Expected tumor_normal_wes or ctdna_umi."
        }

        fastq_validated = FASTQ_VALIDATE(sample_fastqs)
        fastqc_inputs = fastq_validated.map { meta, validation_tsv, r1, r2, bait_bed ->
            tuple(meta, r1, r2, bait_bed)
        }
        fastqc_results = FASTQC(fastqc_inputs)

        align_inputs = fastqc_results.map { meta, r1_html, r1_zip, r2_html, r2_zip, r1, r2, bait_bed ->
            tuple(meta, r1, r2, bait_bed, reference_path.toString(), reference_fasta_path, reference_bwa_index_0123, reference_bwa_index_amb, reference_bwa_index_ann, reference_bwa_index_bwt, reference_bwa_index_pac)
        }
        aligned = ALIGN_BWA_MEM2(align_inputs)

        marked_duplicates = MARK_DUPLICATES(aligned)

        qc_inputs = marked_duplicates.map { meta, bam, bai, dup_metrics, bait_bed, ref_cfg ->
            tuple(meta, bam, bai, dup_metrics, bait_bed, ref_cfg, reference_fasta_path, reference_fasta_index_path, reference_fasta_dict_path, target_intervals_path)
        }
        alignment_qc = ALIGNMENT_QC(qc_inputs)
        sample_qc = MERGE_SAMPLE_QC(alignment_qc.map { meta, qc_tsv -> qc_tsv }.collect())

        paired_bams = marked_duplicates
            .map { meta, bam, bai, dup_metrics, bait_bed, ref_cfg ->
                tuple(meta.patient_id, [meta: meta, bam: bam, bai: bai, bait_bed: bait_bed, ref_cfg: ref_cfg])
            }
            .groupTuple()
            .map { patient_id, records ->
                def tumor = records.find { it.meta.sample_type == 'tumor' }
                def normal = records.find { it.meta.sample_type == 'normal' }
                if (!tumor || !normal) {
                    error "Missing tumor/normal pair for patient ${patient_id}"
                }
                def pairMeta = [
                    patient_id: patient_id,
                    tumor_sample: tumor.meta.sample_id,
                    normal_sample: normal.meta.sample_id,
                    ffpe: tumor.meta.ffpe ?: 'false'
                ]
                tuple(pairMeta, tumor.bam, tumor.bai, normal.bam, normal.bai, tumor.bait_bed, tumor.ref_cfg)
            }

        deepsomatic_inputs = paired_bams.map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, bait_bed, ref_cfg ->
            tuple(meta, tumor_bam, tumor_bai, normal_bam, normal_bai, bait_bed, ref_cfg, reference_fasta_path, reference_fasta_index_path)
        }
        deepsomatic_calls = DEEPSOMATIC(deepsomatic_inputs)

        mutect2_inputs = paired_bams.map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, bait_bed, ref_cfg ->
            tuple(meta, tumor_bam, tumor_bai, normal_bam, normal_bai, bait_bed, ref_cfg, reference_fasta_path, reference_fasta_index_path, reference_fasta_dict_path, germline_resource_path, germline_resource_index_path, common_snps_path, common_snps_index_path, panel_of_normals_path, panel_of_normals_index_path, target_intervals_path)
        }
        mutect2_raw = MUTECT2(mutect2_inputs)
        mutect2_filtered = FILTER_MUTECT(mutect2_raw)

        deepvariant_inputs = marked_duplicates
            .filter { meta, bam, bai, dup_metrics, bait_bed, ref_cfg -> meta.sample_type == 'normal' }
            .map { meta, bam, bai, dup_metrics, bait_bed, ref_cfg ->
                tuple(meta, bam, bai, bait_bed, ref_cfg, reference_fasta_path, reference_fasta_index_path)
            }
        deepvariant_calls = DEEPVARIANT(deepvariant_inputs)

        compare_inputs = deepsomatic_calls
            .join(mutect2_filtered)
            .map { meta, ds_vcf, ds_tbi, m2_vcf, m2_tbi, contamination, segments ->
                tuple(meta, ds_vcf, m2_vcf)
            }
        somatic_comparison = SOMATIC_COMPARE(compare_inputs)

        cnvkit_inputs = paired_bams.map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, bait_bed, ref_cfg ->
            tuple(meta, tumor_bam, tumor_bai, normal_bam, normal_bai, bait_bed, ref_cfg, reference_fasta_path, reference_fasta_index_path)
        }
        cnvkit_calls = CNVKIT(cnvkit_inputs)

        purity_ploidy_inputs = paired_bams.map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, bait_bed, ref_cfg ->
            tuple(meta, tumor_bam, tumor_bai, normal_bam, normal_bai, bait_bed, ref_cfg, pipeline_config_path, common_snps_path, common_snps_index_path, reference_build, facets_cval, facets_ndepth, facets_ndepthmax)
        }
        purity_ploidy = PURITY_PLOIDY(purity_ploidy_inputs)

        arm_level_inputs = cnvkit_calls.map { meta, targetcoverage, antitargetcoverage, cnr, cns, call_cns ->
            tuple(meta, call_cns, reference_path.toString(), chromosome_arms_path)
        }
        arm_level_cnv = ARM_LEVEL_CNV(arm_level_inputs)

        cohort_purity_ploidy = MERGE_PURITY_PLOIDY(purity_ploidy.map { meta, purity_tsv, allele_tsv -> purity_tsv }.collect(), 'purity_ploidy.tsv', false)
        cohort_allele_specific = MERGE_ALLELE_SPECIFIC(purity_ploidy.map { meta, purity_tsv, allele_tsv -> allele_tsv }.collect(), 'allele_specific_cnv.tsv', false)
        arm_level_long = MERGE_ARM_LEVEL_LONG(arm_level_cnv.map { meta, long_tsv, matrix_tsv -> long_tsv }.collect(), 'arm_level_cnv_long.tsv', false)
        arm_level_matrix = MERGE_ARM_MATRICES(arm_level_cnv.map { meta, long_tsv, matrix_tsv -> matrix_tsv }.collect())

        deepsomatic_vep_inputs = deepsomatic_calls.map { meta, vcf, tbi ->
            tuple(meta, vcf, reference_path.toString(), 'deepsomatic', reference_fasta_path, reference_fasta_index_path, vep_cache_path)
        }
        mutect2_vep_inputs = mutect2_filtered.map { meta, vcf, tbi, contamination, segments ->
            tuple(meta, vcf, reference_path.toString(), 'mutect2', reference_fasta_path, reference_fasta_index_path, vep_cache_path)
        }
        deepsomatic_annotated = VEP_DEEPSOMATIC(deepsomatic_vep_inputs)
        mutect2_annotated = VEP_MUTECT2(mutect2_vep_inputs)

        maf_inputs = deepsomatic_annotated
            .join(somatic_comparison)
            .map { meta, annotated_vcf, comparison_tsv, counts_tsv ->
                tuple(meta, annotated_vcf, comparison_tsv)
            }
        maf_rows = MAKE_MAF(maf_inputs)

        driver_inputs = maf_rows.map { meta, maf ->
            tuple(meta, maf, reference_path.toString(), cancer_gene_census_path, hotspots_path)
        }
        drivers = DRIVER_ANNOTATION(driver_inputs)

        cohort_maf = MERGE_MAF(maf_rows.map { meta, maf -> maf }.collect(), 'somatic_mutations.maf', false)
        cohort_drivers_long = MERGE_DRIVERS_LONG(drivers.map { meta, long_tsv, matrix_tsv -> long_tsv }.collect(), 'drivers_long.tsv', false)
        cohort_driver_matrix = MERGE_DRIVER_MATRIX(drivers.map { meta, long_tsv, matrix_tsv -> matrix_tsv }.collect(), 'driver_matrix.tsv', true)

        callable_inputs = paired_bams.map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, bait_bed, ref_cfg ->
            tuple(meta, tumor_bam, tumor_bai, normal_bam, normal_bai, bait_bed, callable_tumor_min_depth, callable_normal_min_depth, callable_min_mapping_quality, callable_min_base_quality)
        }
        callable_regions = CALLABLE_TERRITORY(callable_inputs)

        tmb_inputs = maf_rows
            .join(callable_regions)
            .map { meta, maf, callable_bed, callable_mb ->
                tuple(meta, maf, callable_mb)
            }
        tmb_tables = TMB(tmb_inputs)

        signature_inputs = maf_rows.map { meta, maf ->
            tuple(meta, maf, reference_path.toString(), file(params.pipeline_config, checkIfExists: true).toString())
        }
        signature_tables = SIGNATURES(signature_inputs)

        clonality_inputs = maf_rows
            .join(purity_ploidy)
            .map { meta, maf, purity_tsv, allele_tsv ->
                tuple(meta, maf, purity_tsv, allele_tsv, file(params.pipeline_config, checkIfExists: true).toString())
            }
        clonality_tables = CLONALITY(clonality_inputs)

        cohort_tmb = MERGE_TMB(tmb_tables.map { meta, tmb -> tmb }.collect(), 'tmb.tsv', false)
        cohort_signature_exposure = MERGE_SIGNATURE_EXPOSURES(signature_tables.map { meta, exposure, qc -> exposure }.collect(), 'signature_exposures.tsv', false)
        cohort_clonality = MERGE_CLONALITY_SUMMARY(clonality_tables.map { meta, ccf, summary -> summary }.collect(), 'clonality_summary.tsv', false)

        genomic_features = GENOMIC_FEATURES(
            tmb_tables.map { meta, tmb -> tmb }.collect(),
            purity_ploidy.map { meta, purity_tsv, allele_tsv -> purity_tsv }.collect(),
            signature_tables.map { meta, exposure, qc -> exposure }.collect(),
            clonality_tables.map { meta, ccf, summary -> summary }.collect(),
            drivers.map { meta, long_tsv, matrix_tsv -> matrix_tsv }.collect(),
            arm_level_cnv.map { meta, long_tsv, matrix_tsv -> matrix_tsv }.collect()
        )

        multiqc_report = MULTIQC_STUB(sample_qc)
        execution_report = EXECUTION_REPORT(
            sample_qc,
            cohort_maf,
            cohort_tmb,
            cohort_purity_ploidy,
            cohort_signature_exposure,
            cohort_clonality,
            cohort_driver_matrix,
            arm_level_matrix,
            genomic_features,
            provenance,
            container_manifest,
            validated_reference,
            pipeline_parameters
        )
    }
}
