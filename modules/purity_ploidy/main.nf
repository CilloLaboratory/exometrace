process PURITY_PLOIDY {
    label 'cpu_medium'
    publishDir "${params.results_dir}/cnv/purity_ploidy", mode: 'copy'
    container { params.facets_container }

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(bait_bed), val(reference_config), val(pipeline_config), val(common_snps), val(genome_build), val(facets_cval), val(facets_ndepth), val(facets_ndepthmax)

    output:
    tuple val(meta), path("${meta.patient_id}.purity_ploidy.tsv"), path("${meta.patient_id}.allele_specific_cnv.tsv")

    script:
    """
    snp-pileup-wrapper.R \
      --vcf-file ${common_snps} \
      --normal-bam ${normal_bam} \
      --tumor-bam ${tumor_bam} \
      --output-prefix ${meta.patient_id}

    pileup_rows=\$(gzip -cd ${meta.patient_id}.snp_pileup.gz | wc -l)
    if [[ "\$pileup_rows" -le 1 ]]; then
      printf 'patient_id\\tpurity\\tploidy\\tdiplogr\\n${meta.patient_id}\\tNA\\tNA\\tNA\\n' > ${meta.patient_id}.purity_ploidy.tsv
      printf 'patient_id\\tchromosome\\tstart\\tend\\tlog2\\ttotal_cn\\tmajor_cn\\tminor_cn\\tsegment_cf\\tloh\\n' > ${meta.patient_id}.allele_specific_cnv.tsv
      touch ${meta.patient_id}.facets.rds
    else
      Rscript ${projectDir}/scripts/run_facets.R \
        --counts-file ${meta.patient_id}.snp_pileup.gz \
        --sample-id ${meta.patient_id} \
        --genome-build ${genome_build} \
        --cval ${facets_cval} \
        --ndepth ${facets_ndepth} \
        --ndepthmax ${facets_ndepthmax} \
        --output-purity ${meta.patient_id}.purity_ploidy.tsv \
        --output-allele-specific ${meta.patient_id}.allele_specific_cnv.tsv \
        --output-rds ${meta.patient_id}.facets.rds
    fi
    """

    stub:
    """
    cat > ${meta.patient_id}.purity_ploidy.tsv <<'EOF'
    patient_id	purity	ploidy	diplogr
    ${meta.patient_id}	0.720000	2.800000	-0.180000
    EOF
    cat > ${meta.patient_id}.allele_specific_cnv.tsv <<'EOF'
    patient_id	chromosome	start	end	log2	total_cn	major_cn	minor_cn	segment_cf	loh
    ${meta.patient_id}	chr1	0	6	-0.500000	1	1	0	0.720000	true
    ${meta.patient_id}	chr1	6	12	0.600000	4	3	1	0.720000	false
    EOF
    """
}
