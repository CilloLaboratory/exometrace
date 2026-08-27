process CALLABLE_TERRITORY {
    label 'cpu_small'
    publishDir "${params.results_dir}/callable", mode: 'copy'
    container { params.qc_container }

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(bait_bed), val(tumor_min_depth), val(normal_min_depth), val(min_mapping_quality), val(min_base_quality)

    output:
    tuple val(meta), path("${meta.patient_id}.callable.bed"), path("${meta.patient_id}.callable_mb.txt")

    script:
    """
    samtools depth \
      -a \
      -b ${bait_bed} \
      -q ${min_mapping_quality} \
      -Q ${min_base_quality} \
      ${tumor_bam} ${normal_bam} > ${meta.patient_id}.depth.tsv

    python3 ${projectDir}/scripts/build_callable_bed.py \
      --patient-id ${meta.patient_id} \
      --depth-tsv ${meta.patient_id}.depth.tsv \
      --tumor-min-depth ${tumor_min_depth} \
      --normal-min-depth ${normal_min_depth} \
      --output-bed ${meta.patient_id}.callable.bed \
      --output-mb ${meta.patient_id}.callable_mb.txt
    """

    stub:
    """
    cat > ${meta.patient_id}.depth.tsv <<'EOF'
    chr1	1	25	15
    chr1	2	25	15
    chr1	3	10	15
    chr1	4	25	15
    chr1	5	25	5
    chr1	6	25	15
    EOF
    python3 ${projectDir}/scripts/build_callable_bed.py \
      --patient-id ${meta.patient_id} \
      --depth-tsv ${meta.patient_id}.depth.tsv \
      --tumor-min-depth 20 \
      --normal-min-depth 10 \
      --output-bed ${meta.patient_id}.callable.bed \
      --output-mb ${meta.patient_id}.callable_mb.txt
    """
}
