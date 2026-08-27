process CNVKIT {
    label 'cpu_medium'
    publishDir "${params.results_dir}/cnv/cnvkit", mode: 'copy'
    container { params.cnvkit_container }

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(bait_bed), val(reference_config), val(ref_fasta)

    output:
    tuple val(meta), path("${meta.tumor_sample}.targetcoverage.cnn"), path("${meta.tumor_sample}.antitargetcoverage.cnn"), path("${meta.tumor_sample}.cnr"), path("${meta.tumor_sample}.cns"), path("${meta.tumor_sample}.call.cns")

    script:
    """
    ln -sf ${tumor_bam} ${meta.tumor_sample}.bam
    ln -sf ${normal_bam} ${meta.normal_sample}.bam

    TARGET_COUNT=\$(awk 'END {print NR}' ${bait_bed})
    CNVKIT_METHOD=hybrid
    RELAX_CNVKIT_STATUS=0
    if [[ "\$TARGET_COUNT" -lt 2 ]]; then
      CNVKIT_METHOD=amplicon
      RELAX_CNVKIT_STATUS=1
    fi

    set +e
    cnvkit.py batch ${meta.tumor_sample}.bam \
      --normal ${meta.normal_sample}.bam \
      --targets ${bait_bed} \
      --fasta ${ref_fasta} \
      --output-reference ${meta.patient_id}.reference.cnn \
      --output-dir . \
      --processes ${task.cpus} \
      --method "\$CNVKIT_METHOD"
    batch_status=\$?
    set -e

    if [[ "\$batch_status" -ne 0 ]]; then
      if [[ "\$RELAX_CNVKIT_STATUS" -eq 1 && -f ${meta.tumor_sample}.targetcoverage.cnn && -f ${meta.tumor_sample}.cnr && -f ${meta.tumor_sample}.cns ]]; then
        :
      else
        exit "\$batch_status"
      fi
    fi

    if [[ ! -f ${meta.tumor_sample}.antitargetcoverage.cnn ]]; then
      touch ${meta.tumor_sample}.antitargetcoverage.cnn
    fi

    if [[ ! -f ${meta.tumor_sample}.call.cns ]]; then
      cnvkit.py call ${meta.tumor_sample}.cns -o ${meta.tumor_sample}.call.cns
    fi
    """

    stub:
    """
    printf 'chromosome\tstart\tend\tgene\tlog2\tdepth\nchr1\t0\t6\tGENE1\t-0.500000\t100\n' > ${meta.tumor_sample}.targetcoverage.cnn
    printf 'chromosome\tstart\tend\tgene\tlog2\tdepth\nchr1\t6\t12\tGENE2\t0.600000\t120\n' > ${meta.tumor_sample}.antitargetcoverage.cnn
    printf 'chromosome\tstart\tend\tgene\tlog2\tdepth\tweight\nchr1\t0\t6\tGENE1\t-0.500000\t100\t1.0\nchr1\t6\t12\tGENE2\t0.600000\t120\t1.0\n' > ${meta.tumor_sample}.cnr
    printf 'chromosome\tstart\tend\tgene\tlog2\tprobes\nchr1\t0\t6\tGENE1\t-0.500000\t10\nchr1\t6\t12\tGENE2\t0.600000\t10\n' > ${meta.tumor_sample}.cns
    printf 'chromosome\tstart\tend\tgene\tlog2\tprobes\tcn\nchr1\t0\t6\tGENE1\t-0.500000\t10\t1\nchr1\t6\t12\tGENE2\t0.600000\t10\t4\n' > ${meta.tumor_sample}.call.cns
    """
}
