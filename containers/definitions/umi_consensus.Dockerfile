FROM mambaorg/micromamba:2.3.2

USER root

RUN micromamba create -y -n umi -c conda-forge -c bioconda \
    openjdk=17 \
    bwa-mem2=2.3 \
    samtools \
    fgbio=2.5.1 \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/envs/umi/bin:$PATH
