FROM mambaorg/micromamba:2.3.2

USER root

RUN micromamba install -y -n base -c conda-forge openjdk=8 \
    && micromamba clean --all --yes

RUN micromamba create -y -n cfsnv -c conda-forge -c bioconda \
    r-base=4.0 \
    r-rcpp \
    r-reticulate \
    python=3.8 \
    pip \
    samtools=1.11 \
    bedtools=2.30.0 \
    bwa=0.7.17 \
    tabix \
    flash=1.2.11 \
    curl \
    make \
    boost-cpp \
    cxx-compiler \
    gzip \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/envs/cfsnv/bin:$PATH

RUN python -m pip install --no-cache-dir \
    numpy==1.19.5 \
    pandas==1.1.5 \
    scipy==1.5.4 \
    scikit-learn==0.24.1

RUN mkdir -p /usr/local/share/cfsnv-tools /opt/gatk3 \
    && curl -fsSL -o /usr/local/share/cfsnv-tools/picard.jar \
      https://github.com/broadinstitute/picard/releases/download/2.18.4/picard.jar \
    && /opt/conda/bin/jar tf /usr/local/share/cfsnv-tools/picard.jar \
      | grep -q 'org/broadinstitute/barclay/argparser/CommandLineProgramProperties.class' \
    && curl -fsSL -o /tmp/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 \
      https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 \
    && tar -xjf /tmp/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 -C /opt/gatk3 --strip-components=1 \
    && ln -s /opt/gatk3/GenomeAnalysisTK.jar /usr/local/share/cfsnv-tools/GenomeAnalysisTK.jar \
    && rm -f /tmp/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2

ENV CFSNV_PICARD_JAR=/usr/local/share/cfsnv-tools/picard.jar
ENV CFSNV_GATK_JAR=/usr/local/share/cfsnv-tools/GenomeAnalysisTK.jar
ENV JAVA_HOME=/opt/conda

RUN curl -fsSL -o /tmp/cfSNV_0.99.0.tar.gz \
      https://github.com/jasminezhoulab/cfSNV/raw/main/cfSNV_0.99.0.tar.gz \
    && R CMD INSTALL /tmp/cfSNV_0.99.0.tar.gz \
    && rm -f /tmp/cfSNV_0.99.0.tar.gz
