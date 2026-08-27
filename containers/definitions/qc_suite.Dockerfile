FROM quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0 AS fastqc

FROM quay.io/biocontainers/samtools:1.23--h96c455f_0 AS samtools

FROM quay.io/biocontainers/mosdepth:0.3.11--h0ec343a_1 AS mosdepth

FROM broadinstitute/gatk:4.7.0.0

RUN python3 -m pip install --no-cache-dir PyYAML

COPY --from=fastqc /usr/local/opt/fastqc-0.12.1 /usr/local/opt/fastqc-0.12.1
COPY --from=samtools /usr/local/bin/samtools /usr/local/bin/samtools
COPY --from=samtools /usr/local/lib /usr/local/lib
COPY --from=mosdepth /usr/local/bin/mosdepth /usr/local/bin/mosdepth

RUN ln -sf /usr/local/opt/fastqc-0.12.1/fastqc /usr/local/bin/fastqc
