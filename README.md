# Tumor/Normal WES Pipeline

This repository contains a Nextflow DSL2 workflow for paired tumor/normal whole-exome sequencing analysis with an emphasis on ctDNA.

The workflow currently covers:

- FASTQ validation and FastQC
- BWA-MEM2 alignment
- duplicate marking and alignment QC
- DeepSomatic somatic calling
- GATK Mutect2 somatic calling plus filtering
- DeepVariant germline calling on the matched normal
- CNVkit copy-number calling
- purity/ploidy summary export
- arm-level CNV summaries
- VEP annotation
- MAF-like somatic table generation
- driver table generation
- callable territory generation from depth-filtered tumor/normal BAM coverage
- TMB, signatures, clonality, and genomic feature matrix generation
- provenance manifests and cohort-facing HTML reports

## Scope

This README is written for running the current implementation on real samples.

Important limitations:

- The workflow executes real command lines in normal mode, but several downstream modules remain simplified:
  - the HTML reporting layer is static and file-based rather than an interactive application
- Container handling is profile-dependent:
  - generic `local` and `slurm` runs still assume required executables are available in `PATH`
  - the `pitt_htc` profile attaches cached local Singularity `.sif` paths to the major tool processes and loads `singularity`

If you need strict production isolation outside the Pitt profile, extend the same container mapping approach to your site-specific profile.

## Prerequisites

Required runtime:

- Linux
- Java 17+
- Nextflow 26.04.6 or compatible
- Python 3 with `PyYAML`

This repository already includes a local launcher:

```bash
./bin/nextflow -version
```

Required command-line tools in `PATH` for a real run:

- `fastqc`
- `bwa-mem2`
- `samtools`
- `gatk`
- `mosdepth`
- `bcftools`
- `run_deepsomatic`
- `run_deepvariant`
- `cnvkit.py`
- `vep`
- `SigProfilerAssignment` with the required SigProfiler genome/reference assets
- `snp-pileup-wrapper.R`, `Rscript`, and FACETS/facets-suite for purity/ploidy and clonality

Optional but expected in HPC environments:

- `singularity`
- SLURM if using cluster execution

## Required Inputs

### 1. Sample sheet

Provide a CSV with one row per patient pair.

Required columns:

- `patient_id`
- `tumor_sample`
- `normal_sample`
- `tumor_r1`
- `tumor_r2`
- `normal_r1`
- `normal_r2`
- `bait_bed`

Optional columns currently recognized:

- `sex`
- `FFPE`

Example:

```csv
patient_id,tumor_sample,normal_sample,tumor_r1,tumor_r2,normal_r1,normal_r2,bait_bed,sex,FFPE
P001,P001_T,P001_N,/data/P001_T_R1.fastq.gz,/data/P001_T_R2.fastq.gz,/data/P001_N_R1.fastq.gz,/data/P001_N_R2.fastq.gz,/refs/exome_targets.bed,F,false
```

Validation command:

```bash
python3 scripts/validate_samplesheet.py --samplesheet /path/to/samplesheet.csv --check-files
```

### 2. Reference bundle

For a real run, bootstrap the Broad-style GRCh38 bundle with:

```bash
export TARGETS_BED=/absolute/path/to/assay_targets.bed
export INTERVAL_LIST=/absolute/path/to/assay_targets.interval_list
export CANCER_GENE_CENSUS=/absolute/path/to/cancer_gene_census.tsv
export HOTSPOTS_TSV=/absolute/path/to/hotspots.tsv
export GENE_COORDINATES_TSV=/absolute/path/to/gene_coordinates.tsv
export CHROMOSOME_ARMS_BED=/absolute/path/to/chromosome_arms.bed
bash scripts/bootstrap_references.sh
```

The pinned source manifest is [config/reference_sources.yaml](/home/arc85/Desktop/wes_workflow/config/reference_sources.yaml). The bootstrap script populates the tree described by [config/references.yaml](/home/arc85/Desktop/wes_workflow/config/references.yaml):

```text
references/
└── GRCh38/
    ├── fasta/
    │   ├── GRCh38.fa
    │   ├── GRCh38.fa.fai
    │   └── GRCh38.dict
    ├── intervals/
    │   ├── exome_targets.bed
    │   ├── exome_targets.interval_list
    │   └── callable_regions.bed
    ├── gatk/
    │   ├── af-only-gnomad.vcf.gz
    │   ├── common_biallelic_snps.vcf.gz
    │   └── panel_of_normals.vcf.gz
    ├── cfsnv/
    │   └── blocked_positions.vcf.gz
    ├── vep/
    │   └── cache/
    └── annotations/
        ├── cancer_gene_census.tsv
        ├── hotspots.tsv
        ├── gene_coordinates.tsv
        └── chromosome_arms.bed
```

Reference validation command:

```bash
python3 scripts/validate_reference.py --config config/references.yaml
```

The validator now checks real-run prerequisites, including:

- BWA-MEM2 sidecar indexes for the configured FASTA
- indexed GATK resource VCFs
- derived cfSNV blocked-position blacklist restricted to the configured assay targets
- populated offline VEP cache
- provisioned SigProfiler asset volume
- non-placeholder interval and annotation tables suitable for a real WES run

### 3. Tool configuration

Review these files before a real run:

- [config/default.yaml](/home/arc85/Desktop/wes_workflow/config/default.yaml)
- [config/references.yaml](/home/arc85/Desktop/wes_workflow/config/references.yaml)
- [config/containers.yaml](/home/arc85/Desktop/wes_workflow/config/containers.yaml)
- [nextflow.config](/home/arc85/Desktop/wes_workflow/nextflow.config)

## Pinned Software Versions

Pinned versions were recorded on August 19, 2026 in [config/containers.yaml](/home/arc85/Desktop/wes_workflow/config/containers.yaml):

- `bwa-mem2` 2.3
- `FastQC` 0.12.1
- `samtools` 1.23
- `mosdepth` 0.3.11
- `GATK` 4.7.0.0
- `DeepVariant` 1.10.0
- `DeepSomatic` 1.10.0
- `CNVkit` 0.9.14
- `Ensembl VEP` 116.0
- `SigProfilerAssignment` 1.1.5
- `FACETS` / `facets-suite` 0.6.2

Upstream sources:

- https://github.com/bwa-mem2/bwa-mem2/releases
- https://github.com/s-andrews/FastQC
- https://github.com/samtools/samtools/releases
- https://github.com/brentp/mosdepth
- https://github.com/broadinstitute/gatk/releases
- https://github.com/google/deepvariant/releases
- https://github.com/google/deepsomatic/releases
- https://github.com/etal/cnvkit/releases
- https://github.com/Ensembl/ensembl-vep/releases
- https://github.com/SigProfilerSuite/SigProfilerAssignment
- https://github.com/mskcc/facets
- https://github.com/mskcc/facets-suite

## Environment Preparation For Real Runs

You have two viable ways to run this workflow today.

### Option A: Native tools in `PATH`

Install all required executables in the active environment and run Nextflow directly.

This is the simplest approach with the current codebase.

### Option B: Wrapper-based container execution

If your site requires Singularity/Apptainer, provide shell wrappers or modulefiles so these commands resolve in `PATH`:

- `fastqc`
- `bwa-mem2`
- `samtools`
- `gatk`
- `mosdepth`
- `bcftools`
- `run_deepsomatic`
- `run_deepvariant`
- `cnvkit.py`
- `vep`
- `SigProfilerAssignment`

Those wrappers can call your pinned SIF images with the correct binds and `--nv` where needed.
SigProfiler also needs the provisioned asset volume under `references/GRCh38/sigprofiler` or an equivalent configured path.

### Option C: Pitt HTC profile

If you are running on the University of Pittsburgh CRCD HTC cluster, use the built-in `pitt_htc` profile and the submission wrapper in [scripts/pitt_htc_nextflow.sbatch](/home/arc85/Desktop/wes_workflow/scripts/pitt_htc_nextflow.sbatch). Details are in [docs/pitt_htc.md](/home/arc85/Desktop/wes_workflow/docs/pitt_htc.md).

If you are running on this workstation with cached SIF images under `/home/arc85/Desktop/wes_singularity_cache`, first materialize local sandbox containers:

```bash
bash scripts/prepare_local_singularity_sandboxes.sh
```

If the ctDNA-specific cache entries are missing, build them first:

```bash
bash scripts/build_ctdna_containers.sh
```

Then run with the machine-local profile:

```bash
NXF_OFFLINE=true NXF_DISABLE_CHECK_LATEST=true ./bin/nextflow run main.nf \
  --samplesheet /absolute/path/to/samplesheet.csv \
  --reference_config /absolute/path/to/config/references.yaml \
  -profile arc85_workstation
```

### Rerun and resume

For reruns, prefer `-resume` and keep the launch inputs stable:

```bash
NXF_OFFLINE=true NXF_DISABLE_CHECK_LATEST=true ./bin/nextflow run main.nf \
  --samplesheet /absolute/path/to/samplesheet.csv \
  --reference_config /absolute/path/to/config/references.yaml \
  --pipeline_config /absolute/path/to/config/default.yaml \
  -profile arc85_workstation \
  -resume
```

To maximize cache reuse:

- run from the same repository checkout
- keep the same `work/` directory
- keep the same profile
- keep the same absolute paths for `--samplesheet`, `--reference_config`, and `--pipeline_config`
- keep the same container versions and sandbox/SIF paths

What should resume reliably now:

- `VALIDATE_SAMPLESHEET`
- `VALIDATE_REFERENCE`
- `FASTQ_VALIDATE`
- `FASTQC`
- `UMI_TEMPLATE_TRIM`

Those early immutable steps use `cache 'deep'`, so unchanged input file contents can be reused even if file mtimes changed.

What still invalidates cache by design:

- editing `main.nf` or any module script
- changing container paths or versions
- changing config values that feed a process
- deleting `work/` or `.nextflow/`
- switching profiles or running from a different clone/location

If a single task fails, rerun with `-resume` after patching. Nextflow should regenerate the failed task and any downstream tasks while reusing successful upstream results.

For a real end-to-end smoke test on this machine, use the tiny fixture sheet with real references:

```bash
NXF_OFFLINE=true NXF_DISABLE_CHECK_LATEST=true ./bin/nextflow run main.nf \
  --samplesheet tests/test_data/samplesheet_local.csv \
  --reference_config config/references.yaml \
  --results_dir results/smoke_tiny_real \
  -profile arc85_workstation
```

## Running The Workflow

### Local real run

```bash
NXF_OFFLINE=true NXF_DISABLE_CHECK_LATEST=true ./bin/nextflow run main.nf \
  --samplesheet /absolute/path/to/samplesheet.csv \
  --reference_config /absolute/path/to/config/references.yaml \
  -profile local
```

### SLURM run

The repository contains a basic SLURM profile in [config/profiles/slurm.config](/home/arc85/Desktop/wes_workflow/config/profiles/slurm.config).

Run:

```bash
NXF_OFFLINE=true NXF_DISABLE_CHECK_LATEST=true ./bin/nextflow run main.nf \
  --samplesheet /absolute/path/to/samplesheet.csv \
  --reference_config /absolute/path/to/config/references.yaml \
  -profile slurm
```

### Pitt HTC run

Use the CRCD batch wrapper:

```bash
sbatch scripts/pitt_htc_nextflow.sbatch \
  /absolute/path/to/samplesheet.csv \
  /absolute/path/to/config/references.yaml
```

Optional Pitt-specific flags can be appended after the two required positional arguments:

```bash
sbatch scripts/pitt_htc_nextflow.sbatch \
  /absolute/path/to/samplesheet.csv \
  /absolute/path/to/config/references.yaml \
  --pitt_account my_allocation \
  --pitt_gpu_partition l40s \
  --pitt_constraint amd,genoa
```

What the `pitt_htc` profile does:

- submits all tasks to `--cluster=htc` and partition `htc`
- enables `singularity`
- loads `singularity/3.9.6`
- uses Nextflow scratch execution for task work directories
- keeps CPU-labeled tasks on HTC and routes GPU-labeled DeepSomatic and DeepVariant tasks to the Pitt GPU cluster
- assigns pinned cached local `.sif` paths for the major tool processes and does not download images on demand

### GPU-sensitive stages

The workflow labels DeepSomatic and DeepVariant with the `gpu` resource label. In the `pitt_htc` profile, those tasks are submitted to the Pitt GPU cluster with `--nv` container execution. Outside that profile, GPU use still depends on your environment exposing working `run_deepsomatic` and `run_deepvariant` commands with GPU-capable runtimes.

### Useful custom parameters

Defined in [nextflow.config](/home/arc85/Desktop/wes_workflow/nextflow.config):

- `--samplesheet`
- `--reference_config`

Additional runtime params with defaults:

- `params.copy_number_arm_deep_del = -1.1`
- `params.copy_number_arm_loss = -0.3`
- `params.copy_number_arm_gain = 0.2`
- `params.copy_number_arm_amp = 0.7`
- `params.pitt_account = null`
- `params.pitt_qos = null`
- `params.pitt_constraint = null`
- `params.pitt_gpu_partition = 'l40s'`
- `params.pitt_gpu_gres = 'gpu:1'`

Example:

```bash
NXF_OFFLINE=true NXF_DISABLE_CHECK_LATEST=true ./bin/nextflow run main.nf \
  --samplesheet /absolute/path/to/samplesheet.csv \
  --reference_config /absolute/path/to/config/references.yaml \
  -profile local
```

## Expected Outputs

The workflow writes into [results](/home/arc85/Desktop/wes_workflow/results).

Key outputs:

- `results/qc/fastqc/`
- `results/qc/alignment/`
- `results/qc/coverage/`
- `results/qc/multiqc/multiqc_report.html`
- `results/alignment/`
- `results/somatic/deepsomatic/`
- `results/somatic/mutect2/`
- `results/somatic/comparison/`
- `results/germline/deepvariant/`
- `results/cnv/cnvkit/`
- `results/cnv/purity_ploidy/`
- `results/annotation/`
- `results/callable/`
- `results/signatures/`
- `results/clonality/`
- `results/cohort/sample_qc.tsv`
- `results/cohort/somatic_mutations.maf`
- `results/cohort/tmb.tsv`
- `results/cohort/drivers_long.tsv`
- `results/cohort/driver_matrix.tsv`
- `results/cohort/purity_ploidy.tsv`
- `results/cohort/allele_specific_cnv.tsv`
- `results/cohort/arm_level_cnv_long.tsv`
- `results/cohort/arm_level_cnv_matrix.tsv`
- `results/cohort/signature_exposures.tsv`
- `results/cohort/clonality_summary.tsv`
- `results/cohort/genomic_features.tsv`
- `results/provenance/software_versions.tsv`
- `results/provenance/containers.tsv`
- `results/provenance/reference_manifest.tsv`
- `results/provenance/pipeline_parameters.yaml`
- `results/provenance/execution_report.html`

## Pre-run Checklist For Real Samples

Before running a patient cohort, verify all of the following:

- sample-sheet paths are absolute or resolve correctly relative to the sample sheet
- all FASTQs are gzipped and paired correctly
- `bait_bed` matches the assay used for every sample
- the reference FASTA, `.fai`, and `.dict` are all from the same GRCh38 build
- `exome_targets.interval_list` uses the same contig naming convention as the FASTA
- GATK germline-resource and panel-of-normals files match the same build
- the VEP cache matches GRCh38 and release 116
- all required executables are callable in the same environment used by Nextflow
- GPU-capable tools are available if you intend to use accelerated DeepSomatic or DeepVariant

## Recommended Validation Before Production Use

Run these checks in order:

1. `python3 -m unittest discover -s tests -v`
2. `python3 scripts/validate_samplesheet.py --samplesheet /path/to/samplesheet.csv --check-files`
3. `python3 scripts/validate_reference.py --config config/references.yaml`
4. `NXF_OFFLINE=true NXF_DISABLE_CHECK_LATEST=true ./bin/nextflow run main.nf --samplesheet tests/test_data/samplesheet_local.csv --reference_config tests/test_data/references_test.yaml -stub-run -profile local`
5. a one-patient real run with downsampled or non-critical data

Do not start with a large cohort run.

## Real-Run Caveats

These are the main technical caveats for actual analysis:

- The workflow is structurally end to end, but some analysis modules are simplified placeholders for scientific methodology.
- The reporting layer is static HTML generated from workflow tables rather than a live interactive dashboard.
- Outside the `pitt_htc` profile, containers are still not injected automatically into each process.
- The identity-check stage from the original specification is not yet implemented.
- BAM-input mode and tumor-only mode are not implemented.
- The current callable/TMB/signature/clonality stages are suitable for workflow development and integration testing, but should be reviewed before publication-grade analysis.

## Troubleshooting

### Sample-sheet validation fails

- Check required columns.
- Check for duplicate patient IDs or sample IDs.
- Check that FASTQ and BED paths exist from the sample-sheet location.

### Reference validation fails

- Check every path in `config/references.yaml`.
- Check that relative paths are correct from the config location.
- Check contig naming consistency across FASTA, interval list, BED, and VCF resources.

### Nextflow starts but tools are not found

- The current workflow expects tools in `PATH`.
- Load your modules or activate the correct environment before running `./bin/nextflow`.

### GPU stages do not actually use GPUs

- Ensure your `run_deepsomatic` and `run_deepvariant` installations are GPU-capable.
- Ensure your wrappers or environment expose CUDA-compatible runtimes.
- If using containers, ensure your wrappers call `singularity exec --nv`.

## Repository Files You Will Edit Most Often

- [config/default.yaml](/home/arc85/Desktop/wes_workflow/config/default.yaml)
- [config/references.yaml](/home/arc85/Desktop/wes_workflow/config/references.yaml)
- [config/containers.yaml](/home/arc85/Desktop/wes_workflow/config/containers.yaml)
- [nextflow.config](/home/arc85/Desktop/wes_workflow/nextflow.config)
- [main.nf](/home/arc85/Desktop/wes_workflow/main.nf)

## Current Validation Status

As of Wednesday, August 19, 2026:

- unit tests pass locally
- the full workflow passes `-stub-run` end to end
- real command paths are implemented for the major stages listed above
- full production validation on a representative real tumor/normal pair is still pending
