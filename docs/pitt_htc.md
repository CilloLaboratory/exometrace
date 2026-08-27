# Pitt HTC Adaptation

This workflow now includes a `pitt_htc` Nextflow profile for the University of Pittsburgh CRCD HTC cluster.

## What The Profile Does

- uses the Slurm executor
- targets `--cluster=htc`
- submits to the `htc` partition
- enables `singularity`
- loads `singularity/3.9.6`
- enables Nextflow scratch execution for task work directories
- keeps CPU-labeled tasks on HTC and routes GPU-labeled tasks to the Pitt GPU cluster
- executes against cached local `.sif` images from `containers/sif/` rather than downloading images on demand
- uses a dedicated `qc_suite` SIF path for FastQC, samtools, mosdepth, and Picard-style QC tasks

## Files Added Or Updated

- [nextflow.config](/home/arc85/Desktop/wes_workflow/nextflow.config)
- [scripts/pitt_htc_nextflow.sbatch](/home/arc85/Desktop/wes_workflow/scripts/pitt_htc_nextflow.sbatch)
- process-level container hooks in the workflow modules

## Required Pitt-Specific Preparation

1. Build or pull the required SIF images into `containers/sif/`.
2. Confirm the SIF filenames match those configured in [nextflow.config](/home/arc85/Desktop/wes_workflow/nextflow.config).
3. In particular, provide:
   - `containers/sif/qc_suite_2026.08.sif`
   - `containers/sif/bwa_mem2_2.3.sif`
   - `containers/sif/gatk_4.7.0.0.sif`
   - `containers/sif/deepvariant_1.10.0.sif`
   - `containers/sif/deepsomatic_1.10.0.sif`
   - `containers/sif/cnvkit_0.9.14.sif`
   - `containers/sif/vep_116.0.sif`
4. Validate that the CRCD environment exposes the needed host-side tools for any non-containerized script-only steps:
   - `python3`
   - shell utilities used by Nextflow task scripts
5. Submit the workflow with:

```bash
sbatch scripts/pitt_htc_nextflow.sbatch /path/to/samplesheet.csv /path/to/config/references.yaml
```

## Pitt Parameters

Optional parameters exposed in [nextflow.config](/home/arc85/Desktop/wes_workflow/nextflow.config):

- `--pitt_account`
- `--pitt_qos`
- `--pitt_constraint`
- `--pitt_gpu_partition`
- `--pitt_gpu_gres`

Example:

```bash
sbatch scripts/pitt_htc_nextflow.sbatch \
  /path/to/samplesheet.csv \
  /path/to/config/references.yaml \
  --pitt_account my_allocation \
  --pitt_gpu_partition l40s \
  --pitt_constraint amd,genoa
```

## Notes

- The Pitt HTC cluster is intended for single-node genomics workflows.
- CRCD documentation recommends node-local scratch and Slurm batch submission rather than running real analyses on login nodes.
- GPU-labeled processes such as DeepSomatic and DeepVariant are expected to run from the Pitt GPU cluster with `--nv` container execution.
- The profile expects the referenced `.sif` files to already exist locally under `containers/sif/`.
