# Troubleshooting

## Sample sheet validation fails

- Ensure all required columns are present.
- Ensure tumor and normal sample IDs are distinct.
- Ensure FASTQ paths end with `.fastq.gz`.

## Reference validation fails

- Populate the exact paths from `config/references.yaml`.
- Confirm indices and dictionaries are present before running the workflow.

## Nextflow run exits immediately

- Provide `--samplesheet`.
- Install Nextflow locally or on the target HPC environment.

## Singularity SIF execution fails locally

- If Singularity cannot execute `.sif` images directly on this workstation, prepare sandbox containers with [scripts/prepare_local_singularity_sandboxes.sh](/home/arc85/Desktop/wes_workflow/scripts/prepare_local_singularity_sandboxes.sh).
- Run the workflow with `-profile arc85_workstation` so Nextflow targets `containers/sandbox_local/` instead of raw SIF files.
- The workstation profile uses `singularity.runOptions = '--userns --cleanenv'`, which is intended for unprivileged local execution of those sandbox directories.
