#!/bin/bash
#SBATCH --job-name=wes-nextflow
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=wes-nextflow_%j.out
#SBATCH --error=wes-nextflow_%j.err

set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: sbatch scripts/pitt_htc_nextflow.sbatch /path/to/samplesheet.csv /path/to/references.yaml [extra nextflow args...]"
  exit 1
fi

SAMPLESHEET="$1"
REFERENCE_CONFIG="$2"
shift 2

module load singularity/3.9.6 || true

export NXF_OFFLINE=true
export NXF_DISABLE_CHECK_LATEST=true
export NXF_HOME="${SLURM_SUBMIT_DIR}/.nextflow"
export NXF_SINGULARITY_CACHEDIR="${SLURM_SUBMIT_DIR}/containers/sif"

cd "${SLURM_SUBMIT_DIR}"

./bin/nextflow run main.nf \
  --samplesheet "${SAMPLESHEET}" \
  --reference_config "${REFERENCE_CONFIG}" \
  -profile pitt_htc \
  "$@"
