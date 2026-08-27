#!/usr/bin/env bash
set -euo pipefail

output_path="${1:-software_versions.tsv}"

{
  printf 'tool\tversion\n'
  printf 'python3\t%s\n' "$(python3 --version 2>&1 | awk '{print $2}')"
  printf 'nextflow\t%s\n' "$(nextflow -version 2>/dev/null | awk '/version/ {print $NF; exit}' || printf 'not_installed')"
  printf 'singularity\t%s\n' "$(singularity --version 2>/dev/null || printf 'not_installed')"
  printf 'apptainer\t%s\n' "$(apptainer --version 2>/dev/null || printf 'not_installed')"
} > "${output_path}"
