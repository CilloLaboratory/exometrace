#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SOURCE_DIR=${1:-/home/arc85/Desktop/wes_singularity_cache}
DEST_DIR=${2:-"${PROJECT_DIR}/containers/sandbox_local"}
TMP_DIR=${TMPDIR:-/tmp}

declare -A SANDBOXES=(
  ["bwa_mem2_samtools_2.3.sif"]="bwa_mem2_samtools_2.3"
  ["qc_suite_2026.08.sif"]="qc_suite_2026.08"
  ["bwa_mem2_2.3.sif"]="bwa_mem2_2.3"
  ["cfsnv_1.0.sif"]="cfsnv_1.0"
  ["gatk_4.7.0.0.sif"]="gatk_4.7.0.0"
  ["umi_consensus_1.5.1.sif"]="umi_consensus_1.5.1"
  ["deepvariant_1.10.0.sif"]="deepvariant_1.10.0"
  ["deepsomatic_1.10.0.sif"]="deepsomatic_1.10.0"
  ["cnvkit_0.9.14.sif"]="cnvkit_0.9.14"
  ["vep_116.0.sif"]="vep_116.0"
  ["sigprofiler_1.1.5.sif"]="sigprofiler_1.1.5"
  ["facets_suite_0.6.2.sif"]="facets_suite_0.6.2"
)

log() {
  printf '[prepare_local_singularity_sandboxes] %s\n' "$*"
}

require_tool() {
  local tool=$1
  if ! command -v "${tool}" >/dev/null 2>&1; then
    printf 'ERROR: missing required tool: %s\n' "${tool}" >&2
    exit 1
  fi
}

extract_rootfs_id() {
  local sif=$1
  singularity sif list "${sif}" | awk '/\|FS[[:space:]]*\(/ { gsub(/[[:space:]]+/, "", $1); print $1; exit }'
}

prepare_one() {
  local source_name=$1
  local sandbox_name=$2
  local source_path="${SOURCE_DIR}/${source_name}"
  local sandbox_path="${DEST_DIR}/${sandbox_name}"
  local stamp_path="${sandbox_path}/.source_sif"
  local fs_id
  local tmp_squash

  if [[ ! -e "${source_path}" ]]; then
    printf 'ERROR: missing source container: %s\n' "${source_path}" >&2
    exit 1
  fi

  if [[ -f "${stamp_path}" ]] && [[ "$(cat "${stamp_path}")" == "${source_path}" ]] && [[ "${sandbox_path}" -nt "${source_path}" ]]; then
    log "Reusing ${sandbox_name}"
    return
  fi

  rm -rf "${sandbox_path}"
  mkdir -p "${sandbox_path}"

  if [[ -d "${source_path}" ]]; then
    log "Copying sandbox ${source_name} -> ${sandbox_name}"
    cp -a "${source_path}/." "${sandbox_path}/"
    printf '%s\n' "${source_path}" > "${stamp_path}"
    return
  fi

  fs_id=$(extract_rootfs_id "${source_path}")
  if [[ -z "${fs_id}" ]]; then
    printf 'ERROR: unable to find filesystem object in %s\n' "${source_path}" >&2
    exit 1
  fi

  tmp_squash="${TMP_DIR}/${sandbox_name}.squashfs"
  rm -f "${tmp_squash}"

  log "Extracting ${source_name} -> ${sandbox_name}"
  singularity sif dump "${fs_id}" "${source_path}" > "${tmp_squash}"
  unsquashfs -f -d "${sandbox_path}" "${tmp_squash}" >/dev/null
  rm -f "${tmp_squash}"
  printf '%s\n' "${source_path}" > "${stamp_path}"
}

main() {
  require_tool singularity
  require_tool unsquashfs
  mkdir -p "${DEST_DIR}"

  for sif_name in "${!SANDBOXES[@]}"; do
    prepare_one "${sif_name}" "${SANDBOXES[${sif_name}]}"
  done

  log "Prepared sandbox containers under ${DEST_DIR}"
  log "Run Nextflow with: -profile arc85_workstation"
}

main "$@"
