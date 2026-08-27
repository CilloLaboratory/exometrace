#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SOURCE_DIR=${1:-/home/arc85/Desktop/wes_singularity_cache}
TARGET_DIR=${2:-/home/arc85/Desktop/wes_singularity_cache}
TARGET_NAME=${3:-bwa_mem2_samtools_2.3.sif}
TMP_DIR=${TMPDIR:-/tmp}

BASE_SIF="${SOURCE_DIR}/qc_suite_2026.08.sif"
BWA_SIF="${SOURCE_DIR}/bwa_mem2_2.3.sif"
TARGET_PATH="${TARGET_DIR}/${TARGET_NAME}"
WORK_NAME="${TARGET_NAME%.sif}"
BASE_SANDBOX="${TMP_DIR}/${WORK_NAME}.base.$$"
BWA_SANDBOX="${TMP_DIR}/${WORK_NAME}.bwa.$$"
MERGED_SANDBOX="${TMP_DIR}/${WORK_NAME}.merged.$$"

cleanup() {
  rm -rf "${BASE_SANDBOX}" "${BWA_SANDBOX}" "${MERGED_SANDBOX}"
}

log() {
  printf '[build_align_container] %s\n' "$*"
}

require_tool() {
  local tool=$1
  command -v "${tool}" >/dev/null 2>&1 || {
    printf 'ERROR: missing required tool: %s\n' "${tool}" >&2
    exit 1
  }
}

extract_rootfs_id() {
  local sif=$1
  singularity sif list "${sif}" | awk '/\|FS[[:space:]]*\(/ { gsub(/[[:space:]]+/, "", $1); print $1; exit }'
}

extract_sif() {
  local sif_path=$1
  local dest_dir=$2
  local fs_id
  local squash_path

  fs_id=$(extract_rootfs_id "${sif_path}")
  [[ -n "${fs_id}" ]] || {
    printf 'ERROR: unable to find filesystem object in %s\n' "${sif_path}" >&2
    exit 1
  }

  squash_path="${TMP_DIR}/$(basename "${dest_dir}").squashfs"
  rm -f "${squash_path}"
  mkdir -p "${dest_dir}"
  singularity sif dump "${fs_id}" "${sif_path}" > "${squash_path}"
  unsquashfs -f -d "${dest_dir}" "${squash_path}" >/dev/null
  rm -f "${squash_path}"
}

copy_if_exists() {
  local src=$1
  local dest=$2
  if [[ -e "${src}" ]]; then
    mkdir -p "$(dirname "${dest}")"
    cp -a "${src}" "${dest}"
  fi
}

main() {
  trap cleanup EXIT

  require_tool singularity
  require_tool unsquashfs

  [[ -f "${BASE_SIF}" ]] || {
    printf 'ERROR: missing base SIF: %s\n' "${BASE_SIF}" >&2
    exit 1
  }
  [[ -f "${BWA_SIF}" ]] || {
    printf 'ERROR: missing BWA SIF: %s\n' "${BWA_SIF}" >&2
    exit 1
  }

  log "Extracting ${BASE_SIF}"
  extract_sif "${BASE_SIF}" "${BASE_SANDBOX}"

  log "Extracting ${BWA_SIF}"
  extract_sif "${BWA_SIF}" "${BWA_SANDBOX}"

  log "Copying bwa-mem2 binaries into base sandbox"
  mkdir -p "${BASE_SANDBOX}/usr/local/bin"
  cp -a "${BWA_SANDBOX}/usr/local/bin/." "${BASE_SANDBOX}/usr/local/bin/"

  copy_if_exists "${BWA_SANDBOX}/usr/local/lib" "${BASE_SANDBOX}/usr/local/lib"
  copy_if_exists "${BWA_SANDBOX}/usr/local/share" "${BASE_SANDBOX}/usr/local/share"

  mkdir -p "${TARGET_DIR}"
  rm -rf "${MERGED_SANDBOX}"
  mv "${BASE_SANDBOX}" "${MERGED_SANDBOX}"
  rm -f "${TARGET_PATH}"

  log "Packing combined sandbox into ${TARGET_PATH}"
  singularity build "${TARGET_PATH}" "${MERGED_SANDBOX}" >/dev/null

  log "Built combined alignment SIF at ${TARGET_PATH}"
  log "Contains qc_suite base plus bwa-mem2 binaries"
}

main "$@"
