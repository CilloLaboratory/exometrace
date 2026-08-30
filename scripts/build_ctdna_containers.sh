#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CACHE_DIR=${1:-/home/arc85/Desktop/wes_singularity_cache}
SIF_DIR="${PROJECT_DIR}/containers/sif"
TMP_DIR=${TMPDIR:-/tmp}

log() {
  printf '[build_ctdna_containers] %s\n' "$*"
}

require_tool() {
  local tool=$1
  command -v "${tool}" >/dev/null 2>&1 || {
    printf 'ERROR: missing required tool: %s\n' "${tool}" >&2
    exit 1
  }
}

build_one() {
  local name=$1
  local version=$2
  local dockerfile=$3
  local local_tag="wes/${name}:${version}"
  local sif_name="${name}_${version}.sif"
  local sif_path="${CACHE_DIR}/${sif_name}"
  local link_path="${SIF_DIR}/${sif_name}"
  local archive_path="${TMP_DIR}/${name}_${version}.docker.tar"

  log "Building OCI image ${local_tag}"
  docker build -f "${PROJECT_DIR}/containers/definitions/${dockerfile}" -t "${local_tag}" "${PROJECT_DIR}"

  mkdir -p "${CACHE_DIR}"
  rm -f "${archive_path}"
  log "Saving ${local_tag} as docker archive"
  docker save -o "${archive_path}" "${local_tag}"
  log "Packing ${sif_name}"
  singularity build --force "${sif_path}" "docker-archive://${archive_path}"
  rm -f "${archive_path}"

  mkdir -p "${SIF_DIR}"
  ln -sfn "${sif_path}" "${link_path}"
  log "Linked ${link_path} -> ${sif_path}"
}

main() {
  require_tool docker
  require_tool singularity

  build_one cfsnv 1.0 cfsnv.Dockerfile
  build_one umi_consensus 2.5.1 umi_consensus.Dockerfile

  log "Built ctDNA container set into ${CACHE_DIR}"
}

main "$@"
