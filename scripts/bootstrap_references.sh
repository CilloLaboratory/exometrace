#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
SOURCES_CONFIG="${PROJECT_DIR}/config/reference_sources.yaml"
TARGET_CONFIG="${PROJECT_DIR}/config/references.yaml"
TARGETS_BED=${TARGETS_BED:-}
INTERVAL_LIST=${INTERVAL_LIST:-}
CALLABLE_REGIONS_BED=${CALLABLE_REGIONS_BED:-}
CANCER_GENE_CENSUS=${CANCER_GENE_CENSUS:-}
HOTSPOTS_TSV=${HOTSPOTS_TSV:-}
GENE_COORDINATES_TSV=${GENE_COORDINATES_TSV:-}
CHROMOSOME_ARMS_BED=${CHROMOSOME_ARMS_BED:-}

log() {
  printf '[bootstrap_references] %s\n' "$*"
}

fail() {
  printf '[bootstrap_references] ERROR: %s\n' "$*" >&2
  exit 1
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || fail "Required command not found: $1"
}

copy_required_input() {
  local src=$1
  local dest=$2
  local label=$3
  [[ -n "${src}" ]] || fail "Missing required input ${label}. Export ${label}=/absolute/path/to/file before running bootstrap."
  [[ -f "${src}" ]] || fail "Input file for ${label} does not exist: ${src}"
  cp -f "${src}" "${dest}"
}

yaml_get() {
  local key=$1
  python3 - "$SOURCES_CONFIG" "$key" <<'PY'
from pathlib import Path
import sys
import yaml

config_path = Path(sys.argv[1])
key = sys.argv[2].split(".")
data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
node = data
for part in key:
    node = node[part]
print(node)
PY
}

download_file() {
  local url=$1
  local dest=$2
  local tmp
  tmp=$(mktemp "${dest}.tmp.XXXXXX")
  log "Downloading ${url} -> ${dest}"
  curl -fL --retry 5 --retry-delay 2 -o "${tmp}" "${url}"
  mv "${tmp}" "${dest}"
}

ensure_dir() {
  mkdir -p "$1"
}

ensure_tabix_index() {
  local vcf_gz=$1
  if [[ -f "${vcf_gz}.tbi" ]]; then
    return
  fi
  log "Indexing ${vcf_gz}"
  if command -v bcftools >/dev/null 2>&1; then
    bcftools index -t "${vcf_gz}"
    return
  fi
  if command -v docker >/dev/null 2>&1; then
    docker run --rm \
      -v "${PROJECT_DIR}:${PROJECT_DIR}" \
      -w "${PROJECT_DIR}" \
      broadinstitute/gatk:4.7.0.0 \
      bcftools index -t "${vcf_gz}"
    return
  fi
  fail "Unable to index ${vcf_gz}. Install bcftools or provide Docker."
}

count_vcf_records() {
  local vcf_gz=$1
  python3 - "$vcf_gz" <<'PY'
from pathlib import Path
import gzip
import sys

path = Path(sys.argv[1])
count = 0
with gzip.open(path, "rt", encoding="utf-8") as handle:
    for line in handle:
        if line.strip() and not line.startswith("#"):
            count += 1
print(count)
PY
}

ensure_bwa_index() {
  local fasta=$1
  local prefix=$2
  local required=("${prefix}.0123" "${prefix}.amb" "${prefix}.ann" "${prefix}.bwt.2bit.64" "${prefix}.pac")
  local missing=0
  for path in "${required[@]}"; do
    [[ -f "${path}" ]] || missing=1
  done
  if [[ ${missing} -eq 0 ]]; then
    return
  fi
  log "Building BWA-MEM2 index for ${fasta}"
  if command -v bwa-mem2 >/dev/null 2>&1; then
    bwa-mem2 index "${fasta}"
    return
  fi
  if command -v docker >/dev/null 2>&1; then
    docker run --rm \
      -v "${PROJECT_DIR}:${PROJECT_DIR}" \
      -w "${PROJECT_DIR}" \
      quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0 \
      bwa-mem2 index "${fasta}"
    return
  fi
  fail "Unable to build BWA-MEM2 index. Install bwa-mem2 or provide Docker."
}

ensure_common_biallelic_snps() {
  local source_vcf=$1
  local output_vcf=$2
  local min_af=$3
  if [[ -f "${output_vcf}" && -f "${output_vcf}.tbi" ]]; then
    if [[ $(count_vcf_records "${output_vcf}") -gt 0 ]]; then
      return
    fi
    log "Regenerating empty common SNP resource ${output_vcf}"
    rm -f "${output_vcf}" "${output_vcf}.tbi" "${output_vcf}.idx" "${output_vcf}.csi"
  fi
  log "Deriving common biallelic SNP set from ${source_vcf}"
  if command -v bcftools >/dev/null 2>&1; then
    bcftools view \
      -m2 \
      -M2 \
      -v snps \
      -i "INFO/AF>=${min_af}" \
      -Oz \
      -o "${output_vcf}" \
      "${source_vcf}"
  elif command -v docker >/dev/null 2>&1; then
    docker run --rm \
      -v "${PROJECT_DIR}:${PROJECT_DIR}" \
      -w "${PROJECT_DIR}" \
      broadinstitute/gatk:4.7.0.0 \
      bcftools view \
      -m2 \
      -M2 \
      -v snps \
      -i "INFO/AF>=${min_af}" \
      -Oz \
      -o "${output_vcf}" \
      "${source_vcf}"
  else
    fail "Unable to derive common biallelic SNPs. Install bcftools or provide Docker."
  fi
  ensure_tabix_index "${output_vcf}"
  if [[ $(count_vcf_records "${output_vcf}") -eq 0 ]]; then
    fail "Derived common SNP resource is empty: ${output_vcf}. Check the bcftools filter and source VCF."
  fi
}

ensure_cfsnv_blocked_positions() {
  local source_vcf=$1
  local targets_bed=$2
  local output_vcf=$3
  if [[ -f "${output_vcf}" && -f "${output_vcf}.tbi" ]]; then
    if [[ $(count_vcf_records "${output_vcf}") -gt 0 ]]; then
      return
    fi
    log "Regenerating empty cfSNV blocked-position resource ${output_vcf}"
    rm -f "${output_vcf}" "${output_vcf}.tbi" "${output_vcf}.idx" "${output_vcf}.csi"
  fi
  log "Deriving cfSNV blocked-position blacklist from ${source_vcf} on targets ${targets_bed}"
  if command -v bcftools >/dev/null 2>&1; then
    bcftools view \
      -R "${targets_bed}" \
      -m2 \
      -M2 \
      -v snps \
      -Oz \
      -o "${output_vcf}" \
      "${source_vcf}"
  elif command -v docker >/dev/null 2>&1; then
    docker run --rm \
      -v "${PROJECT_DIR}:${PROJECT_DIR}" \
      -w "${PROJECT_DIR}" \
      broadinstitute/gatk:4.7.0.0 \
      bcftools view \
      -R "${targets_bed}" \
      -m2 \
      -M2 \
      -v snps \
      -Oz \
      -o "${output_vcf}" \
      "${source_vcf}"
  else
    fail "Unable to derive cfSNV blocked positions. Install bcftools or provide Docker."
  fi
  ensure_tabix_index "${output_vcf}"
  if [[ $(count_vcf_records "${output_vcf}") -eq 0 ]]; then
    fail "Derived cfSNV blocked-position resource is empty: ${output_vcf}. Check the target BED and source VCF."
  fi
}

ensure_vep_cache() {
  local cache_dir=$1
  local cache_url=$2
  local archive_name=$3
  if find "${cache_dir}" -mindepth 1 -maxdepth 1 | grep -q .; then
    return
  fi
  local archive_path="${cache_dir}/${archive_name}"
  ensure_dir "${cache_dir}"
  download_file "${cache_url}" "${archive_path}"
  log "Extracting VEP cache ${archive_name}"
  tar -xzf "${archive_path}" -C "${cache_dir}"
  rm -f "${archive_path}"
}

ensure_sigprofiler_assets() {
  local genome_build=$1
  local volume_dir=$2
  ensure_dir "${volume_dir}"
  if [[ -d "${volume_dir}/tsb/${genome_build}" || -d "${volume_dir}/references/chromosomes/tsb/${genome_build}" || -d "${volume_dir}/references/chromosomes/transcripts/${genome_build}" ]]; then
    return
  fi

  if command -v python3 >/dev/null 2>&1; then
    if python3 -c "from SigProfilerMatrixGenerator.install import install" >/dev/null 2>&1; then
      log "Installing SigProfiler assets with local Python"
      python3 -c "from SigProfilerMatrixGenerator.install import install; install('${genome_build}', volume='${volume_dir}')" &&
        return
    fi
  fi

  if command -v docker >/dev/null 2>&1 && docker image inspect sigprofilertoolkit >/dev/null 2>&1; then
    log "Installing SigProfiler assets with Docker image sigprofilertoolkit"
    docker run --rm \
      -v "${PROJECT_DIR}:${PROJECT_DIR}" \
      -w "${PROJECT_DIR}" \
      --entrypoint python \
      sigprofilertoolkit \
      -c "from SigProfilerMatrixGenerator.install import install; install('${genome_build}', volume='${volume_dir}')"
    return
  fi

  fail "Unable to provision SigProfiler assets. Install SigProfilerMatrixGenerator locally or provide the sigprofilertoolkit Docker image."
}

emit_manifest() {
  python3 "${PROJECT_DIR}/scripts/validate_reference.py" \
    --config "${TARGET_CONFIG}" \
    --emit-manifest "${PROJECT_DIR}/references/GRCh38/reference_manifest.tsv"
}

main() {
  need_cmd python3
  need_cmd curl
  need_cmd tar

  local fasta_path fasta_index_path dict_path germline_path pon_path vep_cache_dir vep_cache_url vep_archive
  local ann_census ann_hotspots ann_coords ann_arms targets_bed interval_list callable_regions
  local bwa_prefix common_path common_min_af cfsnv_blocked_path sig_volume sig_genome

  fasta_path=$(yaml_get assets.fasta.path)
  fasta_index_path=$(yaml_get assets.fasta_index.path)
  dict_path=$(yaml_get assets.sequence_dictionary.path)
  germline_path=$(yaml_get assets.germline_resource.path)
  pon_path=$(yaml_get assets.panel_of_normals.path)
  vep_cache_dir=$(yaml_get assets.vep_cache.path)
  vep_cache_url=$(yaml_get assets.vep_cache.url)
  vep_archive=$(yaml_get assets.vep_cache.archive_name)
  targets_bed=$(yaml_get manual_inputs.intervals.targets_bed.path)
  interval_list=$(yaml_get manual_inputs.intervals.interval_list.path)
  callable_regions=$(yaml_get manual_inputs.intervals.callable_regions.path)
  ann_census=$(yaml_get manual_inputs.annotations.cancer_gene_census.path)
  ann_hotspots=$(yaml_get manual_inputs.annotations.hotspots.path)
  ann_coords=$(yaml_get manual_inputs.annotations.gene_coordinates.path)
  ann_arms=$(yaml_get manual_inputs.annotations.chromosome_arms.path)
  bwa_prefix=$(yaml_get derived.bwa_mem2_index.prefix)
  common_path=$(yaml_get derived.common_biallelic_snps.path)
  common_min_af=$(yaml_get derived.common_biallelic_snps.min_af)
  cfsnv_blocked_path=$(yaml_get derived.cfsnv_blocked_positions.path)
  sig_volume=$(yaml_get derived.sigprofiler.volume_dir)
  sig_genome=$(yaml_get derived.sigprofiler.genome_build)

  ensure_dir "${PROJECT_DIR}/references/GRCh38/fasta"
  ensure_dir "${PROJECT_DIR}/references/GRCh38/gatk"
  ensure_dir "${PROJECT_DIR}/references/GRCh38/cfsnv"
  ensure_dir "${PROJECT_DIR}/references/GRCh38/intervals"
  ensure_dir "${PROJECT_DIR}/references/GRCh38/annotations"
  ensure_dir "${PROJECT_DIR}/references/GRCh38/vep/cache"

  [[ -f "${PROJECT_DIR}/${fasta_path}" ]] || download_file "$(yaml_get assets.fasta.url)" "${PROJECT_DIR}/${fasta_path}"
  [[ -f "${PROJECT_DIR}/${fasta_index_path}" ]] || download_file "$(yaml_get assets.fasta_index.url)" "${PROJECT_DIR}/${fasta_index_path}"
  [[ -f "${PROJECT_DIR}/${dict_path}" ]] || download_file "$(yaml_get assets.sequence_dictionary.url)" "${PROJECT_DIR}/${dict_path}"
  [[ -f "${PROJECT_DIR}/${germline_path}" ]] || download_file "$(yaml_get assets.germline_resource.url)" "${PROJECT_DIR}/${germline_path}"
  [[ -f "${PROJECT_DIR}/${pon_path}" ]] || download_file "$(yaml_get assets.panel_of_normals.url)" "${PROJECT_DIR}/${pon_path}"

  copy_required_input "${TARGETS_BED}" "${PROJECT_DIR}/${targets_bed}" "TARGETS_BED"
  copy_required_input "${INTERVAL_LIST}" "${PROJECT_DIR}/${interval_list}" "INTERVAL_LIST"
  if [[ -n "${CALLABLE_REGIONS_BED}" ]]; then
    copy_required_input "${CALLABLE_REGIONS_BED}" "${PROJECT_DIR}/${callable_regions}" "CALLABLE_REGIONS_BED"
  else
    cp -f "${PROJECT_DIR}/${targets_bed}" "${PROJECT_DIR}/${callable_regions}"
  fi
  copy_required_input "${CANCER_GENE_CENSUS}" "${PROJECT_DIR}/${ann_census}" "CANCER_GENE_CENSUS"
  copy_required_input "${HOTSPOTS_TSV}" "${PROJECT_DIR}/${ann_hotspots}" "HOTSPOTS_TSV"
  copy_required_input "${GENE_COORDINATES_TSV}" "${PROJECT_DIR}/${ann_coords}" "GENE_COORDINATES_TSV"
  copy_required_input "${CHROMOSOME_ARMS_BED}" "${PROJECT_DIR}/${ann_arms}" "CHROMOSOME_ARMS_BED"

  ensure_tabix_index "${PROJECT_DIR}/${germline_path}"
  ensure_tabix_index "${PROJECT_DIR}/${pon_path}"
  ensure_bwa_index "${PROJECT_DIR}/${fasta_path}" "${PROJECT_DIR}/${bwa_prefix}"
  ensure_common_biallelic_snps "${PROJECT_DIR}/${germline_path}" "${PROJECT_DIR}/${common_path}" "${common_min_af}"
  ensure_cfsnv_blocked_positions "${PROJECT_DIR}/${common_path}" "${PROJECT_DIR}/${targets_bed}" "${PROJECT_DIR}/${cfsnv_blocked_path}"
  ensure_vep_cache "${PROJECT_DIR}/${vep_cache_dir}" "${vep_cache_url}" "${vep_archive}"
  ensure_sigprofiler_assets "${sig_genome}" "${PROJECT_DIR}/${sig_volume}"
  emit_manifest
  log "Reference bootstrap complete."
}

main "$@"
