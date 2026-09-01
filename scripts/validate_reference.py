#!/usr/bin/env python3
"""Validate required reference bundle files and emit a manifest."""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import yaml


BWA_MEM2_INDEX_SUFFIXES = [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]
BWA_CLASSIC_INDEX_SUFFIXES = [".amb", ".ann", ".bwt", ".pac", ".sa"]
SIGPROFILER_MARKERS = [
    "tsb/GRCh38",
    "references/chromosomes/transcripts/GRCh38",
    "references/chromosomes/tsb/GRCh38",
]
MIN_LINES = {
    "intervals.targets_bed": 10,
    "intervals.interval_list": 10,
    "intervals.callable_regions": 10,
    "annotations.cancer_gene_census": 10,
    "annotations.hotspots": 10,
    "annotations.gene_coordinates": 10,
    "annotations.chromosome_arms": 10,
}


def flatten_paths(config: dict) -> list[tuple[str, Path]]:
    items: list[tuple[str, Path]] = []

    def walk(prefix: str, node: object) -> None:
        if isinstance(node, dict):
            for key, value in node.items():
                walk(f"{prefix}.{key}" if prefix else key, value)
        elif isinstance(node, str) and looks_like_path(node):
            items.append((prefix, Path(node)))

    walk("", config)
    return items


def looks_like_path(value: str) -> bool:
    return "/" in value or value.endswith((".fa", ".fai", ".dict", ".bed", ".gz", ".tsv", ".cnn", ".vcf", ".tbi"))


def dir_has_entries(path: Path) -> bool:
    return any(path.iterdir())


def find_sidecar_index(path: Path) -> Path | None:
    for suffix in (".tbi", ".idx", ".csi"):
        candidate = Path(f"{path}{suffix}")
        if candidate.exists():
            return candidate
    return None


def count_data_lines(path: Path) -> int:
    if path.suffix == ".gz":
        import gzip

        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return sum(1 for line in handle if line.strip() and not line.startswith("#"))
    with path.open(encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def extra_requirements(label: str, path: Path, analysis_mode: str) -> list[tuple[str, Path, str]]:
    checks: list[tuple[str, Path, str]] = []
    if label == "reference.bwa_index_prefix":
        for suffix in BWA_MEM2_INDEX_SUFFIXES:
            checks.append((f"{label}{suffix}", Path(f"{path}{suffix}"), "file"))
        if analysis_mode == "ctdna_umi":
            for suffix in BWA_CLASSIC_INDEX_SUFFIXES:
                checks.append((f"{label}{suffix}", Path(f"{path}{suffix}"), "file"))
    if label.startswith("gatk.") and path.suffixes[-2:] == [".vcf", ".gz"]:
        checks.append((f"{label}.index", Path(f"{path}.tbi"), "file_or_alt_index"))
    if label == "vep.cache_dir":
        checks.append((label, path, "nonempty_dir"))
    if label == "sigprofiler.volume_dir":
        checks.append((label, path, "sigprofiler_dir"))
    if label in MIN_LINES:
        checks.append((label, path, f"min_lines:{MIN_LINES[label]}"))
    return checks


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--pipeline-config", type=Path)
    parser.add_argument("--emit-manifest", type=Path)
    return parser


def resolve_config_path(config_path: Path, raw_path: Path) -> Path:
    if raw_path.is_absolute():
        return raw_path
    base_dir = config_path.parent.parent if config_path.parent.name == "config" else config_path.parent
    return (base_dir / raw_path).resolve()


def get_nested(config: dict, key: str):
    value = config
    for part in key.split("."):
        value = value[part]
    return value


def main() -> int:
    args = build_parser().parse_args()
    with args.config.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    pipeline_config = {}
    if args.pipeline_config:
        with args.pipeline_config.open(encoding="utf-8") as handle:
            pipeline_config = yaml.safe_load(handle)

    analysis_mode = "tumor_normal_wes"
    if pipeline_config:
        analysis_mode = get_nested(pipeline_config, "analysis.mode")

    missing: list[str] = []
    manifest_rows: list[tuple[str, Path, str, int | str]] = []

    for label, raw_path in flatten_paths(config):
        path = resolve_config_path(args.config, raw_path)
        if not path.exists():
            missing.append(f"{label}: {raw_path}")
            continue
        checksum = sha256sum(path) if path.is_file() else "DIRECTORY"
        size = path.stat().st_size if path.is_file() else 0
        manifest_rows.append((label, path, checksum, size))
        for extra_label, extra_path, mode in extra_requirements(label, path, analysis_mode):
            if mode == "file":
                if not extra_path.exists():
                    missing.append(f"{extra_label}: {extra_path}")
                    continue
                manifest_rows.append((extra_label, extra_path, sha256sum(extra_path), extra_path.stat().st_size))
            elif mode == "file_or_alt_index":
                sidecar = find_sidecar_index(path)
                if sidecar is None:
                    missing.append(f"{extra_label}: expected one of {path}.tbi/.idx/.csi")
                    continue
                manifest_rows.append((extra_label, sidecar, sha256sum(sidecar), sidecar.stat().st_size))
            elif mode == "nonempty_dir":
                if not path.is_dir() or not dir_has_entries(path):
                    missing.append(f"{extra_label}: expected populated directory at {raw_path}")
                    continue
            elif mode == "sigprofiler_dir":
                if not path.is_dir():
                    missing.append(f"{extra_label}: expected directory at {raw_path}")
                    continue
                if not any((path / marker).exists() for marker in SIGPROFILER_MARKERS):
                    missing.append(f"{extra_label}: expected SigProfiler assets under {raw_path}")
                    continue
            elif mode.startswith("min_lines:"):
                minimum = int(mode.split(":", 1)[1])
                if count_data_lines(path) < minimum:
                    missing.append(f"{extra_label}: expected at least {minimum} non-header lines in {raw_path}")
                    continue

    if analysis_mode == "ctdna_umi":
        for label in ("cfsnv.blocked_positions_vcf", "cfsnv.blocked_positions_index"):
            raw_path = get_nested(config, label)
            path = resolve_config_path(args.config, Path(raw_path))
            if not path.exists():
                missing.append(f"{label}: {raw_path}")
                continue
            manifest_rows.append((label, path, sha256sum(path), path.stat().st_size))

    if missing:
        for entry in missing:
            print(f"ERROR: missing reference path {entry}", file=sys.stderr)
        return 1

    if args.emit_manifest:
        with args.emit_manifest.open("w", encoding="utf-8") as handle:
            handle.write("label\tpath\tsha256\tsize_bytes\n")
            for label, path, checksum, size in manifest_rows:
                handle.write(f"{label}\t{path}\t{checksum}\t{size}\n")

    print(f"Validated {len(manifest_rows)} reference bundle path(s).", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
