#!/usr/bin/env python3
"""Summarize alignment QC metrics into one TSV row."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def parse_flagstat(path: Path) -> tuple[int, int, float, float]:
    total = 0
    mapped = 0
    paired = 0
    proper = 0

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if " in total " in line:
                total = int(line.split()[0])
            elif " mapped " in line and "primary" not in line:
                mapped = int(line.split()[0])
            elif " paired in sequencing" in line:
                paired = int(line.split()[0])
            elif " properly paired " in line:
                proper = int(line.split()[0])

    mapping_rate = (mapped / total) if total else 0.0
    proper_pair_rate = (proper / paired) if paired else 0.0
    return total, mapped, mapping_rate, proper_pair_rate


def parse_duplication_metrics(path: Path) -> float:
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("LIBRARY"):
                values = next(handle).strip().split("\t")
                if len(values) >= 9:
                    return float(values[8])
    return 0.0


def parse_insert_metrics(path: Path) -> float:
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("MEDIAN_INSERT_SIZE"):
                values = next(handle).strip().split("\t")
                if values and values[0] not in {"", "NA"}:
                    return float(values[0])
    return 0.0


def parse_hs_metrics(path: Path) -> dict[str, float]:
    metrics = {
        "mean_target_coverage": 0.0,
        "median_target_coverage": 0.0,
        "pct_target_10x": 0.0,
        "pct_target_20x": 0.0,
        "pct_target_30x": 0.0,
        "pct_target_50x": 0.0,
        "pct_target_100x": 0.0,
    }
    with path.open(encoding="utf-8") as handle:
        header = None
        for line in handle:
            if line.startswith("BAIT_SET"):
                header = line.strip().split("\t")
                values = next(handle).strip().split("\t")
                data = dict(zip(header, values))
                metrics["mean_target_coverage"] = float(data.get("MEAN_TARGET_COVERAGE", 0.0))
                metrics["median_target_coverage"] = float(data.get("MEDIAN_TARGET_COVERAGE", 0.0))
                metrics["pct_target_10x"] = float(data.get("PCT_TARGET_BASES_10X", 0.0))
                metrics["pct_target_20x"] = float(data.get("PCT_TARGET_BASES_20X", 0.0))
                metrics["pct_target_30x"] = float(data.get("PCT_TARGET_BASES_30X", 0.0))
                metrics["pct_target_50x"] = float(data.get("PCT_TARGET_BASES_50X", 0.0))
                metrics["pct_target_100x"] = float(data.get("PCT_TARGET_BASES_100X", 0.0))
                break
    return metrics


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--sample-type", required=True)
    parser.add_argument("--flagstat", required=True, type=Path)
    parser.add_argument("--duplication-metrics", required=True, type=Path)
    parser.add_argument("--insert-metrics", required=True, type=Path)
    parser.add_argument("--hs-metrics", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    total, mapped, mapping_rate, proper_pair_rate = parse_flagstat(args.flagstat)
    duplicate_rate = parse_duplication_metrics(args.duplication_metrics)
    mean_insert_size = parse_insert_metrics(args.insert_metrics)
    hs_metrics = parse_hs_metrics(args.hs_metrics)

    columns = [
        "patient_id",
        "sample",
        "sample_type",
        "total_reads",
        "mapped_reads",
        "mapping_rate",
        "proper_pair_rate",
        "duplicate_rate",
        "mean_target_coverage",
        "median_target_coverage",
        "pct_target_10x",
        "pct_target_20x",
        "pct_target_30x",
        "pct_target_50x",
        "pct_target_100x",
        "mean_insert_size",
    ]

    values = [
        args.patient_id,
        args.sample,
        args.sample_type,
        str(total),
        str(mapped),
        f"{mapping_rate:.6f}",
        f"{proper_pair_rate:.6f}",
        f"{duplicate_rate:.6f}",
        f"{hs_metrics['mean_target_coverage']:.6f}",
        f"{hs_metrics['median_target_coverage']:.6f}",
        f"{hs_metrics['pct_target_10x']:.6f}",
        f"{hs_metrics['pct_target_20x']:.6f}",
        f"{hs_metrics['pct_target_30x']:.6f}",
        f"{hs_metrics['pct_target_50x']:.6f}",
        f"{hs_metrics['pct_target_100x']:.6f}",
        f"{mean_insert_size:.6f}",
    ]

    with args.output.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        handle.write("\t".join(values) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
