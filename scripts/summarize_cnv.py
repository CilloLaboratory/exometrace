#!/usr/bin/env python3
"""Summarize CNVkit segments into arm-level calls and a matrix row."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def read_segments(path: Path) -> list[dict[str, str]]:
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            break
        except OverflowError:
            limit //= 10
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_arms(path: Path) -> list[dict[str, str]]:
    arms = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, arm = line.rstrip("\n").split("\t")[:4]
            arms.append({"chromosome": chrom, "start": int(start), "end": int(end), "arm": arm})
    return arms


def overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def classify_log2(log2_value: float, thresholds: dict[str, float]) -> str:
    if log2_value <= thresholds["deep_deletion_log2"]:
        return "deep_deletion"
    if log2_value <= thresholds["loss_log2"]:
        return "loss"
    if log2_value >= thresholds["amplification_log2"]:
        return "amplification"
    if log2_value >= thresholds["gain_log2"]:
        return "gain"
    return "neutral"


def summarize(patient_id: str, segments: list[dict[str, str]], arms: list[dict[str, str]], thresholds: dict[str, float]) -> list[dict[str, str]]:
    summaries: list[dict[str, str]] = []
    for arm in arms:
        numerator = 0.0
        denominator = 0
        for segment in segments:
            if segment["chromosome"] != arm["chromosome"]:
                continue
            seg_start = int(segment["start"])
            seg_end = int(segment["end"])
            width = overlap(seg_start, seg_end, arm["start"], arm["end"])
            if width <= 0:
                continue
            numerator += float(segment["log2"]) * width
            denominator += width
        mean_log2 = numerator / denominator if denominator else 0.0
        summaries.append(
            {
                "patient_id": patient_id,
                "arm": arm["arm"],
                "status": classify_log2(mean_log2, thresholds),
                "log2_ratio": f"{mean_log2:.6f}",
            }
        )
    return summaries


def write_long(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["patient_id", "arm", "status", "log2_ratio"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_matrix(path: Path, rows: list[dict[str, str]]) -> None:
    arms = [row["arm"] for row in rows]
    statuses = [row["status"] for row in rows]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["patient_id"] + arms)
        writer.writerow([rows[0]["patient_id"]] + statuses)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--segments", required=True, type=Path)
    parser.add_argument("--arms-bed", required=True, type=Path)
    parser.add_argument("--output-long", required=True, type=Path)
    parser.add_argument("--output-matrix", required=True, type=Path)
    parser.add_argument("--deep-deletion-log2", type=float, default=-1.1)
    parser.add_argument("--loss-log2", type=float, default=-0.3)
    parser.add_argument("--gain-log2", type=float, default=0.2)
    parser.add_argument("--amplification-log2", type=float, default=0.7)
    args = parser.parse_args()

    thresholds = {
        "deep_deletion_log2": args.deep_deletion_log2,
        "loss_log2": args.loss_log2,
        "gain_log2": args.gain_log2,
        "amplification_log2": args.amplification_log2,
    }
    rows = summarize(args.patient_id, read_segments(args.segments), read_arms(args.arms_bed), thresholds)
    write_long(args.output_long, rows)
    write_matrix(args.output_matrix, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
