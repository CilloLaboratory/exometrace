#!/usr/bin/env python3
"""Build callable BED intervals from depth-filtered tumor/normal coverage loci."""

from __future__ import annotations

import argparse
from pathlib import Path


def total_bases(path: Path) -> int:
    total = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *_rest = line.rstrip("\n").split("\t")
            total += int(end) - int(start)
    return total


def depth_to_intervals(depth_path: Path, tumor_min_depth: int, normal_min_depth: int) -> list[tuple[str, int, int]]:
    intervals: list[tuple[str, int, int]] = []
    current_chrom: str | None = None
    current_start: int | None = None
    current_end: int | None = None

    with depth_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, pos, tumor_depth, normal_depth = line.rstrip("\n").split("\t")[:4]
            locus_start = int(pos) - 1
            locus_end = int(pos)
            is_callable = int(tumor_depth) >= tumor_min_depth and int(normal_depth) >= normal_min_depth

            if not is_callable:
                if current_chrom is not None:
                    intervals.append((current_chrom, current_start, current_end))
                    current_chrom = None
                    current_start = None
                    current_end = None
                continue

            if current_chrom == chrom and current_end == locus_start:
                current_end = locus_end
                continue

            if current_chrom is not None:
                intervals.append((current_chrom, current_start, current_end))

            current_chrom = chrom
            current_start = locus_start
            current_end = locus_end

    if current_chrom is not None:
        intervals.append((current_chrom, current_start, current_end))

    return intervals


def write_bed(intervals: list[tuple[str, int, int]], output_bed: Path) -> None:
    with output_bed.open("w", encoding="utf-8") as handle:
        for chrom, start, end in intervals:
            handle.write(f"{chrom}\t{start}\t{end}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--depth-tsv", required=True, type=Path)
    parser.add_argument("--tumor-min-depth", required=True, type=int)
    parser.add_argument("--normal-min-depth", required=True, type=int)
    parser.add_argument("--output-bed", required=True, type=Path)
    parser.add_argument("--output-mb", required=True, type=Path)
    args = parser.parse_args()

    intervals = depth_to_intervals(args.depth_tsv, args.tumor_min_depth, args.normal_min_depth)
    write_bed(intervals, args.output_bed)
    mb = total_bases(args.output_bed) / 1_000_000
    with args.output_mb.open("w", encoding="utf-8") as handle:
        handle.write("patient_id\tcallable_mb\n")
        handle.write(f"{args.patient_id}\t{mb:.6f}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
