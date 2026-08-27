#!/usr/bin/env python3
"""Validate paired FASTQ files and emit a simple QC record."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def count_reads(path: Path) -> int:
    line_count = 0
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_count, _ in enumerate(handle, start=1):
            pass
    if line_count % 4 != 0:
        raise ValueError(f"{path} does not contain complete FASTQ records")
    return line_count // 4


def validate_pair(r1: Path, r2: Path) -> tuple[int, int]:
    if not r1.name.endswith(".fastq.gz") or not r2.name.endswith(".fastq.gz"):
        raise ValueError("FASTQ files must end with .fastq.gz")
    reads_r1 = count_reads(r1)
    reads_r2 = count_reads(r2)
    if reads_r1 != reads_r2:
        raise ValueError(f"read count mismatch: {reads_r1} != {reads_r2}")
    return reads_r1, reads_r2


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    reads_r1, reads_r2 = validate_pair(args.r1, args.r2)

    with args.output.open("w", encoding="utf-8") as handle:
        handle.write("sample\treads_r1\treads_r2\tstatus\n")
        handle.write(f"{args.sample_id}\t{reads_r1}\t{reads_r2}\tPASS\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
