#!/usr/bin/env python3
"""Merge per-patient arm matrices into one cohort matrix."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_matrix(path: Path) -> tuple[list[str], list[str]]:
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    return rows[0], rows[1]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", nargs="+", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    headers = None
    rows: list[list[str]] = []
    for path in args.inputs:
        header, row = read_matrix(path)
        if headers is None:
            headers = header
        elif header != headers:
            raise ValueError(f"inconsistent arm matrix header in {path}")
        rows.append(row)

    if headers is None:
        raise ValueError("no matrices provided")

    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
