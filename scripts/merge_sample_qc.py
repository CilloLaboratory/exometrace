#!/usr/bin/env python3
"""Merge per-sample QC TSVs into a cohort QC table."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", nargs="+", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    header = None
    rows: list[str] = []
    for path in args.inputs:
        with path.open(encoding="utf-8") as handle:
            lines = [line.rstrip("\n") for line in handle]
        if not lines:
            continue
        if header is None:
            header = lines[0]
        elif lines[0] != header:
            raise ValueError(f"inconsistent header in {path}")
        rows.extend(lines[1:])

    if header is None:
        raise ValueError("no QC rows were provided")

    with args.output.open("w", encoding="utf-8") as handle:
        handle.write(header + "\n")
        for row in rows:
            handle.write(row + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
