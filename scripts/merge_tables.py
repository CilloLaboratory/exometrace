#!/usr/bin/env python3
"""Merge TSV files with identical headers."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def merge_identical_headers(inputs: list[Path], output: Path) -> int:
    header = None
    rows: list[str] = []
    for path in inputs:
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
        raise ValueError("no input rows")

    with output.open("w", encoding="utf-8") as handle:
        handle.write(header + "\n")
        for row in rows:
            handle.write(row + "\n")
    return 0


def merge_union_by_first_column(inputs: list[Path], output: Path, fill_value: str) -> int:
    merged_rows: list[dict[str, str]] = []
    union_columns: list[str] = []
    union_seen: set[str] = set()
    key_column: str | None = None

    for path in inputs:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                continue
            if not reader.fieldnames:
                continue
            current_key = reader.fieldnames[0]
            if key_column is None:
                key_column = current_key
            elif current_key != key_column:
                raise ValueError(f"inconsistent key column in {path}: expected {key_column}, found {current_key}")

            for column in reader.fieldnames[1:]:
                if column not in union_seen:
                    union_seen.add(column)
                    union_columns.append(column)

            for row in reader:
                if not row:
                    continue
                row_key = row.get(current_key, "")
                if not row_key:
                    continue
                merged_rows.append({current_key: row_key, **{column: row.get(column, fill_value) or fill_value for column in reader.fieldnames[1:]}})

    if key_column is None:
        raise ValueError("no input rows")

    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([key_column] + union_columns)
        for row in merged_rows:
            writer.writerow([row[key_column]] + [row.get(column, fill_value) for column in union_columns])
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", nargs="+", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--union-by-first-column", action="store_true")
    parser.add_argument("--fill-value", default="0")
    args = parser.parse_args()
    if args.union_by_first_column:
        return merge_union_by_first_column(args.inputs, args.output, args.fill_value)
    return merge_identical_headers(args.inputs, args.output)


if __name__ == "__main__":
    raise SystemExit(main())
