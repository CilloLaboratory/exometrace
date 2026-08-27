#!/usr/bin/env python3
"""Calculate tumor mutational burden from MAF-like and callable-mb inputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def calculate_tmb(qualifying_mutations: int, callable_mb: float) -> float:
    if callable_mb <= 0:
        if qualifying_mutations == 0:
            return 0.0
        raise ValueError("callable_mb must be greater than zero when qualifying mutations are present")
    return qualifying_mutations / callable_mb


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=False)
    parser.add_argument("--maf", type=Path, required=False)
    parser.add_argument("--callable-mb", type=Path, required=False)
    parser.add_argument("--patient-id", required=False)
    parser.add_argument("--tumor-sample", required=False)
    parser.add_argument("--output", type=Path, required=False)
    args = parser.parse_args()
    if args.input:
        with args.input.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                tmb = calculate_tmb(int(row["qualifying_mutations"]), float(row["callable_mb"]))
                print(f'{row["patient_id"]}\t{tmb:.6f}')
        return 0

    if args.maf and args.callable_mb and args.output and args.patient_id and args.tumor_sample:
        with args.maf.open(newline="", encoding="utf-8") as handle:
            maf_rows = list(csv.DictReader(handle, delimiter="\t"))
        qualifying = len(maf_rows)
        with args.callable_mb.open(newline="", encoding="utf-8") as handle:
            callable_rows = list(csv.DictReader(handle, delimiter="\t"))
        callable_mb = float(callable_rows[0]["callable_mb"])
        tmb = calculate_tmb(qualifying, callable_mb)
        with args.output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["patient_id", "tumor_sample", "qualifying_mutations", "callable_mb", "tmb_mut_per_mb"])
            writer.writerow([args.patient_id, args.tumor_sample, qualifying, f"{callable_mb:.6f}", f"{tmb:.6f}"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
