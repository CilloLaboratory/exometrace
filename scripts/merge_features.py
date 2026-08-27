#!/usr/bin/env python3
"""Merge per-patient feature tables into a genomic features matrix."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_single_row_table(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows[0] if rows else {}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tmb", nargs="+", required=True, type=Path)
    parser.add_argument("--purity", nargs="+", required=True, type=Path)
    parser.add_argument("--signatures", nargs="+", required=True, type=Path)
    parser.add_argument("--clonality", nargs="+", required=True, type=Path)
    parser.add_argument("--drivers", nargs="+", required=True, type=Path)
    parser.add_argument("--arm-matrices", nargs="+", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    features: dict[str, dict[str, str]] = {}

    for path in args.tmb:
        row = load_single_row_table(path)
        features.setdefault(row["patient_id"], {}).update({"TMB": row["tmb_mut_per_mb"]})
    for path in args.purity:
        row = load_single_row_table(path)
        features.setdefault(row["patient_id"], {}).update({"purity": row["purity"], "ploidy": row["ploidy"]})
    for path in args.signatures:
        row = load_single_row_table(path)
        features.setdefault(row["patient_id"], {}).update({"UV_signature": row["SBS7_UV"], "APOBEC_signature": row["SBS2_APOBEC"]})
    for path in args.clonality:
        row = load_single_row_table(path)
        features.setdefault(row["patient_id"], {}).update(
            {
                "clonal_mutation_count": row["clonal_mutations"],
                "subclonal_mutation_count": row["subclonal_mutations"],
                "ambiguous_mutation_count": row.get("ambiguous_mutations", "0"),
            }
        )
    for path in args.drivers:
        with path.open(encoding="utf-8") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        header, *data = rows
        for row in data:
            patient_id = row[0]
            for gene, value in zip(header[1:], row[1:]):
                features.setdefault(patient_id, {})[f"driver_{gene}"] = value
    for path in args.arm_matrices:
        with path.open(encoding="utf-8") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        header, *data = rows
        for row in data:
            patient_id = row[0]
            for arm, value in zip(header[1:], row[1:]):
                features.setdefault(patient_id, {})[f"{arm}_status"] = value

    fieldnames = sorted({key for row in features.values() for key in row.keys()})
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["patient_id"] + fieldnames)
        for patient_id in sorted(features):
            writer.writerow([patient_id] + [features[patient_id].get(field, "NA") for field in fieldnames])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
