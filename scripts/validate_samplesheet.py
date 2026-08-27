#!/usr/bin/env python3
"""Validate tumor/normal WES sample sheets."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

REQUIRED_FIELDS = [
    "patient_id",
    "tumor_sample",
    "normal_sample",
    "tumor_r1",
    "tumor_r2",
    "normal_r1",
    "normal_r2",
    "bait_bed",
]

OPTIONAL_FIELDS = [
    "sex",
    "library_id",
    "sequencing_run",
    "capture_kit",
    "batch",
    "FFPE",
    "diagnosis",
]


def read_samplesheet(path: Path) -> tuple[list[dict[str, str]], str]:
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError("sample sheet is missing a header row")
        rows = [normalize_row(row) for row in reader]
    if not rows:
        raise ValueError("sample sheet contains no data rows")
    return rows, delimiter


def normalize_row(row: dict[str, str | None]) -> dict[str, str]:
    return {str(key).strip(): (value or "").strip() for key, value in row.items()}


def validate_rows(
    rows: list[dict[str, str]],
    samplesheet_path: Path,
    check_files: bool = False,
) -> list[str]:
    errors: list[str] = []
    fieldnames = set(rows[0].keys())

    missing_fields = [field for field in REQUIRED_FIELDS if field not in fieldnames]
    if missing_fields:
        errors.append(f"missing required columns: {', '.join(missing_fields)}")
        return errors

    seen_patients: set[str] = set()
    seen_samples: set[str] = set()

    for index, row in enumerate(rows, start=2):
        for field in REQUIRED_FIELDS:
            if not row.get(field):
                errors.append(f"row {index}: empty required field '{field}'")

        patient_id = row["patient_id"]
        if patient_id in seen_patients:
            errors.append(f"row {index}: duplicate patient_id '{patient_id}'")
        seen_patients.add(patient_id)

        tumor_sample = row["tumor_sample"]
        normal_sample = row["normal_sample"]

        if tumor_sample == normal_sample:
            errors.append(f"row {index}: tumor_sample and normal_sample are identical")

        for sample_field, sample_name in (
            ("tumor_sample", tumor_sample),
            ("normal_sample", normal_sample),
        ):
            if sample_name in seen_samples:
                errors.append(f"row {index}: duplicate sample identifier '{sample_name}'")
            seen_samples.add(sample_name)
            if "," in sample_name or "\t" in sample_name:
                errors.append(f"row {index}: invalid delimiter character in {sample_field}")

        for fastq_field in ("tumor_r1", "tumor_r2", "normal_r1", "normal_r2"):
            fastq_path = row[fastq_field]
            if fastq_path and not fastq_path.endswith(".fastq.gz"):
                errors.append(f"row {index}: {fastq_field} must end with .fastq.gz")
            if check_files and fastq_path:
                resolved = resolve_path(samplesheet_path, fastq_path)
                if not resolved.exists():
                    errors.append(f"row {index}: missing FASTQ '{fastq_path}'")

        bait_bed = row["bait_bed"]
        if bait_bed and not bait_bed.endswith((".bed", ".bed.gz")):
            errors.append(f"row {index}: bait_bed must be a BED path")
        if check_files and bait_bed:
            resolved = resolve_path(samplesheet_path, bait_bed)
            if not resolved.exists():
                errors.append(f"row {index}: missing bait BED '{bait_bed}'")

    return errors


def resolve_path(samplesheet_path: Path, raw_path: str) -> Path:
    path = Path(raw_path)
    return path if path.is_absolute() else (samplesheet_path.parent / path).resolve()


def emit_normalized(rows: list[dict[str, str]], output_path: Path, samplesheet_path: Path) -> None:
    path_fields = {"tumor_r1", "tumor_r2", "normal_r1", "normal_r2", "bait_bed"}
    ordered_fields = REQUIRED_FIELDS + [field for field in OPTIONAL_FIELDS if field in rows[0]]
    extra_fields = [field for field in rows[0] if field not in ordered_fields]
    fieldnames = ordered_fields + extra_fields

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            normalized = dict(row)
            for field in path_fields:
                if normalized.get(field):
                    normalized[field] = str(resolve_path(samplesheet_path, normalized[field]))
            writer.writerow(normalized)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", required=True, type=Path)
    parser.add_argument("--check-files", action="store_true")
    parser.add_argument("--emit-normalized", type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    rows, _ = read_samplesheet(args.samplesheet)
    errors = validate_rows(rows, args.samplesheet, check_files=args.check_files)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    if args.emit_normalized:
        emit_normalized(rows, args.emit_normalized, args.samplesheet)
    print(f"Validated {len(rows)} paired tumor/normal record(s).", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
