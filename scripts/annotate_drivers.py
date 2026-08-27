#!/usr/bin/env python3
"""Annotate potential drivers from a MAF-like table and curated resources."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_gene_set(path: Path) -> set[str]:
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["gene"] for row in reader}


def load_hotspots(path: Path) -> set[tuple[str, str]]:
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {(row["gene"], row["protein_change"]) for row in reader}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maf", required=True, type=Path)
    parser.add_argument("--census", required=True, type=Path)
    parser.add_argument("--hotspots", required=True, type=Path)
    parser.add_argument("--output-long", required=True, type=Path)
    parser.add_argument("--output-matrix", required=True, type=Path)
    args = parser.parse_args()

    cancer_genes = load_gene_set(args.census)
    hotspots = load_hotspots(args.hotspots)

    long_rows: list[dict[str, str]] = []
    matrix: dict[str, dict[str, str]] = {}
    all_genes: set[str] = set()

    with args.maf.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene = row["Hugo_Symbol"]
            protein_change = row["HGVSp"]
            driver_reasons = []
            if gene in cancer_genes:
                driver_reasons.append("CancerGeneCensus")
            if (gene, protein_change) in hotspots:
                driver_reasons.append("Hotspot")
            if not driver_reasons:
                continue

            patient_id = row["Tumor_Sample_Barcode"]
            long_rows.append(
                {
                    "patient_id": patient_id,
                    "gene": gene,
                    "variant": f'{row["Chromosome"]}:{row["Start_Position"]}{row["Reference_Allele"]}>{row["Tumor_Seq_Allele2"]}',
                    "protein_change": protein_change,
                    "classification": row["Variant_Classification"],
                    "vaf": row["tumor_vaf"],
                    "driver_evidence": ";".join(driver_reasons),
                }
            )
            matrix.setdefault(patient_id, {})[gene] = "1"
            all_genes.add(gene)

    with args.output_long.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["patient_id", "gene", "variant", "protein_change", "classification", "vaf", "driver_evidence"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(long_rows)

    genes = sorted(all_genes)
    with args.output_matrix.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["patient_id"] + genes)
        for patient_id in sorted(matrix):
            writer.writerow([patient_id] + [matrix[patient_id].get(gene, "0") for gene in genes])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
