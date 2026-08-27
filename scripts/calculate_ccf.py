#!/usr/bin/env python3
"""Estimate CCF and clonality from MAF plus FACETS-derived purity/CN segments."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import yaml


def load_single_row(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows[0] if rows else {}


def load_segments(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalize_chromosome(value: str) -> str:
    return value[3:] if value.startswith("chr") else value


def find_segment(segments: list[dict[str, str]], chromosome: str, position_1based: int) -> dict[str, str] | None:
    chrom = normalize_chromosome(chromosome)
    for segment in segments:
        if normalize_chromosome(segment["chromosome"]) != chrom:
            continue
        start = int(segment["start"])
        end = int(segment["end"])
        if start <= position_1based - 1 < end:
            return segment
    return None


def expected_mutant_copies(vaf: float, total_cn: float, purity: float) -> int:
    adjusted_total_cn = 1.0 if total_cn == 0 else total_cn
    mu = vaf * (1.0 / purity) * (purity * adjusted_total_cn + (1.0 - purity) * 2.0)
    alt_copies = 1.0 if mu < 1.0 else abs(mu)
    return round(alt_copies)


def expected_vaf(ccf: float, purity: float, total_cn: float, mutant_copies: int) -> float:
    return purity * ccf * mutant_copies / (2.0 * (1.0 - purity) + purity * total_cn)


def binomial_pmf(k: int, n: int, p: float) -> float:
    if p <= 0.0:
        return 1.0 if k == 0 else 0.0
    if p >= 1.0:
        return 1.0 if k == n else 0.0
    log_coeff = math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    log_prob = log_coeff + (k * math.log(p)) + ((n - k) * math.log(1.0 - p))
    return math.exp(log_prob)


def infer_ccf(t_alt_count: int, t_depth: int, purity: float, total_cn: float, mutant_copies: int) -> tuple[float, float, float, float]:
    ccfs = [index / 1000.0 for index in range(1, 1001)]
    probs = [
        binomial_pmf(t_alt_count, t_depth, expected_vaf(ccf, purity, total_cn, mutant_copies))
        for ccf in ccfs
    ]
    total_prob = sum(probs)
    if total_prob == 0.0:
        return (float("nan"), float("nan"), float("nan"), 0.0)
    probs = [value / total_prob for value in probs]
    max_prob = max(probs)
    max_index = probs.index(max_prob)
    half_max = [index for index, value in enumerate(probs) if value > max_prob / 2.0]
    lower_index = max(half_max[0] - 1, 0)
    upper_index = min(half_max[-1] + 1, len(ccfs) - 1)
    prob_clonal = sum(prob for ccf, prob in zip(ccfs, probs) if ccf >= 0.9)
    return (ccfs[max_index], ccfs[lower_index], ccfs[upper_index], prob_clonal)


def classify_clonality(ccf: float, ccf_lower: float, ccf_upper: float, clonal_threshold: float) -> str:
    if math.isnan(ccf):
        return "unknown"
    if ccf_lower >= clonal_threshold:
        return "clonal"
    if ccf_upper < clonal_threshold:
        return "subclonal"
    return "ambiguous"


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--maf", required=True, type=Path)
    parser.add_argument("--purity", required=True, type=Path)
    parser.add_argument("--segments", required=True, type=Path)
    parser.add_argument("--pipeline-config", required=True, type=Path)
    parser.add_argument("--output-ccf", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    args = parser.parse_args()

    with args.maf.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    purity_row = load_single_row(args.purity)
    segments = load_segments(args.segments)
    pipeline_config = load_config(args.pipeline_config)
    purity_value = purity_row.get("purity", "NA")
    purity = float(purity_value) if purity_value not in {"NA", "", None} else float("nan")
    clonal_threshold = float(pipeline_config.get("clonality", {}).get("clonal_ccf_threshold", 0.9))

    detailed_rows: list[dict[str, str]] = []
    counts = {"clonal": 0, "subclonal": 0, "ambiguous": 0, "unknown": 0}

    for row in rows:
        depth = int(row["t_depth"] or 0)
        alt_count = int(row["t_alt_count"] or 0)
        vaf = float(row["tumor_vaf"])
        segment = find_segment(segments, row["Chromosome"], int(row["Start_Position"]))

        if segment is None or math.isnan(purity) or purity <= 0.0 or depth <= 0:
            clonality = "unknown"
            total_cn = "NA"
            major_cn = "NA"
            minor_cn = "NA"
            multiplicity = "NA"
            ccf = float("nan")
            ccf_lower = float("nan")
            ccf_upper = float("nan")
            prob_clonal = 0.0
        else:
            total_cn_value = float(segment["total_cn"])
            major_cn_value = int(round(float(segment["major_cn"])))
            minor_cn_value = int(round(float(segment["minor_cn"])))
            max_mutant_copies = max(1, major_cn_value if major_cn_value > 0 else int(round(total_cn_value)) or 1)
            multiplicity_value = min(max(1, expected_mutant_copies(vaf, total_cn_value, purity)), max_mutant_copies)
            ccf, ccf_lower, ccf_upper, prob_clonal = infer_ccf(
                alt_count,
                depth,
                purity,
                total_cn_value if total_cn_value > 0 else 1.0,
                multiplicity_value,
            )
            clonality = classify_clonality(ccf, ccf_lower, ccf_upper, clonal_threshold)
            total_cn = f"{total_cn_value:.6f}"
            major_cn = str(major_cn_value)
            minor_cn = str(minor_cn_value)
            multiplicity = str(multiplicity_value)

        counts[clonality] += 1
        detailed_rows.append(
            {
                "patient_id": args.patient_id,
                "gene": row["Hugo_Symbol"],
                "chromosome": row["Chromosome"],
                "start_position": row["Start_Position"],
                "tumor_vaf": f"{vaf:.6f}",
                "local_total_copy_number": total_cn,
                "major_copy_number": major_cn,
                "minor_copy_number": minor_cn,
                "multiplicity": multiplicity,
                "cancer_cell_fraction": "NA" if math.isnan(ccf) else f"{ccf:.6f}",
                "ccf_lower": "NA" if math.isnan(ccf_lower) else f"{ccf_lower:.6f}",
                "ccf_upper": "NA" if math.isnan(ccf_upper) else f"{ccf_upper:.6f}",
                "prob_clonal": f"{prob_clonal:.6f}",
                "clonality": clonality,
            }
        )

    with args.output_ccf.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "patient_id",
                "gene",
                "chromosome",
                "start_position",
                "tumor_vaf",
                "local_total_copy_number",
                "major_copy_number",
                "minor_copy_number",
                "multiplicity",
                "cancer_cell_fraction",
                "ccf_lower",
                "ccf_upper",
                "prob_clonal",
                "clonality",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(detailed_rows)

    classified_total = counts["clonal"] + counts["subclonal"] + counts["ambiguous"]
    fraction_subclonal = counts["subclonal"] / classified_total if classified_total else 0.0
    fraction_ambiguous = counts["ambiguous"] / classified_total if classified_total else 0.0
    number_of_clusters = len([label for label in ("clonal", "subclonal", "ambiguous") if counts[label] > 0])

    with args.output_summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "patient_id",
                "clonal_mutations",
                "subclonal_mutations",
                "ambiguous_mutations",
                "unknown_mutations",
                "fraction_subclonal",
                "fraction_ambiguous",
                "number_of_clusters",
            ]
        )
        writer.writerow(
            [
                args.patient_id,
                counts["clonal"],
                counts["subclonal"],
                counts["ambiguous"],
                counts["unknown"],
                f"{fraction_subclonal:.6f}",
                f"{fraction_ambiguous:.6f}",
                number_of_clusters,
            ]
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
