#!/usr/bin/env python3
"""Generate a richer static QC dashboard from the cohort sample QC table."""

from __future__ import annotations

import argparse
import csv
import statistics
from html import escape
from pathlib import Path


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def to_float(row: dict[str, str], key: str) -> float:
    return float(row.get(key, "0") or 0.0)


def to_int(row: dict[str, str], key: str) -> int:
    return int(float(row.get(key, "0") or 0.0))


def fmt_float(value: float, digits: int = 2) -> str:
    return f"{value:.{digits}f}"


def metric_card(title: str, value: str, subtitle: str) -> str:
    return f"""
    <div class="card metric">
      <div class="eyebrow">{escape(title)}</div>
      <div class="value">{escape(value)}</div>
      <div class="subtitle">{escape(subtitle)}</div>
    </div>
    """


def threshold_flag(row: dict[str, str]) -> list[str]:
    flags: list[str] = []
    if to_float(row, "mapping_rate") < 0.95:
        flags.append("mapping<95%")
    if to_float(row, "duplicate_rate") > 0.50:
        flags.append("dup>50%")
    if row.get("sample_type") == "tumor" and to_float(row, "mean_target_coverage") < 80:
        flags.append("tumor_cov<80x")
    if row.get("sample_type") == "normal" and to_float(row, "mean_target_coverage") < 40:
        flags.append("normal_cov<40x")
    if to_float(row, "pct_target_20x") < 0.90:
        flags.append("20x<90%")
    return flags


def patient_rows(rows: list[dict[str, str]]) -> list[str]:
    rendered: list[str] = []
    for row in sorted(rows, key=lambda entry: (entry["patient_id"], entry["sample_type"], entry["sample"])):
        flags = threshold_flag(row)
        flag_html = (
            "".join(f'<span class="chip warn">{escape(flag)}</span>' for flag in flags)
            if flags
            else '<span class="chip ok">PASS</span>'
        )
        rendered.append(
            f"""
            <tr>
              <td>{escape(row['patient_id'])}</td>
              <td>{escape(row['sample'])}</td>
              <td>{escape(row['sample_type'])}</td>
              <td>{to_int(row, 'mapped_reads')}</td>
              <td>{fmt_float(100 * to_float(row, 'mapping_rate'))}%</td>
              <td>{fmt_float(to_float(row, 'mean_target_coverage'))}x</td>
              <td>{fmt_float(100 * to_float(row, 'pct_target_20x'))}%</td>
              <td>{fmt_float(100 * to_float(row, 'duplicate_rate'))}%</td>
              <td>{fmt_float(to_float(row, 'mean_insert_size'))}</td>
              <td>{flag_html}</td>
            </tr>
            """
        )
    return rendered


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-qc", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    rows = load_rows(args.sample_qc)
    mapping_rates = [100 * to_float(row, "mapping_rate") for row in rows]
    coverages = [to_float(row, "mean_target_coverage") for row in rows]
    duplicates = [100 * to_float(row, "duplicate_rate") for row in rows]
    flagged = sum(1 for row in rows if threshold_flag(row))
    unique_patients = len({row["patient_id"] for row in rows})

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Tumor/Normal WES QC Report</title>
  <style>
    :root {{
      --bg: #f4efe7;
      --panel: #fffdf8;
      --ink: #182126;
      --muted: #66737d;
      --accent: #0f7b6c;
      --warn: #b85c38;
      --border: #d9d1c6;
    }}
    * {{ box-sizing: border-box; }}
    body {{ margin: 0; font-family: "IBM Plex Sans", "Segoe UI", sans-serif; background: linear-gradient(180deg, #f7f1e8 0%, #f1f6f5 100%); color: var(--ink); }}
    main {{ max-width: 1280px; margin: 0 auto; padding: 32px 24px 56px; }}
    h1, h2 {{ margin: 0; }}
    p {{ color: var(--muted); }}
    .hero {{ display: grid; gap: 10px; margin-bottom: 24px; }}
    .hero h1 {{ font-size: 2.4rem; }}
    .metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 14px; margin: 24px 0; }}
    .card {{ background: var(--panel); border: 1px solid var(--border); border-radius: 16px; padding: 16px 18px; box-shadow: 0 10px 30px rgba(24, 33, 38, 0.05); }}
    .metric .eyebrow {{ color: var(--muted); font-size: 0.78rem; text-transform: uppercase; letter-spacing: 0.08em; }}
    .metric .value {{ font-size: 2rem; font-weight: 700; margin: 8px 0 4px; }}
    .metric .subtitle {{ color: var(--muted); font-size: 0.92rem; }}
    table {{ width: 100%; border-collapse: collapse; }}
    th, td {{ padding: 12px 10px; border-bottom: 1px solid var(--border); text-align: left; vertical-align: top; font-size: 0.94rem; }}
    th {{ color: var(--muted); font-size: 0.78rem; text-transform: uppercase; letter-spacing: 0.06em; }}
    .table-card {{ overflow: hidden; }}
    .chip {{ display: inline-block; border-radius: 999px; padding: 4px 8px; margin: 0 4px 4px 0; font-size: 0.75rem; }}
    .chip.ok {{ background: #e0f3eb; color: #0a6a4c; }}
    .chip.warn {{ background: #f9e4dc; color: #934629; }}
  </style>
</head>
<body>
  <main>
    <section class="hero">
      <h1>Quality Control Dashboard</h1>
      <p>Static cohort summary generated from <code>{escape(args.sample_qc.name)}</code>. This report highlights sample-level alignment, coverage, and duplication metrics together with simple threshold flags.</p>
    </section>
    <section class="metrics">
      {metric_card("Patients", str(unique_patients), "paired tumor/normal cases")}
      {metric_card("Samples", str(len(rows)), "rows in cohort QC table")}
      {metric_card("Median Mapping", f"{fmt_float(statistics.median(mapping_rates))}%", "across all samples")}
      {metric_card("Median Coverage", f"{fmt_float(statistics.median(coverages))}x", "mean target coverage")}
      {metric_card("Median Duplication", f"{fmt_float(statistics.median(duplicates))}%", "Picard duplication rate")}
      {metric_card("Flagged Samples", str(flagged), "samples with threshold flags")}
    </section>
    <section class="card table-card">
      <h2>Sample QC Table</h2>
      <p>Thresholds used here mirror the current workflow defaults: mapping rate ≥95%, duplicate rate ≤50%, target coverage ≥80x for tumors and ≥40x for normals, and ≥90% of target bases at 20x.</p>
      <table>
        <thead>
          <tr>
            <th>Patient</th>
            <th>Sample</th>
            <th>Type</th>
            <th>Mapped Reads</th>
            <th>Mapping Rate</th>
            <th>Mean Target Cov</th>
            <th>Pct Target 20x</th>
            <th>Duplicate Rate</th>
            <th>Insert Size</th>
            <th>Flags</th>
          </tr>
        </thead>
        <tbody>
          {''.join(patient_rows(rows))}
        </tbody>
      </table>
    </section>
  </main>
</body>
</html>
"""
    args.output.write_text(html, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
