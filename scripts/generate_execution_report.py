#!/usr/bin/env python3
"""Generate a fuller static HTML cohort report from workflow outputs."""

from __future__ import annotations

import argparse
import csv
import statistics
from html import escape
from pathlib import Path

import yaml


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_yaml(path: Path):
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def to_float(value: str | None) -> float | None:
    if value in (None, "", "NA", "nan"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def metric_card(title: str, value: str, subtitle: str) -> str:
    return f"""
    <div class="metric-card">
      <div class="metric-title">{escape(title)}</div>
      <div class="metric-value">{escape(value)}</div>
      <div class="metric-subtitle">{escape(subtitle)}</div>
    </div>
    """


def summary_table(rows: list[dict[str, str]], columns: list[tuple[str, str]]) -> str:
    head = "".join(f"<th>{escape(label)}</th>" for label, _key in columns)
    body_rows = []
    for row in rows:
        body_rows.append(
            "<tr>" + "".join(f"<td>{escape(str(row.get(key, 'NA')))}</td>" for _label, key in columns) + "</tr>"
        )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def top_genes(maf_rows: list[dict[str, str]], limit: int = 10) -> list[dict[str, str]]:
    counts: dict[str, int] = {}
    for row in maf_rows:
        gene = row.get("Hugo_Symbol", "NA") or "NA"
        counts[gene] = counts.get(gene, 0) + 1
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[:limit]
    return [{"gene": gene, "variants": str(count)} for gene, count in ranked]


def driver_hits(driver_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    if not driver_rows:
        return []
    row = driver_rows[0]
    patient_key = "patient_id"
    results = []
    for key, value in row.items():
        if key == patient_key or value in {"0", "NA", "", None}:
            continue
        results.append({"feature": key, "value": value})
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-qc", required=True, type=Path)
    parser.add_argument("--maf", required=True, type=Path)
    parser.add_argument("--tmb", required=True, type=Path)
    parser.add_argument("--purity", required=True, type=Path)
    parser.add_argument("--signatures", required=True, type=Path)
    parser.add_argument("--clonality", required=True, type=Path)
    parser.add_argument("--drivers", required=True, type=Path)
    parser.add_argument("--arm-matrix", required=True, type=Path)
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--software-versions", required=True, type=Path)
    parser.add_argument("--containers", required=True, type=Path)
    parser.add_argument("--reference-manifest", required=True, type=Path)
    parser.add_argument("--pipeline-parameters", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    sample_qc = load_tsv(args.sample_qc)
    maf_rows = load_tsv(args.maf)
    tmb_rows = load_tsv(args.tmb)
    purity_rows = load_tsv(args.purity)
    signature_rows = load_tsv(args.signatures)
    clonality_rows = load_tsv(args.clonality)
    driver_rows = load_tsv(args.drivers)
    arm_rows = load_tsv(args.arm_matrix)
    feature_rows = load_tsv(args.features)
    software_rows = load_tsv(args.software_versions)
    container_rows = load_tsv(args.containers)
    reference_rows = load_tsv(args.reference_manifest)
    pipeline_parameters = load_yaml(args.pipeline_parameters)

    patient_count = len({row["patient_id"] for row in feature_rows}) if feature_rows else len({row["patient_id"] for row in sample_qc})
    mean_tmb_values = [to_float(row.get("tmb_mut_per_mb")) for row in tmb_rows]
    mean_tmb_values = [value for value in mean_tmb_values if value is not None]
    purity_values = [to_float(row.get("purity")) for row in purity_rows]
    purity_values = [value for value in purity_values if value is not None]
    subclonal_fracs = [to_float(row.get("fraction_subclonal")) for row in clonality_rows]
    subclonal_fracs = [value for value in subclonal_fracs if value is not None]

    metrics = [
        metric_card("Patients", str(patient_count), "cohort cases in feature matrix"),
        metric_card("Somatic Variants", str(len(maf_rows)), "rows in cohort MAF"),
        metric_card("Median TMB", f"{statistics.median(mean_tmb_values):.2f}" if mean_tmb_values else "NA", "mut/Mb across patients"),
        metric_card("Median Purity", f"{statistics.median(purity_values):.2f}" if purity_values else "NA", "FACETS estimate"),
        metric_card("Median Subclonal Fraction", f"{statistics.median(subclonal_fracs):.2f}" if subclonal_fracs else "NA", "from clonality summary"),
        metric_card("Tracked Software", str(len(software_rows)), "software manifest entries"),
    ]

    pipeline_table = "".join(
        f"<tr><th>{escape(str(key))}</th><td>{escape(str(value))}</td></tr>"
        for key, value in pipeline_parameters.items()
    )
    top_gene_rows = summary_table(top_genes(maf_rows), [("Gene", "gene"), ("Variants", "variants")])
    tmb_table = summary_table(tmb_rows, [("Patient", "patient_id"), ("Tumor", "tumor_sample"), ("Callable Mb", "callable_mb"), ("TMB", "tmb_mut_per_mb")])
    purity_table = summary_table(purity_rows, [("Patient", "patient_id"), ("Purity", "purity"), ("Ploidy", "ploidy"), ("DipLogR", "diplogr")])
    signature_table = summary_table(signature_rows, [("Patient", "patient_id"), ("Top Signature", "top_signature"), ("SBS7 UV", "SBS7_UV"), ("SBS2 APOBEC", "SBS2_APOBEC"), ("Mutation Count", "mutation_count")])
    clonality_table = summary_table(clonality_rows, [("Patient", "patient_id"), ("Clonal", "clonal_mutations"), ("Subclonal", "subclonal_mutations"), ("Ambiguous", "ambiguous_mutations"), ("Fraction Subclonal", "fraction_subclonal")])
    feature_table = summary_table(feature_rows[: min(len(feature_rows), 12)], [(key, key) for key in feature_rows[0].keys()] if feature_rows else [])
    driver_table = summary_table(driver_hits(driver_rows), [("Feature", "feature"), ("Value", "value")]) if driver_rows else "<p>No driver matrix rows available.</p>"
    arm_table = summary_table(arm_rows[: min(len(arm_rows), 10)], [(key, key) for key in arm_rows[0].keys()] if arm_rows else []) if arm_rows else "<p>No arm-level CNV matrix rows available.</p>"
    sample_qc_table = summary_table(sample_qc, [("Patient", "patient_id"), ("Sample", "sample"), ("Type", "sample_type"), ("Mean Cov", "mean_target_coverage"), ("Pct 20x", "pct_target_20x"), ("Dup Rate", "duplicate_rate")])
    software_table = summary_table(software_rows[: min(len(software_rows), 12)], [("Tool", "tool"), ("Version", "version"), ("Source", "source")])
    container_table = summary_table(container_rows[: min(len(container_rows), 12)], [("Tool", "tool"), ("Image", "image"), ("SIF", "sif")])
    reference_table = summary_table(reference_rows[: min(len(reference_rows), 20)], [(key, key) for key in reference_rows[0].keys()] if reference_rows else [])

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Tumor/Normal WES Cohort Report</title>
  <style>
    :root {{
      --bg: #f5f2eb;
      --paper: #fffdfa;
      --ink: #162229;
      --muted: #5e6d77;
      --accent: #0d7c74;
      --accent-soft: #d9f1ec;
      --border: #ded6cb;
    }}
    * {{ box-sizing: border-box; }}
    body {{ margin: 0; font-family: "IBM Plex Sans", "Segoe UI", sans-serif; color: var(--ink); background:
      radial-gradient(circle at top left, #f1ebe1 0%, transparent 32%),
      linear-gradient(180deg, #f7f3ec 0%, #eef6f5 100%); }}
    main {{ max-width: 1360px; margin: 0 auto; padding: 32px 24px 64px; }}
    h1, h2, h3 {{ margin: 0; }}
    p {{ color: var(--muted); line-height: 1.5; }}
    .hero {{ display: grid; gap: 10px; margin-bottom: 26px; }}
    .hero h1 {{ font-size: 2.7rem; }}
    .hero-meta {{ display: flex; flex-wrap: wrap; gap: 10px; }}
    .pill {{ background: var(--accent-soft); color: var(--accent); border-radius: 999px; padding: 6px 10px; font-size: 0.82rem; font-weight: 600; }}
    .metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap: 14px; margin-bottom: 24px; }}
    .metric-card, .panel {{ background: var(--paper); border: 1px solid var(--border); border-radius: 18px; box-shadow: 0 12px 30px rgba(20, 30, 36, 0.05); }}
    .metric-card {{ padding: 16px 18px; }}
    .metric-title {{ color: var(--muted); text-transform: uppercase; letter-spacing: 0.07em; font-size: 0.78rem; }}
    .metric-value {{ font-size: 2rem; font-weight: 700; margin: 8px 0 4px; }}
    .metric-subtitle {{ color: var(--muted); font-size: 0.9rem; }}
    .grid-2 {{ display: grid; grid-template-columns: 1.1fr 0.9fr; gap: 16px; margin-bottom: 16px; }}
    .grid-3 {{ display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 16px; margin-bottom: 16px; }}
    .panel {{ padding: 18px 20px; overflow: hidden; }}
    .panel h2 {{ margin-bottom: 6px; }}
    table {{ width: 100%; border-collapse: collapse; margin-top: 12px; }}
    th, td {{ padding: 10px 8px; border-bottom: 1px solid var(--border); text-align: left; vertical-align: top; font-size: 0.92rem; }}
    th {{ color: var(--muted); font-size: 0.76rem; text-transform: uppercase; letter-spacing: 0.06em; }}
    code {{ background: #f0ebe3; padding: 2px 6px; border-radius: 6px; }}
    .small {{ font-size: 0.85rem; color: var(--muted); }}
    @media (max-width: 1100px) {{
      .grid-2, .grid-3 {{ grid-template-columns: 1fr; }}
    }}
  </style>
</head>
<body>
  <main>
    <section class="hero">
      <h1>Tumor/Normal WES Cohort Report</h1>
      <p>Comprehensive static summary of QC, somatic burden, FACETS purity/ploidy, SigProfiler signatures, clonality, drivers, arm-level CNV states, and provenance for this workflow run.</p>
      <div class="hero-meta">
        <span class="pill">{escape(str(patient_count))} patients</span>
        <span class="pill">{escape(str(len(maf_rows)))} somatic variants</span>
        <span class="pill">{escape(args.pipeline_parameters.name)}</span>
        <span class="pill">{escape(args.reference_manifest.name)}</span>
      </div>
    </section>

    <section class="metrics">
      {''.join(metrics)}
    </section>

    <section class="grid-2">
      <div class="panel">
        <h2>Run Parameters</h2>
        <p>Key run metadata captured for this analysis.</p>
        <table>{pipeline_table}</table>
      </div>
      <div class="panel">
        <h2>Somatic Landscape</h2>
        <p>Most frequently observed genes in the cohort MAF.</p>
        {top_gene_rows}
      </div>
    </section>

    <section class="grid-3">
      <div class="panel">
        <h2>TMB</h2>
        <p>Per-patient mutational burden using callable-megabase denominators.</p>
        {tmb_table}
      </div>
      <div class="panel">
        <h2>Purity & Ploidy</h2>
        <p>FACETS-derived sample-level estimates used downstream for CCF modeling.</p>
        {purity_table}
      </div>
      <div class="panel">
        <h2>Signatures</h2>
        <p>Condensed SigProfiler assignment summary.</p>
        {signature_table}
      </div>
    </section>

    <section class="grid-2">
      <div class="panel">
        <h2>Clonality Summary</h2>
        <p>Counts of clonal, subclonal, ambiguous, and unknown mutation assignments from the CCF model.</p>
        {clonality_table}
      </div>
      <div class="panel">
        <h2>Sample QC</h2>
        <p>Compact alignment/coverage table for all tumor and normal samples.</p>
        {sample_qc_table}
      </div>
    </section>

    <section class="grid-2">
      <div class="panel">
        <h2>Driver Matrix Snapshot</h2>
        <p>First cohort row from the binary driver matrix.</p>
        {driver_table}
      </div>
      <div class="panel">
        <h2>Arm-Level CNV Snapshot</h2>
        <p>First cohort rows from the arm-level status matrix.</p>
        {arm_table}
      </div>
    </section>

    <section class="panel">
      <h2>Feature Matrix Snapshot</h2>
      <p>First rows from the assembled genomic feature matrix used for downstream modeling.</p>
      {feature_table if feature_rows else '<p>No feature rows available.</p>'}
    </section>

    <section class="grid-2" style="margin-top: 16px;">
      <div class="panel">
        <h2>Software Manifest</h2>
        <p>Tools observed in the software version manifest.</p>
        {software_table}
      </div>
      <div class="panel">
        <h2>Container Manifest</h2>
        <p>Pinned container references for workflow processes.</p>
        {container_table}
      </div>
    </section>

    <section class="panel" style="margin-top: 16px;">
      <h2>Reference Manifest</h2>
      <p>Validated reference resources used by the workflow.</p>
      {reference_table if reference_rows else '<p>No reference manifest rows available.</p>'}
      <p class="small">Generated from <code>{escape(args.reference_manifest.name)}</code>, <code>{escape(args.software_versions.name)}</code>, and <code>{escape(args.containers.name)}</code>.</p>
    </section>
  </main>
</body>
</html>
"""
    args.output.write_text(html, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
