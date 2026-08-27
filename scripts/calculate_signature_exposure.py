#!/usr/bin/env python3
"""Run SigProfilerAssignment and emit condensed signature exposure/QC tables."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import yaml


SIGPROFILER_MARKERS = (
    "tsb/GRCh38",
    "references/chromosomes/transcripts/GRCh38",
    "references/chromosomes/tsb/GRCh38",
)


def load_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def resolve_config_path(config_path: Path, raw_path: str | None) -> Path | None:
    if not raw_path:
        return None
    path = Path(raw_path)
    if path.is_absolute():
        return path
    base_dir = config_path.parent.parent if config_path.parent.name == "config" else config_path.parent
    return (base_dir / path).resolve()


def load_maf_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalize_contig(chromosome: str, contig_style: str) -> str:
    if contig_style == "chr" and not chromosome.startswith("chr"):
        return f"chr{chromosome}"
    if contig_style != "chr" and chromosome.startswith("chr"):
        return chromosome[3:]
    return chromosome


def write_vcf_from_maf(maf_rows: list[dict[str, str]], output_vcf: Path, sample_name: str, contig_style: str) -> None:
    with output_vcf.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write(f"##tumor_sample={sample_name}\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in maf_rows:
            handle.write(
                "\t".join(
                    [
                        normalize_contig(row["Chromosome"], contig_style),
                        row["Start_Position"],
                        ".",
                        row["Reference_Allele"],
                        row["Tumor_Seq_Allele2"],
                        ".",
                        "PASS",
                        ".",
                    ]
                )
                + "\n"
            )


def parse_activity_table(path: Path, patient_id: str) -> dict[str, float]:
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        return {}

    header = rows[0]
    if len(header) > 1 and any(cell.startswith("SBS") for cell in header[1:]):
        for row in rows[1:]:
            if row and row[0] == patient_id:
                return {
                    signature: float(value or 0.0)
                    for signature, value in zip(header[1:], row[1:])
                    if signature.startswith("SBS")
                }

    if patient_id in header:
        sample_index = header.index(patient_id)
        exposures: dict[str, float] = {}
        for row in rows[1:]:
            if len(row) > sample_index and row and row[0].startswith("SBS"):
                exposures[row[0]] = float(row[sample_index] or 0.0)
        return exposures

    return {}


def read_signature_activities(output_dir: Path, patient_id: str) -> dict[str, float]:
    candidate_paths = sorted(
        [
            path
            for path in output_dir.rglob("*")
            if path.is_file()
            and path.suffix in {".txt", ".tsv"}
            and ("activities" in path.name.lower() or "assignment_solution" in path.name.lower())
        ],
        key=lambda path: (0 if "activities" in path.name.lower() else 1, len(path.parts), path.as_posix()),
    )
    for path in candidate_paths:
        exposures = parse_activity_table(path, patient_id)
        if exposures:
            return exposures
    raise FileNotFoundError(f"Unable to locate SigProfiler activities table for {patient_id} under {output_dir}")


def summarize_exposures(patient_id: str, mutation_count: int, exposures: dict[str, float], min_mutations_for_qc_pass: int) -> tuple[dict[str, str], dict[str, str]]:
    total_assigned = sum(exposures.values())
    top_signature = max(exposures, key=exposures.get) if exposures else "NA"
    uv = sum(value for name, value in exposures.items() if name.startswith("SBS7"))
    apobec = sum(value for name, value in exposures.items() if name in {"SBS2", "SBS13"})

    if not exposures or total_assigned <= 0:
        status = "NO_SIGNATURES"
    elif mutation_count < min_mutations_for_qc_pass:
        status = "LOW_COUNT"
    else:
        status = "PASS"

    exposure_row: dict[str, str] = {
        "patient_id": patient_id,
        "SBS7_UV": f"{uv:.6f}",
        "SBS2_APOBEC": f"{apobec:.6f}",
        "mutation_count": str(mutation_count),
        "top_signature": top_signature,
        "total_assigned": f"{total_assigned:.6f}",
    }
    for signature in sorted(exposures):
        exposure_row[signature] = f"{exposures[signature]:.6f}"

    qc_row = {
        "patient_id": patient_id,
        "mutation_count": str(mutation_count),
        "top_signature": top_signature,
        "total_assigned": f"{total_assigned:.6f}",
        "status": status,
    }
    return exposure_row, qc_row


def write_single_row_tsv(path: Path, row: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def write_no_signatures_outputs(
    patient_id: str,
    mutation_count: int,
    pipeline_config: dict,
    reference_config: dict,
    output_exposure: Path,
    output_qc: Path,
) -> None:
    exposure_row, qc_row = summarize_exposures(
        patient_id,
        mutation_count,
        {},
        int(pipeline_config.get("signatures", {}).get("min_mutations_for_qc_pass", 50)),
    )
    qc_row.update(
        {
            "context_type": str(pipeline_config.get("signatures", {}).get("context_type", "96")),
            "cosmic_version": str(pipeline_config.get("signatures", {}).get("cosmic_version", 3.6)),
            "genome_build": str(reference_config["reference"]["build"]),
            "exome": str(bool(pipeline_config.get("signatures", {}).get("exome", True))).lower(),
        }
    )
    write_single_row_tsv(output_exposure, exposure_row)
    write_single_row_tsv(output_qc, qc_row)


def run_sigprofiler(samples_dir: Path, output_dir: Path, genome_build: str, cosmic_version: float, exome: bool, context_type: str, make_plots: bool, sample_reconstruction_plots: str, cpu: int, volume: str | None) -> None:
    try:
        from SigProfilerAssignment import Analyzer as Analyze
    except ImportError as exc:
        raise RuntimeError("SigProfilerAssignment is not installed in the active environment") from exc

    kwargs = {
        "samples": str(samples_dir),
        "output": str(output_dir),
        "input_type": "vcf",
        "context_type": context_type,
        "collapse_to_SBS96": context_type == "96",
        "cosmic_version": cosmic_version,
        "exome": exome,
        "genome_build": genome_build,
        "export_probabilities": False,
        "export_probabilities_per_mutation": False,
        "make_plots": make_plots,
        "sample_reconstruction_plots": sample_reconstruction_plots,
        "verbose": False,
        "cpu": cpu,
    }
    if volume:
        kwargs["volume"] = volume
    Analyze.cosmic_fit(**kwargs)


def validate_sigprofiler_volume(path: Path) -> None:
    if not path.is_dir():
        raise FileNotFoundError(f"Configured SigProfiler volume does not exist or is not a directory: {path}")
    if not any((path / marker).exists() for marker in SIGPROFILER_MARKERS):
        raise FileNotFoundError(
            "Configured SigProfiler volume is missing required assets under "
            f"{path}. Expected one of: {', '.join(SIGPROFILER_MARKERS)}"
        )


def resolve_sigprofiler_volume(
    pipeline_config: dict,
    pipeline_config_path: Path,
    reference_config: dict,
    reference_config_path: Path,
) -> str | None:
    signatures_cfg = pipeline_config.get("signatures", {})
    reference_cfg = reference_config.get("sigprofiler", {})
    configured_volume = signatures_cfg.get("volume")
    if configured_volume:
        path = resolve_config_path(pipeline_config_path, str(configured_volume))
        validate_sigprofiler_volume(path)
        return str(path)

    reference_volume = reference_cfg.get("volume_dir")
    if reference_volume:
        path = resolve_config_path(reference_config_path, str(reference_volume))
        validate_sigprofiler_volume(path)
        return str(path)

    environment_volume = os.environ.get("SIGPROFILERASSIGNMENT_VOLUME")
    if environment_volume:
        path = Path(environment_volume).resolve()
        validate_sigprofiler_volume(path)
        return str(path)

    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--maf", required=True, type=Path)
    parser.add_argument("--reference-config", required=True, type=Path)
    parser.add_argument("--pipeline-config", required=True, type=Path)
    parser.add_argument("--output-exposure", required=True, type=Path)
    parser.add_argument("--output-qc", required=True, type=Path)
    parser.add_argument("--cpu", type=int, default=1)
    parser.add_argument("--stub", action="store_true")
    args = parser.parse_args()

    maf_rows = load_maf_rows(args.maf)
    mutation_count = len(maf_rows)
    pipeline_config = load_yaml(args.pipeline_config)
    reference_config = load_yaml(args.reference_config)

    if mutation_count == 0:
        write_no_signatures_outputs(
            args.patient_id,
            mutation_count,
            pipeline_config,
            reference_config,
            args.output_exposure,
            args.output_qc,
        )
        return 0

    if args.stub:
        exposures = {"SBS1": float(mutation_count), "SBS7a": 0.0, "SBS13": 0.0}
        exposure_row, qc_row = summarize_exposures(
            args.patient_id,
            mutation_count,
            exposures,
            int(pipeline_config.get("signatures", {}).get("min_mutations_for_qc_pass", 50)),
        )
        qc_row.update(
            {
                "context_type": str(pipeline_config.get("signatures", {}).get("context_type", "96")),
                "cosmic_version": str(pipeline_config.get("signatures", {}).get("cosmic_version", 3.6)),
                "genome_build": str(reference_config["reference"]["build"]),
                "exome": str(bool(pipeline_config.get("signatures", {}).get("exome", True))).lower(),
            }
        )
        write_single_row_tsv(args.output_exposure, exposure_row)
        write_single_row_tsv(args.output_qc, qc_row)
        return 0

    signatures_cfg = pipeline_config.get("signatures", {})
    sigprofiler_work = args.output_exposure.parent / f"{args.patient_id}.sigprofiler"
    input_dir = sigprofiler_work / "input"
    output_dir = sigprofiler_work / "output"
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_vcf = input_dir / f"{args.patient_id}.vcf"
    write_vcf_from_maf(maf_rows, sample_vcf, args.patient_id, str(reference_config["reference"].get("contig_style", "chr")))

    volume = resolve_sigprofiler_volume(
        pipeline_config,
        args.pipeline_config,
        reference_config,
        args.reference_config,
    )
    run_sigprofiler(
        input_dir,
        output_dir,
        genome_build=str(reference_config["reference"]["build"]),
        cosmic_version=float(signatures_cfg.get("cosmic_version", 3.6)),
        exome=bool(signatures_cfg.get("exome", True)),
        context_type=str(signatures_cfg.get("context_type", "96")),
        make_plots=bool(signatures_cfg.get("make_plots", False)),
        sample_reconstruction_plots=str(signatures_cfg.get("sample_reconstruction_plots", "none")),
        cpu=args.cpu,
        volume=str(volume) if volume else None,
    )

    exposures = read_signature_activities(output_dir, args.patient_id)
    exposure_row, qc_row = summarize_exposures(
        args.patient_id,
        mutation_count,
        exposures,
        int(signatures_cfg.get("min_mutations_for_qc_pass", 50)),
    )
    qc_row.update(
        {
            "context_type": str(signatures_cfg.get("context_type", "96")),
            "cosmic_version": str(signatures_cfg.get("cosmic_version", 3.6)),
            "genome_build": str(reference_config["reference"]["build"]),
            "exome": str(bool(signatures_cfg.get("exome", True))).lower(),
        }
    )

    write_single_row_tsv(args.output_exposure, exposure_row)
    write_single_row_tsv(args.output_qc, qc_row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
