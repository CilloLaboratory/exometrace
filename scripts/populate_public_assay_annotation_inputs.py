#!/usr/bin/env python3
"""Download public assay/annotation inputs and materialize workflow-ready tables."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import stat
import shutil
import urllib.request
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


KAPA_HYPEREXOME_BB_URL = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/exomeProbesets/KAPA_HyperExome_hg38_capture_targets.bb"
UCSC_BIGBEDTOBED_URL = "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigBedToBed"
ONCOKB_CANCER_GENES_URL = "https://www.oncokb.org/api/v1/utils/allCuratedGenes"
CANCER_HOTSPOTS_URL = "https://www.cancerhotspots.org/api/hotspots/single?version=v3"
GENCODE_GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz"
UCSC_CYTOBAND_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"

AA3 = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "*": "Ter",
}


@dataclass
class Targets:
    chrom: str
    start0: int
    end: int
    name: str


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0:
        return
    request = urllib.request.Request(url, headers={"User-Agent": "wes-workflow/1.0"})
    with urllib.request.urlopen(request) as response, destination.open("wb") as handle:
        shutil.copyfileobj(response, handle)


def ensure_executable(path: Path) -> None:
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def convert_bigbed_to_bed(bigbed_to_bed: Path, bigbed_path: Path, bed_path: Path) -> None:
    import subprocess

    bed_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run([str(bigbed_to_bed), str(bigbed_path), str(bed_path)], check=True)


def parse_dict_headers(dict_path: Path) -> tuple[list[str], dict[str, int]]:
    headers: list[str] = []
    contig_order: dict[str, int] = {}
    with dict_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("@"):
                continue
            headers.append(line.rstrip("\n"))
            if not line.startswith("@SQ"):
                continue
            fields = dict(item.split(":", 1) for item in line.rstrip("\n").split("\t")[1:])
            if "SN" in fields:
                contig_order[fields["SN"]] = len(contig_order)
    return headers, contig_order


def load_targets(bed_path: Path, contig_order: dict[str, int]) -> list[Targets]:
    rows: list[Targets] = []
    with bed_path.open(encoding="utf-8") as handle:
        for index, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            if chrom not in contig_order:
                continue
            start0 = int(fields[1])
            end = int(fields[2])
            name = fields[3] if len(fields) > 3 and fields[3] else f"target_{index}"
            rows.append(Targets(chrom=chrom, start0=start0, end=end, name=name))
    rows.sort(key=lambda row: (contig_order[row.chrom], row.start0, row.end, row.name))
    return rows


def write_interval_list(interval_list_path: Path, dict_headers: list[str], targets: list[Targets]) -> None:
    interval_list_path.parent.mkdir(parents=True, exist_ok=True)
    with interval_list_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("@HD\tVN:1.6\tSO:coordinate\n")
        for header in dict_headers:
            if header.startswith("@SQ"):
                handle.write(f"{header}\n")
        for row in targets:
            handle.write(f"{row.chrom}\t{row.start0 + 1}\t{row.end}\t+\t{row.name}\n")


def write_targets_bed(target_bed_path: Path, targets: list[Targets]) -> None:
    with target_bed_path.open("w", encoding="utf-8", newline="") as handle:
        for row in targets:
            handle.write(f"{row.chrom}\t{row.start0}\t{row.end}\t{row.name}\n")


def write_callable_regions(callable_path: Path, target_bed_path: Path) -> None:
    callable_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(target_bed_path, callable_path)


def write_cancer_gene_census(output_path: Path, raw_json_path: Path) -> None:
    genes = json.loads(raw_json_path.read_text(encoding="utf-8"))
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["gene", "source", "gene_type", "oncokb_annotated", "sanger_cgc", "vogelstein"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in genes:
            writer.writerow(
                {
                    "gene": row["hugoSymbol"],
                    "source": "OncoKB_allCuratedGenes",
                    "gene_type": row.get("geneType") or "",
                    "oncokb_annotated": "1" if row.get("oncokbAnnotated") else "0",
                    "sanger_cgc": "1" if row.get("sangerCGC") else "0",
                    "vogelstein": "1" if row.get("vogelstein") else "0",
                }
            )


def aa1_to_aa3(residue: str) -> str:
    if residue not in AA3:
        raise ValueError(f"Unsupported amino acid code: {residue}")
    return AA3[residue]


def write_hotspots(output_path: Path, raw_json_path: Path) -> None:
    entries = json.loads(raw_json_path.read_text(encoding="utf-8"))
    rows: list[dict[str, str]] = []
    for entry in entries:
        if entry.get("type") != "single residue":
            continue
        gene = entry["hugoSymbol"]
        residue = entry["residue"]
        if len(residue) < 2 or not residue[1:].isdigit():
            continue
        transcript_id = entry.get("transcriptId") or ""
        if residue[0] not in AA3:
            continue
        ref_aa = aa1_to_aa3(residue[0])
        position = residue[1:]
        variant_aas = entry.get("variantAminoAcid") or {}
        for alt_aa in sorted(variant_aas):
            if alt_aa not in AA3:
                continue
            rows.append(
                {
                    "gene": gene,
                    "protein_change": f"p.{ref_aa}{position}{aa1_to_aa3(alt_aa)}",
                    "transcript_id": transcript_id,
                    "source": "CancerHotspots_v3",
                }
            )
    rows.sort(key=lambda row: (row["gene"], row["protein_change"]))
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene", "protein_change", "transcript_id", "source"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_gtf_attributes(field: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for item in field.split(";"):
        item = item.strip()
        if not item:
            continue
        key, value = item.split(" ", 1)
        values[key] = value.strip().strip('"')
    return values


def write_gene_coordinates(output_path: Path, gtf_path: Path) -> None:
    genes: dict[str, tuple[str, int, int, str, str]] = {}
    with gzip.open(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = parse_gtf_attributes(fields[8])
            gene_name = attrs.get("gene_name")
            if not gene_name:
                continue
            genes[gene_name] = (
                fields[0],
                int(fields[3]),
                int(fields[4]),
                attrs.get("gene_id", ""),
                attrs.get("gene_type", ""),
            )
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene", "chrom", "start", "end", "gene_id", "gene_type"], delimiter="\t")
        writer.writeheader()
        for gene in sorted(genes):
            chrom, start, end, gene_id, gene_type = genes[gene]
            writer.writerow(
                {
                    "gene": gene,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "gene_id": gene_id,
                    "gene_type": gene_type,
                }
            )


def write_chromosome_arms(output_path: Path, cytoband_path: Path) -> None:
    grouped: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    with gzip.open(cytoband_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            chrom, start, end, band, _stain = line.rstrip("\n").split("\t")
            if chrom in {"chrM", "chrUn"} or "_" in chrom:
                continue
            arm = band[0]
            if arm not in {"p", "q"}:
                continue
            chrom_label = chrom[3:] if chrom.startswith("chr") else chrom
            grouped[(chrom, f"{chrom_label}{arm}")].append((int(start), int(end)))
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        for (chrom, arm), spans in sorted(grouped.items(), key=lambda item: (item[0][0], item[0][1])):
            start = min(span[0] for span in spans)
            end = max(span[1] for span in spans)
            handle.write(f"{chrom}\t{start}\t{end}\t{arm}\n")


def write_provenance(output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["label", "url", "notes"], delimiter="\t")
        writer.writeheader()
        writer.writerows(
            [
                {
                    "label": "targets_bed",
                    "url": KAPA_HYPEREXOME_BB_URL,
                    "notes": "UCSC mirrored Roche KAPA HyperExome hg38 capture-target BigBed converted locally to BED.",
                },
                {
                    "label": "targets_bed_converter",
                    "url": UCSC_BIGBEDTOBED_URL,
                    "notes": "UCSC bigBedToBed utility used to materialize the mirrored KAPA assay track as BED.",
                },
                {
                    "label": "cancer_gene_census",
                    "url": ONCOKB_CANCER_GENES_URL,
                    "notes": "Public OncoKB curated cancer gene JSON transformed into workflow TSV.",
                },
                {
                    "label": "hotspots",
                    "url": CANCER_HOTSPOTS_URL,
                    "notes": "Cancer Hotspots single-residue v3 API transformed into HGVS protein-change TSV.",
                },
                {
                    "label": "gene_coordinates",
                    "url": GENCODE_GTF_URL,
                    "notes": "GENCODE gene features collapsed to one row per gene symbol.",
                },
                {
                    "label": "chromosome_arms",
                    "url": UCSC_CYTOBAND_URL,
                    "notes": "UCSC hg38 cytobands aggregated to p/q arm spans.",
                },
            ]
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dict", type=Path, default=Path("references/GRCh38/fasta/GRCh38.dict"))
    parser.add_argument("--targets-bed", type=Path, default=Path("references/GRCh38/intervals/exome_targets.bed"))
    parser.add_argument("--interval-list", type=Path, default=Path("references/GRCh38/intervals/exome_targets.interval_list"))
    parser.add_argument("--callable-regions", type=Path, default=Path("references/GRCh38/intervals/callable_regions.bed"))
    parser.add_argument("--cancer-gene-census", type=Path, default=Path("references/GRCh38/annotations/cancer_gene_census.tsv"))
    parser.add_argument("--hotspots", type=Path, default=Path("references/GRCh38/annotations/hotspots.tsv"))
    parser.add_argument("--gene-coordinates", type=Path, default=Path("references/GRCh38/annotations/gene_coordinates.tsv"))
    parser.add_argument("--chromosome-arms", type=Path, default=Path("references/GRCh38/annotations/chromosome_arms.bed"))
    parser.add_argument("--provenance", type=Path, default=Path("references/GRCh38/annotations/public_annotation_sources.tsv"))
    parser.add_argument("--scratch-dir", type=Path, default=Path("references/GRCh38/.downloads"))
    parser.add_argument("--bigbed-to-bed", type=Path, default=Path("references/GRCh38/.downloads/bigBedToBed"))
    return parser


def main() -> int:
    args = build_parser().parse_args()
    args.scratch_dir.mkdir(parents=True, exist_ok=True)

    kapa_bigbed = args.scratch_dir / "KAPA_HyperExome_hg38_capture_targets.bb"
    bed_download = args.scratch_dir / "KAPA_HyperExome_hg38_capture_targets.bed"
    oncokb_json = args.scratch_dir / "oncokb_allCuratedGenes.json"
    hotspots_json = args.scratch_dir / "cancerhotspots_single_v3.json"
    gtf_gz = args.scratch_dir / "gencode.v50.annotation.gtf.gz"
    cytoband_gz = args.scratch_dir / "cytoBand.txt.gz"

    download(UCSC_BIGBEDTOBED_URL, args.bigbed_to_bed)
    ensure_executable(args.bigbed_to_bed)
    download(KAPA_HYPEREXOME_BB_URL, kapa_bigbed)
    convert_bigbed_to_bed(args.bigbed_to_bed, kapa_bigbed, bed_download)
    download(ONCOKB_CANCER_GENES_URL, oncokb_json)
    download(CANCER_HOTSPOTS_URL, hotspots_json)
    download(GENCODE_GTF_URL, gtf_gz)
    download(UCSC_CYTOBAND_URL, cytoband_gz)

    dict_headers, contig_order = parse_dict_headers(args.reference_dict)
    targets = load_targets(bed_download, contig_order)
    write_targets_bed(args.targets_bed, targets)
    write_interval_list(args.interval_list, dict_headers, targets)
    write_callable_regions(args.callable_regions, args.targets_bed)
    write_cancer_gene_census(args.cancer_gene_census, oncokb_json)
    write_hotspots(args.hotspots, hotspots_json)
    write_gene_coordinates(args.gene_coordinates, gtf_gz)
    write_chromosome_arms(args.chromosome_arms, cytoband_gz)
    write_provenance(args.provenance)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
