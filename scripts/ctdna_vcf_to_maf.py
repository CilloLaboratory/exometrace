#!/usr/bin/env python3
"""Build ctDNA MAF-like tables from a VEP-annotated cfSNV VCF and concordance data."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


CLASS_MAP = {
    "missense_variant": "Missense_Mutation",
    "stop_gained": "Nonsense_Mutation",
    "frameshift_variant": "Frame_Shift_Ins",
    "splice_acceptor_variant": "Splice_Site",
    "splice_donor_variant": "Splice_Site",
    "inframe_deletion": "In_Frame_Del",
    "inframe_insertion": "In_Frame_Ins",
}


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def parse_info(info: str) -> dict[str, str]:
    values = {}
    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            values[key] = value
    return values


def parse_csq_header(line: str) -> list[str]:
    marker = "Format: "
    start = line.index(marker) + len(marker)
    end = line.rindex('"')
    return line[start:end].split("|")


def parse_sample_depths(fmt: str, sample: str) -> tuple[int, int, int]:
    mapping = dict(zip(fmt.split(":"), sample.split(":")))
    dp = int(mapping.get("DP", "0") or 0)
    ad = mapping.get("AD", "")
    if not ad:
        return dp, 0, 0
    counts = [int(value or 0) for value in ad.split(",")]
    ref = counts[0]
    alt = counts[1] if len(counts) > 1 else 0
    return dp, ref, alt


def load_comparison(path: Path) -> dict[tuple[str, str, str, str], dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {(row["CHROM"], row["POS"], row["REF"], row["ALT"]): row for row in reader}


def load_sample_umi_qc(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        row = next(reader, {})
    return row


def build_rows(annotated_vcf: Path, comparison: dict[tuple[str, str, str, str], dict[str, str]], sample_qc: dict[str, str], min_family_support: int) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    sensitivity_rows: list[dict[str, str]] = []
    confidence_rows: list[dict[str, str]] = []
    csq_fields: list[str] = []
    tumor_sample = "PLASMA"
    normal_sample = "WBC"
    with open_text(annotated_vcf) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                csq_fields = parse_csq_header(line.rstrip("\n"))
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                tumor_sample = header[9]
                normal_sample = header[10] if len(header) > 10 else "WBC"
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, _vid, ref, alt, _qual, filt, info, fmt = fields[:9]
            if filt != "PASS":
                continue
            tumor_value = fields[9]
            normal_value = fields[10] if len(fields) > 10 else ""
            info_map = parse_info(info)
            csq_entry = info_map.get("CSQ", "").split(",")[0]
            csq = dict(zip(csq_fields, csq_entry.split("|"))) if csq_fields and csq_entry else {}
            tumor_dp, tumor_ref, tumor_alt = parse_sample_depths(fmt, tumor_value)
            normal_dp, normal_ref, normal_alt = parse_sample_depths(fmt, normal_value) if normal_value else (0, 0, int(info_map.get("WBC_ALT_COUNT", "0") or 0))
            key = (chrom, pos, ref, alt)
            concordance = comparison.get(key, {})
            umi_family_count = int(concordance.get("umi_family_count") or info_map.get("UMI_FAMILY_COUNT", "0") or 0)
            umi_max_family_size = concordance.get("umi_max_family_size") or info_map.get("UMI_MAX_FAMILY_SIZE", "0")
            consensus_t_depth = concordance.get("consensus_t_depth", "0")
            consensus_alt_count = concordance.get("consensus_alt_count", "0")
            mutect2_support = concordance.get("mutect2", "0")
            tumor_vaf = tumor_alt / (tumor_ref + tumor_alt) if (tumor_ref + tumor_alt) else 0.0
            row = {
                "Hugo_Symbol": csq.get("SYMBOL", "NA"),
                "Chromosome": chrom.replace("chr", ""),
                "Start_Position": pos,
                "End_Position": str(int(pos) + len(ref) - 1),
                "Reference_Allele": ref,
                "Tumor_Seq_Allele2": alt,
                "Tumor_Sample_Barcode": tumor_sample,
                "Matched_Norm_Sample_Barcode": normal_sample,
                "Variant_Classification": CLASS_MAP.get(csq.get("Consequence", "").split("&")[0], "Targeted_Region"),
                "Variant_Type": "SNP" if len(ref) == 1 and len(alt) == 1 else "INDEL",
                "HGVSc": csq.get("HGVSc", ""),
                "HGVSp": csq.get("HGVSp", ""),
                "t_depth": str(tumor_dp),
                "t_ref_count": str(tumor_ref),
                "t_alt_count": str(tumor_alt),
                "n_depth": str(normal_dp),
                "n_ref_count": str(normal_ref),
                "n_alt_count": str(normal_alt),
                "tumor_vaf": f"{tumor_vaf:.6f}",
                "primary_caller": "cfSNV",
                "mutect2_support": mutect2_support,
                "call_tier": "high_sensitivity",
                "consensus_t_depth": consensus_t_depth,
                "consensus_alt_count": consensus_alt_count,
                "umi_family_count": str(umi_family_count or int(sample_qc.get("umi_family_count", "0") or 0)),
                "umi_max_family_size": str(umi_max_family_size or sample_qc.get("umi_max_family_size", "0")),
                "wbc_alt_count": concordance.get("wbc_alt_count", str(normal_alt)),
            }
            sensitivity_rows.append(row)
            if mutect2_support == "1" or umi_family_count >= min_family_support:
                confident = dict(row)
                confident["call_tier"] = "high_confidence"
                confidence_rows.append(confident)
    return sensitivity_rows, confidence_rows


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
        "Tumor_Sample_Barcode",
        "Matched_Norm_Sample_Barcode",
        "Variant_Classification",
        "Variant_Type",
        "HGVSc",
        "HGVSp",
        "t_depth",
        "t_ref_count",
        "t_alt_count",
        "n_depth",
        "n_ref_count",
        "n_alt_count",
        "tumor_vaf",
        "primary_caller",
        "mutect2_support",
        "call_tier",
        "consensus_t_depth",
        "consensus_alt_count",
        "umi_family_count",
        "umi_max_family_size",
        "wbc_alt_count",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--comparison", required=True, type=Path)
    parser.add_argument("--umi-qc", required=True, type=Path)
    parser.add_argument("--high-sensitivity-output", required=True, type=Path)
    parser.add_argument("--high-confidence-output", required=True, type=Path)
    parser.add_argument("--min-family-support", required=True, type=int)
    args = parser.parse_args()

    comparison = load_comparison(args.comparison)
    sample_qc = load_sample_umi_qc(args.umi_qc)
    sensitivity_rows, confidence_rows = build_rows(args.input, comparison, sample_qc, args.min_family_support)
    write_rows(args.high_sensitivity_output, sensitivity_rows)
    write_rows(args.high_confidence_output, confidence_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
