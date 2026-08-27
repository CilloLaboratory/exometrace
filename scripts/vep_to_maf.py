#!/usr/bin/env python3
"""Convert a VEP-annotated somatic VCF into a simple MAF-like table."""

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
    marker = 'Format: '
    start = line.index(marker) + len(marker)
    end = line.rindex('"')
    return line[start:end].split("|")


def parse_sample_ad(sample_field: str, sample_value: str) -> tuple[int, int, int]:
    fmt = sample_field.split(":")
    values = sample_value.split(":")
    mapping = dict(zip(fmt, values))
    dp = int(mapping.get("DP", "0") or 0)
    ad = mapping.get("AD", "")
    if ad:
        counts = [int(value or 0) for value in ad.split(",")]
        ref = counts[0]
        alt = counts[1] if len(counts) > 1 else 0
    else:
        ref = 0
        alt = 0
    return dp, ref, alt


def parse_normal_counts_from_tumor_field(sample_field: str, sample_value: str) -> tuple[int, int, int]:
    fmt = sample_field.split(":")
    values = sample_value.split(":")
    mapping = dict(zip(fmt, values))
    dp = int(mapping.get("NDP", "0") or 0)
    ad = mapping.get("NAD", "")
    if ad:
        counts = [int(value or 0) for value in ad.split(",")]
        ref = counts[0]
        alt = counts[1] if len(counts) > 1 else 0
    else:
        ref = 0
        alt = 0
    return dp, ref, alt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--comparison", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--caller", default="deepsomatic")
    args = parser.parse_args()

    comparison_support: dict[tuple[str, str, str, str], str] = {}
    with args.comparison.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            comparison_support[(row["CHROM"], row["POS"], row["REF"], row["ALT"])] = row["mutect2"]

    maf_rows: list[dict[str, str]] = []
    csq_fields: list[str] = []
    tumor_sample = "TUMOR"
    normal_sample = "NORMAL"
    has_normal_sample_column = False
    with open_text(args.input) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                csq_fields = parse_csq_header(line.rstrip("\n"))
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                if len(header) < 10:
                    raise ValueError("VCF header is missing sample columns")
                tumor_sample = header[9]
                if len(header) > 10:
                    normal_sample = header[10]
                    has_normal_sample_column = True
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, _vid, ref, alt, _qual, filt, info, fmt = fields[:9]
            tumor_val = fields[9]
            normal_val = fields[10] if len(fields) > 10 else ""
            if filt != "PASS":
                continue
            info_map = parse_info(info)
            csq_entries = info_map.get("CSQ", "").split(",")
            csq = dict(zip(csq_fields, csq_entries[0].split("|"))) if csq_fields and csq_entries and csq_entries[0] else {}
            consequence = csq.get("Consequence", "Sequence_Ontology_variant").split("&")[0]
            tumor_dp, tumor_ref, tumor_alt = parse_sample_ad(fmt, tumor_val)
            if has_normal_sample_column:
                normal_dp, normal_ref, normal_alt = parse_sample_ad(fmt, normal_val)
            else:
                normal_dp, normal_ref, normal_alt = parse_normal_counts_from_tumor_field(fmt, tumor_val)
            tumor_vaf = tumor_alt / (tumor_ref + tumor_alt) if (tumor_ref + tumor_alt) else 0.0

            maf_rows.append(
                {
                    "Hugo_Symbol": csq.get("SYMBOL", "NA"),
                    "Chromosome": chrom.replace("chr", ""),
                    "Start_Position": pos,
                    "End_Position": str(int(pos) + len(ref) - 1),
                    "Reference_Allele": ref,
                    "Tumor_Seq_Allele2": alt,
                    "Tumor_Sample_Barcode": tumor_sample,
                    "Matched_Norm_Sample_Barcode": normal_sample,
                    "Variant_Classification": CLASS_MAP.get(consequence, "Targeted_Region"),
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
                    "caller": args.caller,
                    "mutect2_support": comparison_support.get((chrom, pos, ref, alt), "0"),
                    "caller_count": "2" if comparison_support.get((chrom, pos, ref, alt), "0") == "1" else "1",
                }
            )

    fieldnames = [
        "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele",
        "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
        "Variant_Classification", "Variant_Type", "HGVSc", "HGVSp", "t_depth", "t_ref_count",
        "t_alt_count", "n_depth", "n_ref_count", "n_alt_count", "tumor_vaf", "caller",
        "mutect2_support", "caller_count",
    ]
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(maf_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
