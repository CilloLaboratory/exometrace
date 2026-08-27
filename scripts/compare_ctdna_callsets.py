#!/usr/bin/env python3
"""Compare cfSNV and Mutect2 ctDNA callsets by exact allele key."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


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


def parse_depths(fmt: str, sample: str) -> tuple[int, int]:
    fields = fmt.split(":")
    values = sample.split(":")
    mapping = dict(zip(fields, values))
    dp = int(mapping.get("DP", "0") or 0)
    ad = mapping.get("AD", "")
    if not ad:
        return dp, 0
    counts = [int(value or 0) for value in ad.split(",")]
    return dp, counts[1] if len(counts) > 1 else 0


def read_vcf(path: Path, caller: str) -> dict[tuple[str, str, str, str], dict[str, str]]:
    calls: dict[tuple[str, str, str, str], dict[str, str]] = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, _vid, ref, alt, _qual, filt, info, fmt = fields[:9]
            tumor_sample = fields[9]
            normal_sample = fields[10] if len(fields) > 10 else ""
            info_map = parse_info(info)
            depth, alt_count = parse_depths(fmt, tumor_sample)
            _normal_depth, wbc_alt_count = parse_depths(fmt, normal_sample) if normal_sample else (0, int(info_map.get("WBC_ALT_COUNT", "0") or 0))
            key = (chrom, pos, ref, alt)
            calls[key] = {
                f"{caller}_filter": filt,
                f"{caller}_depth": str(depth),
                f"{caller}_alt_count": str(alt_count),
                "wbc_alt_count": str(wbc_alt_count),
                "umi_family_count": info_map.get("UMI_FAMILY_COUNT", info_map.get("FAMILY_COUNT", "0")),
                "umi_max_family_size": info_map.get("UMI_MAX_FAMILY_SIZE", info_map.get("MAX_FAMILY_SIZE", "0")),
            }
    return calls


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cfsnv", required=True, type=Path)
    parser.add_argument("--mutect2", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    cfsnv_calls = read_vcf(args.cfsnv, "cfsnv")
    mutect_calls = read_vcf(args.mutect2, "mutect2")
    all_keys = sorted(set(cfsnv_calls) | set(mutect_calls))

    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "CHROM",
                "POS",
                "REF",
                "ALT",
                "cfsnv",
                "mutect2",
                "support_class",
                "consensus_t_depth",
                "consensus_alt_count",
                "umi_family_count",
                "umi_max_family_size",
                "wbc_alt_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for chrom, pos, ref, alt in all_keys:
            cfsnv = cfsnv_calls.get((chrom, pos, ref, alt))
            mutect = mutect_calls.get((chrom, pos, ref, alt))
            writer.writerow(
                {
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "cfsnv": "1" if cfsnv and cfsnv["cfsnv_filter"] == "PASS" else "0",
                    "mutect2": "1" if mutect and mutect["mutect2_filter"] == "PASS" else "0",
                    "support_class": "shared" if cfsnv and mutect else ("cfSNV-only" if cfsnv else "Mutect2-only"),
                    "consensus_t_depth": mutect["mutect2_depth"] if mutect else "0",
                    "consensus_alt_count": mutect["mutect2_alt_count"] if mutect else "0",
                    "umi_family_count": cfsnv["umi_family_count"] if cfsnv else (mutect["umi_family_count"] if mutect else "0"),
                    "umi_max_family_size": cfsnv["umi_max_family_size"] if cfsnv else (mutect["umi_max_family_size"] if mutect else "0"),
                    "wbc_alt_count": cfsnv["wbc_alt_count"] if cfsnv else (mutect["wbc_alt_count"] if mutect else "0"),
                }
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
