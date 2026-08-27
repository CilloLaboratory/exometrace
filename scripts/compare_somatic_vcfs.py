#!/usr/bin/env python3
"""Compare DeepSomatic and Mutect2 somatic callsets."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def load_variants(path: Path, caller_name: str) -> dict[tuple[str, str, str, str], dict[str, str]]:
    variants: dict[tuple[str, str, str, str], dict[str, str]] = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, pos, _vid, ref, alt, _qual, filt, _info, *_rest = line.rstrip("\n").split("\t")
            key = (chrom, pos, ref, alt)
            variants[key] = {"caller": caller_name, "filter": filt}
    return variants


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deepsomatic", required=True, type=Path)
    parser.add_argument("--mutect2", required=True, type=Path)
    parser.add_argument("--output-table", required=True, type=Path)
    parser.add_argument("--output-counts", required=True, type=Path)
    args = parser.parse_args()

    ds = load_variants(args.deepsomatic, "deepsomatic")
    m2 = load_variants(args.mutect2, "mutect2")

    all_keys = sorted(set(ds) | set(m2), key=lambda item: (item[0], int(item[1]), item[2], item[3]))

    with args.output_table.open("w", encoding="utf-8") as handle:
        handle.write("CHROM\tPOS\tREF\tALT\tdeepsomatic\tmutect2\tdeepsomatic_filter\tmutect2_filter\n")
        for chrom, pos, ref, alt in all_keys:
            ds_call = "1" if (chrom, pos, ref, alt) in ds else "0"
            m2_call = "1" if (chrom, pos, ref, alt) in m2 else "0"
            ds_filter = ds.get((chrom, pos, ref, alt), {}).get("filter", ".")
            m2_filter = m2.get((chrom, pos, ref, alt), {}).get("filter", ".")
            handle.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{ds_call}\t{m2_call}\t{ds_filter}\t{m2_filter}\n")

    ds_pass = sum(1 for value in ds.values() if value["filter"] == "PASS")
    m2_pass = sum(1 for value in m2.values() if value["filter"] == "PASS")
    intersection = sum(1 for key in set(ds) & set(m2) if ds[key]["filter"] == "PASS" and m2[key]["filter"] == "PASS")
    ds_only = sum(1 for key in set(ds) - set(m2) if ds[key]["filter"] == "PASS")
    m2_only = sum(1 for key in set(m2) - set(ds) if m2[key]["filter"] == "PASS")

    with args.output_counts.open("w", encoding="utf-8") as handle:
        handle.write("deepsomatic_pass\tmutect2_pass\tintersection\tdeepsomatic_only\tmutect2_only\n")
        handle.write(f"{ds_pass}\t{m2_pass}\t{intersection}\t{ds_only}\t{m2_only}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
