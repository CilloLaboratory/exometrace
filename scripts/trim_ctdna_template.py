#!/usr/bin/env python3
"""Trim leading non-template cycles from paired ctDNA FASTQs using read structures."""

from __future__ import annotations

import argparse
import gzip
import re
from dataclasses import dataclass
from pathlib import Path


READ_STRUCTURE_PATTERN = re.compile(r"(\d+)([A-Z])")
TRIMMABLE_SEGMENTS = {"M", "S"}


@dataclass(frozen=True)
class ReadStructure:
    raw: str
    template_trim_bases: int


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def parse_read_structure(value: str) -> ReadStructure:
    if not value:
        raise ValueError("read structure must not be empty")
    if not value.endswith("+T"):
        raise ValueError(f"unsupported read structure '{value}'")

    prefix = value[:-2]
    if not prefix:
        return ReadStructure(raw=value, template_trim_bases=0)

    position = 0
    trim_bases = 0

    for match in READ_STRUCTURE_PATTERN.finditer(prefix):
        if match.start() != position:
            raise ValueError(f"unsupported read structure '{value}'")
        length = int(match.group(1))
        segment = match.group(2)
        if length <= 0 or segment not in TRIMMABLE_SEGMENTS:
            raise ValueError(f"unsupported read structure segment '{segment}' in '{value}'")
        trim_bases += length
        position = match.end()

    if position != len(prefix):
        raise ValueError(f"unsupported read structure '{value}'")
    return ReadStructure(raw=value, template_trim_bases=trim_bases)


def trim_read(header: str, sequence: str, plus: str, quality: str, trim_bases: int) -> tuple[str, str, str, str]:
    trimmed_sequence = sequence.rstrip("\n")[trim_bases:] + "\n"
    trimmed_quality = quality.rstrip("\n")[trim_bases:] + "\n"
    return header, trimmed_sequence, plus, trimmed_quality


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--read-structure-r1", required=True)
    parser.add_argument("--read-structure-r2", required=True)
    parser.add_argument("--output-r1", required=True, type=Path)
    parser.add_argument("--output-r2", required=True, type=Path)
    parser.add_argument("--qc-output", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    args = parser.parse_args()

    read_structure_r1 = parse_read_structure(args.read_structure_r1)
    read_structure_r2 = parse_read_structure(args.read_structure_r2)

    read_pairs = 0
    with open_text(args.r1, "r") as r1_in, open_text(args.r2, "r") as r2_in, open_text(args.output_r1, "w") as r1_out, open_text(args.output_r2, "w") as r2_out:
        while True:
            r1_block = [r1_in.readline() for _ in range(4)]
            r2_block = [r2_in.readline() for _ in range(4)]
            if not r1_block[0] and not r2_block[0]:
                break
            if any(not line for line in r1_block + r2_block):
                raise ValueError("encountered truncated FASTQ record while trimming ctDNA template")
            trimmed_r1 = trim_read(*r1_block, trim_bases=read_structure_r1.template_trim_bases)
            trimmed_r2 = trim_read(*r2_block, trim_bases=read_structure_r2.template_trim_bases)
            r1_out.writelines(trimmed_r1)
            r2_out.writelines(trimmed_r2)
            read_pairs += 1

    with args.qc_output.open("w", encoding="utf-8") as handle:
        handle.write(
            "sample_id\tread_pairs\tread_structure_r1\tread_structure_r2\t"
            "template_trim_bases_r1\ttemplate_trim_bases_r2\tstatus\n"
        )
        handle.write(
            f"{args.sample_id}\t{read_pairs}\t{read_structure_r1.raw}\t{read_structure_r2.raw}\t"
            f"{read_structure_r1.template_trim_bases}\t{read_structure_r2.template_trim_bases}\tPASS\n"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
