#!/usr/bin/env python3
"""Calculate variant allele fraction."""

from __future__ import annotations

import argparse


def calculate_vaf(ref_count: int, alt_count: int) -> float:
    total = ref_count + alt_count
    if total <= 0:
        raise ValueError("ref_count + alt_count must be greater than zero")
    return alt_count / total


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-count", type=int, required=True)
    parser.add_argument("--alt-count", type=int, required=True)
    args = parser.parse_args()
    print(f"{calculate_vaf(args.ref_count, args.alt_count):.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
