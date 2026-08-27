#!/usr/bin/env python3
"""Dump selected pipeline parameters to YAML."""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--reference-config", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    payload = {
        "run_date": "2026-08-19",
        "samplesheet": args.samplesheet,
        "reference_config": args.reference_config,
        "results_dir": args.results_dir,
    }
    with args.output.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
