#!/usr/bin/env python3
"""Read a scalar value from a YAML file using dotted-key syntax."""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml


def resolve_value(config_path: Path, value):
    if isinstance(value, str):
        candidate = Path(value)
        if candidate.is_absolute():
            return str(candidate)
        base_dir = config_path.parent.parent if config_path.parent.name == "config" else config_path.parent
        if "/" in value or value.endswith((".fa", ".fai", ".dict", ".bed", ".gz", ".tsv", ".vcf", ".tbi")):
            return str((base_dir / candidate).resolve())
    return value


def get_value(document: dict, key: str):
    value = document
    for part in key.split("."):
        value = value[part]
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--key", required=True)
    args = parser.parse_args()

    with args.config.open(encoding="utf-8") as handle:
        document = yaml.safe_load(handle)

    value = get_value(document, args.key)
    print(resolve_value(args.config, value))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
