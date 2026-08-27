#!/usr/bin/env python3
"""Generate a container manifest TSV from containers.yaml."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import yaml


def split_image_reference(image: str) -> tuple[str, str]:
    if "://" in image:
        reference = image.split("://", 1)[1]
    else:
        reference = image
    if ":" in reference:
        return tuple(reference.rsplit(":", 1))
    return reference, "NA"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    with args.config.open(encoding="utf-8") as handle:
        document = yaml.safe_load(handle)

    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["tool", "container_path", "source_registry", "tag", "digest", "build_date"])
        for tool, entry in document["tools"].items():
            image = entry["image"]
            source_registry, tag = split_image_reference(image)
            writer.writerow([tool, entry["sif"], source_registry, tag, "NA", document["metadata"]["pinned_on"]])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
