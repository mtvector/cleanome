#!/usr/bin/env python3
"""Normalize whitespace in the GTF source column."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_gtf", type=Path)
    parser.add_argument("output_gtf", type=Path)
    return parser.parse_args()


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def sanitize_source(value: str) -> str:
    value = value.strip()
    if not value:
        return "unknown"
    return re.sub(r"\s+", "_", value)


def main() -> int:
    args = parse_args()
    args.output_gtf.parent.mkdir(parents=True, exist_ok=True)

    with open_text(args.input_gtf, "r") as src, open_text(args.output_gtf, "w") as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                dst.write(line)
                continue

            parts[1] = sanitize_source(parts[1])
            dst.write("\t".join(parts) + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
