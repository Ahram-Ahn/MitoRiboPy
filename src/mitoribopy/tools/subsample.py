#!/usr/bin/env python3
"""Create a reproducible random subsample of a BED file."""

from __future__ import annotations

import argparse
import random
from pathlib import Path

from ..console import log_info


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Randomly subsample lines from a BED file.")
    parser.add_argument("--input", required=True, help="Input BED file path.")
    parser.add_argument("--output", required=True, help="Output BED file path.")
    parser.add_argument(
        "--n",
        type=int,
        required=True,
        help="Number of reads (lines) to keep.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility.",
    )
    return parser.parse_args()


def reservoir_sample_lines(path: Path, n: int, seed: int) -> list[str]:
    if n <= 0:
        raise ValueError("--n must be > 0")

    random.seed(seed)
    reservoir: list[str] = []

    with path.open("r", encoding="utf-8") as fh:
        for idx, line in enumerate(fh):
            if idx < n:
                reservoir.append(line)
                continue

            j = random.randint(0, idx)
            if j < n:
                reservoir[j] = line

    return reservoir


def main() -> None:
    args = parse_args()
    in_path = Path(args.input)
    out_path = Path(args.output)

    if not in_path.exists():
        raise FileNotFoundError(f"Input BED not found: {in_path}")

    sampled = reservoir_sample_lines(in_path, args.n, args.seed)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out:
        out.writelines(sampled)

    log_info("SUBSAMPLE", f"Input:  {in_path}")
    log_info("SUBSAMPLE", f"Output: {out_path}")
    log_info(
        "SUBSAMPLE",
        f"Requested n={args.n}, wrote n={len(sampled)}, seed={args.seed}.",
    )


if __name__ == "__main__":
    main()

