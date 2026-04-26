#!/usr/bin/env python3
"""Create a reproducible random subsample of a FASTQ file.

Mirrors the BED subsampler in :mod:`mitoribopy.tools.subsample`: reads
the input as 4-line records (header / seq / plus / quality), reservoir
samples N records via Algorithm R, and writes them to ``--output``.
Both ``.fastq`` and ``.fastq.gz`` (gzip auto-detected by suffix) are
supported on either end.
"""

from __future__ import annotations

import argparse
import gzip
import random
from pathlib import Path

from ..console import log_info


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Randomly subsample reads from a FASTQ file (gzip auto-detected)."
    )
    parser.add_argument("--input", required=True, help="Input FASTQ (.fastq or .fastq.gz).")
    parser.add_argument("--output", required=True, help="Output FASTQ (.fastq or .fastq.gz).")
    parser.add_argument(
        "--n",
        type=int,
        required=True,
        help="Number of reads (4-line records) to keep.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility.",
    )
    return parser.parse_args()


def _open_text(path: Path, mode: str):
    """Open ``path`` in text mode, transparently handling .gz suffix."""
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def reservoir_sample_fastq(path: Path, n: int, seed: int) -> list[tuple[str, str, str, str]]:
    """Reservoir-sample exactly ``n`` 4-line FASTQ records.

    The first record's header line must start with ``@``; otherwise we
    raise so a misnamed input is caught immediately rather than
    producing a silently truncated output.
    """
    if n <= 0:
        raise ValueError("--n must be > 0")

    random.seed(seed)
    reservoir: list[tuple[str, str, str, str]] = []

    with _open_text(path, "r") as fh:
        idx = 0
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not (seq and plus and qual):
                raise ValueError(
                    f"Truncated FASTQ record at index {idx} in {path}"
                )
            if idx == 0 and not header.startswith("@"):
                raise ValueError(
                    f"{path}: first record header does not start with '@'; "
                    "is this really a FASTQ file?"
                )
            record = (header, seq, plus, qual)
            if idx < n:
                reservoir.append(record)
            else:
                j = random.randint(0, idx)
                if j < n:
                    reservoir[j] = record
            idx += 1

    return reservoir


def main() -> None:
    args = parse_args()
    in_path = Path(args.input)
    out_path = Path(args.output)

    if not in_path.exists():
        raise FileNotFoundError(f"Input FASTQ not found: {in_path}")

    sampled = reservoir_sample_fastq(in_path, args.n, args.seed)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with _open_text(out_path, "w") as out:
        for record in sampled:
            out.writelines(record)

    log_info("SUBSAMPLE_FASTQ", f"Input:  {in_path}")
    log_info("SUBSAMPLE_FASTQ", f"Output: {out_path}")
    log_info(
        "SUBSAMPLE_FASTQ",
        f"Requested n={args.n}, wrote n={len(sampled)}, seed={args.seed}.",
    )


if __name__ == "__main__":
    main()
