"""Stream-split a FASTQ into two pseudo-replicates by record parity.

Used by from-FASTQ mode when a condition has only one biological sample
on the RNA-seq or Ribo-seq side: pyDESeq2 needs at least two replicates
per condition to estimate dispersion, and the existing delta-TE math
falls back to a "single_replicate_no_statistics" note when there is no
internal Ribo log2FC available.

Splitting by record parity (even-indexed records → rep1, odd-indexed →
rep2) keeps the position / fragment-length distribution identical
across the two halves, which is what we want when the goal is "let the
math run" rather than "estimate biological variance" — the latter
genuinely needs more libraries, not more software.

Both pseudo-replicates are obviously NOT biological replicates: they
share the exact same read population minus the order. DESeq2's
dispersion estimates will be artificially low. The orchestrator emits
a stderr WARNING every time a split happens so the user is never
surprised.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import IO

from .fastq_pairing import FastqSample


def _open_text_in(path: Path) -> IO[str]:
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    return open(p, "rt")


def _open_text_out(path: Path) -> IO[str]:
    path.parent.mkdir(parents=True, exist_ok=True)
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "wt")
    return open(p, "wt")


def stream_split_by_parity(src: Path, even_dst: Path, odd_dst: Path) -> int:
    """Split a FASTQ into two halves by record parity.

    Even-indexed records (0, 2, 4, ...) go to ``even_dst``; odd-indexed
    records go to ``odd_dst``. Returns the total number of records
    written across both halves.
    """
    src = Path(src)
    even_dst = Path(even_dst)
    odd_dst = Path(odd_dst)

    n = 0
    with _open_text_in(src) as fh, \
            _open_text_out(even_dst) as a, \
            _open_text_out(odd_dst) as b:
        record: list[str] = []
        for line in fh:
            record.append(line)
            if len(record) == 4:
                target = a if (n % 2 == 0) else b
                target.writelines(record)
                record = []
                n += 1
        if record:
            raise ValueError(
                f"Truncated FASTQ record at end of {src} "
                f"(got {len(record)} lines after the last complete record)"
            )
    return n


def split_sample_into_pseudo_replicates(
    sample: FastqSample, workdir: Path
) -> tuple[FastqSample, FastqSample]:
    """Stream-split ``sample`` into two pseudo-replicates by record parity.

    For paired-end samples both R1 and R2 are split; record N's R1
    always lands on the same parity as record N's R2 so the pairs stay
    matched.

    Output FASTQs are written under ``workdir/<rep_name>/`` preserving
    the original filename so users can recognise the source.
    """
    workdir = Path(workdir)
    rep1_name = f"{sample.sample}_rep1"
    rep2_name = f"{sample.sample}_rep2"
    rep1_dir = workdir / rep1_name
    rep2_dir = workdir / rep2_name

    rep1_r1 = rep1_dir / sample.r1.name
    rep2_r1 = rep2_dir / sample.r1.name
    stream_split_by_parity(sample.r1, rep1_r1, rep2_r1)

    rep1_r2: Path | None = None
    rep2_r2: Path | None = None
    if sample.r2 is not None:
        rep1_r2 = rep1_dir / sample.r2.name
        rep2_r2 = rep2_dir / sample.r2.name
        stream_split_by_parity(sample.r2, rep1_r2, rep2_r2)

    return (
        FastqSample(sample=rep1_name, r1=rep1_r1, r2=rep1_r2),
        FastqSample(sample=rep2_name, r1=rep2_r1, r2=rep2_r2),
    )
