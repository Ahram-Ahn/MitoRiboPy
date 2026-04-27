"""Per-sample UMI auto-detection from FASTQ entropy.

A UMI region is, by design, a uniform random N-mer at one end of every
read. Per-position Shannon entropy across many reads is therefore close
to 2.0 bits inside a UMI and well below 2.0 bits in the biological
insert (which inherits codon / 5'UTR composition biases).

We scan up to ``n_reads`` reads, compute per-position entropy on the
first 16 nt and last 16 nt, and call the longest run of high-entropy
positions at either end as the UMI. The run length is snapped DOWN to
the nearest of ``{4, 6, 8, 10, 12}`` because real UMI kits use those
exact lengths and snapping conservatively avoids over-trimming.

This module is intentionally conservative: when in doubt, return
``length=0``. A wrong UMI length silently corrupts dedup; a missed UMI
just means dedup is skipped, which is recoverable.
"""

from __future__ import annotations

import gzip
import math
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Literal


_WINDOW = 16
_BASES = ("A", "C", "G", "T")
_MAX_ENTROPY = 2.0
_FLAT_THRESHOLD = 1.85
_SNAP_LENGTHS = (12, 10, 8, 6, 4)


@dataclass(frozen=True)
class UmiDetectionResult:
    """Result of an entropy-based UMI scan over a FASTQ.

    ``length`` is one of 0, 4, 6, 8, 10, 12 (zero means "no UMI
    detected"). ``position`` indicates which end the UMI sits at; it
    is meaningless when ``length == 0`` and defaults to ``"5p"`` then.
    ``entropy_5p`` / ``entropy_3p`` are per-position Shannon entropies
    (bits) across the first / last :data:`_WINDOW` positions, exposed
    for diagnostics and tests.
    """

    length: int
    position: Literal["5p", "3p"]
    entropy_5p: tuple[float, ...]
    entropy_3p: tuple[float, ...]
    n_reads_scanned: int


def _open_fastq(path: Path) -> IO[str]:
    if str(path).endswith(".gz"):
        return gzip.open(str(path), "rt")
    return open(str(path), "rt")


def _iter_sequences(fh: IO[str], n_reads: int) -> Iterator[str]:
    count = 0
    while count < n_reads:
        header = fh.readline()
        if not header:
            return
        seq = fh.readline()
        if not seq:
            return
        fh.readline()
        fh.readline()
        yield seq.rstrip("\n").rstrip("\r")
        count += 1


def _shannon(counts: list[int]) -> float:
    total = sum(counts)
    if total == 0:
        return 0.0
    h = 0.0
    for c in counts:
        if c == 0:
            continue
        p = c / total
        h -= p * math.log2(p)
    return h


def _flat_run_length(entropy: tuple[float, ...]) -> int:
    """Return the longest contiguous run of positions with entropy >= threshold."""
    best = 0
    current = 0
    for e in entropy:
        if e >= _FLAT_THRESHOLD:
            current += 1
            if current > best:
                best = current
        else:
            current = 0
    return best


def _snap_down(run_length: int) -> int:
    for snap in _SNAP_LENGTHS:
        if run_length >= snap:
            return snap
    return 0


def detect_umi(fastq: Path, *, n_reads: int = 5000) -> UmiDetectionResult:
    """Estimate UMI length and position from per-position entropy.

    Reads up to ``n_reads`` records (gzip-aware via ``.gz`` suffix).
    Skips reads shorter than :data:`_WINDOW` nt so a stray short read
    cannot disturb the entropy buckets.

    Returns ``length=0`` whenever:

    * fewer than 4 contiguous high-entropy positions are found at either
      end (no detectable UMI), or
    * the FASTQ is empty / only contained sub-window reads.
    """
    counts_5p: list[list[int]] = [[0, 0, 0, 0] for _ in range(_WINDOW)]
    counts_3p: list[list[int]] = [[0, 0, 0, 0] for _ in range(_WINDOW)]
    base_index = {b: i for i, b in enumerate(_BASES)}
    n_scanned = 0

    with _open_fastq(Path(fastq)) as handle:
        for seq in _iter_sequences(handle, n_reads):
            if len(seq) < _WINDOW:
                continue
            n_scanned += 1
            head = seq[:_WINDOW]
            tail = seq[-_WINDOW:]
            for i, b in enumerate(head):
                idx = base_index.get(b)
                if idx is not None:
                    counts_5p[i][idx] += 1
            for i, b in enumerate(tail):
                idx = base_index.get(b)
                if idx is not None:
                    counts_3p[i][idx] += 1

    entropy_5p = tuple(_shannon(c) for c in counts_5p)
    # 3' window is reported in 5'->3' order along the read tail, so
    # entropy_3p[-1] is the very last base of the read.
    entropy_3p = tuple(_shannon(c) for c in counts_3p)

    if n_scanned == 0:
        return UmiDetectionResult(
            length=0,
            position="5p",
            entropy_5p=entropy_5p,
            entropy_3p=entropy_3p,
            n_reads_scanned=0,
        )

    # 5' run: longest streak starting from position 0 of the head window.
    head_run = 0
    for e in entropy_5p:
        if e >= _FLAT_THRESHOLD:
            head_run += 1
        else:
            break
    # 3' run: longest streak ending at the last position of the tail window.
    tail_run = 0
    for e in reversed(entropy_3p):
        if e >= _FLAT_THRESHOLD:
            tail_run += 1
        else:
            break

    # If the flat run covers the ENTIRE window we cannot distinguish a
    # UMI from a uniformly random biological insert (typical of normal
    # RNA-seq / Ribo-seq libraries that span many genes and start
    # positions). Return length=0 in that case — false-negative is
    # recoverable, false-positive silently corrupts dedup or trims real
    # bases off every read.
    if head_run >= _WINDOW:
        head_run = 0
    if tail_run >= _WINDOW:
        tail_run = 0

    head_snap = _snap_down(head_run)
    tail_snap = _snap_down(tail_run)

    if head_snap == 0 and tail_snap == 0:
        return UmiDetectionResult(
            length=0,
            position="5p",
            entropy_5p=entropy_5p,
            entropy_3p=entropy_3p,
            n_reads_scanned=n_scanned,
        )

    if head_snap >= tail_snap:
        return UmiDetectionResult(
            length=head_snap,
            position="5p",
            entropy_5p=entropy_5p,
            entropy_3p=entropy_3p,
            n_reads_scanned=n_scanned,
        )
    return UmiDetectionResult(
        length=tail_snap,
        position="3p",
        entropy_5p=entropy_5p,
        entropy_3p=entropy_3p,
        n_reads_scanned=n_scanned,
    )
