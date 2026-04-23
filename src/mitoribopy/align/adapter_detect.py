"""Adapter auto-detection from the head of a FASTQ.

The single failure mode this module exists to catch: a user supplies the
wrong ``--kit-preset``, cutadapt finds the named adapter in ~0% of reads,
and the pipeline then silently filters the untrimmed 150 nt reads as
"too long" without raising any error. The user sees ``post_trim`` drop by
99%+ in the read-counts table only after the run finishes.

This module reads the first ``n_reads`` records of a FASTQ, looks for the
leading window of every known kit's adapter sequence, and reports per-kit
match rates. The orchestrator uses the result to either auto-select the
right preset (when the user left ``--kit-preset custom``) or hard-fail
with a message that names the detected preset and the per-kit rates.

Detection is deliberately simple: a fixed-length substring search of the
adapter prefix. Adapter trimmers like cutadapt do something far more
sophisticated; we are not trying to replace them, only to sanity-check
the user's preset choice before the pipeline burns CPU on data that
cannot survive trimming.
"""

from __future__ import annotations

import gzip
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import IO

from ._types import KIT_PRESETS


@dataclass(frozen=True)
class DetectionResult:
    """Outcome of an adapter scan over the head of a FASTQ.

    ``preset_name`` names the best-matching preset, or ``None`` when no
    preset's leading adapter window appears in at least
    ``min_match_rate`` of the scanned reads. ``per_kit_rates`` exposes
    the full per-preset match rate for diagnostics. ``ambiguous`` is
    true when a second preset's rate is within ``ambiguity_window`` of
    the winner, which happens for kits that share an adapter sequence
    (notably ``nebnext_smallrna`` vs ``nebnext_ultra_umi``).
    """

    preset_name: str | None
    match_rate: float
    per_kit_rates: dict[str, float]
    n_reads_scanned: int
    ambiguous: bool


def _open_fastq(path: Path) -> IO[str]:
    """Open a FASTQ file (gzipped or plain) in text mode."""
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    return open(p, "rt")


def _iter_sequences(fh: IO[str], n_reads: int) -> Iterator[str]:
    """Yield up to ``n_reads`` sequence lines from a FASTQ-formatted stream."""
    count = 0
    while count < n_reads:
        header = fh.readline()
        if not header:
            return
        seq = fh.readline()
        if not seq:
            return
        # Skip the +/quality lines without validating them; we only need
        # the sequence and a rough record boundary.
        fh.readline()
        fh.readline()
        yield seq.rstrip("\n").rstrip("\r")
        count += 1


def detect_adapter(
    fastq_path: Path,
    *,
    n_reads: int = 5000,
    min_match_len: int = 12,
    min_match_rate: float = 0.30,
    ambiguity_window: float = 0.05,
    opener: Callable[[Path], IO[str]] = _open_fastq,
) -> DetectionResult:
    """Scan the head of ``fastq_path`` for known adapter signatures.

    For every preset in :data:`KIT_PRESETS` whose ``adapter`` is not
    ``None``, count how many of the first ``n_reads`` reads contain the
    adapter's leading ``min_match_len`` nt as a substring. The preset
    with the highest match rate wins; ties are broken alphabetically so
    the result is deterministic across runs and platforms.

    Parameters
    ----------
    fastq_path:
        Path to the FASTQ to scan. ``.gz`` is auto-detected.
    n_reads:
        Maximum number of reads to inspect. 5000 is enough to drive the
        match rate to ~0% or ~100% on every real library we have seen
        and keeps the scan well under one second on typical hardware.
    min_match_len:
        Length of the adapter prefix used as the search needle. 12 nt is
        long enough that random matches are negligible (4^12 ~= 1.7e7)
        while short enough to survive an occasional sequencing error
        in the adapter region.
    min_match_rate:
        Minimum fraction of reads that must contain the winning needle
        for the result to be considered a positive call. Below this
        threshold ``preset_name`` is ``None``.
    ambiguity_window:
        Two kits whose match rates are within this window are reported
        as ``ambiguous``. The default of 5 percentage points cleanly
        separates the NEB family pair (which always tie at ~100%) from
        unrelated kits.
    opener:
        Injection point for tests. Must return an object that supports
        the file protocol (``readline`` + context-manager) for text
        mode. Defaults to a gzip-aware opener over real paths.
    """
    candidates = {
        name: preset.adapter[:min_match_len]
        for name, preset in KIT_PRESETS.items()
        if preset.adapter is not None and len(preset.adapter) >= min_match_len
    }
    counts = {name: 0 for name in candidates}
    n_scanned = 0

    with opener(fastq_path) as fh:
        for seq in _iter_sequences(fh, n_reads):
            n_scanned += 1
            for name, needle in candidates.items():
                if needle in seq:
                    counts[name] += 1

    rates = {
        name: (counts[name] / n_scanned if n_scanned else 0.0)
        for name in candidates
    }

    if n_scanned == 0 or not rates:
        return DetectionResult(
            preset_name=None,
            match_rate=0.0,
            per_kit_rates=rates,
            n_reads_scanned=n_scanned,
            ambiguous=False,
        )

    best_name = min(rates.keys(), key=lambda name: (-rates[name], name))
    best_rate = rates[best_name]

    if best_rate < min_match_rate:
        return DetectionResult(
            preset_name=None,
            match_rate=best_rate,
            per_kit_rates=rates,
            n_reads_scanned=n_scanned,
            ambiguous=False,
        )

    ambiguous = any(
        other != best_name and abs(rates[other] - best_rate) <= ambiguity_window
        for other in rates
    )

    return DetectionResult(
        preset_name=best_name,
        match_rate=best_rate,
        per_kit_rates=rates,
        n_reads_scanned=n_scanned,
        ambiguous=ambiguous,
    )


def format_per_kit_rates(rates: dict[str, float]) -> str:
    """Format per-kit match rates for human-readable error / log lines."""
    return ", ".join(
        f"{name}={rate * 100:.1f}%" for name, rate in sorted(rates.items())
    )
