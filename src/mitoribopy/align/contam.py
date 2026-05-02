"""bowtie2 contaminant subtraction (rRNA / tRNA / nuclear spike-in).

Step D of the ``mitoribopy align`` pipeline. Runs bowtie2 against a
user-supplied contaminant index and emits only the UNALIGNED reads,
which are the non-contaminant fraction that proceeds to Step E
(mt-transcriptome alignment).

If no ``--contam-index`` is supplied the caller is expected to stop with
instructions to build one; this module refuses to silently skip.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path
from typing import Iterable

from ._types import ContamResult, Strandedness
from .tool_check import bowtie2_strand_flag


_INDEX_SUFFIXES = (".1.bt2", ".1.bt2l")


# bowtie2 stderr summary lines (unpaired reads).
_TOTAL_RE = re.compile(r"^\s*(\d+)\s+reads;\s+of\s+these", re.MULTILINE)
_UNALIGNED_RE = re.compile(r"^\s*(\d+)\s+\([\d.]+%\)\s+aligned\s+0\s+times", re.MULTILINE)


def _validate_index_prefix(contam_index: Path) -> None:
    """Raise with an actionable message if the bowtie2 index is not there."""
    contam_index = Path(contam_index)
    for suffix in _INDEX_SUFFIXES:
        if Path(str(contam_index) + suffix).exists():
            return
    raise FileNotFoundError(
        f"bowtie2 contaminant index not found at prefix {contam_index!s}. "
        "Build one first via:\n"
        f"  bowtie2-build contaminants.fa {contam_index}\n"
        "Typical contaminant references: ribosomal RNAs "
        "(e.g. NR_003286, NR_003287), tRNAs, spliced-leader sequences."
    )


def _build_bowtie2_contam_command(
    *,
    contam_index: Path,
    fastq_in: Path,
    fastq_out_unaligned: Path,
    strandedness: Strandedness,
    threads: int,
    seed: int,
    extra_flags: Iterable[str] = (),
) -> list[str]:
    """Assemble the bowtie2 command for contam subtraction.

    Aligned reads (contaminants) are discarded by routing the SAM stream
    to ``/dev/null``; the interesting output is ``--un-gz``, which writes
    the unaligned (non-contaminant) reads back out as gzipped FASTQ for
    Step E.
    """
    cmd: list[str] = [
        "bowtie2",
        "-x",
        str(contam_index),
        "-U",
        str(fastq_in),
        "--end-to-end",
        "--very-sensitive",
        "-L",
        "18",
        "--no-unal",
        "--un-gz",
        str(fastq_out_unaligned),
        "-p",
        str(max(int(threads), 1)),
        "--seed",
        str(seed),
    ]
    cmd += bowtie2_strand_flag(strandedness)
    cmd += list(extra_flags)
    cmd += ["-S", "/dev/null"]
    return cmd


def parse_bowtie2_stderr(stderr: str) -> tuple[int, int]:
    """Return ``(total_reads, unaligned_reads)`` from a bowtie2 stderr block.

    Raises ValueError when the expected lines are absent; this happens
    only on bowtie2 failures or unexpected versions, and we want the
    error to surface instead of defaulting to zeros.
    """
    total_match = _TOTAL_RE.search(stderr)
    unaligned_match = _UNALIGNED_RE.search(stderr)
    if total_match is None or unaligned_match is None:
        raise ValueError(
            "Could not parse bowtie2 stderr summary. "
            "Expected '<N> reads; of these' and '<M> (...) aligned 0 times' "
            f"lines. Got:\n{stderr}"
        )
    total = int(total_match.group(1))
    unaligned = int(unaligned_match.group(1))
    return total, unaligned


def subtract_contaminants(
    *,
    fastq_in: Path,
    contam_index: Path,
    fastq_out_unaligned: Path,
    strandedness: Strandedness = "forward",
    threads: int = 1,
    seed: int = 42,
    extra_flags: Iterable[str] = (),
    runner=subprocess.run,
) -> ContamResult:
    """Run bowtie2 to remove contaminant-aligning reads.

    Parameters
    ----------
    fastq_in:
        Trimmed FASTQ from Step B (cutadapt).
    contam_index:
        Path prefix of a bowtie2 index covering the contaminant panel
        (rRNA + tRNA + any nuclear / spike-in sequences the user wants
        subtracted). ``bowtie2-build`` creates six files with a shared
        prefix; pass the shared prefix here, not one of the sidecar
        files.
    fastq_out_unaligned:
        Destination for unaligned reads (these are the non-contaminant
        fraction passed to Step E).
    strandedness:
        Library strandedness; drives the bowtie2 ``--norc`` / ``--nofw``
        flag via :func:`tool_check.bowtie2_strand_flag`.
    threads:
        Passed through as ``-p``.
    seed:
        bowtie2 ``--seed`` for reproducibility (recorded in the run manifest).
    extra_flags:
        Optional bowtie2 flags appended verbatim; use for power-user
        tuning after the default recipe.
    runner:
        Injection point for tests; defaults to :func:`subprocess.run`.
    """
    contam_index = Path(contam_index)
    fastq_in = Path(fastq_in)
    fastq_out_unaligned = Path(fastq_out_unaligned)

    _validate_index_prefix(contam_index)
    fastq_out_unaligned.parent.mkdir(parents=True, exist_ok=True)

    cmd = _build_bowtie2_contam_command(
        contam_index=contam_index,
        fastq_in=fastq_in,
        fastq_out_unaligned=fastq_out_unaligned,
        strandedness=strandedness,
        threads=threads,
        seed=seed,
        extra_flags=extra_flags,
    )

    completed = runner(cmd, check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        stderr = (getattr(completed, "stderr", "") or "").strip()
        raise RuntimeError(
            "bowtie2 (contam subtraction) failed with exit code "
            f"{completed.returncode}: {stderr or '<no stderr captured>'}"
        )

    stderr = getattr(completed, "stderr", "") or ""
    total, unaligned = parse_bowtie2_stderr(stderr)
    return ContamResult(
        total_reads=total,
        aligned_to_contam=total - unaligned,
        unaligned_reads=unaligned,
    )
