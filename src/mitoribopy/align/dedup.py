"""Deduplication step for the ``mitoribopy align`` pipeline (Phase 3.5).

Implements the Approach-D strategy decided in Phase 3.5:

* ``umi-tools``       UMI-aware collapse; the only correct choice when
                      UMIs are present. Uses ``--method=unique`` by
                      default so reads are collapsed only when both the
                      coordinate AND the exact UMI match - directional
                      clustering over-collapses in the low-complexity
                      mt-mRNA regime.

* ``skip``            no deduplication. The required default for
                      UMI-less mt-Ribo-seq libraries. Coordinate-only
                      dedup destroys codon-occupancy signal because
                      ribosomes translating the same codon naturally
                      produce reads at the same transcript coordinate.

* ``mark-duplicates`` picard MarkDuplicates. Dangerous on mt-Ribo-seq
                      data. Gated behind the intentionally long
                      confirmation flag
                      ``--i-understand-mark-duplicates-destroys-mt-ribo-seq-signal``
                      so nobody uses it by accident.

* ``auto``            resolves to ``umi-tools`` if ``umi_length > 0``
                      else ``skip``. This is the default.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pysam

from ._types import DedupResult, DedupStrategy
from .bam_utils import count_mapped_reads as _count_mapped_reads
from .tool_check import ToolNotFoundError


CONFIRM_MARK_DUPLICATES_FLAG = "--i-understand-mark-duplicates-destroys-mt-ribo-seq-signal"

_MARK_DUPLICATES_DANGER_MESSAGE = (
    "mark-duplicates uses coordinate-only deduplication. mt-Ribo-seq "
    "libraries are low-complexity: ribosomes translating the same codon "
    "produce reads at the same transcript coordinate. Coordinate-only "
    "dedup cannot distinguish genuine ribosomes from PCR duplicates and "
    "will flatten codon-occupancy signal (polyproline stalls, TACO1-KO "
    "phenotype). Only use when you have an independent reason to believe "
    "PCR duplication dominates your library (carrier nuclear RNA, >22 "
    "PCR cycles, very low input). Confirm by passing "
    f"{CONFIRM_MARK_DUPLICATES_FLAG} on the command line."
)


# ---------------------------------------------------------------------------
# Strategy resolution
# ---------------------------------------------------------------------------


def resolve_dedup_strategy(
    strategy: DedupStrategy,
    umi_length: int,
    confirm_mark_duplicates: bool = False,
) -> DedupStrategy:
    """Resolve ``auto`` and validate explicit strategies against the config.

    Parameters
    ----------
    strategy:
        One of ``auto``, ``umi-tools``, ``skip``, ``mark-duplicates``.
    umi_length:
        The ResolvedKit umi_length. Drives the ``auto`` branch and
        validates that ``umi-tools`` has UMIs to work with.
    confirm_mark_duplicates:
        True when the user passed the long confirmation flag.

    Returns
    -------
    The effective strategy, never ``auto``. The caller can then dispatch
    directly.

    Raises
    ------
    ValueError
        on incompatible combinations (umi-tools without UMIs; mark-
        duplicates without confirmation).
    """
    if strategy == "auto":
        return "umi-tools" if umi_length > 0 else "skip"

    if strategy == "umi-tools":
        if umi_length <= 0:
            raise ValueError(
                "--dedup-strategy umi-tools requires --umi-length > 0; "
                "cutadapt must have extracted UMIs into the read name. "
                "Either set a kit preset / --umi-length > 0, or pick "
                "--dedup-strategy skip."
            )
        return "umi-tools"

    if strategy == "skip":
        return "skip"

    if strategy == "mark-duplicates":
        if not confirm_mark_duplicates:
            raise ValueError(_MARK_DUPLICATES_DANGER_MESSAGE)
        return "mark-duplicates"

    raise ValueError(  # pragma: no cover - literal guarded
        f"Unknown --dedup-strategy: {strategy!r}"
    )


# ---------------------------------------------------------------------------
# Per-strategy implementations
# ---------------------------------------------------------------------------


def run_umi_tools_dedup(
    *,
    bam_in: Path,
    bam_out: Path,
    log_path: Path,
    method: str = "unique",
    umi_separator: str = "_",
    runner=subprocess.run,
    counter=_count_mapped_reads,
    indexer=pysam.index,
) -> DedupResult:
    """Run ``umi_tools dedup`` on a sorted, indexed BAM.

    We default to ``--method=unique`` because directional clustering
    over-collapses in low-complexity mt-mRNA regions where neighboring
    codons carry correlated UMIs by chance. Users can override via the
    CLI's ``--umi-dedup-method`` flag if they have a specific library
    quirk that justifies it.
    """
    if shutil.which("umi_tools") is None:
        raise ToolNotFoundError(
            "'umi_tools' is required for --dedup-strategy umi-tools but "
            "was not found on PATH. Install with:\n"
            "  conda install -c bioconda umi_tools"
        )
    bam_in = Path(bam_in)
    bam_out = Path(bam_out)
    log_path = Path(log_path)
    bam_out.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "umi_tools",
        "dedup",
        f"--stdin={bam_in}",
        f"--stdout={bam_out}",
        f"--method={method}",
        f"--umi-separator={umi_separator}",
        f"--log={log_path}",
    ]
    completed = runner(cmd, check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        stderr = (getattr(completed, "stderr", "") or "").strip()
        raise RuntimeError(
            f"umi_tools dedup failed with exit code {completed.returncode}: "
            f"{stderr or '<no stderr captured>'}"
        )

    indexer(str(bam_out))

    input_reads = counter(bam_in)
    output_reads = counter(bam_out)
    return DedupResult(
        strategy="umi-tools",
        input_reads=input_reads,
        output_reads=output_reads,
        bam_path=bam_out,
    )


def run_mark_duplicates(
    *,
    bam_in: Path,
    bam_out: Path,
    log_path: Path,
    runner=subprocess.run,
    counter=_count_mapped_reads,
    indexer=pysam.index,
) -> DedupResult:
    """Run picard MarkDuplicates with REMOVE_DUPLICATES=true.

    Pre-gated by :func:`resolve_dedup_strategy` requiring the
    confirmation flag. Still emits a prominent warning banner at the
    start of the picard log so any reader of the log sees the risk.
    """
    if shutil.which("picard") is None:
        raise ToolNotFoundError(
            "'picard' is required for --dedup-strategy mark-duplicates "
            "but was not found on PATH. Install with:\n"
            "  conda install -c bioconda picard"
        )
    bam_in = Path(bam_in)
    bam_out = Path(bam_out)
    log_path = Path(log_path)
    bam_out.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("w", encoding="utf-8") as logfile:
        logfile.write(
            "# WARNING: MarkDuplicates uses coordinate-only deduplication.\n"
            "# This destroys codon-occupancy signal on mt-Ribo-seq data.\n"
            "# User confirmed with "
            f"{CONFIRM_MARK_DUPLICATES_FLAG}.\n\n"
        )

    cmd = [
        "picard",
        "MarkDuplicates",
        f"I={bam_in}",
        f"O={bam_out}",
        f"M={log_path.with_suffix('.metrics.txt')}",
        "REMOVE_DUPLICATES=true",
    ]
    completed = runner(cmd, check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        stderr = (getattr(completed, "stderr", "") or "").strip()
        raise RuntimeError(
            "picard MarkDuplicates failed with exit code "
            f"{completed.returncode}: {stderr or '<no stderr captured>'}"
        )

    indexer(str(bam_out))

    input_reads = counter(bam_in)
    output_reads = counter(bam_out)
    return DedupResult(
        strategy="mark-duplicates",
        input_reads=input_reads,
        output_reads=output_reads,
        bam_path=bam_out,
    )


def skip_dedup_in_place(
    bam_path: Path,
    *,
    counter=_count_mapped_reads,
) -> DedupResult:
    """Record a 'skip' dedup result without touching the filesystem.

    The default code path in :func:`run_dedup` (when strategy='skip')
    creates a hardlink/copy of the input BAM under deduped/. That
    duplicated artefact is wasted disk for the common no-UMI case
    because the downstream :func:`mitoribopy.align.bam_utils.bam_to_bed6`
    reads BAM content, not a path. The CLI orchestrator calls THIS
    function instead, then feeds the original (mapq.bam) BAM straight
    into bam_to_bed6.

    Returns a :class:`DedupResult` with ``strategy='skip'`` and equal
    input/output counts so the read-counts table still records the
    stage clearly.
    """
    bam_path = Path(bam_path)
    count = counter(bam_path)
    return DedupResult(
        strategy="skip",
        input_reads=count,
        output_reads=count,
        bam_path=bam_path,
    )


def skip_dedup(
    *,
    bam_in: Path,
    bam_out: Path,
    counter=_count_mapped_reads,
    indexer=pysam.index,
) -> DedupResult:
    """Pass-through 'dedup': copy/symlink the input BAM to the output path.

    Used when:

    * ``--umi-length 0`` and ``--dedup-strategy auto`` (the default for
      no-UMI libraries), OR
    * the user explicitly passed ``--dedup-strategy skip``.

    Returns a :class:`DedupResult` with ``strategy='skip'`` and equal
    input/output counts so the Phase H read-counts table records the
    stage clearly.

    The CLI orchestrator now prefers :func:`skip_dedup_in_place` for the
    common no-UMI case (no on-disk duplication); this function is kept
    for callers that explicitly want a separate output path (e.g. a
    test harness that asserts the file exists).
    """
    bam_in = Path(bam_in)
    bam_out = Path(bam_out)
    bam_out.parent.mkdir(parents=True, exist_ok=True)

    if bam_out.exists() or bam_out.is_symlink():
        bam_out.unlink()
    # Hardlink first (cheapest, atomic); copy as a fallback for cross-fs.
    try:
        bam_out.hardlink_to(bam_in)
    except (AttributeError, OSError):
        import shutil as _shutil

        _shutil.copyfile(bam_in, bam_out)

    indexer(str(bam_out))

    count = counter(bam_in)
    return DedupResult(
        strategy="skip",
        input_reads=count,
        output_reads=count,
        bam_path=bam_out,
    )


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------


def run_dedup(
    *,
    strategy: DedupStrategy,
    umi_length: int,
    confirm_mark_duplicates: bool,
    bam_in: Path,
    bam_out: Path,
    log_path: Path | None = None,
    umi_tools_method: str = "unique",
    umi_separator: str = "_",
    runner=subprocess.run,
    counter=_count_mapped_reads,
    indexer=pysam.index,
) -> DedupResult:
    """Resolve the strategy and dispatch to the right implementation."""
    resolved = resolve_dedup_strategy(
        strategy, umi_length, confirm_mark_duplicates
    )

    if resolved == "skip":
        return skip_dedup(
            bam_in=bam_in,
            bam_out=bam_out,
            counter=counter,
            indexer=indexer,
        )

    if log_path is None:
        log_path = Path(bam_out).with_suffix(".dedup.log")

    if resolved == "umi-tools":
        return run_umi_tools_dedup(
            bam_in=bam_in,
            bam_out=bam_out,
            log_path=log_path,
            method=umi_tools_method,
            umi_separator=umi_separator,
            runner=runner,
            counter=counter,
            indexer=indexer,
        )

    # resolved == "mark-duplicates"
    return run_mark_duplicates(
        bam_in=bam_in,
        bam_out=bam_out,
        log_path=log_path,
        runner=runner,
        counter=counter,
        indexer=indexer,
    )
