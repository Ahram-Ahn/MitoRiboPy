"""Deduplication step for the ``mitoribopy align`` pipeline.

Pipeline contract
-----------------

UMI handling is a two-step process spread across the align pipeline:

1. **UMI extraction** happens BEFORE alignment in
   :mod:`mitoribopy.align.trim`. The cutadapt invocation removes the
   UMI bases from the biological insert and writes the UMI into the
   read name (``QNAME_UMI`` separator). This is mandatory: without
   extraction, the UMI bases would be aligned as if they were
   biological sequence and the per-read coordinate would be wrong.

2. **UMI deduplication** happens AFTER alignment, here. Collapsing PCR
   duplicates requires both the UMI sequence AND the aligned
   coordinate, because two reads with the same UMI but different
   alignment coordinates are genuinely distinct molecules. Coordinate-
   only deduplication is biologically wrong for ribo-seq: ribosomes
   translating the same codon naturally produce reads at the same
   transcript coordinate, and there is no way to distinguish those
   genuine biological duplicates from PCR artefacts without UMI
   information. The TACO1-KO regression in
   ``docs/validation/taco1_ko_regression.md`` is the empirical evidence
   that motivated removal of the legacy ``mark-duplicates`` strategy —
   coordinate-only dedup flattens the polyproline stall signal to
   baseline.

Strategies
----------

* ``umi-tools``  UMI-aware coordinate+UMI collapse. The only correct
                 choice when UMIs are present. Uses ``--method=unique``
                 by default so reads collapse only on exact coordinate
                 AND UMI match — directional clustering over-collapses
                 in the low-complexity mt-mRNA regime.

* ``skip``       no deduplication. The required default for UMI-less
                 mt-Ribo-seq libraries.

* ``auto``       resolves to ``umi-tools`` if ``umi_length > 0`` else
                 ``skip``. Default.

QC output
---------

After every align run :func:`write_umi_qc_tsv` aggregates the per-
sample :class:`DedupResult`s into ``umi_qc.tsv`` so a reviewer can see
at a glance how many reads each sample had pre-/post-dedup and how
much PCR duplication was inferred.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pysam

from ._types import DedupResult, DedupStrategy
from .bam_utils import count_mapped_reads as _count_mapped_reads
from .tool_check import ToolNotFoundError


# UMI-tools method recommendation thresholds. These pick a default
# method based on UMI length and observed duplication, mirroring the
# guidance in the publication-readiness review:
#   * UMI length >= 8 nt + duplication > 0.20  -> 'directional'
#   * UMI length < 8 nt OR low duplication     -> 'unique'
_UMI_METHOD_LENGTH_THRESHOLD = 8
_UMI_METHOD_DUPLICATION_THRESHOLD = 0.20


def recommend_umi_method(umi_length: int, duplicate_fraction: float) -> tuple[str, str]:
    """Return ``(method, warning_code)`` for an observed sample.

    A pure helper exposed so the CLI / sample-resolver can use the same
    rule when emitting per-sample warnings. ``warning_code`` is a stable
    short identifier suitable for inclusion in ``warnings.tsv`` (e.g.
    ``"short_umi_collision_risk"``); ``"none"`` signals no warning.
    """
    umi_length = int(umi_length or 0)
    if umi_length <= 0:
        return ("skip", "no_umi")
    if umi_length < _UMI_METHOD_LENGTH_THRESHOLD:
        return ("unique", "short_umi_collision_risk")
    if duplicate_fraction > _UMI_METHOD_DUPLICATION_THRESHOLD:
        return ("directional", "none")
    return ("unique", "none")


@dataclass(frozen=True)
class UmiQCRow:
    """One row of ``umi_qc.tsv``."""

    sample_id: str
    umi_present: bool
    umi_length: int
    umi_position: str | None
    n_reads_pre_dedup: int
    n_reads_post_dedup: int
    duplicate_fraction: float
    dedup_method: str
    dedup_strategy: str
    warning_code: str

    def as_row(self) -> dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "umi_present": "true" if self.umi_present else "false",
            "umi_length": int(self.umi_length),
            "umi_position": self.umi_position or "",
            "n_reads_pre_dedup": int(self.n_reads_pre_dedup),
            "n_reads_post_dedup": int(self.n_reads_post_dedup),
            "duplicate_fraction": f"{float(self.duplicate_fraction):.6g}",
            "dedup_strategy": self.dedup_strategy,
            "dedup_method": self.dedup_method,
            "warning_code": self.warning_code,
        }


_UMI_QC_COLUMNS = (
    "sample_id",
    "umi_present",
    "umi_length",
    "umi_position",
    "n_reads_pre_dedup",
    "n_reads_post_dedup",
    "duplicate_fraction",
    "dedup_strategy",
    "dedup_method",
    "warning_code",
)


def build_umi_qc_row(
    *,
    sample_id: str,
    umi_length: int,
    umi_position: str | None,
    dedup_strategy: str,
    dedup_method: str,
    pre_count: int,
    post_count: int,
) -> UmiQCRow:
    """Assemble one :class:`UmiQCRow` from observed counts.

    Returns the recommended-method warning when the chosen method is
    not the recommended one; the row is still written so a reviewer
    sees the discrepancy.
    """
    umi_present = bool(umi_length and int(umi_length) > 0)
    pre_count = int(pre_count or 0)
    post_count = int(post_count or 0)
    duplicate_fraction = (
        0.0 if pre_count <= 0 else max(0.0, 1.0 - (post_count / pre_count))
    )
    if umi_present:
        recommended, recommend_warn = recommend_umi_method(
            umi_length, duplicate_fraction
        )
        if dedup_strategy == "skip":
            warning = "umi_present_but_skipped"
        elif dedup_method != recommended and recommend_warn == "none":
            warning = "method_off_recommended"
        else:
            warning = recommend_warn
    else:
        warning = "no_umi"
    return UmiQCRow(
        sample_id=str(sample_id),
        umi_present=umi_present,
        umi_length=int(umi_length or 0),
        umi_position=umi_position,
        n_reads_pre_dedup=pre_count,
        n_reads_post_dedup=post_count,
        duplicate_fraction=duplicate_fraction,
        dedup_method=str(dedup_method),
        dedup_strategy=str(dedup_strategy),
        warning_code=warning,
    )


def write_umi_qc_tsv(rows: list[UmiQCRow], path: Path) -> Path:
    """Write a publication-grade UMI QC table.

    Columns match the names used in the publication-readiness review
    (one row per sample). Always overwrites; callers are expected to
    invoke this once per align run after every sample has finished.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as h:
        h.write("\t".join(_UMI_QC_COLUMNS) + "\n")
        for row in rows:
            d = row.as_row()
            h.write("\t".join(str(d[col]) for col in _UMI_QC_COLUMNS) + "\n")
    return path


# ---------------------------------------------------------------------------
# Strategy resolution
# ---------------------------------------------------------------------------


# Legacy aliases for the canonical strategy names. Accepted on input
# (CLI / YAML / TSV sample-overrides) and rewritten to the canonical
# value as the very first step inside :func:`resolve_dedup_strategy`.
# Keep these in sync with the ``--dedup-strategy`` argparse choices
# and with :data:`mitoribopy.config.migrate.DEDUP_STRATEGY_VALUE_REWRITES`.
_DEDUP_STRATEGY_ALIASES: dict[str, str] = {
    "umi-tools": "umi_coordinate",
    "umi_tools": "umi_coordinate",
}


def canonicalize_dedup_strategy(strategy: str) -> str:
    """Return the canonical name for a dedup-strategy string.

    Lowercases the input and rewrites the legacy ``umi-tools`` /
    ``umi_tools`` aliases to ``umi_coordinate``. Returns the unmodified
    input when no alias matches; downstream validation still rejects
    truly unknown values.
    """
    if not isinstance(strategy, str):
        return strategy  # type: ignore[return-value]
    norm = strategy.strip().lower()
    return _DEDUP_STRATEGY_ALIASES.get(norm, norm)


def resolve_dedup_strategy(
    strategy: DedupStrategy,
    umi_length: int,
) -> DedupStrategy:
    """Resolve ``auto`` and validate explicit strategies against the config.

    Parameters
    ----------
    strategy:
        One of ``auto``, ``umi_coordinate`` (canonical) / ``umi-tools``
        / ``umi_tools`` (deprecated aliases), ``skip``.
    umi_length:
        The ResolvedKit umi_length. Drives the ``auto`` branch and
        validates that ``umi_coordinate`` has UMIs to work with.

    Returns
    -------
    The effective strategy, never ``auto``. The caller can then dispatch
    directly. The returned value is one of the IMPLEMENTATION-level
    strings ``umi-tools`` and ``skip`` so existing call sites keep
    working; the canonical public name (``umi_coordinate``) is used in
    config + manifest output, not in the runtime dispatcher.

    Raises
    ------
    ValueError
        when ``umi_coordinate`` is selected without UMIs.
    """
    canonical = canonicalize_dedup_strategy(strategy) if isinstance(strategy, str) else strategy

    if canonical == "auto":
        return "umi-tools" if umi_length > 0 else "skip"

    if canonical == "umi_coordinate":
        if umi_length <= 0:
            raise ValueError(
                "--dedup-strategy umi_coordinate (or legacy umi-tools) "
                "requires --umi-length > 0; cutadapt must have "
                "extracted UMIs into the read name. Either set a kit "
                "preset / --umi-length > 0, or pick "
                "--dedup-strategy skip."
            )
        return "umi-tools"

    if canonical == "skip":
        return "skip"

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
    resolved = resolve_dedup_strategy(strategy, umi_length)

    if resolved == "skip":
        return skip_dedup(
            bam_in=bam_in,
            bam_out=bam_out,
            counter=counter,
            indexer=indexer,
        )

    if log_path is None:
        log_path = Path(bam_out).with_suffix(".dedup.log")

    # resolved == "umi-tools"
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
