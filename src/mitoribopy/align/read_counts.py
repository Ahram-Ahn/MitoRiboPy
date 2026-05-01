"""Per-sample read-counts provenance table for ``mitoribopy align``.

Step H of the pipeline. After every sample has been trimmed, filtered,
aligned, MAPQ-filtered, and deduplicated, we assemble a
:class:`~mitoribopy.align._types.SampleCounts` from the per-stage result
dataclasses and write a tab-delimited ``read_counts.tsv`` in the run
root.

This file is the provenance spine: it documents how many reads survived
each step so downstream Phase-6 manifest consumers and human reviewers
can reconstruct where the read attrition happened.

P5.6: the invariants documented in :func:`assemble_sample_counts` are
now enforced at run time by :func:`check_count_invariants`. By default
a violation raises :class:`CountInvariantError` (mapped to
``E_OUTPUT_COUNT_INVARIANT`` in :mod:`mitoribopy.errors`); the align
CLI exposes ``--allow-count-invariant-warning`` to demote them to
warnings for legacy / debugging use.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from pathlib import Path
from typing import Iterable

from ..io.schema_versions import schema_header_line
from ._types import (
    AlignResult,
    ContamResult,
    CutadaptResult,
    DedupResult,
    SampleCounts,
)


def assemble_sample_counts(
    *,
    sample: str,
    cutadapt: CutadaptResult,
    contam: ContamResult,
    align: AlignResult,
    mapq_count: int,
    dedup: DedupResult,
) -> SampleCounts:
    """Combine per-stage results into a single :class:`SampleCounts`.

    Invariants the caller relies on:

    * ``total_reads == cutadapt.input_reads``
    * ``post_trim == cutadapt.reads_passing_filters == contam.total_reads``
    * ``rrna_aligned + post_rrna_filter == post_trim``
    * ``mt_aligned + unaligned_to_mt == post_rrna_filter``
    * ``mt_aligned_after_dedup <= mt_aligned_after_mapq <= mt_aligned``
    """
    total_reads = int(cutadapt.input_reads)
    post_trim = int(cutadapt.reads_passing_filters)
    rrna_aligned = int(contam.aligned_to_contam)
    post_rrna_filter = int(contam.unaligned_reads)
    mt_aligned = int(align.aligned)
    unaligned_to_mt = max(post_rrna_filter - mt_aligned, 0)
    mt_aligned_after_mapq = int(mapq_count)
    mt_aligned_after_dedup = int(dedup.output_reads)

    return SampleCounts(
        sample=sample,
        total_reads=total_reads,
        post_trim=post_trim,
        rrna_aligned=rrna_aligned,
        post_rrna_filter=post_rrna_filter,
        mt_aligned=mt_aligned,
        unaligned_to_mt=unaligned_to_mt,
        mt_aligned_after_mapq=mt_aligned_after_mapq,
        mt_aligned_after_dedup=mt_aligned_after_dedup,
    )


def read_counts_columns() -> tuple[str, ...]:
    """Return the TSV column order in the same sequence as the dataclass.

    Surfaced as a helper so tests (and downstream docs) have one source
    of truth for the column layout.
    """
    return tuple(f.name for f in fields(SampleCounts))


def format_row(counts: SampleCounts) -> str:
    """Format a single :class:`SampleCounts` as a tab-delimited row."""
    return "\t".join(str(getattr(counts, name)) for name in read_counts_columns())


@dataclass(frozen=True)
class CountInvariantViolation:
    """One failed read-count invariant for a single sample."""

    sample: str
    rule: str
    message: str


class CountInvariantError(ValueError):
    """Raised when one or more read-count invariants are violated."""

    def __init__(self, violations: list[CountInvariantViolation]) -> None:
        self.violations = violations
        super().__init__(
            "read_counts invariant violation(s): "
            + "; ".join(f"{v.sample}/{v.rule}: {v.message}" for v in violations)
        )


def check_count_invariants(
    rows: Iterable[SampleCounts],
) -> list[CountInvariantViolation]:
    """Return the list of invariant violations across *rows*.

    The four invariants enforced here are the ones documented at the
    top of :func:`assemble_sample_counts`. Each rule is keyed by a
    short stable name so a downstream report can index by it:

    * ``rrna_split``  — rrna_aligned + post_rrna_filter == post_trim
    * ``mt_split``    — mt_aligned + unaligned_to_mt == post_rrna_filter
    * ``mapq_le_mt``  — mt_aligned_after_mapq <= mt_aligned
    * ``dedup_le_mapq`` — mt_aligned_after_dedup <= mt_aligned_after_mapq
    """
    violations: list[CountInvariantViolation] = []
    for row in rows:
        if row.rrna_aligned + row.post_rrna_filter != row.post_trim:
            violations.append(
                CountInvariantViolation(
                    sample=row.sample,
                    rule="rrna_split",
                    message=(
                        f"rrna_aligned ({row.rrna_aligned}) + "
                        f"post_rrna_filter ({row.post_rrna_filter}) "
                        f"!= post_trim ({row.post_trim})"
                    ),
                )
            )
        if row.mt_aligned + row.unaligned_to_mt != row.post_rrna_filter:
            violations.append(
                CountInvariantViolation(
                    sample=row.sample,
                    rule="mt_split",
                    message=(
                        f"mt_aligned ({row.mt_aligned}) + "
                        f"unaligned_to_mt ({row.unaligned_to_mt}) "
                        f"!= post_rrna_filter ({row.post_rrna_filter})"
                    ),
                )
            )
        if row.mt_aligned_after_mapq > row.mt_aligned:
            violations.append(
                CountInvariantViolation(
                    sample=row.sample,
                    rule="mapq_le_mt",
                    message=(
                        f"mt_aligned_after_mapq ({row.mt_aligned_after_mapq}) "
                        f"> mt_aligned ({row.mt_aligned})"
                    ),
                )
            )
        if row.mt_aligned_after_dedup > row.mt_aligned_after_mapq:
            violations.append(
                CountInvariantViolation(
                    sample=row.sample,
                    rule="dedup_le_mapq",
                    message=(
                        f"mt_aligned_after_dedup ({row.mt_aligned_after_dedup}) "
                        f"> mt_aligned_after_mapq ({row.mt_aligned_after_mapq})"
                    ),
                )
            )
    return violations


def enforce_count_invariants(
    rows: Iterable[SampleCounts],
    *,
    allow_warning: bool = False,
) -> list[CountInvariantViolation]:
    """Run :func:`check_count_invariants` and act on the result.

    When *allow_warning* is True (CLI ``--allow-count-invariant-warning``),
    every violation is recorded via :mod:`mitoribopy.io.warnings_log` and
    the function returns the (possibly empty) list. Otherwise the first
    violation list raises :class:`CountInvariantError` so the caller can
    abort the run with a clear exit code.
    """
    rows_list = list(rows)
    violations = check_count_invariants(rows_list)
    if not violations:
        return violations
    if allow_warning:
        from ..errors import E_OUTPUT_COUNT_INVARIANT
        from ..io.warnings_log import record as _record_warning

        for v in violations:
            _record_warning(
                "align",
                f"{v.rule}: {v.message}",
                sample_id=v.sample,
                code=E_OUTPUT_COUNT_INVARIANT,
                suggested_action=(
                    "Investigate the upstream stage for this sample. "
                    "Run with the default (no --allow-count-invariant-warning) "
                    "to abort on these violations."
                ),
            )
        return violations
    raise CountInvariantError(violations)


def write_read_counts_table(
    rows: Iterable[SampleCounts],
    output_path: Path,
) -> Path:
    """Write ``read_counts.tsv`` at *output_path* and return the path.

    The output is a tab-delimited file with a single header row matching
    :func:`read_counts_columns` and one data row per :class:`SampleCounts`.
    Sample rows are written in iteration order - the caller decides
    stability (for align, rows are in the order samples are processed;
    the orchestrator sorts by sample name before writing so the file is
    deterministic).
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cols = read_counts_columns()
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(schema_header_line("read_counts.tsv"))
        handle.write("\t".join(cols) + "\n")
        for row in rows:
            handle.write(format_row(row) + "\n")
    return output_path
