"""Per-sample read-counts provenance table for ``mitoribopy align``.

Step H of the pipeline. After every sample has been trimmed, filtered,
aligned, MAPQ-filtered, and deduplicated, we assemble a
:class:`~mitoribopy.align._types.SampleCounts` from the per-stage result
dataclasses and write a tab-delimited ``read_counts.tsv`` in the run
root.

This file is the provenance spine: it documents how many reads survived
each step so downstream Phase-6 manifest consumers and human reviewers
can reconstruct where the read attrition happened.
"""

from __future__ import annotations

from dataclasses import fields
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
