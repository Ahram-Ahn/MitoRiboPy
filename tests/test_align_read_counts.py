"""Unit tests for ``mitoribopy.align.read_counts``."""

from __future__ import annotations

from pathlib import Path

import pytest

from mitoribopy.align._types import (
    AlignResult,
    ContamResult,
    CutadaptResult,
    DedupResult,
    SampleCounts,
)
from mitoribopy.align.read_counts import (
    assemble_sample_counts,
    format_row,
    read_counts_columns,
    write_read_counts_table,
)


# ---------- column order -----------------------------------------------------


def test_read_counts_columns_are_in_dataclass_order() -> None:
    # The column order matches SampleCounts dataclass field order so
    # positional indexing stays stable across releases.
    assert read_counts_columns() == (
        "sample",
        "total_reads",
        "post_trim",
        "rrna_aligned",
        "post_rrna_filter",
        "mt_aligned",
        "unaligned_to_mt",
        "mt_aligned_after_mapq",
        "mt_aligned_after_dedup",
    )


# ---------- assemble_sample_counts -------------------------------------------


def _stub_stage_results(tmp_path: Path) -> dict:
    return dict(
        cutadapt=CutadaptResult(
            input_reads=1_000_000,
            reads_with_adapter=950_000,
            reads_passing_filters=900_000,
            log_json_path=tmp_path / "cut.json",
        ),
        contam=ContamResult(
            total_reads=900_000,
            aligned_to_contam=400_000,  # rRNA bound
            unaligned_reads=500_000,  # enters mt-align
        ),
        align=AlignResult(
            total_reads=500_000,
            aligned=480_000,
            bam_path=tmp_path / "mt.bam",
        ),
        mapq_count=470_000,
        dedup=DedupResult(
            strategy="umi-tools",
            input_reads=470_000,
            output_reads=310_000,
            bam_path=tmp_path / "dedup.bam",
        ),
    )


def test_assemble_sample_counts_populates_every_stage(tmp_path) -> None:
    stages = _stub_stage_results(tmp_path)

    counts = assemble_sample_counts(sample="sampleA", **stages)

    assert isinstance(counts, SampleCounts)
    assert counts.sample == "sampleA"
    assert counts.total_reads == 1_000_000
    assert counts.post_trim == 900_000
    assert counts.rrna_aligned == 400_000
    assert counts.post_rrna_filter == 500_000
    assert counts.mt_aligned == 480_000
    assert counts.unaligned_to_mt == 500_000 - 480_000
    assert counts.mt_aligned_after_mapq == 470_000
    assert counts.mt_aligned_after_dedup == 310_000


def test_assemble_sample_counts_invariants_hold(tmp_path) -> None:
    stages = _stub_stage_results(tmp_path)
    counts = assemble_sample_counts(sample="sampleA", **stages)

    # Key invariants documented in the module docstring:
    assert counts.total_reads == stages["cutadapt"].input_reads
    assert counts.post_trim == stages["cutadapt"].reads_passing_filters
    assert counts.rrna_aligned + counts.post_rrna_filter == counts.post_trim
    assert counts.mt_aligned + counts.unaligned_to_mt == counts.post_rrna_filter
    assert counts.mt_aligned_after_dedup <= counts.mt_aligned_after_mapq
    assert counts.mt_aligned_after_mapq <= counts.mt_aligned


def test_assemble_sample_counts_unaligned_to_mt_clamps_at_zero(tmp_path) -> None:
    # Defensive: if mt_aligned exceeds post_rrna_filter (should never
    # happen under valid bowtie2 output) we clamp the derived field to
    # zero rather than emit a negative count.
    stages = _stub_stage_results(tmp_path)
    stages["align"] = AlignResult(
        total_reads=500_000,
        aligned=600_000,  # nonsense on purpose
        bam_path=tmp_path / "mt.bam",
    )

    counts = assemble_sample_counts(sample="sampleA", **stages)
    assert counts.unaligned_to_mt == 0


# ---------- write_read_counts_table -----------------------------------------


def test_write_read_counts_table_writes_header_and_rows_in_order(tmp_path) -> None:
    output = tmp_path / "read_counts.tsv"
    rows = [
        SampleCounts(
            sample="A",
            total_reads=10,
            post_trim=9,
            rrna_aligned=3,
            post_rrna_filter=6,
            mt_aligned=5,
            unaligned_to_mt=1,
            mt_aligned_after_mapq=4,
            mt_aligned_after_dedup=3,
        ),
        SampleCounts(
            sample="B",
            total_reads=20,
            post_trim=18,
            rrna_aligned=8,
            post_rrna_filter=10,
            mt_aligned=9,
            unaligned_to_mt=1,
            mt_aligned_after_mapq=8,
            mt_aligned_after_dedup=7,
        ),
    ]

    result = write_read_counts_table(rows, output)

    assert result == output
    lines = output.read_text().splitlines()
    assert lines[0].split("\t") == list(read_counts_columns())
    assert lines[1].split("\t")[0] == "A"
    assert lines[2].split("\t")[0] == "B"


def test_write_read_counts_table_creates_parent_directory(tmp_path) -> None:
    output = tmp_path / "nested" / "a" / "b" / "read_counts.tsv"

    write_read_counts_table([], output)

    assert output.exists()
    # Empty rows list still writes a header-only file.
    assert output.read_text().splitlines() == ["\t".join(read_counts_columns())]


def test_format_row_produces_tab_separated_values() -> None:
    row = SampleCounts(
        sample="sampleX",
        total_reads=100,
        post_trim=95,
        rrna_aligned=30,
        post_rrna_filter=65,
        mt_aligned=60,
        unaligned_to_mt=5,
        mt_aligned_after_mapq=58,
        mt_aligned_after_dedup=50,
    )

    formatted = format_row(row)

    assert formatted == "sampleX\t100\t95\t30\t65\t60\t5\t58\t50"
