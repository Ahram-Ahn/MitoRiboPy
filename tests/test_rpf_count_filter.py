"""Unit tests for the read-length auto-filter (``_apply_rpf_count_filter``).

The filter is data-driven: it drops read-length bins whose total count
across all samples is below ``--rpf_min_count_frac`` x the most-enriched
length's count. This shrinks ``context.rpf_range`` and re-filters
``context.filtered_bed_df`` so the downstream offset-enrichment step
only sees signal-bearing bins.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from mitoribopy.pipeline.context import PipelineContext
from mitoribopy.pipeline.steps import _apply_rpf_count_filter


def _make_context(read_lengths: list[int], counts: list[int], frac: float) -> PipelineContext:
    """Build a minimal PipelineContext with a synthetic filtered_bed_df."""
    rows = []
    for rl, n in zip(read_lengths, counts):
        for i in range(n):
            rows.append(
                {
                    "chrom": "tx",
                    "start": i,
                    "end": i + rl,
                    "read_length": rl,
                    "sample_name": "S1",
                }
            )
    df = pd.DataFrame(rows)
    args = SimpleNamespace(rpf_min_count_frac=frac)
    return PipelineContext(
        args=args,
        base_output_dir=Path("/tmp/unused"),
        plot_output_dir=Path("/tmp/unused/plots"),
        csv_dir=Path("/tmp/unused/plots/csv"),
        plot_subdir=Path("/tmp/unused/plots/plots"),
        annotation_df=pd.DataFrame(),
        resolved_codon_table={},
        resolved_start_codons=[],
        rpf_range=list(read_lengths),
        unfiltered_read_length_range=(min(read_lengths), max(read_lengths)),
        filtered_bed_df=df,
        sample_dirs=["S1"],
    )


def _capture_status() -> tuple[list[str], callable]:
    captured: list[str] = []
    return captured, captured.append


def test_filter_drops_bins_below_threshold() -> None:
    # Most enriched length (30) has 1000 reads; threshold at 20% = 200.
    # Read length 27 has 100 (drop), 28 has 250 (keep), 29 has 800 (keep),
    # 30 has 1000 (keep), 31 has 50 (drop).
    ctx = _make_context(
        read_lengths=[27, 28, 29, 30, 31],
        counts=[100, 250, 800, 1000, 50],
        frac=0.20,
    )
    log, emit = _capture_status()

    _apply_rpf_count_filter(ctx, emit)

    assert ctx.rpf_range == [28, 29, 30]
    assert sorted(ctx.filtered_bed_df["read_length"].unique().tolist()) == [28, 29, 30]
    assert any("kept [28, 29, 30]" in line for line in log)
    assert any("dropped [27, 31]" in line for line in log)
    assert any("20%" in line for line in log)


def test_filter_disabled_when_frac_zero() -> None:
    ctx = _make_context(
        read_lengths=[27, 30],
        counts=[10, 1000],
        frac=0.0,
    )
    log, emit = _capture_status()

    _apply_rpf_count_filter(ctx, emit)

    assert ctx.rpf_range == [27, 30]
    assert log == []
    # The DataFrame is untouched.
    assert sorted(ctx.filtered_bed_df["read_length"].unique().tolist()) == [27, 30]


def test_filter_quiet_when_nothing_pruned() -> None:
    # All bins are equally populated => nothing to drop.
    ctx = _make_context(
        read_lengths=[28, 29, 30],
        counts=[500, 500, 500],
        frac=0.20,
    )
    log, emit = _capture_status()

    _apply_rpf_count_filter(ctx, emit)

    assert ctx.rpf_range == [28, 29, 30]
    assert log == []  # quiet when there is nothing to report


def test_filter_handles_empty_dataframe() -> None:
    df = pd.DataFrame(columns=["chrom", "start", "end", "read_length", "sample_name"])
    args = SimpleNamespace(rpf_min_count_frac=0.20)
    ctx = PipelineContext(
        args=args,
        base_output_dir=Path("/tmp/unused"),
        plot_output_dir=Path("/tmp/unused/plots"),
        csv_dir=Path("/tmp/unused/plots/csv"),
        plot_subdir=Path("/tmp/unused/plots/plots"),
        annotation_df=pd.DataFrame(),
        resolved_codon_table={},
        resolved_start_codons=[],
        rpf_range=[28, 29, 30],
        unfiltered_read_length_range=(28, 30),
        filtered_bed_df=df,
        sample_dirs=[],
    )
    log, emit = _capture_status()

    _apply_rpf_count_filter(ctx, emit)

    assert ctx.rpf_range == [28, 29, 30]
    assert log == []


def test_filter_threshold_is_inclusive() -> None:
    # A bin exactly at the threshold (20% of max) is kept.
    ctx = _make_context(
        read_lengths=[28, 30],
        counts=[200, 1000],  # 200 == 0.20 * 1000
        frac=0.20,
    )
    log, emit = _capture_status()

    _apply_rpf_count_filter(ctx, emit)

    assert ctx.rpf_range == [28, 30]  # both retained
    assert log == []  # nothing pruned, quiet
