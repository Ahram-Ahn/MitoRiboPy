"""Shared runtime state for package-native pipeline execution.

Convention (Task 5c): after ``build_pipeline_context`` returns, the
``args`` namespace is treated as **read-only input**. Every piece of
derived state lives on ``PipelineContext`` itself (e.g.
``total_counts_map``, ``filtered_bed_df``, ``sample_dirs``,
``selected_offsets_df``). This keeps the "input parameters" record
clean so a reviewer reading ``run_settings.json`` can trust that what
they see is what was passed in, not a mutated snapshot taken mid-run.

The one remaining shim - ``context.args.total_mrna_map =
context.total_counts_map`` - is there only to keep older consumers
working through one deprecation cycle; ``run_coverage_profile_plots``
now accepts ``total_mrna_map`` as an explicit keyword argument and
emits a deprecation warning when it has to fall back to the args copy.
The shim will be removed in v0.4.0.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass
class PipelineContext:
    """Mutable state shared across standalone pipeline steps."""

    args: argparse.Namespace
    base_output_dir: Path
    # ``plot_output_dir`` is the parent of ``csv_dir`` and ``plot_subdir``;
    # it is kept for backward-compatible callers that just need the
    # offset-diagnostics root. New code writes CSVs to ``csv_dir`` and
    # plots to ``plot_subdir`` so the two file types stay separated.
    plot_output_dir: Path
    csv_dir: Path
    plot_subdir: Path
    annotation_df: pd.DataFrame
    resolved_codon_table: dict[str, str]
    resolved_start_codons: list[str]
    rpf_range: list[int]
    unfiltered_read_length_range: tuple[int, int]
    total_counts_map: dict[str, int | float] = field(default_factory=dict)
    filtered_bed_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    sample_dirs: list[str] = field(default_factory=list)
    offset_summary_df: pd.DataFrame | None = None
    offset_details_df: pd.DataFrame | None = None
    # Combined-across-samples selected offsets. Always populated when
    # offsets were selected, and surfaced as a diagnostic regardless of
    # `--offset_mode`.
    selected_offsets_df: pd.DataFrame | None = None
    # Per-sample offset summaries (input to the per-sample offset
    # picker). One DataFrame per sample, keyed by sample name.
    per_sample_offset_summaries: dict[str, pd.DataFrame] = field(default_factory=dict)
    # Per-sample selected offsets. Keys are sample names. Empty when
    # `--offset_mode combined` is active (downstream then uses
    # `selected_offsets_df` for every sample).
    selected_offsets_by_sample: dict[str, pd.DataFrame] = field(default_factory=dict)
    # Read-length filter bookkeeping populated by `_apply_rpf_count_filter`.
    # These fields drive the rpf_counts.metadata.json sidecar so a
    # reviewer can audit which read lengths actually fed the counts and
    # why the others got dropped.
    requested_rpf_range: list[int] = field(default_factory=list)
    observed_counts_by_length: dict[int, int] = field(default_factory=dict)
    dropped_lengths: list[int] = field(default_factory=list)
    read_length_filter_threshold_rule: str | None = None
    read_length_filter_threshold_value: float | None = None
    extra: dict[str, Any] = field(default_factory=dict)
