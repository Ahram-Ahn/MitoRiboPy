"""Shared runtime state for package-native pipeline execution."""

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
    plot_output_dir: Path
    annotation_df: pd.DataFrame
    rpf_range: list[int]
    total_counts_map: dict[str, int | float] = field(default_factory=dict)
    filtered_bed_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    sample_dirs: list[str] = field(default_factory=list)
    offset_summary_df: pd.DataFrame | None = None
    offset_details_df: pd.DataFrame | None = None
    selected_offsets_df: pd.DataFrame | None = None
    extra: dict[str, Any] = field(default_factory=dict)
