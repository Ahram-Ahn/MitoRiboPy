"""Legacy import wrapper for Phase II BED I/O module migration."""

from __future__ import annotations

import warnings


warnings.warn(
    "mitoribopy._legacy.bed_processing is deprecated; use mitoribopy.io modules.",
    DeprecationWarning,
    stacklevel=2,
)

from mitoribopy.io.bed_reader import (  # noqa: E402,F401
    compute_unfiltered_read_length_summary,
    plot_unfiltered_read_length_heatmap,
    process_bed_files,
)
from mitoribopy.io.read_counts import compute_total_counts  # noqa: E402,F401

__all__ = [
    "compute_total_counts",
    "process_bed_files",
    "compute_unfiltered_read_length_summary",
    "plot_unfiltered_read_length_heatmap",
]
