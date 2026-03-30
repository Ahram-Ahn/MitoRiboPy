"""Backward-compatible wrapper for legacy `bed_processing` imports.

Phase II note:
The implementation moved to:
- ``mitoribopy.io.bed_reader``
- ``mitoribopy.io.read_counts``
"""

from __future__ import annotations

from pathlib import Path
import sys
import warnings


# Keep `python main.py ...` working from repository root without requiring install.
_src_dir = Path(__file__).resolve().parent / "src"
if _src_dir.exists() and str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))


warnings.warn(
    "bed_processing is deprecated; import from mitoribopy.io instead.",
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
