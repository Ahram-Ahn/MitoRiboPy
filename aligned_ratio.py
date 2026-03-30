"""Backward-compatible wrapper for legacy `aligned_ratio` script."""

from __future__ import annotations

from pathlib import Path
import sys
import warnings


_src_dir = Path(__file__).resolve().parent / "src"
if _src_dir.exists() and str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))


warnings.warn(
    "aligned_ratio is deprecated; use mitoribopy.tools.read_composition instead.",
    DeprecationWarning,
    stacklevel=2,
)

from mitoribopy.tools.read_composition import (  # noqa: E402,F401
    EnrichmentConfig,
    compute_enrichment_table,
    main,
    plot_barh,
    plot_pie,
    summarize_and_plot,
)

__all__ = [
    "EnrichmentConfig",
    "compute_enrichment_table",
    "plot_barh",
    "plot_pie",
    "summarize_and_plot",
    "main",
]


if __name__ == "__main__":
    main()
