"""Backward-compatible wrapper for legacy `igv_style_plot` imports."""

from __future__ import annotations

from pathlib import Path
import sys
import warnings


_src_dir = Path(__file__).resolve().parent / "src"
if _src_dir.exists() and str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))


warnings.warn(
    "igv_style_plot is deprecated; import from mitoribopy.plotting instead.",
    DeprecationWarning,
    stacklevel=2,
)

from mitoribopy.plotting.igv_plots import run_igv_style_plot  # noqa: E402,F401

__all__ = ["run_igv_style_plot"]
