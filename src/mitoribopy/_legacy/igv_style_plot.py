"""Legacy import wrapper for Phase II IGV plotting migration."""

from __future__ import annotations

import warnings


warnings.warn(
    "mitoribopy._legacy.igv_style_plot is deprecated; use mitoribopy.plotting modules.",
    DeprecationWarning,
    stacklevel=2,
)

from mitoribopy.plotting.igv_plots import run_igv_style_plot  # noqa: E402,F401

__all__ = ["run_igv_style_plot"]
