"""Legacy import wrapper for Phase II VARNA export migration."""

from __future__ import annotations

import warnings


warnings.warn(
    "mitoribopy._legacy.varna_plot is deprecated; use mitoribopy.plotting modules.",
    DeprecationWarning,
    stacklevel=2,
)

from mitoribopy.plotting.varna_export import run_varna_plot  # noqa: E402,F401

__all__ = ["run_varna_plot"]
