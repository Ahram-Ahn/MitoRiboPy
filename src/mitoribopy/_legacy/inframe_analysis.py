"""Legacy import wrapper for Phase II in-frame analysis migration."""

from __future__ import annotations

import warnings


warnings.warn(
    "mitoribopy._legacy.inframe_analysis is deprecated; use mitoribopy.analysis modules.",
    DeprecationWarning,
    stacklevel=2,
)

from mitoribopy.analysis.inframe_analysis import run_inframe_codon_analysis  # noqa: E402,F401

__all__ = ["run_inframe_codon_analysis"]
