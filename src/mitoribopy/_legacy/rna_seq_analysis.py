"""Legacy import wrapper for Phase II RNA-seq analysis migration."""

from __future__ import annotations

import warnings


warnings.warn(
    "mitoribopy._legacy.rna_seq_analysis is deprecated; use mitoribopy.analysis modules.",
    DeprecationWarning,
    stacklevel=2,
)

from mitoribopy.analysis.rna_seq import run_rna_seq_analysis  # noqa: E402,F401

__all__ = ["run_rna_seq_analysis"]
