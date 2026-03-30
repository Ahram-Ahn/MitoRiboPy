"""Legacy import wrapper for Phase II codon-correlation migration."""

from __future__ import annotations

import warnings


warnings.warn(
    "mitoribopy._legacy.codon_correlation is deprecated; use mitoribopy.analysis modules.",
    DeprecationWarning,
    stacklevel=2,
)

from mitoribopy.analysis.codon_correlation import run_codon_correlation  # noqa: E402,F401

__all__ = ["run_codon_correlation"]
