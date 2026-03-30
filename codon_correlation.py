"""Backward-compatible wrapper for legacy `codon_correlation` imports."""

from __future__ import annotations

from pathlib import Path
import sys
import warnings


_src_dir = Path(__file__).resolve().parent / "src"
if _src_dir.exists() and str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))


warnings.warn(
    "codon_correlation is deprecated; import from mitoribopy.analysis instead.",
    DeprecationWarning,
    stacklevel=2,
)

from mitoribopy.analysis.codon_correlation import run_codon_correlation  # noqa: E402,F401

__all__ = ["run_codon_correlation"]
