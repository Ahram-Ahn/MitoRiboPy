"""Backward-compatible wrapper for legacy `inframe_analysis` imports."""

from __future__ import annotations

from pathlib import Path
import sys
import warnings


_src_dir = Path(__file__).resolve().parent / "src"
if _src_dir.exists() and str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))


warnings.warn(
    "inframe_analysis is deprecated; import from mitoribopy.analysis instead.",
    DeprecationWarning,
    stacklevel=2,
)

from mitoribopy.analysis.inframe_analysis import run_inframe_codon_analysis  # noqa: E402,F401

__all__ = ["run_inframe_codon_analysis"]
