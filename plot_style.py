"""Backward-compatible wrapper for legacy `plot_style` imports."""

from __future__ import annotations

from pathlib import Path
import sys
import warnings


_src_dir = Path(__file__).resolve().parent / "src"
if _src_dir.exists() and str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))


warnings.warn(
    "plot_style is deprecated; import from mitoribopy.plotting.style instead.",
    DeprecationWarning,
    stacklevel=2,
)

from mitoribopy.plotting.style import (  # noqa: E402,F401
    PUBLICATION_STYLE,
    apply_publication_style,
)

__all__ = ["PUBLICATION_STYLE", "apply_publication_style"]
