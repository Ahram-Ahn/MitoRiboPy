#!/usr/bin/env python3
"""Backward-compatible wrapper for legacy `subsample_bed` script."""

from __future__ import annotations

from pathlib import Path
import sys
import warnings


_src_dir = Path(__file__).resolve().parent / "src"
if _src_dir.exists() and str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))


warnings.warn(
    "subsample_bed is deprecated; use mitoribopy.tools.subsample instead.",
    DeprecationWarning,
    stacklevel=2,
)

from mitoribopy.tools.subsample import (  # noqa: E402,F401
    main,
    parse_args,
    reservoir_sample_lines,
)

__all__ = ["parse_args", "reservoir_sample_lines", "main"]


if __name__ == "__main__":
    main()
