#!/usr/bin/env python3
"""Repository-root compatibility wrapper for the package-native pipeline."""

from __future__ import annotations

from pathlib import Path
import sys


_SRC_DIR = Path(__file__).resolve().parent / "src"
if _SRC_DIR.exists() and str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

from mitoribopy.pipeline.runner import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
