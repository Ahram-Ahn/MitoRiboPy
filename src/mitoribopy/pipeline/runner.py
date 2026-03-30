"""Phase I pipeline runner wrappers.

This module keeps package entrypoints working while the legacy script-based
implementation is being migrated into package modules.
"""

from __future__ import annotations

import importlib
from pathlib import Path
import sys


def run_legacy_pipeline(argv: list[str]) -> int:
    """Run the current legacy pipeline entrypoint with forwarded CLI args."""
    legacy_dir = Path(__file__).resolve().parents[1] / "_legacy"
    if not legacy_dir.exists():
        raise RuntimeError(
            f"Packaged legacy directory missing: {legacy_dir}"
        )

    original_argv = sys.argv
    original_sys_path = list(sys.path)
    try:
        # Keep compatibility with legacy-style imports by prioritizing packaged _legacy.
        sys.path = [str(legacy_dir)] + original_sys_path
        legacy_main = importlib.import_module("main").main
        # Keep legacy argparse behavior identical to `python main.py ...`.
        sys.argv = ["mitoribopy"] + argv
        legacy_main()
    finally:
        sys.path = original_sys_path
        sys.argv = original_argv
    return 0
