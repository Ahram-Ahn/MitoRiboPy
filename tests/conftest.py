"""Pytest configuration for local source-tree imports and shared markers."""

from __future__ import annotations

import shutil
from pathlib import Path
import sys

import pytest


SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


_REAL_ALIGN_TOOLS = ("cutadapt", "bowtie2", "bowtie2-build", "samtools")


def pytest_collection_modifyitems(config, items):
    """Auto-skip @pytest.mark.requires_tools when any required tool is missing."""
    missing = [name for name in _REAL_ALIGN_TOOLS if shutil.which(name) is None]
    if not missing:
        return
    skip_marker = pytest.mark.skip(
        reason=f"requires external tools on PATH: missing {', '.join(missing)}"
    )
    for item in items:
        if "requires_tools" in item.keywords:
            item.add_marker(skip_marker)

