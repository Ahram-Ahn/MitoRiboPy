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
    """Auto-skip optional / opt-in test groups.

    * ``@pytest.mark.requires_tools`` is skipped when any of the
      external bioinformatics tools is missing from ``PATH``.
    * ``@pytest.mark.smoke`` is opt-in (per the marker description in
      ``pyproject.toml``): it only runs when the user explicitly
      selects it with ``pytest -m smoke``. The smoke fixture under
      ``examples/smoke/`` is still labelled "Planned" in CHANGELOG
      and is not part of the default suite.
    """
    missing = [name for name in _REAL_ALIGN_TOOLS if shutil.which(name) is None]
    if missing:
        skip_tools = pytest.mark.skip(
            reason=f"requires external tools on PATH: missing {', '.join(missing)}"
        )
        for item in items:
            if "requires_tools" in item.keywords:
                item.add_marker(skip_tools)

    marker_expr = config.getoption("-m", default="") or ""
    if "smoke" not in marker_expr:
        skip_smoke = pytest.mark.skip(
            reason="smoke fixture is opt-in; run with `pytest -m smoke`"
        )
        for item in items:
            if "smoke" in item.keywords:
                item.add_marker(skip_smoke)

