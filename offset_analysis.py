"""Backward-compatible wrapper for legacy `offset_analysis` imports.

Phase II note:
The implementation moved to:
- ``mitoribopy.analysis.offset_enrichment``
- ``mitoribopy.analysis.offset_selection``
"""

from __future__ import annotations

import warnings
from pathlib import Path
import sys


# Keep `python main.py ...` working from repository root without requiring install.
_src_dir = Path(__file__).resolve().parent / "src"
if _src_dir.exists() and str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))


warnings.warn(
    "offset_analysis is deprecated; import from mitoribopy.analysis instead.",
    DeprecationWarning,
    stacklevel=2,
)

from mitoribopy.analysis.offset_enrichment import (  # noqa: E402,F401
    compute_offsets,
    create_csv_for_offset_enrichment,
    plot_offset_enrichment,
)
from mitoribopy.analysis.offset_selection import (  # noqa: E402,F401
    determine_p_site_offsets,
)

__all__ = [
    "compute_offsets",
    "create_csv_for_offset_enrichment",
    "plot_offset_enrichment",
    "determine_p_site_offsets",
]
