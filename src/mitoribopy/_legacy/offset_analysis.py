"""Legacy import wrapper for Phase II offset module migration."""

from __future__ import annotations

import warnings


warnings.warn(
    "mitoribopy._legacy.offset_analysis is deprecated; use mitoribopy.analysis modules.",
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
