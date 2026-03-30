"""Backward-compatible wrapper for legacy `codon_data` imports.

Phase II note:
The implementation moved to:
- ``mitoribopy.data.codon_tables``
- ``mitoribopy.data.transcript_annotations``
"""

from __future__ import annotations

from pathlib import Path
import sys
import warnings


# Keep `python main.py ...` working from repository root without requiring install.
_src_dir = Path(__file__).resolve().parent / "src"
if _src_dir.exists() and str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))


warnings.warn(
    "codon_data is deprecated; import from mitoribopy.data instead.",
    DeprecationWarning,
    stacklevel=2,
)

from mitoribopy.data.codon_tables import (  # noqa: E402,F401
    human_mitochondrial_codon_table,
    standard_codon_table,
    yeast_mitochondrial_codon_table,
)
from mitoribopy.data.transcript_annotations import (  # noqa: E402,F401
    human_annotation_df,
    yeast_annotation_df,
)

__all__ = [
    "standard_codon_table",
    "yeast_mitochondrial_codon_table",
    "human_mitochondrial_codon_table",
    "yeast_annotation_df",
    "human_annotation_df",
]
