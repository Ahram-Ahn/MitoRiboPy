"""Legacy import wrapper for Phase II codon-data module migration."""

from __future__ import annotations

import warnings


warnings.warn(
    "mitoribopy._legacy.codon_data is deprecated; use mitoribopy.data modules.",
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
