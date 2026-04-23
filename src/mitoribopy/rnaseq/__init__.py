"""Phase 5 RNA-seq integration for MitoRiboPy.

Consumes a pre-computed differential-expression table (DESeq2 / Xtail /
Anota2Seq) plus a prior ``mitoribopy rpf`` run and produces TE and
delta-TE tables + plots. A SHA256 reference-consistency gate ensures
Ribo-seq and RNA-seq were aligned to the SAME transcript set - a hard
requirement for TE comparisons.
"""

from ._types import (
    DE_COLUMN_ALIASES,
    GENE_ID_CONVENTIONS,
    DeColumnMap,
    DeFormat,
    DeTable,
    DTeRow,
    GeneIdConvention,
    TeRow,
)
from .de_loader import detect_de_format, load_de_table
from .gene_ids import (
    HUMAN_MT_MRNAS,
    YEAST_MT_MRNAS,
    MtGene,
    match_mt_mrnas,
    normalize_gene_id,
)
from .reference_gate import (
    ReferenceMismatchError,
    compute_reference_checksum,
    verify_reference_consistency,
)
from .te import compute_delta_te, compute_te
from .counts import load_ribo_counts

__all__ = [
    "DE_COLUMN_ALIASES",
    "DTeRow",
    "DeColumnMap",
    "DeFormat",
    "DeTable",
    "GENE_ID_CONVENTIONS",
    "GeneIdConvention",
    "HUMAN_MT_MRNAS",
    "MtGene",
    "ReferenceMismatchError",
    "TeRow",
    "YEAST_MT_MRNAS",
    "compute_delta_te",
    "compute_reference_checksum",
    "compute_te",
    "detect_de_format",
    "load_de_table",
    "load_ribo_counts",
    "match_mt_mrnas",
    "normalize_gene_id",
    "verify_reference_consistency",
]
