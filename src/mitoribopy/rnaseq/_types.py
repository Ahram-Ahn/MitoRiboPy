"""Shared dataclasses + literals for ``mitoribopy.rnaseq``."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


# ---------------------------------------------------------------------------
# Literals
# ---------------------------------------------------------------------------


DeFormat = Literal["deseq2", "xtail", "anota2seq", "custom"]
"""Supported differential-expression output schemas.

``deseq2``    DESeq2::results() CSV/TSV.
``xtail``     Xtail::xtail() results table.
``anota2seq`` Anota2Seq output table.
``custom``    Caller supplies explicit column mapping via DeColumnMap.
"""

GeneIdConvention = Literal["ensembl", "refseq", "hgnc", "mt_prefixed", "bare"]
"""Gene identifier conventions used in DE tables.

``ensembl``       e.g. ``ENSG00000198888`` (human mt-ND1)
``refseq``        e.g. ``NM_001385440`` (mt-mRNA NCBI accession)
``hgnc``          HGNC symbol, e.g. ``MT-ND1``
``mt_prefixed``   ``MT-`` prefix, identical to HGNC for mt-mRNAs
``bare``          no prefix, e.g. ``ND1``
"""


GENE_ID_CONVENTIONS: tuple[str, ...] = (
    "ensembl",
    "refseq",
    "hgnc",
    "mt_prefixed",
    "bare",
)


# ---------------------------------------------------------------------------
# DE table abstraction
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DeColumnMap:
    """Explicit mapping of DE-table column names to canonical fields.

    Filled in by ``detect_de_format(header)`` for recognized schemas, or
    by the user for ``--de-format custom``. Missing fields (``None``)
    mean that column is not present in the source; downstream code must
    handle that gracefully.
    """

    gene_id: str
    log2fc: str | None
    padj: str | None
    basemean: str | None


# Column name aliases for each recognized format. First match wins, so
# more specific names (like 'log2FC_TE_final') appear before generic
# ('log2FoldChange'). Every value MUST appear in a real output from the
# named tool; vendor docs cross-checked.
DE_COLUMN_ALIASES: dict[str, DeColumnMap] = {
    "deseq2": DeColumnMap(
        gene_id="gene_id",  # fallback to first column on load
        log2fc="log2FoldChange",
        padj="padj",
        basemean="baseMean",
    ),
    "xtail": DeColumnMap(
        gene_id="gene_id",  # first column
        log2fc="log2FC_TE_final",
        padj="pvalue.adjust",
        basemean="mRNA_log2FC_A",  # Xtail does not emit a baseMean;
        # we reuse mRNA_log2FC_A as the 'mRNA' reference column.
    ),
    "anota2seq": DeColumnMap(
        gene_id="gene_id",
        log2fc="slopeTranslation",
        padj="pAdjustedRvmRvm",
        basemean="TotalExpression",  # anota2seq's baseline mRNA abundance
    ),
}


@dataclass(frozen=True)
class DeTable:
    """Normalized DE-table representation used across the rnaseq subcommand.

    ``rows`` is a list of dicts keyed by canonical names:
    ``gene_id``, ``log2fc``, ``padj``, ``basemean`` (any of the latter
    three may be ``None`` if absent in the source).
    ``format`` records which auto-detected schema was used so the run
    manifest can reproduce the decision.
    """

    format: DeFormat
    column_map: DeColumnMap
    rows: list[dict]


# ---------------------------------------------------------------------------
# TE / delta-TE result rows
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TeRow:
    """One per sample per gene, emitted into ``te.tsv``."""

    sample: str
    gene: str
    rpf_count: int
    mrna_abundance: float
    te: float


@dataclass(frozen=True)
class DTeRow:
    """One per gene, emitted into ``delta_te.tsv``."""

    gene: str
    mrna_log2fc: float | None
    rpf_log2fc: float | None
    delta_te_log2: float | None
    padj: float | None
    note: str = ""
