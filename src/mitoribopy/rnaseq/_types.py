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
# Spec-aligned note codes for TE / dTE rows (§5).
#
# These are part of the public TSV contract: a downstream filter that
# greps for ``insufficient_replicates`` should keep working across
# minor releases. Add new codes additively; do not silently rename.
# ---------------------------------------------------------------------------


NOTE_PUBLICATION_GRADE = "publication_grade"
NOTE_EXPLORATORY_MT_ONLY = "exploratory_mt_only"
NOTE_PSEUDO_REPLICATE_NO_STATISTICS = "pseudo_replicate_no_statistics"
NOTE_INSUFFICIENT_REPLICATES = "insufficient_replicates"
NOTE_MISSING_RNA_GENE = "missing_rna_gene"
NOTE_MISSING_RPF_GENE = "missing_rpf_gene"
NOTE_ZERO_COUNT_IN_CONDITION = "zero_count_in_condition"
NOTE_GENE_ID_UNMATCHED = "gene_id_unmatched"

CANONICAL_NOTE_CODES: tuple[str, ...] = (
    NOTE_PUBLICATION_GRADE,
    NOTE_EXPLORATORY_MT_ONLY,
    NOTE_PSEUDO_REPLICATE_NO_STATISTICS,
    NOTE_INSUFFICIENT_REPLICATES,
    NOTE_MISSING_RNA_GENE,
    NOTE_MISSING_RPF_GENE,
    NOTE_ZERO_COUNT_IN_CONDITION,
    NOTE_GENE_ID_UNMATCHED,
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
# TE / delta-TE result rows (§5 spec schema)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TeRow:
    """One per (sample, gene) pair, emitted into ``te.tsv``.

    §5 spec columns: ``sample_id, condition, assay, gene, rpf_count,
    rna_abundance, te, log2_te, note``. The dataclass field names match
    those columns 1:1; the writer in ``cli/rnaseq.py`` is the only
    place that converts to TSV.

    Backward-compatible properties ``sample`` and ``mrna_abundance``
    keep older callers / plotting helpers working without an immediate
    rewrite.
    """

    sample_id: str
    gene: str
    rpf_count: int
    rna_abundance: float
    te: float
    log2_te: float | None = None
    condition: str | None = None
    assay: str = "ribo"
    note: str = ""

    # ----- backward-compat aliases -----------------------------------------
    @property
    def sample(self) -> str:  # noqa: D401
        """Deprecated alias for :attr:`sample_id`."""
        return self.sample_id

    @property
    def mrna_abundance(self) -> float:  # noqa: D401
        """Deprecated alias for :attr:`rna_abundance`."""
        return self.rna_abundance


@dataclass(frozen=True)
class DTeRow:
    """One per gene, emitted into ``delta_te.tsv``.

    §5 spec columns: ``gene, base_condition, compare_condition,
    mrna_log2fc, rpf_log2fc, delta_te_log2, padj_mrna, padj_rpf,
    padj_delta_te, method, note``.

    The previous schema had a single ``padj`` field that carried the
    DE table's mRNA padj; that value now lives in ``padj_mrna`` and is
    exposed under the legacy attribute name as a back-compat property
    so existing plotting code keeps working.
    """

    gene: str
    base_condition: str | None = None
    compare_condition: str | None = None
    mrna_log2fc: float | None = None
    rpf_log2fc: float | None = None
    delta_te_log2: float | None = None
    padj_mrna: float | None = None
    padj_rpf: float | None = None
    padj_delta_te: float | None = None
    method: str = ""
    note: str = ""

    # ----- backward-compat aliases -----------------------------------------
    @property
    def padj(self) -> float | None:  # noqa: D401
        """Deprecated alias: the legacy single ``padj`` column carried
        the DE table's mRNA padj for the gene. New code should pick
        between :attr:`padj_mrna`, :attr:`padj_rpf`, and
        :attr:`padj_delta_te` explicitly."""
        return self.padj_mrna
