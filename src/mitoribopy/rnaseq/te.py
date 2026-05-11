"""Translation efficiency (TE) and delta-TE computation.

TE per sample per gene::

    TE_i = RPF_count_i / mRNA_abundance_i

For single-replicate data where no Ribo-seq group comparison is
possible, the delta-TE log2 is computed as
``log2(RPF_fc) - log2(mRNA_fc)`` using the DE table's mRNA log2FC and
an internally-computed Ribo-seq log2FC. When only one sample exists
per condition we emit the point estimate with a note and no p-value.

This module intentionally does NOT run dispersion-based statistics
(DESeq2 / Xtail) internally over the 13 mt-mRNAs; that universe is too
small for the shrinkage assumptions to hold. The user runs DE on the
full transcriptome and hands us the result.

§5 spec: every row carries a ``note`` field drawn from
:data:`mitoribopy.rnaseq._types.CANONICAL_NOTE_CODES` so a reviewer
can filter by ``publication_grade`` / ``exploratory_mt_only`` /
``insufficient_replicates`` without reading the YAML config.
"""

from __future__ import annotations

import math

from ._types import (
    DTeRow,
    DeTable,
    NOTE_EXPLORATORY_MT_ONLY,
    NOTE_INSUFFICIENT_REPLICATES,
    NOTE_MISSING_RNA_GENE,
    NOTE_PSEUDO_REPLICATE_NO_STATISTICS,
    NOTE_PUBLICATION_GRADE,
    NOTE_ZERO_COUNT_IN_CONDITION,
    TeRow,
)


_PSEUDOCOUNT = 0.5


# Method labels written into the ``method`` column of delta_te.tsv.
METHOD_DE_TABLE = "external_de_table"
METHOD_INTERNAL_MEAN_LOG2FC = "internal_mean_log2fc"
METHOD_PSEUDO_REPLICATE = "pseudo_replicate_log2fc"
METHOD_NO_STATISTICS = "single_sample_log2fc"


def _safe_log2(value: float | None) -> float | None:
    if value is None:
        return None
    if value <= 0:
        return None
    return math.log2(value)


def _row_note_for_mode(mode: str | None, *, pseudo_replicate: bool) -> str:
    """Choose the per-row ``note`` value for a TE / dTE row.

    Rules:

    * Pseudo-replicate runs always tag every row
      ``pseudo_replicate_no_statistics``.
    * ``de_table`` mode tags rows ``publication_grade``.
    * ``from_fastq`` tags rows ``exploratory_mt_only``.
    * Any other mode (or unset) leaves the note blank so the writer
      does not assert a status it can't justify.
    """
    if pseudo_replicate:
        return NOTE_PSEUDO_REPLICATE_NO_STATISTICS
    if mode == "de_table":
        return NOTE_PUBLICATION_GRADE
    if mode == "from_fastq":
        return NOTE_EXPLORATORY_MT_ONLY
    return ""


def compute_te(
    ribo_counts: dict[str, dict[str, int]],
    mrna_abundances: dict[str, float],
    *,
    pseudocount: float = _PSEUDOCOUNT,
    condition_map: dict[str, str] | None = None,
    assay: str = "ribo",
    mode: str | None = None,
    pseudo_replicate: bool = False,
) -> list[TeRow]:
    """Return one TeRow per (sample, gene) present in *ribo_counts*.

    Parameters
    ----------
    ribo_counts:
        ``{gene: {sample: count}}`` from
        :func:`mitoribopy.rnaseq.counts.load_ribo_counts`.
    mrna_abundances:
        ``{gene: baseMean}``. When a gene's mRNA abundance is missing
        we skip it: TE is undefined without an mRNA denominator.
    pseudocount:
        Added to numerator and denominator to avoid div-by-zero and
        log(0). Default 0.5 is a conservative Laplace-style prior.
    condition_map:
        Optional ``{sample_id: condition}``. Populates the
        ``condition`` column on each TeRow.
    assay:
        Recorded verbatim in the ``assay`` column. RPF-side TE rows
        use ``"ribo"``; the from_fastq orchestrator currently only
        produces ribo TE rows.
    mode:
        Run mode (``"de_table"`` / ``"from_fastq"`` / ``"rna_only"``);
        chooses the row note via :func:`_row_note_for_mode`.
    pseudo_replicate:
        When ``True``, override mode-driven notes with
        ``pseudo_replicate_no_statistics``.
    """
    rows: list[TeRow] = []
    base_note = _row_note_for_mode(mode, pseudo_replicate=pseudo_replicate)
    cmap = condition_map or {}
    for gene, per_sample in ribo_counts.items():
        mrna = mrna_abundances.get(gene)
        if mrna is None:
            continue
        mrna_den = float(mrna) + pseudocount
        for sample, count in per_sample.items():
            numerator = int(count) + pseudocount
            te = numerator / mrna_den
            rows.append(
                TeRow(
                    sample_id=sample,
                    gene=gene,
                    rpf_count=int(count),
                    rna_abundance=float(mrna),
                    te=te,
                    log2_te=_safe_log2(te),
                    condition=cmap.get(sample),
                    assay=assay,
                    note=base_note,
                )
            )
    return rows


def _group_counts_by_condition(
    ribo_counts: dict[str, dict[str, int]],
    condition_map: dict[str, str],
) -> dict[str, dict[str, list[int]]]:
    """Group ``{gene: {sample: count}}`` by condition.

    Returns ``{gene: {condition: [count, count, ...]}}`` with samples
    whose condition is not in *condition_map* silently dropped.
    """
    grouped: dict[str, dict[str, list[int]]] = {}
    for gene, per_sample in ribo_counts.items():
        for sample, count in per_sample.items():
            cond = condition_map.get(sample)
            if cond is None:
                continue
            grouped.setdefault(gene, {}).setdefault(cond, []).append(int(count))
    return grouped


def compute_delta_te(
    ribo_counts: dict[str, dict[str, int]],
    de_table: DeTable,
    *,
    condition_map: dict[str, str] | None = None,
    base_condition: str | None = None,
    compare_condition: str | None = None,
    method: str | None = None,
    mode: str | None = None,
    pseudo_replicate: bool = False,
    # Legacy aliases (kept for older test/caller code):
    condition_a: str | None = None,
    condition_b: str | None = None,
) -> list[DTeRow]:
    """Compute per-gene delta-TE log2 using the DE table's mRNA log2FC.

    The §5 spec changes the schema of the result (no field rename in
    Python — the dataclass already carries the new names — but new
    fields are populated):

    * ``base_condition`` / ``compare_condition``  — explicit labels
    * ``padj_mrna`` (was ``padj``)                — DE table's padj
    * ``padj_rpf`` / ``padj_delta_te``            — populated only when
      a real test was run; left ``None`` for single-replicate /
      pseudo-replicate runs
    * ``method``                                  — string code
      describing how the row was computed (see ``METHOD_*`` constants)
    * ``note``                                    — one of
      :data:`CANONICAL_NOTE_CODES`

    Parameters
    ----------
    base_condition / compare_condition:
        Replace the legacy ``condition_a`` / ``condition_b`` (still
        accepted for back-compat). ``compare_condition`` is the
        numerator of the RPF log2FC.
    method:
        Optional override; defaults are picked from the resolved
        replicate / mode state.
    mode:
        Run mode for tagging row notes when no per-row override fires.
    pseudo_replicate:
        Forces ``pseudo_replicate_no_statistics`` on every row.
    """
    base = base_condition or condition_a
    compare = compare_condition or condition_b

    # Back-compat: previous callers passed ``condition_a`` and ``condition_b``
    # only; the names are preserved through the local variables ``base`` /
    # ``compare`` below.

    mrna_log2fc_by_gene = {row["gene_id"]: row["log2fc"] for row in de_table.rows}
    padj_by_gene = {row["gene_id"]: row["padj"] for row in de_table.rows}

    have_replicates = (
        condition_map is not None
        and base is not None
        and compare is not None
    )

    grouped: dict[str, dict[str, list[int]]] = (
        _group_counts_by_condition(ribo_counts, condition_map)
        if have_replicates
        else {}
    )

    default_note = _row_note_for_mode(mode, pseudo_replicate=pseudo_replicate)

    rows: list[DTeRow] = []
    for gene in sorted(ribo_counts):
        mrna_log2fc = mrna_log2fc_by_gene.get(gene)
        if mrna_log2fc is None:
            rows.append(
                DTeRow(
                    gene=gene,
                    base_condition=base,
                    compare_condition=compare,
                    mrna_log2fc=None,
                    rpf_log2fc=None,
                    delta_te_log2=None,
                    padj_mrna=None,
                    padj_rpf=None,
                    padj_delta_te=None,
                    method=(method or METHOD_DE_TABLE),
                    note=NOTE_MISSING_RNA_GENE,
                )
            )
            continue

        rpf_log2fc: float | None = None
        row_note = default_note
        row_method = method
        if have_replicates:
            cond_counts = grouped.get(gene, {})
            a_counts = cond_counts.get(base, [])
            b_counts = cond_counts.get(compare, [])
            if a_counts and b_counts:
                if not any(a_counts) or not any(b_counts):
                    row_note = NOTE_ZERO_COUNT_IN_CONDITION
                mean_a = sum(a_counts) / len(a_counts) + _PSEUDOCOUNT
                mean_b = sum(b_counts) / len(b_counts) + _PSEUDOCOUNT
                ratio = mean_b / mean_a
                rpf_log2fc = _safe_log2(ratio)
                if row_method is None:
                    row_method = (
                        METHOD_PSEUDO_REPLICATE
                        if pseudo_replicate
                        else METHOD_INTERNAL_MEAN_LOG2FC
                    )
            else:
                row_note = NOTE_INSUFFICIENT_REPLICATES
                if row_method is None:
                    row_method = METHOD_INTERNAL_MEAN_LOG2FC
        else:
            row_note = NOTE_INSUFFICIENT_REPLICATES
            if row_method is None:
                row_method = METHOD_NO_STATISTICS

        delta_te = (
            (rpf_log2fc - mrna_log2fc) if rpf_log2fc is not None else None
        )
        rows.append(
            DTeRow(
                gene=gene,
                base_condition=base,
                compare_condition=compare,
                mrna_log2fc=mrna_log2fc,
                rpf_log2fc=rpf_log2fc,
                delta_te_log2=delta_te,
                padj_mrna=padj_by_gene.get(gene),
                padj_rpf=None,
                padj_delta_te=None,
                method=row_method or METHOD_INTERNAL_MEAN_LOG2FC,
                note=row_note,
            )
        )
    return rows
