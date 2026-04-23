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
"""

from __future__ import annotations

import math
from dataclasses import replace
from typing import Iterable

from ._types import DTeRow, DeTable, TeRow


_PSEUDOCOUNT = 0.5


def _safe_log2(value: float | None) -> float | None:
    if value is None:
        return None
    if value <= 0:
        return None
    return math.log2(value)


def compute_te(
    ribo_counts: dict[str, dict[str, int]],
    mrna_abundances: dict[str, float],
    *,
    pseudocount: float = _PSEUDOCOUNT,
) -> list[TeRow]:
    """Return one TeRow per (sample, gene) present in *ribo_counts*.

    ``mrna_abundances`` maps gene -> baseMean (or equivalent). When a
    gene's mRNA abundance is missing we skip it: TE is undefined
    without an mRNA denominator.

    ``pseudocount`` is added to both numerator and denominator to avoid
    div-by-zero and log(0). Default 0.5 is the conservative Laplace-style
    prior that does not dominate real counts at mt-Ribo-seq depths
    (~10^4-10^6 per sample per gene is typical).
    """
    rows: list[TeRow] = []
    for gene, per_sample in ribo_counts.items():
        mrna = mrna_abundances.get(gene)
        if mrna is None:
            continue
        mrna_den = float(mrna) + pseudocount
        for sample, count in per_sample.items():
            numerator = int(count) + pseudocount
            rows.append(
                TeRow(
                    sample=sample,
                    gene=gene,
                    rpf_count=int(count),
                    mrna_abundance=float(mrna),
                    te=numerator / mrna_den,
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
    condition_a: str | None = None,
    condition_b: str | None = None,
) -> list[DTeRow]:
    """Compute per-gene delta-TE log2 using the DE table's mRNA log2FC.

    Parameters
    ----------
    ribo_counts:
        ``{gene: {sample: rpf_count}}`` from
        :func:`mitoribopy.rnaseq.counts.load_ribo_counts`.
    de_table:
        Parsed DE table (:class:`DeTable`). We use its ``log2fc`` as
        the mRNA log2FC between the same two conditions that
        ``condition_a`` / ``condition_b`` name on the Ribo-seq side.
    condition_map, condition_a, condition_b:
        Optional Ribo-seq-side condition assignments. Required for
        computing an INTERNAL Ribo-seq log2FC from replicate counts.
        When omitted (single-replicate / no condition labels), the
        delta-TE rows carry only the mRNA log2FC and a note.

    Returns
    -------
    list[DTeRow]
        One row per gene that appears in BOTH ``ribo_counts`` and the
        DE table.
    """
    mrna_log2fc_by_gene = {row["gene_id"]: row["log2fc"] for row in de_table.rows}
    padj_by_gene = {row["gene_id"]: row["padj"] for row in de_table.rows}

    have_replicates = (
        condition_map is not None
        and condition_a is not None
        and condition_b is not None
    )

    grouped: dict[str, dict[str, list[int]]] = (
        _group_counts_by_condition(ribo_counts, condition_map)
        if have_replicates
        else {}
    )

    rows: list[DTeRow] = []
    for gene in sorted(ribo_counts):
        mrna_log2fc = mrna_log2fc_by_gene.get(gene)
        if mrna_log2fc is None:
            rows.append(
                DTeRow(
                    gene=gene,
                    mrna_log2fc=None,
                    rpf_log2fc=None,
                    delta_te_log2=None,
                    padj=None,
                    note="missing_from_de_table",
                )
            )
            continue

        rpf_log2fc: float | None = None
        note = ""
        if have_replicates:
            cond_counts = grouped.get(gene, {})
            a_counts = cond_counts.get(condition_a, [])
            b_counts = cond_counts.get(condition_b, [])
            if a_counts and b_counts:
                mean_a = sum(a_counts) / len(a_counts) + _PSEUDOCOUNT
                mean_b = sum(b_counts) / len(b_counts) + _PSEUDOCOUNT
                ratio = mean_b / mean_a
                rpf_log2fc = _safe_log2(ratio)
            else:
                note = "insufficient_ribo_replicates"
        else:
            note = "single_replicate_no_statistics"

        delta_te = (
            (rpf_log2fc - mrna_log2fc) if rpf_log2fc is not None else None
        )
        rows.append(
            DTeRow(
                gene=gene,
                mrna_log2fc=mrna_log2fc,
                rpf_log2fc=rpf_log2fc,
                delta_te_log2=delta_te,
                padj=padj_by_gene.get(gene),
                note=note,
            )
        )
    return rows
