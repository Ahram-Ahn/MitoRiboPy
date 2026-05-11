"""pyDESeq2 wrapper for the from-FASTQ rnaseq mode.

Builds a sample-sheet from per-sample alignment results, runs pyDESeq2
on the RNA-seq side, and serialises the result back into the canonical
:class:`~mitoribopy.rnaseq._types.DeTable` DESeq2 schema so the existing
TE / delta-TE / plot path can consume it via
:func:`mitoribopy.rnaseq.de_loader.load_de_table`.

pyDESeq2 is soft-imported. The package's existing pre-computed-DE flow
must keep working for users who do not install the ``[fastq]`` extra.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

from ._types import DE_COLUMN_ALIASES, DeTable

if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd

    from .alignment import SampleAlignmentResult


_FASTQ_EXTRA_HINT = (
    "pyDESeq2 is required for from-FASTQ mode. "
    "Install with `pip install 'mitoribopy[fastq]'`."
)


def _import_pandas():
    try:
        import pandas as pd
    except ImportError as exc:  # pragma: no cover - pandas is a hard dep
        raise RuntimeError(
            "pandas is required for from-FASTQ DE analysis."
        ) from exc
    return pd


def build_sample_sheet(
    rna_results: "Sequence[SampleAlignmentResult]",
    ribo_results: "Sequence[SampleAlignmentResult]",
    condition_map: dict[str, str],
) -> "tuple[list[str], pd.DataFrame, pd.DataFrame]":
    """Combine per-sample counts + a condition map into a DESeq2-ready sheet.

    Returns
    -------
    samples:
        Ordered sample names (RNA first, then Ribo).
    counts_df:
        rows = samples, columns = union of genes, integer counts,
        zero-filled.
    metadata_df:
        index = samples, columns ``assay`` (``"rna"`` or ``"ribo"``)
        and ``condition``. Both columns are pandas categorical-friendly
        strings — pyDESeq2 will make them categorical itself.

    Raises
    ------
    ValueError
        if any sample is missing from ``condition_map`` (the message
        names every missing sample), or if the resulting metadata has
        fewer than two distinct condition levels.
    """
    pd = _import_pandas()

    rna_results = list(rna_results)
    ribo_results = list(ribo_results)

    samples: list[str] = [r.sample for r in rna_results] + [r.sample for r in ribo_results]
    if len(set(samples)) != len(samples):
        dups = [s for s in samples if samples.count(s) > 1]
        raise ValueError(
            "Sample names collide between RNA-seq and Ribo-seq sides: "
            + ", ".join(sorted(set(dups)))
        )

    missing = [s for s in samples if s not in condition_map]
    if missing:
        raise ValueError(
            "Condition map is missing entries for sample(s): "
            + ", ".join(missing)
        )

    genes: set[str] = set()
    for r in rna_results + ribo_results:
        genes.update(r.counts.keys())
    sorted_genes = sorted(genes)

    counts_data: dict[str, list[int]] = {gene: [] for gene in sorted_genes}
    for r in rna_results + ribo_results:
        for gene in sorted_genes:
            counts_data[gene].append(int(r.counts.get(gene, 0)))

    counts_df = pd.DataFrame(counts_data, index=samples, dtype=int)

    assays: list[str] = ["rna"] * len(rna_results) + ["ribo"] * len(ribo_results)
    conditions: list[str] = [condition_map[s] for s in samples]
    metadata_df = pd.DataFrame(
        {"assay": assays, "condition": conditions},
        index=samples,
    )

    if metadata_df["condition"].nunique() < 2:
        raise ValueError(
            "Condition map yields fewer than two condition levels; "
            "pyDESeq2 needs at least two levels for a contrast. "
            f"Found: {sorted(metadata_df['condition'].unique().tolist())}"
        )

    return samples, counts_df, metadata_df


def run_deseq2(
    counts_df: "pd.DataFrame",
    metadata_df: "pd.DataFrame",
    *,
    contrast_factor: str,
    contrast_a: str,
    contrast_b: str,
    additional_factors: tuple[str, ...] = (),
    assay: str = "rna",
) -> "pd.DataFrame":
    """Fit pyDESeq2 on one assay subset and return the contrast results.

    When ``metadata_df`` has an ``assay`` column, rows are filtered to
    ``assay == <assay>`` before fitting (default ``"rna"``); pass
    ``assay="ribo"`` to fit the Ribo-seq side. The resulting DataFrame
    is the verbatim ``DeseqStats.results_df`` (DESeq2 column names:
    ``baseMean``, ``log2FoldChange``, ``padj``, ...).
    """
    _import_pandas()

    try:
        from pydeseq2.dds import DeseqDataSet  # type: ignore
        from pydeseq2.ds import DeseqStats  # type: ignore
    except ImportError as exc:
        raise RuntimeError(_FASTQ_EXTRA_HINT) from exc

    if "assay" in metadata_df.columns:
        mask = metadata_df["assay"] == assay
        sub_metadata = metadata_df.loc[mask].drop(columns=["assay"])
        sub_counts = counts_df.loc[mask]
    else:
        sub_metadata = metadata_df
        sub_counts = counts_df

    if sub_counts.empty:
        raise ValueError(f"No {assay!r} samples available for DESeq2.")

    if sub_metadata[contrast_factor].nunique() < 2:
        raise ValueError(
            f"{assay!r} subset has fewer than two levels of {contrast_factor!r}: "
            f"{sorted(sub_metadata[contrast_factor].unique().tolist())}"
        )

    design_factors: list[str] = list(additional_factors) + [contrast_factor]

    dds = DeseqDataSet(
        counts=sub_counts,
        metadata=sub_metadata,
        design_factors=design_factors,
    )
    dds.deseq2()
    ds = DeseqStats(dds, contrast=[contrast_factor, contrast_a, contrast_b])
    ds.summary()
    results_df = ds.results_df
    # pyDESeq2 returns the gene index as the row index; preserve it.
    return results_df


def _safe_value(value: Any) -> float | None:
    if value is None:
        return None
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(f) or math.isinf(f):
        return None
    return f


def deseq2_to_de_table(
    results_df: "pd.DataFrame",
    *,
    nullify_padj: bool = False,
) -> DeTable:
    """Convert a pyDESeq2 ``results_df`` to a canonical :class:`DeTable`.

    NaN / ±Inf are mapped to ``None`` (mirrors
    :func:`mitoribopy.rnaseq.de_loader._safe_float`).

    Parameters
    ----------
    nullify_padj:
        When ``True``, the ``padj`` column on every output row is set
        to ``None`` regardless of the input value. The intended caller
        is the ``rnaseq_mode=from_fastq`` orchestrator, which fits
        DESeq2 on the 13-mt-mRNA subset — too few genes for the
        dispersion-shrinkage estimator to be reliable. We keep the
        log2FoldChange (a reasonable point estimate even with low gene
        count) but suppress the Wald p-value so a downstream user
        cannot accidentally cite "padj < 0.05" from a 13-gene fit. The
        publication-grade flow uses ``rnaseq_mode=de_table`` with an
        external full-transcriptome DE table where padj is meaningful.
    """
    column_map = DE_COLUMN_ALIASES["deseq2"]
    rows: list[dict] = []
    for gene_id, row in results_df.iterrows():
        padj_value = None if nullify_padj else _safe_value(row.get("padj"))
        rows.append(
            {
                "gene_id": str(gene_id),
                "log2fc": _safe_value(row.get("log2FoldChange")),
                "padj": padj_value,
                "basemean": _safe_value(row.get("baseMean")),
            }
        )
    return DeTable(format="deseq2", column_map=column_map, rows=rows)


def write_de_table_tsv(de_table: DeTable, path: Path) -> None:
    """Write ``de_table`` to a DESeq2-shaped TSV.

    The header is exactly ``gene_id\\tbaseMean\\tlog2FoldChange\\tpadj``,
    which :func:`mitoribopy.rnaseq.de_loader.detect_de_format` recognises
    as DESeq2 because it sees the ``log2FoldChange`` column. ``None`` is
    serialised as ``NA``.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    def _fmt(v: float | None) -> str:
        if v is None:
            return "NA"
        return f"{v:.6g}"

    with path.open("w", encoding="utf-8") as handle:
        handle.write("gene_id\tbaseMean\tlog2FoldChange\tpadj\n")
        for row in de_table.rows:
            handle.write(
                f"{row['gene_id']}\t"
                f"{_fmt(row.get('basemean'))}\t"
                f"{_fmt(row.get('log2fc'))}\t"
                f"{_fmt(row.get('padj'))}\n"
            )
