"""Tests for the v0.9.0 from_fastq padj nullification.

Wald inference on the 13-mt-mRNA subset is unreliable; the from_fastq
flow now always sets padj=None on the de_table.tsv it generates so a
downstream user cannot accidentally cite "padj < 0.05" from a 13-gene
DESeq2 fit. The publication-grade flow (rnaseq_mode=de_table) is
unaffected because it consumes a user-supplied external DE table.
"""

from __future__ import annotations

import pandas as pd

from mitoribopy.rnaseq.de_analysis import deseq2_to_de_table


def _mock_results_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "log2FoldChange": [-1.5, 0.5, 2.0],
            "padj": [0.001, 0.04, 0.5],
            "baseMean": [100.0, 50.0, 200.0],
        },
        index=["MT-CO1", "MT-ND1", "MT-ATP6"],
    )


def test_default_keeps_padj() -> None:
    """The publication-grade external-table path must keep padj as-is."""
    df = _mock_results_df()
    table = deseq2_to_de_table(df)
    by_gene = {row["gene_id"]: row for row in table.rows}
    assert by_gene["MT-CO1"]["padj"] == 0.001
    assert by_gene["MT-ND1"]["padj"] == 0.04
    assert by_gene["MT-ATP6"]["padj"] == 0.5


def test_nullify_padj_drops_pvalues_but_keeps_log2fc() -> None:
    df = _mock_results_df()
    table = deseq2_to_de_table(df, nullify_padj=True)
    by_gene = {row["gene_id"]: row for row in table.rows}
    # padj is gone for every row.
    for row in table.rows:
        assert row["padj"] is None, f"padj not nulled for {row['gene_id']}"
    # log2FoldChange and baseMean are preserved verbatim.
    assert by_gene["MT-CO1"]["log2fc"] == -1.5
    assert by_gene["MT-CO1"]["basemean"] == 100.0
    assert by_gene["MT-ND1"]["log2fc"] == 0.5
    assert by_gene["MT-ATP6"]["log2fc"] == 2.0


def test_nullify_padj_works_when_input_padj_is_already_nan() -> None:
    df = pd.DataFrame(
        {
            "log2FoldChange": [1.0],
            "padj": [float("nan")],
            "baseMean": [10.0],
        },
        index=["MT-ND6"],
    )
    table = deseq2_to_de_table(df, nullify_padj=True)
    assert table.rows[0]["padj"] is None
    assert table.rows[0]["log2fc"] == 1.0
