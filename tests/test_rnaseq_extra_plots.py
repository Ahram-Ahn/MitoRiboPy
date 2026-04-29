"""Smoke tests for the additional TE / DE / QC plots.

These tests confirm each new plot function emits a non-empty PNG with
plausible dimensions and recognisable PNG headers. Visual correctness
is impractical to assert in unit tests; the goal here is catching
regressions in input-shape handling, label rendering, and matplotlib
backend issues.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from mitoribopy.rnaseq._types import (
    DE_COLUMN_ALIASES,
    DTeRow,
    DeTable,
    TeRow,
)
from mitoribopy.rnaseq.plots import (
    plot_de_volcano,
    plot_ma,
    plot_pca,
    plot_te_bar_grouped,
    plot_te_compare_scatter,
    plot_te_heatmap,
    plot_te_log2fc_bar,
)


_PNG_HEADER = b"\x89PNG\r\n\x1a\n"


def _is_png(path: Path) -> bool:
    return path.is_file() and path.read_bytes()[:8] == _PNG_HEADER


def _make_de_table() -> DeTable:
    return DeTable(
        format="deseq2",
        column_map=DE_COLUMN_ALIASES["deseq2"],
        rows=[
            {"gene_id": "MT-ND1", "log2fc": 0.4, "padj": 0.01, "basemean": 500.0},
            {"gene_id": "MT-CO1", "log2fc": -1.2, "padj": 1e-30, "basemean": 5000.0},
            {"gene_id": "MT-ND6", "log2fc": 1.7, "padj": 0.5, "basemean": 50.0},
            {"gene_id": "MT-CYB", "log2fc": None, "padj": None, "basemean": None},
        ],
    )


def test_plot_ma_emits_png(tmp_path: Path) -> None:
    out = tmp_path / "ma.png"
    plot_ma(_make_de_table(), out)
    assert _is_png(out)
    assert out.stat().st_size > 1000  # not just a stub


def _make_te_rows() -> list[TeRow]:
    rows: list[TeRow] = []
    for sample, te_factor in (
        ("WT_rep1", 1.0), ("WT_rep2", 1.05),
        ("KO_rep1", 0.4), ("KO_rep2", 0.35),
    ):
        for gene, base in (("ND1", 1.5), ("CO1", 6.0), ("ND6", 0.8)):
            rows.append(
                TeRow(
                    sample=sample,
                    gene=gene,
                    rpf_count=int(base * 100),
                    mrna_abundance=base * 50.0,
                    te=base * te_factor,
                )
            )
    return rows


def _make_condition_map() -> dict[str, str]:
    return {
        "WT_rep1": "WT", "WT_rep2": "WT",
        "KO_rep1": "KO", "KO_rep2": "KO",
    }


def test_plot_te_bar_grouped_emits_png(tmp_path: Path) -> None:
    out = tmp_path / "bar.png"
    plot_te_bar_grouped(_make_te_rows(), _make_condition_map(), out)
    assert _is_png(out)


def test_plot_te_bar_grouped_handles_empty(tmp_path: Path) -> None:
    out = tmp_path / "empty.png"
    plot_te_bar_grouped([], {}, out)
    assert _is_png(out)


def test_plot_te_heatmap_emits_png(tmp_path: Path) -> None:
    out = tmp_path / "heat.png"
    plot_te_heatmap(_make_te_rows(), _make_condition_map(), out)
    assert _is_png(out)


def test_plot_te_heatmap_handles_empty(tmp_path: Path) -> None:
    out = tmp_path / "empty_heat.png"
    plot_te_heatmap([], None, out)
    assert _is_png(out)


def _make_pca_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.DataFrame(
        {
            "ND1": [10, 12, 80, 75, 8, 11],
            "CO1": [200, 220, 180, 175, 210, 205],
            "ND6": [5, 6, 4, 5, 50, 55],
        },
        index=["WT_rep1", "WT_rep2", "KO_rep1", "KO_rep2", "RESC_rep1", "RESC_rep2"],
    )
    meta = pd.DataFrame(
        {
            "assay": ["rna"] * 6,
            "condition": ["WT", "WT", "KO", "KO", "RESC", "RESC"],
        },
        index=counts.index,
    )
    return counts, meta


def test_plot_pca_emits_png(tmp_path: Path) -> None:
    counts, meta = _make_pca_inputs()
    out = tmp_path / "pca.png"
    plot_pca(counts, meta, out)
    assert _is_png(out)


def test_plot_pca_single_sample_falls_back(tmp_path: Path) -> None:
    counts = pd.DataFrame({"ND1": [10], "CO1": [200]}, index=["only"])
    meta = pd.DataFrame({"assay": ["rna"], "condition": ["WT"]}, index=["only"])
    out = tmp_path / "pca_one.png"
    plot_pca(counts, meta, out)
    # Must still produce a real PNG (placeholder, not a crash).
    assert _is_png(out)


def test_plot_de_volcano_emits_png(tmp_path: Path) -> None:
    out = tmp_path / "de_volcano.png"
    plot_de_volcano(_make_de_table(), out, contrast_label="WT vs KO")
    assert _is_png(out)
    assert out.stat().st_size > 1000


def test_plot_de_volcano_handles_padj_zero(tmp_path: Path) -> None:
    """A row with ``padj == 0`` (numerical underflow) must be drawn at the
    largest finite -log10(padj) on the figure, not at infinity."""
    de = DeTable(
        format="deseq2",
        column_map=DE_COLUMN_ALIASES["deseq2"],
        rows=[
            {"gene_id": "ND1", "log2fc": 0.4, "padj": 0.01, "basemean": 500.0},
            {"gene_id": "CO1", "log2fc": -1.2, "padj": 0.0, "basemean": 5000.0},
            {"gene_id": "ND6", "log2fc": 1.7, "padj": 0.5, "basemean": 50.0},
        ],
    )
    out = tmp_path / "de_volcano_zero.png"
    plot_de_volcano(de, out, contrast_label="WT vs KO")
    assert _is_png(out)


def test_plot_de_volcano_handles_empty(tmp_path: Path) -> None:
    de = DeTable(
        format="deseq2",
        column_map=DE_COLUMN_ALIASES["deseq2"],
        rows=[],
    )
    out = tmp_path / "de_volcano_empty.png"
    plot_de_volcano(de, out)
    assert _is_png(out)


def test_plot_te_compare_scatter_emits_png(tmp_path: Path) -> None:
    out = tmp_path / "te_compare.png"
    plot_te_compare_scatter(
        _make_te_rows(), _make_condition_map(), "WT", "KO", out
    )
    assert _is_png(out)


def test_plot_te_compare_scatter_handles_no_shared_genes(tmp_path: Path) -> None:
    """Empty cells from a condition with no replicates must not crash."""
    out = tmp_path / "te_compare_empty.png"
    plot_te_compare_scatter([], {}, "WT", "KO", out)
    assert _is_png(out)


def test_plot_te_log2fc_bar_emits_png(tmp_path: Path) -> None:
    out = tmp_path / "te_log2fc.png"
    plot_te_log2fc_bar(
        _make_te_rows(), _make_condition_map(), "WT", "KO", out
    )
    assert _is_png(out)


def test_plot_te_log2fc_bar_handles_empty(tmp_path: Path) -> None:
    out = tmp_path / "te_log2fc_empty.png"
    plot_te_log2fc_bar([], {}, "WT", "KO", out)
    assert _is_png(out)
