"""TE / delta-TE diagnostic plots."""

from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING, Iterable, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from ._types import DeTable, DTeRow, TeRow

if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd


def plot_mrna_vs_rpf_scatter(
    rows: Iterable[DTeRow],
    output_path: Path,
    *,
    title: str = "log2FC RPF vs mRNA",
) -> Path:
    """Scatter log2FC_mRNA (x) vs log2FC_RPF (y); label the 4 quadrants.

    Q1 upper-right: translation AND transcription up (co-regulated up)
    Q2 upper-left:  translation up, transcription down (buffered-up)
    Q3 lower-left:  both down (co-regulated down)
    Q4 lower-right: translation down, transcription up (buffered-down)
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    xs: list[float] = []
    ys: list[float] = []
    labels: list[str] = []
    for row in rows:
        if row.mrna_log2fc is None or row.rpf_log2fc is None:
            continue
        xs.append(row.mrna_log2fc)
        ys.append(row.rpf_log2fc)
        labels.append(row.gene)

    fig, ax = plt.subplots(figsize=(5, 5))
    if xs:
        ax.scatter(xs, ys, s=40, alpha=0.85, color="#2a7ab0")
        for x, y, label in zip(xs, ys, labels):
            ax.annotate(label, (x, y), fontsize=7, xytext=(3, 3), textcoords="offset points")

    ax.axhline(0, color="grey", linewidth=0.6)
    ax.axvline(0, color="grey", linewidth=0.6)
    ax.set_xlabel("log2FC mRNA (RNA-seq)")
    ax.set_ylabel("log2FC RPF (Ribo-seq)")
    ax.set_title(title)
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return output_path


def plot_delta_te_volcano(
    rows: Iterable[DTeRow],
    output_path: Path,
    *,
    title: str = "Delta-TE (log2) volcano",
) -> Path:
    """Volcano-style scatter of delta_te_log2 (x) vs -log10(padj) (y).

    When the DE table has no padj column we fall back to plotting
    delta-TE along x with y fixed at 0 so the file still exists and
    reviewers see the message in the saved legend.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    xs: list[float] = []
    ys: list[float] = []
    labels: list[str] = []
    for row in rows:
        if row.delta_te_log2 is None:
            continue
        xs.append(row.delta_te_log2)
        padj = row.padj
        if padj is not None and padj > 0:
            ys.append(-math.log10(padj))
        else:
            ys.append(0.0)
        labels.append(row.gene)

    fig, ax = plt.subplots(figsize=(5, 5))
    if xs:
        ax.scatter(xs, ys, s=40, alpha=0.85, color="#b85c00")
        for x, y, label in zip(xs, ys, labels):
            ax.annotate(label, (x, y), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.axvline(0, color="grey", linewidth=0.6)
    ax.set_xlabel("Delta-TE (log2)")
    ax.set_ylabel("-log10(padj)")
    ax.set_title(title)
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return output_path


def plot_ma(
    de_table: DeTable,
    output_path: Path,
    *,
    padj_threshold: float = 0.05,
    title: str = "MA plot (RNA-seq DE)",
) -> Path:
    """Classic DESeq2 MA plot: mean expression (x) vs log2FoldChange (y).

    Points with ``padj < padj_threshold`` are highlighted in red; the
    rest are grey. A horizontal line at y=0 marks no change.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    xs_sig: list[float] = []
    ys_sig: list[float] = []
    xs_ns: list[float] = []
    ys_ns: list[float] = []
    labels_sig: list[tuple[float, float, str]] = []

    for row in de_table.rows:
        bm = row.get("basemean")
        l2 = row.get("log2fc")
        if bm is None or l2 is None or bm <= 0:
            continue
        x = math.log10(bm)
        padj = row.get("padj")
        is_sig = padj is not None and padj < padj_threshold
        if is_sig:
            xs_sig.append(x)
            ys_sig.append(l2)
            labels_sig.append((x, l2, str(row.get("gene_id", ""))))
        else:
            xs_ns.append(x)
            ys_ns.append(l2)

    fig, ax = plt.subplots(figsize=(5, 5))
    if xs_ns:
        ax.scatter(xs_ns, ys_ns, s=24, alpha=0.65, color="#888888",
                   label=f"padj ≥ {padj_threshold}")
    if xs_sig:
        ax.scatter(xs_sig, ys_sig, s=40, alpha=0.9, color="#c0392b",
                   label=f"padj < {padj_threshold}")
        for x, y, label in labels_sig:
            ax.annotate(label, (x, y), fontsize=7, xytext=(3, 3),
                        textcoords="offset points")
    ax.axhline(0, color="grey", linewidth=0.6)
    ax.set_xlabel("log10(baseMean)")
    ax.set_ylabel("log2FoldChange")
    ax.set_title(title)
    ax.grid(True, alpha=0.2)
    if xs_ns or xs_sig:
        ax.legend(loc="best", fontsize=8, framealpha=0.85)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return output_path


def plot_te_bar_grouped(
    te_rows: Iterable[TeRow],
    condition_map: dict[str, str],
    output_path: Path,
    *,
    title: str = "Translation efficiency per gene by condition",
    log2: bool = True,
) -> Path:
    """Per-gene TE bar plot, bars grouped by condition with SE error bars.

    Replicates within each (gene, condition) cell are averaged; error
    bars show the standard error of the mean. ``log2=True`` plots
    ``log2(TE)`` so the y-axis is centred on 0 = no change; pass
    ``log2=False`` for raw TE.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # gene -> condition -> list[TE]
    grouped: dict[str, dict[str, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for row in te_rows:
        cond = condition_map.get(row.sample)
        if cond is None or row.te <= 0:
            continue
        value = math.log2(row.te) if log2 else row.te
        grouped[row.gene][cond].append(value)

    if not grouped:
        # Degenerate; emit an empty figure so the caller still gets a path.
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, "no TE rows with a condition assignment",
                ha="center", va="center")
        ax.axis("off")
        fig.savefig(output_path, dpi=150)
        plt.close(fig)
        return output_path

    genes = sorted(grouped.keys())
    conditions = sorted({c for d in grouped.values() for c in d.keys()})
    n_genes = len(genes)
    n_conds = len(conditions)
    width = 0.8 / max(n_conds, 1)
    x_indices = np.arange(n_genes)

    cmap = plt.get_cmap("tab10")
    fig, ax = plt.subplots(figsize=(max(6, 0.7 * n_genes + 2), 4.5))
    for i, cond in enumerate(conditions):
        means: list[float] = []
        ses: list[float] = []
        for gene in genes:
            values = grouped[gene].get(cond, [])
            if values:
                arr = np.array(values, dtype=float)
                means.append(float(arr.mean()))
                ses.append(
                    float(arr.std(ddof=1) / math.sqrt(len(arr)))
                    if len(arr) > 1
                    else 0.0
                )
            else:
                means.append(0.0)
                ses.append(0.0)
        offset = (i - (n_conds - 1) / 2) * width
        ax.bar(
            x_indices + offset, means, width=width,
            yerr=ses, capsize=3, label=cond,
            color=cmap(i % 10), alpha=0.9, edgecolor="black", linewidth=0.4,
        )

    ax.axhline(0, color="grey", linewidth=0.6)
    ax.set_xticks(x_indices)
    ax.set_xticklabels(genes, rotation=45, ha="right")
    ax.set_ylabel("log2(TE)" if log2 else "TE")
    ax.set_title(title)
    ax.grid(True, axis="y", alpha=0.2)
    ax.legend(loc="best", fontsize=8, framealpha=0.85)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return output_path


def plot_te_heatmap(
    te_rows: Iterable[TeRow],
    condition_map: dict[str, str] | None,
    output_path: Path,
    *,
    title: str = "log2(TE) heatmap (gene × sample)",
) -> Path:
    """Gene × sample heatmap of log2(TE).

    Sample columns are sorted by condition (from ``condition_map``) then
    by name so replicates within a condition cluster together. Values
    are annotated when the matrix is small (≤ 200 cells).
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = [r for r in te_rows if r.te > 0]
    if not rows:
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, "no TE rows to plot", ha="center", va="center")
        ax.axis("off")
        fig.savefig(output_path, dpi=150)
        plt.close(fig)
        return output_path

    samples = sorted(
        {r.sample for r in rows},
        key=lambda s: ((condition_map or {}).get(s, ""), s),
    )
    genes = sorted({r.gene for r in rows})
    matrix = np.full((len(genes), len(samples)), np.nan, dtype=float)
    gene_idx = {g: i for i, g in enumerate(genes)}
    sample_idx = {s: j for j, s in enumerate(samples)}
    for r in rows:
        matrix[gene_idx[r.gene], sample_idx[r.sample]] = math.log2(r.te)

    abs_max = float(np.nanmax(np.abs(matrix))) if np.any(~np.isnan(matrix)) else 1.0
    abs_max = max(abs_max, 0.5)

    fig, ax = plt.subplots(
        figsize=(max(5, 0.55 * len(samples) + 2), 0.45 * len(genes) + 2)
    )
    im = ax.imshow(
        matrix, aspect="auto", cmap="RdBu_r",
        vmin=-abs_max, vmax=abs_max,
    )

    ax.set_xticks(np.arange(len(samples)))
    ax.set_yticks(np.arange(len(genes)))
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(genes, fontsize=9)
    ax.set_title(title)

    # Cell annotations only when not too dense.
    if matrix.size <= 200:
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                v = matrix[i, j]
                if math.isnan(v):
                    continue
                color = "white" if abs(v) > 0.6 * abs_max else "black"
                ax.text(j, i, f"{v:.2f}",
                        ha="center", va="center", fontsize=7, color=color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("log2(TE)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return output_path


def plot_pca(
    counts_df: "pd.DataFrame",
    metadata_df: "pd.DataFrame",
    output_path: Path,
    *,
    title: str = "Sample PCA",
) -> Path:
    """PC1 vs PC2 of samples on log1p-transformed counts.

    Color = ``condition``, marker shape = ``assay`` (when present).
    Falls back to a placeholder figure when fewer than 2 samples are
    available (PCA is undefined). PCA is computed via numpy SVD so we
    do not pull in scikit-learn.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    samples = list(counts_df.index)
    n = len(samples)
    if n < 2:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, f"PCA needs ≥2 samples (got {n})",
                ha="center", va="center")
        ax.axis("off")
        fig.savefig(output_path, dpi=150)
        plt.close(fig)
        return output_path

    X = np.log1p(counts_df.to_numpy(dtype=float))
    # Drop genes that are zero in every sample (no variance contribution).
    keep = X.sum(axis=0) > 0
    X = X[:, keep]
    if X.shape[1] == 0:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "PCA: no genes with non-zero counts",
                ha="center", va="center")
        ax.axis("off")
        fig.savefig(output_path, dpi=150)
        plt.close(fig)
        return output_path

    X_centered = X - X.mean(axis=0, keepdims=True)
    # Truncated SVD of the centered matrix.
    U, S, _Vt = np.linalg.svd(X_centered, full_matrices=False)
    scores = U * S  # (n_samples, n_components)
    var = (S ** 2) / max(n - 1, 1)
    var_ratio = var / var.sum() if var.sum() > 0 else np.zeros_like(var)
    pc1 = scores[:, 0] if scores.shape[1] >= 1 else np.zeros(n)
    pc2 = scores[:, 1] if scores.shape[1] >= 2 else np.zeros(n)

    conditions = (
        list(metadata_df.loc[samples, "condition"])
        if "condition" in metadata_df.columns
        else ["?"] * n
    )
    assays = (
        list(metadata_df.loc[samples, "assay"])
        if "assay" in metadata_df.columns
        else ["?"] * n
    )
    unique_conds = sorted(set(conditions))
    unique_assays = sorted(set(assays))
    cmap = plt.get_cmap("tab10")
    cond_color = {c: cmap(i % 10) for i, c in enumerate(unique_conds)}
    assay_marker = {a: m for a, m in zip(unique_assays, ["o", "s", "^", "D", "P"])}

    fig, ax = plt.subplots(figsize=(6, 5))
    for s, x, y, c, a in zip(samples, pc1, pc2, conditions, assays):
        ax.scatter(
            x, y, s=70, alpha=0.9,
            color=cond_color[c],
            marker=assay_marker.get(a, "o"),
            edgecolor="black", linewidth=0.4,
        )
        ax.annotate(s, (x, y), fontsize=7, xytext=(4, 4),
                    textcoords="offset points")

    # Legend with condition colours and assay markers (deduped).
    cond_handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=cond_color[c], markersize=8,
                   label=f"condition: {c}", markeredgecolor="black")
        for c in unique_conds
    ]
    assay_handles = [
        plt.Line2D([0], [0], marker=assay_marker.get(a, "o"), color="w",
                   markerfacecolor="#aaaaaa", markersize=8,
                   label=f"assay: {a}", markeredgecolor="black")
        for a in unique_assays
    ]
    if len(unique_assays) > 1:
        ax.legend(handles=cond_handles + assay_handles,
                  loc="best", fontsize=8, framealpha=0.85)
    else:
        ax.legend(handles=cond_handles, loc="best", fontsize=8, framealpha=0.85)

    ax.axhline(0, color="grey", linewidth=0.4)
    ax.axvline(0, color="grey", linewidth=0.4)
    ax.set_xlabel(f"PC1 ({var_ratio[0] * 100:.1f}%)" if len(var_ratio) >= 1 else "PC1")
    ax.set_ylabel(f"PC2 ({var_ratio[1] * 100:.1f}%)" if len(var_ratio) >= 2 else "PC2")
    ax.set_title(title)
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return output_path
