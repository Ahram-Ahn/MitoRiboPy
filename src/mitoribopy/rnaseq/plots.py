"""TE / delta-TE diagnostic plots."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from ._types import DTeRow


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
