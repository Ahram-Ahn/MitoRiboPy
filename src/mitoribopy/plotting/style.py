"""Shared plotting style used by package analysis modules."""

from __future__ import annotations

import matplotlib as mpl


PUBLICATION_STYLE = {
    "font.family": "Arial",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 1.0,
    "axes.labelsize": 12,
    "axes.titlesize": 14,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 300,
}


def apply_publication_style() -> None:
    """Set consistent matplotlib defaults across all analysis modules."""
    mpl.rcParams.update(PUBLICATION_STYLE)
