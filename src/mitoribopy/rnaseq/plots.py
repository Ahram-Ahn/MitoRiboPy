"""TE / delta-TE diagnostic plots, publication-styled.

Every figure in this module is drawn under a shared :data:`_PUB_RC`
matplotlib rc-context (sans-serif fonts, top/right spines off,
SVG-as-text, 300 dpi) and saved through :func:`_save_figure`, which
emits both a PNG (for slides / Markdown) and an SVG sidecar (for
Illustrator / publication figures), mirroring the rpf
codon-correlation plot pattern.

Colour palette is the Okabe-Ito colour-blind-safe set (see
https://jfly.uni-koeln.de/color/), with consistent semantics:

* ``_OKABE_ITO['vermillion']`` = up-regulated / above identity
* ``_OKABE_ITO['blue']``       = down-regulated / below identity
* ``_OKABE_ITO['grey']``       = non-significant / reference
"""

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


# ---------------------------------------------------------------------------
# Publication style
# ---------------------------------------------------------------------------


_OKABE_ITO: dict[str, str] = {
    # Colour-blind-safe palette (Okabe & Ito 2002). 8 hues; mapped to
    # the meanings the rnaseq plots use repeatedly.
    "black":      "#000000",
    "orange":     "#E69F00",
    "skyblue":    "#56B4E9",
    "green":      "#009E73",
    "yellow":     "#F0E442",
    "blue":       "#0072B2",
    "vermillion": "#D55E00",
    "purple":     "#CC79A7",
    "grey":       "#999999",
}

# Semantic aliases. Keep the volcano / MA / TE-direction axis consistent:
# up moves in vermillion, down moves in blue, n.s. is light grey.
_C_UP = _OKABE_ITO["vermillion"]
_C_DOWN = _OKABE_ITO["blue"]
_C_NS = "#bdbdbd"
_C_NEUTRAL = _OKABE_ITO["grey"]
_C_GUIDE = "#666666"

_LABEL_BBOX = dict(
    boxstyle="round,pad=0.18",
    facecolor="white",
    edgecolor="none",
    alpha=0.78,
)
_LEADER = dict(arrowstyle="-", color=_C_GUIDE, lw=0.4, shrinkA=0, shrinkB=2)

_PUB_RC: dict[str, object] = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "axes.titleweight": "regular",
    "axes.linewidth": 1.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "axes.axisbelow": True,
    "xtick.major.size": 4.5,
    "ytick.major.size": 4.5,
    "xtick.minor.size": 2.5,
    "ytick.minor.size": 2.5,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.frameon": True,
    "legend.framealpha": 0.92,
    "legend.fontsize": 9,
    "legend.edgecolor": "#cccccc",
    "svg.fonttype": "none",   # keep text editable in Illustrator
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.12,
    "savefig.transparent": False,
    "figure.dpi": 100,        # screen draft; savefig.dpi controls files
}


def _save_figure(
    fig: "plt.Figure",
    output_path: Path,
    *,
    formats: Sequence[str] = ("png", "svg"),
    metadata: dict | None = None,
) -> Path:
    """Save *fig* to *output_path* (whose suffix decides the primary
    format) and a sibling for every other entry in *formats*.

    The path returned is always *output_path* — the primary file the
    caller asked for. SVG sidecars (when emitted) get the same stem
    with the ``.svg`` suffix so downstream LaTeX / Markdown can pick
    whichever format they prefer without renaming.

    When *metadata* is supplied, also write a ``.metadata.json`` sidecar
    next to the plot so :mod:`mitoribopy.plotting.figure_validator` can
    mechanically validate the contract (point counts, label policy,
    palette, etc.).
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Always write the file the caller named first. Use 300 dpi for
    # PNG so slides / Markdown render crisply; SVG is vector.
    # bbox_inches="tight" keeps any out-of-axes legend in the saved
    # figure instead of clipping it.
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    primary_suffix = output_path.suffix.lower().lstrip(".")
    for fmt in formats:
        fmt = fmt.lower()
        if fmt == primary_suffix:
            continue
        sidecar = output_path.with_suffix(f".{fmt}")
        fig.savefig(sidecar, format=fmt, bbox_inches="tight")
    plt.close(fig)

    if metadata is not None:
        # Local import keeps figure_validator out of the matplotlib
        # render-time import graph.
        from ..plotting.figure_validator import write_plot_metadata

        write_plot_metadata(
            output_path,
            stage=metadata.get("stage", "rnaseq"),
            plot_type=metadata.get("plot_type", ""),
            source_data=metadata.get("source_data"),
            n_points_expected=metadata.get("n_points_expected"),
            n_points_drawn=metadata.get("n_points_drawn"),
            n_labels=metadata.get("n_labels"),
            labels_drawn=metadata.get("labels_drawn"),
            label_policy=metadata.get("label_policy"),
            label_overlap_count=metadata.get("label_overlap_count"),
            label_outside_axes_count=metadata.get("label_outside_axes_count"),
            palette=metadata.get("palette", "Okabe-Ito"),
            formats=metadata.get("formats", list(formats)),
            dpi=metadata.get("dpi", 300),
            min_font_size=metadata.get("min_font_size"),
            extra=metadata.get("extra"),
        )
    return output_path


def _annotate_with_bbox(
    ax: "plt.Axes",
    points: Iterable[tuple[float, float, str]],
    *,
    fontsize: int = 8,
    offset: tuple[float, float] = (5.0, 5.0),
    leader: bool = True,
    smart: bool = False,
    distance: float = 12.0,
) -> None:
    """Annotate *points* with white-bbox labels and (optional) leader lines.

    *points* is an iterable of ``(x, y, label)`` tuples. The bbox
    background keeps labels legible when they cross other elements;
    the short leader line makes the point-to-label correspondence
    unambiguous on dense plots.

    When ``smart=True``, each label is placed in the offset direction
    (one of 8 cardinal/diagonal positions ``distance`` pixels from the
    point) that minimises overlap with already-placed labels — a cheap
    greedy alternative to ``adjustText`` that handles the small
    mt-mRNA universe (~13 points) without dragging in a dependency.
    Unsuitable for whole-transcriptome inputs where points and labels
    number in the thousands; pass ``smart=False`` (default) there.
    """
    arrow = _LEADER if leader else None
    pts = list(points)

    if not smart:
        for x, y, label in pts:
            ax.annotate(
                label,
                xy=(x, y),
                xytext=offset,
                textcoords="offset points",
                fontsize=fontsize,
                ha="left",
                va="bottom",
                bbox=_LABEL_BBOX,
                arrowprops=arrow,
            )
        return

    # Smart placement: try 8 candidate directions and pick the one
    # whose estimated bbox overlaps the least with previously placed
    # labels. Display-coordinate maths so we get pixel-accurate
    # overlap regardless of axis aspect / scale.
    fig = ax.figure
    fig.canvas.draw()  # force layout so transData is final

    # Candidates ordered top-first so the greedy `if score == 0: break`
    # early-exit prefers placing labels ABOVE their points. Below-the-
    # point positions are only chosen when every above / side
    # alternative would overlap an already-placed label.
    candidates: list[tuple[float, float, str, str]] = []
    for angle_deg in (90, 45, 135, 0, 180, -45, -135, -90):
        rad = math.radians(angle_deg)
        dx = distance * math.cos(rad)
        dy = distance * math.sin(rad)
        ha = "left" if dx > 0.5 else ("right" if dx < -0.5 else "center")
        va = "bottom" if dy > 0.5 else ("top" if dy < -0.5 else "center")
        candidates.append((dx, dy, ha, va))

    placed: list[tuple[float, float, float, float]] = []
    for x, y, label in pts:
        # Pixel-rough estimate of the rendered bbox (matplotlib's true
        # bbox needs a draw round-trip per candidate; this is plenty
        # for greedy disambiguation on a ~13-gene set).
        text_w = max(1, len(str(label))) * fontsize * 0.62 + 4.0
        text_h = fontsize * 1.45 + 4.0
        try:
            disp_x, disp_y = ax.transData.transform((x, y))
        except Exception:
            disp_x, disp_y = (0.0, 0.0)

        best = candidates[0]
        best_score = float("inf")
        best_box = (disp_x, disp_y, disp_x + text_w, disp_y + text_h)
        for dx, dy, ha, va in candidates:
            cx = disp_x + dx
            cy = disp_y + dy
            if ha == "right":
                bx0, bx1 = cx - text_w, cx
            elif ha == "center":
                bx0, bx1 = cx - text_w / 2, cx + text_w / 2
            else:
                bx0, bx1 = cx, cx + text_w
            if va == "top":
                by0, by1 = cy - text_h, cy
            elif va == "center":
                by0, by1 = cy - text_h / 2, cy + text_h / 2
            else:
                by0, by1 = cy, cy + text_h
            score = 0.0
            for px0, py0, px1, py1 in placed:
                ow = max(0.0, min(bx1, px1) - max(bx0, px0))
                oh = max(0.0, min(by1, py1) - max(by0, py0))
                score += ow * oh
            if score < best_score:
                best_score = score
                best = (dx, dy, ha, va)
                best_box = (bx0, by0, bx1, by1)
                if score == 0.0:
                    break

        dx, dy, ha, va = best
        ax.annotate(
            label,
            xy=(x, y),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=fontsize,
            ha=ha, va=va,
            bbox=_LABEL_BBOX,
            arrowprops=arrow,
        )
        placed.append(best_box)


def _stat_box(ax: "plt.Axes", lines: Sequence[str], *, loc: str = "upper right") -> None:
    """Draw a small text box of summary stats in *loc* corner of *ax*.

    Used by the volcano + scatter plots to show n_up / n_down / total
    or Pearson r without crowding the data area.
    """
    text = "\n".join(lines)
    ha = "right" if "right" in loc else "left"
    va = "top" if "upper" in loc else "bottom"
    x = 0.97 if ha == "right" else 0.03
    y = 0.97 if va == "top" else 0.03
    ax.text(
        x, y, text,
        transform=ax.transAxes,
        ha=ha, va=va,
        fontsize=9,
        bbox=dict(
            boxstyle="round,pad=0.35",
            facecolor="white",
            edgecolor="#cccccc",
            alpha=0.92,
        ),
        zorder=10,
    )


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

    Every gene is plotted (rows missing either log2FC are dropped).
    Axes are square and symmetric about the origin so the four
    quadrants read at the same scale; faint quadrant labels in the
    corners spell out the biological interpretation. The identity
    line ``y = x`` (translation tracks transcription) is drawn dashed.
    """
    rows = [r for r in rows
            if r.mrna_log2fc is not None and r.rpf_log2fc is not None]

    with plt.rc_context(_PUB_RC):
        fig, ax = plt.subplots(figsize=(5.4, 5.4))

        if not rows:
            ax.text(0.5, 0.5, "no genes with both mRNA and RPF log2FC",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            return _save_figure(fig, output_path)

        xs = [r.mrna_log2fc for r in rows]
        ys = [r.rpf_log2fc for r in rows]
        labels = [r.gene for r in rows]

        # Symmetric, square axes — the four quadrants must read at the
        # same scale so a "buffered" gene is visually distinguishable
        # from a "co-regulated" one without measuring tick distances.
        bound = max(1.0, max(abs(v) for v in xs + ys) * 1.15)
        ax.set_xlim(-bound, bound)
        ax.set_ylim(-bound, bound)
        ax.set_aspect("equal", adjustable="box")

        # Quadrant guides + identity line.
        ax.axhline(0, color=_C_GUIDE, linewidth=0.6, zorder=1)
        ax.axvline(0, color=_C_GUIDE, linewidth=0.6, zorder=1)
        ax.plot([-bound, bound], [-bound, bound],
                color=_C_GUIDE, linewidth=0.6, linestyle="--",
                alpha=0.55, label="y = x (TE unchanged)", zorder=1)

        # Faint quadrant captions (fixed-position so they never collide
        # with data even when a quadrant is dense).
        cap = bound * 0.93
        for x, y, ha, va, txt in (
            (cap,  cap,  "right", "top",
                "Q1: co-regulated up\nmRNA + RPF up"),
            (-cap, cap,  "left",  "top",
                "Q2: buffered up\nmRNA down, RPF up"),
            (-cap, -cap, "left",  "bottom",
                "Q3: co-regulated down\nmRNA + RPF down"),
            (cap,  -cap, "right", "bottom",
                "Q4: buffered down\nmRNA up, RPF down"),
        ):
            ax.text(x, y, txt, ha=ha, va=va,
                    fontsize=7.5, color=_C_NEUTRAL, alpha=0.85,
                    style="italic")

        # Points, coloured by quadrant for at-a-glance reading.
        colors = [
            _C_UP if (rx > 0 and ry > 0)
            else _C_DOWN if (rx < 0 and ry < 0)
            else _C_NEUTRAL
            for rx, ry in zip(xs, ys)
        ]
        ax.scatter(xs, ys, s=58, alpha=0.92,
                   c=colors, edgecolor="black", linewidth=0.45, zorder=3)
        _annotate_with_bbox(
            ax, zip(xs, ys, labels), fontsize=8, smart=True,
        )

        ax.set_xlabel("log2FC mRNA  (RNA-seq)")
        ax.set_ylabel("log2FC RPF  (Ribo-seq)")
        ax.set_title(title)
        ax.legend(
            loc="upper left", bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0.0, fontsize=8,
        )

        return _save_figure(fig, output_path)


def plot_delta_te_volcano(
    rows: Iterable[DTeRow],
    output_path: Path,
    *,
    padj_threshold: float = 0.05,
    log2fc_threshold: float = 1.0,
    title: str = "Delta-TE (log2) volcano",
) -> Path:
    """Volcano-style scatter of delta_te_log2 (x) vs -log10(padj) (y).

    Every gene with a finite delta-TE is plotted. Points are coloured
    by combined significance (padj + |delta-TE| thresholds): red for
    TE-up genes, blue for TE-down, grey for n.s. Threshold guides at
    ``±log2fc_threshold`` (vertical) and ``-log10(padj_threshold)``
    (horizontal) make the cut-offs explicit. ``padj == 0`` is rendered
    at the largest finite ``-log10(padj)`` on the figure (a 10%
    headroom above the highest real point) so reviewers see the
    extreme dots without an infinity excursion.

    When the DE table has no padj column the y-axis collapses to 0
    for every point; the file still exists and the title carries the
    note so reviewers see the limitation.
    """
    rows = [r for r in rows if r.delta_te_log2 is not None]

    with plt.rc_context(_PUB_RC):
        fig, ax = plt.subplots(figsize=(5.6, 5.0))

        if not rows:
            ax.text(0.5, 0.5, "no delta-TE values to plot",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            return _save_figure(fig, output_path)

        xs_up: list[float] = []
        ys_up: list[float] = []
        labels_up: list[str] = []
        xs_dn: list[float] = []
        ys_dn: list[float] = []
        labels_dn: list[str] = []
        xs_ns: list[float] = []
        ys_ns: list[float] = []
        labels_ns: list[str] = []
        n_no_padj = 0

        for r in rows:
            x = r.delta_te_log2
            padj = r.padj
            if padj is None:
                n_no_padj += 1
                y = 0.0
                is_sig = False
            elif padj <= 0:
                y = math.inf
                is_sig = True
            else:
                y = -math.log10(padj)
                is_sig = padj < padj_threshold
            is_up = is_sig and x >= log2fc_threshold
            is_dn = is_sig and x <= -log2fc_threshold

            if is_up:
                xs_up.append(x)
                ys_up.append(y)
                labels_up.append(r.gene)
            elif is_dn:
                xs_dn.append(x)
                ys_dn.append(y)
                labels_dn.append(r.gene)
            else:
                xs_ns.append(x)
                ys_ns.append(y)
                labels_ns.append(r.gene)

        finite_ys = [y for y in ys_up + ys_dn + ys_ns if math.isfinite(y)]
        y_cap = (max(finite_ys) if finite_ys else 1.0) * 1.10

        def _cap(ys: list[float]) -> list[float]:
            return [y if math.isfinite(y) else y_cap for y in ys]

        ys_up = _cap(ys_up)
        ys_dn = _cap(ys_dn)
        ys_ns = _cap(ys_ns)

        # Threshold guides — drawn first so points sit on top.
        if any(math.isfinite(y) and y > 0 for y in finite_ys):
            ax.axhline(-math.log10(padj_threshold),
                       color=_C_GUIDE, linewidth=0.7, linestyle="--",
                       alpha=0.7, zorder=1)
        ax.axvline(log2fc_threshold,  color=_C_GUIDE, linewidth=0.7,
                   linestyle="--", alpha=0.7, zorder=1)
        ax.axvline(-log2fc_threshold, color=_C_GUIDE, linewidth=0.7,
                   linestyle="--", alpha=0.7, zorder=1)
        ax.axvline(0, color=_C_GUIDE, linewidth=0.6, zorder=1)

        if xs_ns:
            ax.scatter(xs_ns, ys_ns, s=44, alpha=0.78, c=_C_NS,
                       edgecolor="white", linewidth=0.4,
                       label="n.s.", zorder=2)
        if xs_up:
            ax.scatter(xs_up, ys_up, s=58, alpha=0.95, c=_C_UP,
                       edgecolor="black", linewidth=0.45,
                       label=f"TE up (padj<{padj_threshold}, ΔTE≥{log2fc_threshold})",
                       zorder=3)
        if xs_dn:
            ax.scatter(xs_dn, ys_dn, s=58, alpha=0.95, c=_C_DOWN,
                       edgecolor="black", linewidth=0.45,
                       label=f"TE down (padj<{padj_threshold}, ΔTE≤-{log2fc_threshold})",
                       zorder=3)

        # Symmetric x-axis + generous y headroom so labels at y_cap
        # render fully (no top-edge clipping by savefig.bbox=tight).
        x_bound = max(
            1.2,
            max(abs(r.delta_te_log2) for r in rows) * 1.18,
        )
        ax.set_xlim(-x_bound, x_bound)
        ax.set_ylim(-0.08 * y_cap, y_cap * 1.22)

        # Label every gene; bbox-bg + leader keeps it legible.
        all_labels = (
            list(zip(xs_up, ys_up, labels_up))
            + list(zip(xs_dn, ys_dn, labels_dn))
            + list(zip(xs_ns, ys_ns, labels_ns))
        )
        _annotate_with_bbox(ax, all_labels, fontsize=8, smart=True)

        _stat_box(ax, [
            f"genes plotted: {len(rows)}",
            f"TE up:   {len(xs_up)}",
            f"TE down: {len(xs_dn)}",
            f"n.s.:    {len(xs_ns)}",
        ], loc="lower right")

        ax.set_xlabel("Delta-TE  (log2)")
        ax.set_ylabel("-log10(padj)")
        ax.set_title(
            title + ("  (no padj column — y axis is zero)"
                     if n_no_padj == len(rows) else "")
        )
        legend_handles, legend_labels = ax.get_legend_handles_labels()
        if legend_handles:
            ax.legend(
                legend_handles, legend_labels,
                loc="upper left", bbox_to_anchor=(1.02, 1.0),
                borderaxespad=0.0, fontsize=8,
            )
        return _save_figure(
            fig,
            output_path,
            metadata={
                "stage": "rnaseq",
                "plot_type": "delta_te_volcano",
                "source_data": "rnaseq/delta_te.tsv",
                "n_points_expected": len(rows),
                "n_points_drawn": len(xs_up) + len(xs_dn) + len(xs_ns),
                "n_labels": len(all_labels),
                "labels_drawn": len(all_labels),
                "label_policy": "all_genes" if len(rows) <= 20 else "top_by_significance",
                "min_font_size": 8,
            },
        )


def plot_de_volcano(
    de_table: DeTable,
    output_path: Path,
    *,
    padj_threshold: float = 0.05,
    log2fc_threshold: float = 1.0,
    contrast_label: str | None = None,
    title: str | None = None,
    annotate: bool = True,
) -> Path:
    """DE volcano: log2FoldChange (x) vs -log10(padj) (y).

    Every gene with a finite log2FC is plotted. Points are coloured by
    significance:

    * **vermillion** — padj < ``padj_threshold`` AND log2FC ≥  ``log2fc_threshold`` (up)
    * **blue**       — padj < ``padj_threshold`` AND log2FC ≤ -``log2fc_threshold`` (down)
    * **grey**       — everything else

    Threshold guides are dashed at ±``log2fc_threshold`` (vertical) and
    ``-log10(padj_threshold)`` (horizontal). ``padj == 0`` is rendered
    at 1.10 × the largest finite ``-log10(padj)`` on the figure so
    reviewers see the extreme dots without an infinity excursion. A
    summary stat box (top-right) reports the counts. When ``annotate``
    is true, every gene is labelled with a white-bbox tag and a thin
    leader line — fine for the ~13 mt-mRNA universe; pass
    ``annotate=False`` for whole-transcriptome inputs.
    """
    rows = [
        r for r in de_table.rows
        if r.get("log2fc") is not None
    ]

    with plt.rc_context(_PUB_RC):
        fig, ax = plt.subplots(figsize=(6.0, 5.2))

        if not rows:
            ax.text(0.5, 0.5, "no DE rows with a log2FoldChange",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            return _save_figure(fig, output_path)

        xs_up: list[float] = []
        ys_up: list[float] = []
        xs_dn: list[float] = []
        ys_dn: list[float] = []
        xs_ns: list[float] = []
        ys_ns: list[float] = []
        labels: list[tuple[float, float, str]] = []
        n_no_padj = 0
        for row in rows:
            l2 = row["log2fc"]
            padj = row.get("padj")
            gene = str(row.get("gene_id", ""))
            if padj is None:
                n_no_padj += 1
                y = 0.0
                is_sig = False
            elif padj <= 0:
                y = math.inf
                is_sig = True
            else:
                y = -math.log10(padj)
                is_sig = padj < padj_threshold
            is_up = is_sig and l2 >= log2fc_threshold
            is_dn = is_sig and l2 <= -log2fc_threshold
            if is_up:
                xs_up.append(l2)
                ys_up.append(y)
            elif is_dn:
                xs_dn.append(l2)
                ys_dn.append(y)
            else:
                xs_ns.append(l2)
                ys_ns.append(y)
            labels.append((l2, y, gene))

        finite_ys = [y for y in ys_up + ys_dn + ys_ns if math.isfinite(y)]
        y_cap = (max(finite_ys) if finite_ys else 1.0) * 1.10

        def _cap(ys: list[float]) -> list[float]:
            return [y if math.isfinite(y) else y_cap for y in ys]

        ys_up = _cap(ys_up)
        ys_dn = _cap(ys_dn)
        ys_ns = _cap(ys_ns)
        labels = [(x, y_cap if not math.isfinite(y) else y, g)
                  for x, y, g in labels]

        # Symmetric x-axis for visual balance — the volcano reads at
        # the same scale on both sides of zero.
        x_finite = [r["log2fc"] for r in rows]
        x_bound = max(1.0, max(abs(v) for v in x_finite) * 1.18)
        ax.set_xlim(-x_bound, x_bound)
        # Generous y headroom: top labels on points at y_cap need
        # room to render above their dot without spilling out of the
        # axes (which savefig.bbox=tight then crops). Slight bottom
        # headroom for the same reason on the smart-placer's
        # below-point fallback positions.
        ax.set_ylim(-0.08 * y_cap, y_cap * 1.22)

        # Threshold guides — under data so points / labels read on top.
        if finite_ys and any(y > 0 for y in finite_ys):
            ax.axhline(-math.log10(padj_threshold),
                       color=_C_GUIDE, linewidth=0.7,
                       linestyle="--", alpha=0.7, zorder=1)
        ax.axvline(log2fc_threshold,  color=_C_GUIDE, linewidth=0.7,
                   linestyle="--", alpha=0.7, zorder=1)
        ax.axvline(-log2fc_threshold, color=_C_GUIDE, linewidth=0.7,
                   linestyle="--", alpha=0.7, zorder=1)
        ax.axvline(0, color=_C_GUIDE, linewidth=0.6, zorder=1)

        # Points. Sig dots are larger + black-edged so they pop.
        if xs_ns:
            ax.scatter(xs_ns, ys_ns, s=42, alpha=0.78, c=_C_NS,
                       edgecolor="white", linewidth=0.4,
                       label="n.s.", zorder=2)
        if xs_up:
            ax.scatter(xs_up, ys_up, s=58, alpha=0.95, c=_C_UP,
                       edgecolor="black", linewidth=0.45,
                       label=f"up (padj<{padj_threshold}, L2FC≥{log2fc_threshold})",
                       zorder=3)
        if xs_dn:
            ax.scatter(xs_dn, ys_dn, s=58, alpha=0.95, c=_C_DOWN,
                       edgecolor="black", linewidth=0.45,
                       label=f"down (padj<{padj_threshold}, L2FC≤-{log2fc_threshold})",
                       zorder=3)

        if annotate:
            _annotate_with_bbox(ax, labels, fontsize=8, smart=True)

        _stat_box(ax, [
            f"genes plotted: {len(rows)}",
            f"up:   {len(xs_up)}",
            f"down: {len(xs_dn)}",
            f"n.s.: {len(xs_ns)}",
        ], loc="lower right")

        ax.set_xlabel(
            "log2FoldChange"
            + (f"  ({contrast_label})" if contrast_label else "")
        )
        ax.set_ylabel("-log10(padj)")
        contrast_suffix = f" — {contrast_label}" if contrast_label else ""
        ax.set_title(
            (title or f"DE volcano{contrast_suffix}")
            + ("  (no padj column — y axis is zero)"
               if n_no_padj == len(rows) else "")
        )
        # Legend is anchored OUTSIDE the data axes (right of the
        # volcano cloud) so it never collides with gene labels or the
        # significant-hit cluster. Skip empty bins so the legend
        # doesn't promise dot colours that aren't on the plot.
        legend_handles, legend_labels = ax.get_legend_handles_labels()
        if legend_handles:
            ax.legend(
                legend_handles, legend_labels,
                loc="upper left", bbox_to_anchor=(1.02, 1.0),
                borderaxespad=0.0, fontsize=8,
            )
        n_total = len(xs_up) + len(xs_dn) + len(xs_ns)
        return _save_figure(
            fig,
            output_path,
            metadata={
                "stage": "rnaseq",
                "plot_type": "de_volcano",
                "source_data": "rnaseq/de_table.tsv",
                "n_points_expected": len(rows),
                "n_points_drawn": n_total,
                "n_labels": len(labels) if annotate else 0,
                "labels_drawn": len(labels) if annotate else 0,
                "label_policy": (
                    "all_genes" if (annotate and len(rows) <= 20)
                    else "top_by_significance" if annotate
                    else "none"
                ),
                "min_font_size": 8,
            },
        )


def plot_te_compare_scatter(
    te_rows: Iterable[TeRow],
    condition_map: dict[str, str],
    base: str,
    compare: str,
    output_path: Path,
    *,
    title: str | None = None,
) -> Path:
    """Per-gene mean ``log2(TE)`` scatter: ``base`` (x) vs ``compare`` (y).

    Each gene contributes one point — the **mean** of its replicate
    ``log2(TE)`` values within each condition. Individual replicates
    are also drawn as small dots (jittered by ±0.05 along the gene's
    coordinate to make overlapping replicates visible) so reviewers
    see the variance behind every mean. The identity ``y = x`` line
    is drawn dashed for reference. Pearson r over gene-mean points
    appears in the stat box.

    Points are coloured by direction:

    * **vermillion** — ``compare`` mean above identity (TE up in compare)
    * **blue**       — ``compare`` mean below identity (TE down in compare)
    * **grey**       — within ±0.05 log2 of identity (≈10% TE change)
    """
    grouped: dict[str, dict[str, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for row in te_rows:
        cond = condition_map.get(row.sample)
        if cond not in (base, compare) or row.te <= 0:
            continue
        grouped[row.gene][cond].append(math.log2(row.te))

    genes_xy: list[tuple[str, float, float, list[float], list[float]]] = []
    for gene in sorted(grouped):
        cond_vals = grouped[gene]
        if base not in cond_vals or compare not in cond_vals:
            continue
        b_vals = cond_vals[base]
        c_vals = cond_vals[compare]
        b_mean = sum(b_vals) / len(b_vals)
        c_mean = sum(c_vals) / len(c_vals)
        genes_xy.append((gene, b_mean, c_mean, b_vals, c_vals))

    with plt.rc_context(_PUB_RC):
        fig, ax = plt.subplots(figsize=(5.6, 5.6))

        if not genes_xy:
            ax.text(0.5, 0.5,
                    f"no genes shared between {base!r} and {compare!r}",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            return _save_figure(fig, output_path)

        xs = [t[1] for t in genes_xy]
        ys = [t[2] for t in genes_xy]
        labels = [t[0] for t in genes_xy]

        # Square axes covering both samples + replicates.
        all_pts = list(xs) + list(ys)
        for _, _, _, b_vals, c_vals in genes_xy:
            all_pts.extend(b_vals)
            all_pts.extend(c_vals)
        lo = min(all_pts)
        hi = max(all_pts)
        span = hi - lo
        margin = max(0.15, 0.06 * span)
        bound_lo, bound_hi = lo - margin, hi + margin
        ax.set_xlim(bound_lo, bound_hi)
        ax.set_ylim(bound_lo, bound_hi)
        ax.set_aspect("equal", adjustable="box")

        # Identity line + axes guides.
        ax.plot([bound_lo, bound_hi], [bound_lo, bound_hi],
                color=_C_GUIDE, linewidth=0.7, linestyle="--",
                alpha=0.7, label="y = x (identity)", zorder=1)
        ax.axhline(0, color=_C_GUIDE, linewidth=0.4, alpha=0.5, zorder=1)
        ax.axvline(0, color=_C_GUIDE, linewidth=0.4, alpha=0.5, zorder=1)

        # Individual replicate cloud (light, behind the means). Jitter
        # small so the variance is visible without distorting position.
        rng = np.random.default_rng(0)
        for _, b_mean, c_mean, b_vals, c_vals in genes_xy:
            n = max(len(b_vals), len(c_vals))
            for i in range(n):
                bx = b_vals[i % len(b_vals)] + rng.uniform(-0.04, 0.04)
                cy = c_vals[i % len(c_vals)] + rng.uniform(-0.04, 0.04)
                ax.scatter(bx, cy, s=14, alpha=0.42, c=_C_NS,
                           edgecolor="none", zorder=2)

        # Per-gene means coloured by direction.
        diag_eps = 0.05
        colors = [
            _C_UP if (cy - bx) > diag_eps
            else _C_DOWN if (cy - bx) < -diag_eps
            else _C_NEUTRAL
            for bx, cy in zip(xs, ys)
        ]
        ax.scatter(xs, ys, s=64, alpha=0.95, c=colors,
                   edgecolor="black", linewidth=0.5, zorder=3)
        _annotate_with_bbox(ax, zip(xs, ys, labels),
                            fontsize=8, smart=True)

        # Pearson r over gene-mean points (n is small but still useful).
        if len(xs) >= 3:
            r = float(np.corrcoef(xs, ys)[0, 1])
            stats = [f"genes: {len(xs)}", f"Pearson r = {r:.3f}"]
        else:
            stats = [f"genes: {len(xs)}"]
        # Counts above / below identity.
        n_up = sum(1 for c in colors if c == _C_UP)
        n_dn = sum(1 for c in colors if c == _C_DOWN)
        stats += [f"above y=x: {n_up}", f"below y=x: {n_dn}"]
        _stat_box(ax, stats, loc="upper left")

        ax.set_xlabel(f"log2(TE) — {base}")
        ax.set_ylabel(f"log2(TE) — {compare}")
        ax.set_title(title or f"Per-gene log2(TE): {base} vs {compare}")
        ax.legend(
            loc="upper left", bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0.0, fontsize=8,
        )

        return _save_figure(fig, output_path)


def plot_te_log2fc_bar(
    te_rows: Iterable[TeRow],
    condition_map: dict[str, str],
    base: str,
    compare: str,
    output_path: Path,
    *,
    title: str | None = None,
) -> Path:
    """Sorted bar plot of ``log2(mean_TE_compare / mean_TE_base)`` per gene.

    Bars above zero: TE up in ``compare`` relative to ``base``; below
    zero: TE down. Sorted by value so the strongest condition-specific
    TE movers appear at the extremes. Each bar carries its numeric
    log2FC at the tip and its replicate counts in the legend so the
    figure stays self-contained without a companion table.
    """
    grouped: dict[str, dict[str, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for row in te_rows:
        cond = condition_map.get(row.sample)
        if cond not in (base, compare) or row.te <= 0:
            continue
        grouped[row.gene][cond].append(row.te)

    pairs: list[tuple[str, float, int, int]] = []
    for gene in grouped:
        cond_vals = grouped[gene]
        if base not in cond_vals or compare not in cond_vals:
            continue
        b_vals = cond_vals[base]
        c_vals = cond_vals[compare]
        m_base = sum(b_vals) / len(b_vals)
        m_comp = sum(c_vals) / len(c_vals)
        if m_base <= 0 or m_comp <= 0:
            continue
        pairs.append((gene, math.log2(m_comp / m_base),
                      len(b_vals), len(c_vals)))
    pairs.sort(key=lambda gp: gp[1])

    with plt.rc_context(_PUB_RC):
        fig, ax = plt.subplots(
            figsize=(max(5.4, 0.55 * len(pairs) + 2.2), 4.4)
        )

        if not pairs:
            ax.text(0.5, 0.5,
                    f"no genes shared between {base!r} and {compare!r}",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            return _save_figure(fig, output_path)

        genes = [p[0] for p in pairs]
        values = [p[1] for p in pairs]
        rep_b = pairs[0][2] if pairs else 0
        rep_c = pairs[0][3] if pairs else 0
        colors = [_C_UP if v > 0 else _C_DOWN for v in values]
        bars = ax.bar(
            range(len(genes)), values,
            color=colors, alpha=0.92,
            edgecolor="black", linewidth=0.45,
            label=(f"{compare} (n={rep_c})  vs  {base} (n={rep_b})"
                   if rep_b and rep_c else None),
        )
        ax.axhline(0, color=_C_GUIDE, linewidth=0.7, zorder=1)

        # Numeric log2FC tip labels.
        v_max = max(abs(v) for v in values)
        pad = max(0.04, 0.04 * v_max)
        for bar, v in zip(bars, values):
            x = bar.get_x() + bar.get_width() / 2
            if v >= 0:
                ax.text(x, v + pad, f"{v:+.2f}",
                        ha="center", va="bottom", fontsize=8)
            else:
                ax.text(x, v - pad, f"{v:+.2f}",
                        ha="center", va="top", fontsize=8)

        # Headroom so the tip labels don't get clipped.
        y_low = min(min(values), 0) - 1.6 * pad
        y_high = max(max(values), 0) + 1.6 * pad
        ax.set_ylim(y_low, y_high)

        ax.set_xticks(range(len(genes)))
        ax.set_xticklabels(genes, rotation=40, ha="right")
        ax.set_ylabel(f"log2(TE_{compare} / TE_{base})")
        ax.set_title(
            title or f"Per-gene log2(TE) fold change: {compare} vs {base}"
        )
        # Tiny direction legend so a reviewer who lands on this plot
        # alone reads colour semantics from the figure itself.
        legend_handles = [
            plt.Rectangle((0, 0), 1, 1, color=_C_UP,
                          label=f"TE up in {compare}"),
            plt.Rectangle((0, 0), 1, 1, color=_C_DOWN,
                          label=f"TE down in {compare}"),
        ]
        ax.legend(
            handles=legend_handles,
            loc="upper left", bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0.0, fontsize=8,
        )

        return _save_figure(fig, output_path)


def plot_ma(
    de_table: DeTable,
    output_path: Path,
    *,
    padj_threshold: float = 0.05,
    title: str = "MA plot (RNA-seq DE)",
) -> Path:
    """Classic DESeq2 MA plot: ``log10(baseMean)`` (x) vs log2FoldChange (y).

    Every gene with a positive ``baseMean`` and a finite log2FC is
    plotted. Points are coloured by direction-of-significance:
    vermillion for sig up, blue for sig down, grey for n.s. — the
    same semantics used elsewhere in the rnaseq plot set so the
    palette reads consistently across figures. All sig genes are
    labelled with white-bbox tags + leader lines.
    """
    rows = [
        r for r in de_table.rows
        if r.get("basemean") is not None
        and r.get("log2fc") is not None
        and r["basemean"] > 0
    ]

    with plt.rc_context(_PUB_RC):
        fig, ax = plt.subplots(figsize=(5.6, 5.0))
        if not rows:
            ax.text(0.5, 0.5, "no rows with baseMean > 0",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            return _save_figure(fig, output_path)

        xs_up: list[float] = []
        ys_up: list[float] = []
        xs_dn: list[float] = []
        ys_dn: list[float] = []
        xs_ns: list[float] = []
        ys_ns: list[float] = []
        sig_labels: list[tuple[float, float, str]] = []
        for r in rows:
            x = math.log10(r["basemean"])
            l2 = r["log2fc"]
            padj = r.get("padj")
            is_sig = padj is not None and padj < padj_threshold
            if is_sig and l2 > 0:
                xs_up.append(x)
                ys_up.append(l2)
                sig_labels.append((x, l2, str(r.get("gene_id", ""))))
            elif is_sig and l2 < 0:
                xs_dn.append(x)
                ys_dn.append(l2)
                sig_labels.append((x, l2, str(r.get("gene_id", ""))))
            else:
                xs_ns.append(x)
                ys_ns.append(l2)

        ax.axhline(0, color=_C_GUIDE, linewidth=0.7, zorder=1)
        if xs_ns:
            ax.scatter(xs_ns, ys_ns, s=42, alpha=0.78, c=_C_NS,
                       edgecolor="white", linewidth=0.4,
                       label=f"padj ≥ {padj_threshold}", zorder=2)
        if xs_up:
            ax.scatter(xs_up, ys_up, s=58, alpha=0.95, c=_C_UP,
                       edgecolor="black", linewidth=0.45,
                       label=f"sig up (padj<{padj_threshold})",
                       zorder=3)
        if xs_dn:
            ax.scatter(xs_dn, ys_dn, s=58, alpha=0.95, c=_C_DOWN,
                       edgecolor="black", linewidth=0.45,
                       label=f"sig down (padj<{padj_threshold})",
                       zorder=3)
        _annotate_with_bbox(ax, sig_labels, fontsize=8, smart=True)

        ax.set_xlabel("log10(baseMean)")
        ax.set_ylabel("log2FoldChange")
        ax.set_title(title)
        if xs_ns or xs_up or xs_dn:
            ax.legend(
                loc="upper left", bbox_to_anchor=(1.02, 1.0),
                borderaxespad=0.0, fontsize=8,
            )
        return _save_figure(fig, output_path)


def plot_te_bar_grouped(
    te_rows: Iterable[TeRow],
    condition_map: dict[str, str],
    output_path: Path,
    *,
    title: str = "Translation efficiency per gene by condition",
    log2: bool = True,
) -> Path:
    """Per-gene TE bar plot, bars grouped by condition with SE error bars
    AND every replicate's value overlaid as a small dot.

    Replicates within each (gene, condition) cell are averaged for the
    bar height; error bars show the standard error of the mean. Each
    replicate's individual TE value is then drawn as a black dot
    jittered along x within the bar's footprint so the within-group
    spread is visible by eye — no implied SE from a 2-rep N. The
    publication-grade rule of thumb: bars without dots hide the
    sample size.

    ``log2=True`` plots ``log2(TE)`` (y centred on 0 = no change);
    ``log2=False`` plots raw TE.
    """
    grouped: dict[str, dict[str, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for row in te_rows:
        cond = condition_map.get(row.sample)
        if cond is None or row.te <= 0:
            continue
        value = math.log2(row.te) if log2 else row.te
        grouped[row.gene][cond].append(value)

    with plt.rc_context(_PUB_RC):
        if not grouped:
            fig, ax = plt.subplots(figsize=(5, 3))
            ax.text(0.5, 0.5, "no TE rows with a condition assignment",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            return _save_figure(fig, output_path)

        genes = sorted(grouped.keys())
        conditions = sorted({c for d in grouped.values() for c in d.keys()})
        n_genes = len(genes)
        n_conds = len(conditions)
        width = 0.84 / max(n_conds, 1)
        x_indices = np.arange(n_genes)

        # Use the Okabe-Ito palette in a fixed order so condition
        # colours stay consistent across figures.
        palette = [
            _OKABE_ITO["blue"], _OKABE_ITO["vermillion"],
            _OKABE_ITO["green"], _OKABE_ITO["orange"],
            _OKABE_ITO["purple"], _OKABE_ITO["skyblue"],
            _OKABE_ITO["yellow"], _OKABE_ITO["black"],
        ]
        cond_color = {c: palette[i % len(palette)]
                      for i, c in enumerate(conditions)}

        fig, ax = plt.subplots(
            figsize=(max(6.4, 0.78 * n_genes + 2.2), 4.6)
        )

        rng = np.random.default_rng(0)
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
                color=cond_color[cond], alpha=0.85,
                edgecolor="black", linewidth=0.45, zorder=2,
                error_kw=dict(elinewidth=0.9, ecolor="black"),
            )
            # Overlay every replicate as a dot — the publication-grade
            # signal that the bar is a real summary of >=2 points.
            for gene_idx, gene in enumerate(genes):
                values = grouped[gene].get(cond, [])
                if not values:
                    continue
                center = x_indices[gene_idx] + offset
                jitter_w = 0.55 * width
                for v in values:
                    jx = center + rng.uniform(-jitter_w / 2, jitter_w / 2)
                    ax.scatter(jx, v, s=18, c="black", alpha=0.85,
                               edgecolor="white", linewidth=0.5, zorder=3)

        ax.axhline(0, color=_C_GUIDE, linewidth=0.7, zorder=1)
        ax.set_xticks(x_indices)
        ax.set_xticklabels(genes, rotation=40, ha="right")
        ax.set_ylabel("log2(TE)" if log2 else "TE")
        ax.set_title(title)
        ax.legend(
            loc="upper left", bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0.0, fontsize=9,
            title="condition", title_fontsize=9,
        )
        return _save_figure(fig, output_path)


def plot_te_heatmap(
    te_rows: Iterable[TeRow],
    condition_map: dict[str, str] | None,
    output_path: Path,
    *,
    title: str = "log2(TE) heatmap (gene × sample)",
) -> Path:
    """Gene × sample heatmap of log2(TE).

    Sample columns are sorted by condition (from ``condition_map``) then
    by name so replicates within a condition cluster together. A
    coloured strip above the columns spells the condition assignment
    so reviewers do not have to read sample names. Cells are
    annotated with their numeric log2(TE) value (auto-contrast: white
    on saturated, black on faint) when the matrix is small enough
    that annotations stay legible (≤ 200 cells).
    """
    rows = [r for r in te_rows if r.te > 0]

    with plt.rc_context(_PUB_RC):
        if not rows:
            fig, ax = plt.subplots(figsize=(5, 3))
            ax.text(0.5, 0.5, "no TE rows to plot",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            return _save_figure(fig, output_path)

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

        abs_max = (
            float(np.nanmax(np.abs(matrix)))
            if np.any(~np.isnan(matrix)) else 1.0
        )
        abs_max = max(abs_max, 0.5)

        # Layout: a thin condition strip above the heatmap, then the
        # main heatmap. height_ratios keeps the strip narrow.
        has_strip = bool(condition_map)
        if has_strip:
            fig, (ax_strip, ax) = plt.subplots(
                2, 1,
                figsize=(max(5.6, 0.6 * len(samples) + 2.4),
                         0.46 * len(genes) + 2.4),
                gridspec_kw={"height_ratios": [0.6, max(6, len(genes))],
                             "hspace": 0.04},
            )
        else:
            fig, ax = plt.subplots(
                figsize=(max(5.6, 0.6 * len(samples) + 2.4),
                         0.46 * len(genes) + 2.0)
            )
            ax_strip = None

        im = ax.imshow(
            matrix, aspect="auto", cmap="RdBu_r",
            vmin=-abs_max, vmax=abs_max,
        )

        ax.set_xticks(np.arange(len(samples)))
        ax.set_yticks(np.arange(len(genes)))
        ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(genes, fontsize=9)
        ax.tick_params(top=False, bottom=True, left=True, right=False)

        # Condition strip with one tick per condition group.
        if ax_strip is not None:
            cond_labels = [(condition_map or {}).get(s, "") for s in samples]
            unique_conds = sorted(set(cond_labels))
            palette = [
                _OKABE_ITO["blue"], _OKABE_ITO["vermillion"],
                _OKABE_ITO["green"], _OKABE_ITO["orange"],
                _OKABE_ITO["purple"], _OKABE_ITO["skyblue"],
                _OKABE_ITO["yellow"], _OKABE_ITO["black"],
            ]
            cond_color = {
                c: palette[i % len(palette)]
                for i, c in enumerate(unique_conds)
            }
            colors = np.array([
                tuple(int(cond_color[c].lstrip("#")[i:i+2], 16) / 255
                      for i in (0, 2, 4)) + (1.0,)
                for c in cond_labels
            ]).reshape(1, len(samples), 4)
            ax_strip.imshow(colors, aspect="auto",
                            extent=(-0.5, len(samples) - 0.5, 0, 1))
            ax_strip.set_yticks([])
            ax_strip.set_xticks([])
            for spine in ax_strip.spines.values():
                spine.set_visible(False)
            # Condition labels above each block.
            run_start = 0
            for j in range(1, len(samples) + 1):
                if j == len(samples) or cond_labels[j] != cond_labels[run_start]:
                    mid = (run_start + j - 1) / 2
                    ax_strip.text(mid, 0.5, cond_labels[run_start] or "",
                                  ha="center", va="center",
                                  fontsize=9, color="white",
                                  fontweight="bold")
                    run_start = j

        # Cell annotations (small enough to stay legible).
        if matrix.size <= 200:
            for i in range(matrix.shape[0]):
                for j in range(matrix.shape[1]):
                    v = matrix[i, j]
                    if math.isnan(v):
                        continue
                    color = "white" if abs(v) > 0.55 * abs_max else "black"
                    ax.text(j, i, f"{v:.2f}",
                            ha="center", va="center",
                            fontsize=7, color=color)

        cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
        cbar.set_label("log2(TE)", rotation=270, labelpad=14)
        cbar.outline.set_linewidth(0.6)
        # Use the strip-axes for the title so it sits above the strip
        # (when present); otherwise put it on the heatmap.
        (ax_strip or ax).set_title(title)
        return _save_figure(fig, output_path)


def plot_pca(
    counts_df: "pd.DataFrame",
    metadata_df: "pd.DataFrame",
    output_path: Path,
    *,
    title: str = "Sample PCA",
) -> Path:
    """PC1 vs PC2 of samples on log1p-transformed counts.

    Colour = ``condition``, marker shape = ``assay`` (when present).
    Per-sample dots get a white-bbox label + leader line; samples are
    drawn at a generous size (s=92) so the marker shape is unambiguous
    in printed figures. PCA is computed via numpy SVD so we do not
    pull in scikit-learn. Falls back to a placeholder figure when
    fewer than 2 samples are available (PCA is undefined).
    """
    samples = list(counts_df.index)
    n = len(samples)

    with plt.rc_context(_PUB_RC):
        if n < 2:
            fig, ax = plt.subplots(figsize=(5, 4))
            ax.text(0.5, 0.5, f"PCA needs ≥2 samples (got {n})",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            return _save_figure(fig, output_path)

        X = np.log1p(counts_df.to_numpy(dtype=float))
        keep = X.sum(axis=0) > 0
        X = X[:, keep]
        if X.shape[1] == 0:
            fig, ax = plt.subplots(figsize=(5, 4))
            ax.text(0.5, 0.5, "PCA: no genes with non-zero counts",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            return _save_figure(fig, output_path)

        X_centered = X - X.mean(axis=0, keepdims=True)
        U, S, _Vt = np.linalg.svd(X_centered, full_matrices=False)
        scores = U * S
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
        palette = [
            _OKABE_ITO["blue"], _OKABE_ITO["vermillion"],
            _OKABE_ITO["green"], _OKABE_ITO["orange"],
            _OKABE_ITO["purple"], _OKABE_ITO["skyblue"],
            _OKABE_ITO["yellow"], _OKABE_ITO["black"],
        ]
        cond_color = {c: palette[i % len(palette)]
                      for i, c in enumerate(unique_conds)}
        assay_marker = {
            a: m for a, m in zip(unique_assays,
                                 ["o", "s", "^", "D", "P", "X", "v", "<"])
        }

        fig, ax = plt.subplots(figsize=(6.4, 5.2))

        ax.axhline(0, color=_C_GUIDE, linewidth=0.5, alpha=0.5, zorder=1)
        ax.axvline(0, color=_C_GUIDE, linewidth=0.5, alpha=0.5, zorder=1)

        for s, x, y, c, a in zip(samples, pc1, pc2, conditions, assays):
            ax.scatter(
                x, y, s=92, alpha=0.95,
                color=cond_color[c],
                marker=assay_marker.get(a, "o"),
                edgecolor="black", linewidth=0.55, zorder=3,
            )
        # White-bbox labels with leader lines, smart-placed so dense
        # condition clusters do not produce overlapping sample IDs.
        _annotate_with_bbox(
            ax, zip(pc1, pc2, samples), fontsize=8, smart=True,
        )

        cond_handles = [
            plt.Line2D([0], [0], marker="o", color="w",
                       markerfacecolor=cond_color[c], markersize=9,
                       label=f"{c}", markeredgecolor="black",
                       markeredgewidth=0.55)
            for c in unique_conds
        ]
        assay_handles = [
            plt.Line2D([0], [0], marker=assay_marker.get(a, "o"),
                       color="w", markerfacecolor="#aaaaaa",
                       markersize=9,
                       label=f"{a}", markeredgecolor="black",
                       markeredgewidth=0.55)
            for a in unique_assays
        ]
        # Anchor both legends OUTSIDE the data axes (right of the
        # scatter) so they cannot overlap PCA points or sample labels.
        # Two legends stacked vertically: condition on top, assay below.
        cond_legend = ax.legend(
            handles=cond_handles,
            loc="upper left", bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0.0,
            title="condition", fontsize=8, title_fontsize=9,
        )
        ax.add_artist(cond_legend)
        if len(unique_assays) > 1:
            ax.legend(
                handles=assay_handles,
                loc="upper left", bbox_to_anchor=(1.02, 0.55),
                borderaxespad=0.0,
                title="assay", fontsize=8, title_fontsize=9,
            )

        ax.set_xlabel(
            f"PC1 ({var_ratio[0] * 100:.1f}%)"
            if len(var_ratio) >= 1 else "PC1"
        )
        ax.set_ylabel(
            f"PC2 ({var_ratio[1] * 100:.1f}%)"
            if len(var_ratio) >= 2 else "PC2"
        )
        ax.set_title(title)
        return _save_figure(fig, output_path)
