"""Render the four MitoRiboPy pipeline diagrams.

Produces four publication-quality PNGs under ``docs/diagrams/``:

    01_pipeline_overview.png     -- horizontal LR view of the three
                                    stages (align / rpf / rnaseq)
    02_align_stage.png           -- per-step internals of the align
                                    stage (cutadapt + bowtie2 + dedup)
    03_rpf_stage.png             -- per-step internals of the rpf
                                    stage (offsets + outputs)
    04_rnaseq_stage.png          -- per-step internals of the optional
                                    rnaseq stage (DE + TE + plots)

Run:

    python docs/diagrams/render_diagrams.py

Dependencies: matplotlib (already a project dependency). No Node, no
mermaid-cli required. Each PNG is rendered at 300 dpi and >= 1920 px
wide so it stays readable when embedded in the README.

Design notes
------------
Arrows route from one box's *edge* to the next box's *edge*, never
through a box centre, so a downstream arrow never visually crosses or
covers its source box. ``Box.right_at()`` / ``Box.left_at()`` /
``Box.top_at()`` / ``Box.bottom_at()`` accept an optional offset so
several arrows can share an edge without stacking.

Content tracks the v0.4.5 layout:

* the ``align`` stage drops mark-duplicates from its dedup options
* the ``rpf`` stage uses the flat ``translation_profile/<sample>/`` and
  ``coverage_profile_plots/{p_site,a_site}_density_*`` layout (no
  legacy ``p/`` / ``a/`` subfolders)
* offsets land under ``offset_diagnostics/{csv,plots}/`` (renamed from
  ``plots_and_csv/``)
* per-site ``codon_correlation/{p_site,a_site}/`` and ``igv_tracks/``
  appear as v0.4.4-and-later optional modules
* the ``--rpf_min_count_frac`` auto-prune is shown between filter-BED
  and offset enrichment
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


# Okabe-Ito colour-blind safe palette.
COL_INPUT = "#FFE5B4"     # peach -- raw input
COL_ALIGN = "#9BC6E4"     # blue  -- align stage
COL_RPF = "#A6D9A1"       # green -- rpf stage
COL_RNASEQ = "#F5A89B"    # pink  -- rnaseq stage
COL_OUTPUT = "#E0E0E0"    # grey  -- output artefact
COL_OPTIONAL = "#FFFFFF"  # white -- optional / dashed
EDGE_COL = "#222222"
TEXT_COL = "#111111"
SUBTITLE_COL = "#555555"


HERE = Path(__file__).resolve().parent

# Version label that appears in every figure's footer so a reader can
# tell at a glance which release the diagram describes.
VERSION_LABEL = "MitoRiboPy 0.4.5+"


@dataclass(frozen=True)
class Box:
    """Geometry of a drawn box. Provides edge-anchor helpers for arrows."""
    x: float
    y: float
    w: float
    h: float

    @property
    def cx(self) -> float:
        return self.x + self.w / 2

    @property
    def cy(self) -> float:
        return self.y + self.h / 2

    def left_at(self, frac: float = 0.5) -> tuple[float, float]:
        return self.x, self.y + self.h * frac

    def right_at(self, frac: float = 0.5) -> tuple[float, float]:
        return self.x + self.w, self.y + self.h * frac

    def top_at(self, frac: float = 0.5) -> tuple[float, float]:
        return self.x + self.w * frac, self.y + self.h

    def bottom_at(self, frac: float = 0.5) -> tuple[float, float]:
        return self.x + self.w * frac, self.y


def draw_box(
    ax,
    *,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    body: str = "",
    fill: str = COL_OUTPUT,
    title_size: int = 13,
    body_size: int = 10,
    dashed: bool = False,
    title_pad: float = 0.22,
    body_pad: float = 0.62,
) -> Box:
    """Draw a rounded rectangle with a bold title and an optional body.

    Returns a :class:`Box` capturing the geometry so callers can wire
    arrows into the right anchor edge instead of the centre.
    """
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.08",
        linewidth=1.5,
        linestyle="--" if dashed else "-",
        edgecolor=EDGE_COL,
        facecolor=fill,
    )
    ax.add_patch(box)
    cx, cy = x + w / 2, y + h / 2
    if body:
        ax.text(
            cx, y + h - title_pad,
            title,
            ha="center", va="top",
            fontsize=title_size, fontweight="bold", color=TEXT_COL,
        )
        ax.text(
            cx, y + h - body_pad,
            body,
            ha="center", va="top",
            fontsize=body_size, color=TEXT_COL,
        )
    else:
        ax.text(
            cx, cy,
            title,
            ha="center", va="center",
            fontsize=title_size, fontweight="bold", color=TEXT_COL,
        )
    return Box(x=x, y=y, w=w, h=h)


def draw_arrow(
    ax,
    src: tuple[float, float],
    dst: tuple[float, float],
    *,
    label: str = "",
    label_offset: tuple[float, float] = (0.0, 0.18),
    style: str = "-|>",
    dashed: bool = False,
    color: str = EDGE_COL,
    connection: str = "arc3,rad=0",
    linewidth: float = 1.5,
) -> None:
    arrow = FancyArrowPatch(
        src, dst,
        arrowstyle=style,
        mutation_scale=16,
        linewidth=linewidth,
        linestyle="--" if dashed else "-",
        color=color,
        connectionstyle=connection,
        shrinkA=0, shrinkB=0,
    )
    ax.add_patch(arrow)
    if label:
        mx = (src[0] + dst[0]) / 2 + label_offset[0]
        my = (src[1] + dst[1]) / 2 + label_offset[1]
        ax.text(
            mx, my, label,
            ha="center", va="center",
            fontsize=9, color=TEXT_COL,
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                      edgecolor="none", alpha=0.9),
        )


def draw_footer(ax, x: float, y: float, label: str = VERSION_LABEL) -> None:
    """Tiny version label so readers can match figure to release.

    Anchored to the LEFT of the given (x, y) and with a tiny top
    padding so ``bbox_inches="tight"`` does not clip it on the right
    margin. Place it at e.g. ``(0.1, 0.1)`` for a bottom-left anchor.
    """
    ax.text(
        x, y, label,
        ha="left", va="bottom",
        fontsize=9, style="italic", color=SUBTITLE_COL,
    )


def setup(figsize, xlim, ylim):
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    return fig, ax


def save(fig, name: str, dpi: int = 300) -> Path:
    out = HERE / name
    fig.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# 01 -- Pipeline overview (horizontal, 3 stages)
# ---------------------------------------------------------------------------


def render_pipeline_overview() -> Path:
    fig, ax = setup(figsize=(18, 7), xlim=(0, 18), ylim=(0, 7))

    ax.text(
        9, 6.55,
        "MitoRiboPy pipeline overview",
        ha="center", va="center",
        fontsize=20, fontweight="bold", color=TEXT_COL,
    )
    ax.text(
        9, 6.05,
        "FASTQ → BED + counts → offsets + translation/coverage profiles → (optional) TE / ΔTE",
        ha="center", va="center",
        fontsize=11, style="italic", color=SUBTITLE_COL,
    )

    fastq = draw_box(
        ax, x=0.2, y=2.7, w=2.0, h=1.6,
        title="FASTQ",
        body="(per sample;\nmixed kits OK)",
        fill=COL_INPUT, title_size=14, body_size=9,
    )

    align = draw_box(
        ax, x=2.7, y=1.7, w=4.4, h=3.6,
        title="align",
        body=(
            "cutadapt trim + UMI → QNAME\n"
            "bowtie2 contam subtract\n"
            "bowtie2 mt-transcriptome align\n"
            "MAPQ filter (NUMT suppression)\n"
            "umi-tools dedup or skip\n"
            "BAM → BED6"
        ),
        fill=COL_ALIGN, title_size=15, body_size=10,
        title_pad=0.32, body_pad=0.85,
    )

    align_out = draw_box(
        ax, x=7.55, y=2.5, w=2.05, h=2.0,
        title="align outputs",
        body=(
            "bed/, deduped/,\n"
            "kit_resolution.tsv,\n"
            "read_counts.tsv,\n"
            "run_settings.json"
        ),
        fill=COL_OUTPUT, title_size=11, body_size=8,
        title_pad=0.28, body_pad=0.65,
    )

    rpf = draw_box(
        ax, x=10.05, y=1.7, w=4.4, h=3.6,
        title="rpf",
        body=(
            "filter BED on RPF window\n"
            "auto-prune low-count lengths\n"
            "  (--rpf_min_count_frac)\n"
            "offset enrichment + selection\n"
            "  (combined + per-sample)\n"
            "translation_profile/<sample>/\n"
            "coverage_profile_plots/\n"
            "(opt) cor_plot, igv_export,\n"
            "  structure_density"
        ),
        fill=COL_RPF, title_size=15, body_size=9,
        title_pad=0.32, body_pad=0.85,
    )

    rpf_out = draw_box(
        ax, x=14.9, y=2.5, w=2.6, h=2.0,
        title="rpf outputs",
        body=(
            "translation_profile/, coverage_\n"
            "profile_plots/, offset_\n"
            "diagnostics/, codon_correlation/\n"
            "{p_site,a_site}/, igv_tracks/,\n"
            "rpf_counts.tsv, run_settings.json"
        ),
        fill=COL_OUTPUT, title_size=11, body_size=8,
        title_pad=0.28, body_pad=0.65,
    )

    rnaseq = draw_box(
        ax, x=14.9, y=0.15, w=2.6, h=1.95,
        title="rnaseq (optional)",
        body=(
            "DE table + rpf_counts →\n"
            "SHA256 ref-gate →\n"
            "TE / ΔTE + scatter, volcano"
        ),
        fill=COL_RNASEQ, title_size=11, body_size=9,
        dashed=True,
        title_pad=0.28, body_pad=0.65,
    )

    # Edge-anchored horizontal arrows.
    draw_arrow(ax, fastq.right_at(),     align.left_at())
    draw_arrow(ax, align.right_at(),     align_out.left_at())
    draw_arrow(ax, align_out.right_at(), rpf.left_at())
    draw_arrow(ax, rpf.right_at(),       rpf_out.left_at())
    # rpf outputs feed the optional rnaseq stage. Drop straight down
    # along the same x-column so the dashed link is unambiguous.
    draw_arrow(
        ax,
        rpf_out.bottom_at(),
        rnaseq.top_at(),
        dashed=True, label="optional",
        label_offset=(0.0, 0.0),
    )

    legend_handles = [
        mpatches.Patch(color=COL_INPUT, label="Input"),
        mpatches.Patch(color=COL_ALIGN, label="align stage"),
        mpatches.Patch(color=COL_RPF, label="rpf stage"),
        mpatches.Patch(color=COL_RNASEQ, label="rnaseq stage (optional)"),
        mpatches.Patch(color=COL_OUTPUT, label="output artefact"),
    ]
    ax.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.04),
        ncol=5, frameon=False, fontsize=10,
    )

    draw_footer(ax, x=0.2, y=0.05)
    return save(fig, "01_pipeline_overview.png")


# ---------------------------------------------------------------------------
# 02 -- Align stage internals (top-down)
# ---------------------------------------------------------------------------


def render_align_stage() -> Path:
    fig, ax = setup(figsize=(13, 14), xlim=(0, 13), ylim=(0, 14))

    ax.text(
        6.5, 13.55,
        "align stage  ·  mitoribopy.align",
        ha="center", va="center",
        fontsize=18, fontweight="bold", color=TEXT_COL,
    )
    ax.text(
        6.5, 13.05,
        "Per-sample resolution + cutadapt + bowtie2 + dedup; UMI extracted into the read QNAME before alignment.",
        ha="center", va="center",
        fontsize=11, style="italic", color=SUBTITLE_COL,
    )

    steps = [
        ("FASTQ input",                    "(per sample)",                                                          COL_INPUT,  "user input"),
        ("per-sample resolution",          "kit + UMI + adapter + dedup strategy\nresolved independently per sample",  COL_ALIGN,  "align.sample_resolve"),
        ("cutadapt trim",                  "kit-aware adapter; UMI extracted into\nQNAME (5' single-pass or 3' two-pass)",  COL_ALIGN,  "align.trim"),
        ("contam subtract",                "bowtie2 vs rRNA / tRNA index;\nunaligned reads pass through",            COL_ALIGN,  "align.contam"),
        ("bowtie2 mt-align",               "bowtie2 vs mt-transcriptome\n--end-to-end --very-sensitive -L 18",       COL_ALIGN,  "align.align"),
        ("MAPQ filter",                    "drop reads below --mapq;\nsuppresses NUMT cross-talk",                   COL_ALIGN,  "align.bam_utils"),
        ("dedup",                          "umi-tools (coord + UMI) for UMI samples,\nskip otherwise. mark-duplicates removed in v0.4.5.",  COL_ALIGN,  "align.dedup"),
        ("BAM → BED6",                "strand-aware export\nready for the rpf stage",                            COL_ALIGN,  "align.bam_utils"),
        ("align outputs",                  "bed/<sample>.bed\nread_counts.tsv\nkit_resolution.tsv\nrun_settings.json\n.sample_done/<sample>.json",  COL_OUTPUT, "on disk"),
    ]

    box_w = 7.2
    box_h = 1.10
    spacing = 1.30
    centre_x = (13 - box_w) / 2

    boxes: list[Box] = []
    for i, (title, body, fill, module) in enumerate(steps):
        # Outputs box gets a little more vertical space for its 5 lines.
        h = 1.45 if i == len(steps) - 1 else box_h
        y = 12.0 - i * spacing
        b = draw_box(
            ax, x=centre_x, y=y - h / 2, w=box_w, h=h,
            title=title, body=body, fill=fill,
            title_size=12, body_size=9,
            title_pad=0.22, body_pad=0.55,
        )
        ax.text(
            centre_x + box_w + 0.2, b.cy,
            module,
            ha="left", va="center",
            fontsize=9, style="italic", color=SUBTITLE_COL,
        )
        boxes.append(b)

    for prev, nxt in zip(boxes, boxes[1:]):
        draw_arrow(ax, prev.bottom_at(), nxt.top_at())

    draw_footer(ax, x=0.1, y=0.05)
    return save(fig, "02_align_stage.png")


# ---------------------------------------------------------------------------
# 03 -- rpf stage internals (top-down)
# ---------------------------------------------------------------------------


def render_rpf_stage() -> Path:
    fig, ax = setup(figsize=(15, 19), xlim=(0, 15), ylim=(0, 19))

    ax.text(
        7.5, 18.55,
        "rpf stage  ·  mitoribopy.pipeline + analysis + plotting",
        ha="center", va="center",
        fontsize=18, fontweight="bold", color=TEXT_COL,
    )
    ax.text(
        7.5, 18.05,
        "Offset selection + per-site translation profile + coverage plots; v0.4.5 flat layout with site-prefixed filenames.",
        ha="center", va="center",
        fontsize=11, style="italic", color=SUBTITLE_COL,
    )

    # ----- Backbone (top half) ------------------------------------------
    backbone = [
        ("BED6 input",
         "from align/bed/ or user-supplied",
         COL_INPUT, ""),
        ("filter BED by RPF window",
         "drop reads outside [min, max] nt\n(short / monosome / disome)",
         COL_RPF, "pipeline.steps"),
        ("auto-prune low-count lengths",
         "drop bins below --rpf_min_count_frac × max\n(default 0.20; set 0 to disable)",
         COL_RPF, "pipeline.steps._apply_rpf_count_filter"),
        ("offset enrichment",
         "combined across samples +\nper-sample tables",
         COL_RPF, "analysis.offset_enrichment"),
        ("offset selection",
         "argmax + tie-break per read length\nper sample (drift surfaced)",
         COL_RPF, "analysis.offset_selection"),
    ]

    box_w = 6.5
    box_h = 1.20
    backbone_x = 1.0
    spacing = 1.45
    backbone_boxes: list[Box] = []
    for i, (title, body, fill, module) in enumerate(backbone):
        y = 16.7 - i * spacing
        b = draw_box(
            ax, x=backbone_x, y=y - box_h / 2, w=box_w, h=box_h,
            title=title, body=body, fill=fill,
            title_size=12, body_size=9,
            title_pad=0.24, body_pad=0.58,
        )
        if module:
            ax.text(
                backbone_x + box_w + 0.2, b.cy, module,
                ha="left", va="center",
                fontsize=9, style="italic", color=SUBTITLE_COL,
            )
        backbone_boxes.append(b)

    for prev, nxt in zip(backbone_boxes, backbone_boxes[1:]):
        draw_arrow(ax, prev.bottom_at(), nxt.top_at())

    last_backbone = backbone_boxes[-1]

    # ----- Required outputs row (3 boxes) -------------------------------
    # Below the backbone with clear vertical separation; route through a
    # horizontal "rail" so each branch gets its own straight feeder.
    req_y_top = 7.50
    req_h = 2.20
    tp = draw_box(
        ax, x=0.4, y=req_y_top, w=4.5, h=req_h,
        title="translation_profile/<sample>/",
        body=(
            "footprint_density/\n"
            "  <tr>_footprint_density.csv,\n"
            "  <tr>_{p,a}_site_depth.png\n"
            "translating_frame/\n"
            "  {p,a}_site_frame_usage_*.csv\n"
            "codon_usage/\n"
            "  {p,a}_site_codon_usage_*.csv"
        ),
        fill=COL_RPF, title_size=10, body_size=8,
        title_pad=0.24, body_pad=0.58,
    )
    cov = draw_box(
        ax, x=5.25, y=req_y_top, w=4.5, h=req_h,
        title="coverage_profile_plots/",
        body=(
            "{p_site,a_site}_density_{rpm,raw}/\n"
            "{p_site,a_site}_density_{rpm,raw}_codon/\n"
            "{p_site,a_site}_density_{rpm,raw}_frame/\n"
            "read_coverage_{rpm,raw}{,_codon}/\n"
            "(read_coverage gated by\n"
            "  --read_coverage_{rpm,raw})"
        ),
        fill=COL_RPF, title_size=10, body_size=8,
        title_pad=0.24, body_pad=0.58,
    )
    diag = draw_box(
        ax, x=10.10, y=req_y_top, w=4.5, h=req_h,
        title="offset_diagnostics/",
        body=(
            "csv/\n"
            "  p_site_offsets_*.csv,\n"
            "  offset_*.csv,\n"
            "  per_sample_offset/<sample>/\n"
            "    offset_applied.csv\n"
            "plots/\n"
            "  offset_drift_*.svg, offset_*.svg"
        ),
        fill=COL_RPF, title_size=10, body_size=8,
        title_pad=0.24, body_pad=0.58,
    )

    # Trunk from offset selection down to a horizontal rail.
    fork_y = last_backbone.bottom_at()[1] - 0.55
    rail_y = req_y_top + req_h + 0.55
    draw_arrow(ax, last_backbone.bottom_at(), (last_backbone.cx, fork_y), style="-")
    ax.plot(
        [last_backbone.cx, last_backbone.cx], [fork_y, rail_y],
        color=EDGE_COL, linewidth=1.5, zorder=1,
    )
    ax.plot(
        [tp.cx, diag.cx], [rail_y, rail_y],
        color=EDGE_COL, linewidth=1.5, zorder=1,
    )
    for branch in (tp, cov, diag):
        draw_arrow(ax, (branch.cx, rail_y), branch.top_at())

    # ----- Optional modules row (3 dashed boxes) -----------------------
    opt_y_top = 3.40
    opt_h = 2.20
    cor = draw_box(
        ax, x=0.4, y=opt_y_top, w=4.5, h=opt_h,
        title="codon_correlation/   (--cor_plot)",
        body=(
            "{p_site,a_site}/\n"
            "  <base>_vs_<sample>_{all,masked}.{csv,svg,png}\n"
            "publication-quality scatter;\n"
            "per-site folders produced when\n"
            "analysis_sites=both. Reads from\n"
            "translation_profile codon_usage CSVs."
        ),
        fill=COL_OPTIONAL, title_size=10, body_size=8,
        dashed=True, title_pad=0.24, body_pad=0.58,
    )
    igv = draw_box(
        ax, x=5.25, y=opt_y_top, w=4.5, h=opt_h,
        title="igv_tracks/   (--igv_export)",
        body=(
            "<sample>/\n"
            "  <sample>_p_site.bedgraph\n"
            "  <sample>_a_site.bedgraph\n"
            "BedGraph tracks for IGV;\n"
            "P-site forest green, A-site\n"
            "dark orange. Sourced from\n"
            "translation_profile footprint CSVs."
        ),
        fill=COL_OPTIONAL, title_size=10, body_size=8,
        dashed=True, title_pad=0.24, body_pad=0.58,
    )
    sd = draw_box(
        ax, x=10.10, y=opt_y_top, w=4.5, h=opt_h,
        title="structure_density/   (--structure_density)",
        body=(
            "log2 + percentile-scaled density\n"
            "from translation_profile footprint\n"
            "CSVs; cap controlled by\n"
            "--structure_density_norm_perc.\n"
            "Useful for comparing translation\n"
            "rate vs RNA-structure tracks."
        ),
        fill=COL_OPTIONAL, title_size=10, body_size=8,
        dashed=True, title_pad=0.24, body_pad=0.58,
    )

    # Each optional module reads from translation_profile (the parent
    # required output). Connect each via its own dashed feeder so the
    # data dependency is explicit.
    draw_arrow(ax, tp.bottom_at(0.5), cor.top_at(0.5), dashed=True)
    draw_arrow(ax, tp.bottom_at(0.5), igv.top_at(0.5), dashed=True,
               connection="arc3,rad=-0.10")
    draw_arrow(ax, tp.bottom_at(0.5), sd.top_at(0.5), dashed=True,
               connection="arc3,rad=-0.20")

    # ----- Bottom summary ----------------------------------------------
    out = draw_box(
        ax, x=2.5, y=0.6, w=10.0, h=1.6,
        title="rpf-stage outputs (on disk)",
        body=(
            "rpf_counts.tsv  +  run_settings.json (with reference_checksum)\n"
            "plus every directory in the green/dashed boxes above"
        ),
        fill=COL_OUTPUT, title_size=11, body_size=9,
        title_pad=0.30, body_pad=0.78,
    )
    draw_arrow(ax, cor.bottom_at(),  (out.cx - 3.0, out.y + out.h), dashed=True)
    draw_arrow(ax, igv.bottom_at(),  (out.cx,         out.y + out.h), dashed=True)
    draw_arrow(ax, sd.bottom_at(),   (out.cx + 3.0, out.y + out.h), dashed=True)

    draw_footer(ax, x=0.1, y=0.05)
    return save(fig, "03_rpf_stage.png")


# ---------------------------------------------------------------------------
# 04 -- rnaseq stage internals (top-down)
# ---------------------------------------------------------------------------


def render_rnaseq_stage() -> Path:
    fig, ax = setup(figsize=(13, 12), xlim=(0, 13), ylim=(0, 12))

    ax.text(
        6.5, 11.55,
        "rnaseq stage  ·  mitoribopy.rnaseq  (optional)",
        ha="center", va="center",
        fontsize=18, fontweight="bold", color=TEXT_COL,
    )
    ax.text(
        6.5, 11.05,
        "Translation efficiency + ΔTE from a paired DE table; SHA256 reference-consistency gate prevents mismatched runs.",
        ha="center", va="center",
        fontsize=11, style="italic", color=SUBTITLE_COL,
    )

    de_in = draw_box(
        ax, x=0.5, y=8.7, w=5.0, h=1.7,
        title="DE table (external)",
        body=(
            "DESeq2 / Xtail / Anota2Seq output;\n"
            "format auto-detected by column headers."
        ),
        fill=COL_INPUT, title_size=11, body_size=9,
        title_pad=0.30, body_pad=0.75,
    )
    rpf_in = draw_box(
        ax, x=7.5, y=8.7, w=5.0, h=1.7,
        title="rpf outputs",
        body=(
            "rpf_counts.tsv  +  run_settings.json\n"
            "(reference_checksum recorded)"
        ),
        fill=COL_OUTPUT, title_size=11, body_size=9,
        title_pad=0.30, body_pad=0.75,
    )

    gate = draw_box(
        ax, x=4.0, y=6.5, w=5.0, h=1.5,
        title="SHA256 reference-consistency gate",
        body="hard-fail if Ribo-seq and DE references mismatch",
        fill=COL_RNASEQ, title_size=11, body_size=9,
        title_pad=0.32, body_pad=0.85,
    )
    draw_arrow(ax, de_in.bottom_at(),  gate.top_at(0.25))
    draw_arrow(ax, rpf_in.bottom_at(), gate.top_at(0.75))

    te = draw_box(
        ax, x=4.0, y=4.55, w=5.0, h=1.25,
        title="TE per (sample, gene)",
        body="te.tsv  --  rpf_count, mrna_abundance, te",
        fill=COL_RNASEQ, title_size=11, body_size=9,
        title_pad=0.28, body_pad=0.72,
    )
    draw_arrow(ax, gate.bottom_at(), te.top_at())

    dte = draw_box(
        ax, x=4.0, y=2.6, w=5.0, h=1.25,
        title="ΔTE per gene",
        body="delta_te.tsv  --  mrna_log2fc, rpf_log2fc, delta_te_log2, padj",
        fill=COL_RNASEQ, title_size=11, body_size=9,
        title_pad=0.28, body_pad=0.72,
    )
    draw_arrow(ax, te.bottom_at(), dte.top_at())

    plots = draw_box(
        ax, x=2.0, y=0.6, w=9.0, h=1.4,
        title="diagnostic plots",
        body=(
            "mrna_vs_rpf scatter (4-quadrant log2FC)   |   "
            "delta_te_volcano (ΔTE vs -log10 padj)"
        ),
        fill=COL_OUTPUT, title_size=11, body_size=9,
        title_pad=0.25, body_pad=0.65,
    )
    draw_arrow(ax, dte.bottom_at(), plots.top_at())

    draw_footer(ax, x=0.1, y=0.05)
    return save(fig, "04_rnaseq_stage.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    print("Rendering MitoRiboPy diagrams to", HERE)
    for func in (
        render_pipeline_overview,
        render_align_stage,
        render_rpf_stage,
        render_rnaseq_stage,
    ):
        out = func()
        size_kb = out.stat().st_size // 1024
        print(f"  {out.name}  ({size_kb} KB)")


if __name__ == "__main__":
    main()
