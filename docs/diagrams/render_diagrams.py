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
"""

from __future__ import annotations

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


HERE = Path(__file__).resolve().parent


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
) -> tuple[float, float]:
    """Draw a rounded rectangle with a bold title and an optional body.

    Returns the (x, y) of the box centre so callers can wire arrows
    between centres without recomputing.
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
            cx, y + h - 0.18,
            title,
            ha="center", va="top",
            fontsize=title_size, fontweight="bold", color=TEXT_COL,
        )
        ax.text(
            cx, y + h - 0.55,
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
    return cx, cy


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
) -> None:
    arrow = FancyArrowPatch(
        src, dst,
        arrowstyle=style,
        mutation_scale=14,
        linewidth=1.4,
        linestyle="--" if dashed else "-",
        color=color,
    )
    ax.add_patch(arrow)
    if label:
        mx = (src[0] + dst[0]) / 2 + label_offset[0]
        my = (src[1] + dst[1]) / 2 + label_offset[1]
        ax.text(
            mx, my, label,
            ha="center", va="center",
            fontsize=9, color=TEXT_COL,
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                      edgecolor="none", alpha=0.85),
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
# 01 -- Pipeline overview (horizontal, 3 stages, ~12 nodes max)
# ---------------------------------------------------------------------------


def render_pipeline_overview() -> Path:
    fig, ax = setup(figsize=(18, 6), xlim=(0, 18), ylim=(0, 6))

    # Title
    ax.text(
        9, 5.6,
        "MitoRiboPy pipeline overview",
        ha="center", va="center",
        fontsize=18, fontweight="bold", color=TEXT_COL,
    )

    # Inputs
    fastq = draw_box(
        ax, x=0.2, y=2.3, w=2.0, h=1.5,
        title="FASTQ", body="(one or more samples;\nmixed kits OK)",
        fill=COL_INPUT, title_size=14, body_size=9,
    )

    # Align stage
    align = draw_box(
        ax, x=2.8, y=1.5, w=4.2, h=3.0,
        title="align",
        body=(
            "cutadapt trim + UMI -> QNAME\n"
            "bowtie2 contam subtract\n"
            "bowtie2 mt-transcriptome align\n"
            "MAPQ filter (NUMT suppression)\n"
            "umi-tools dedup or skip\n"
            "BAM -> BED6"
        ),
        fill=COL_ALIGN, title_size=15, body_size=10,
    )

    # Align outputs (intermediate artefact box)
    align_out = draw_box(
        ax, x=7.4, y=2.3, w=2.0, h=1.5,
        title="align output",
        body="bed/, kit_resolution.tsv,\nread_counts.tsv",
        fill=COL_OUTPUT, title_size=11, body_size=8,
    )

    # RPF stage
    rpf = draw_box(
        ax, x=10.0, y=1.5, w=4.2, h=3.0,
        title="rpf",
        body=(
            "filter BED on RPF window\n"
            "offset enrichment (combined +\n"
            "  per-sample) -> selection\n"
            "translation_profile/{p,a}/...\n"
            "coverage_profile_plots/...\n"
            "(optional) structure_density,\n"
            "  codon_correlation"
        ),
        fill=COL_RPF, title_size=15, body_size=10,
    )

    # RPF outputs
    rpf_out = draw_box(
        ax, x=14.6, y=2.3, w=1.7, h=1.5,
        title="rpf output",
        body="rpf_counts.tsv,\noffsets, plots",
        fill=COL_OUTPUT, title_size=11, body_size=8,
    )

    # rnaseq (optional, dashed)
    rnaseq = draw_box(
        ax, x=14.6, y=0.1, w=3.2, h=1.7,
        title="rnaseq (optional)",
        body=(
            "DE table + rpf_counts ->\n"
            "SHA256 ref-gate -> TE / dTE\n"
            "scatter + volcano"
        ),
        fill=COL_RNASEQ, title_size=11, body_size=9,
        dashed=True,
    )

    # Arrows
    draw_arrow(ax, fastq, (align[0] - 2.1, align[1]))
    draw_arrow(ax, (align[0] + 2.1, align[1]), (align_out[0] - 1.0, align_out[1]))
    draw_arrow(ax, (align_out[0] + 1.0, align_out[1]), (rpf[0] - 2.1, rpf[1]))
    draw_arrow(ax, (rpf[0] + 2.1, rpf[1]), (rpf_out[0] - 0.85, rpf_out[1]))
    # rpf_out -> rnaseq (dashed = optional)
    draw_arrow(
        ax,
        (rpf_out[0], rpf_out[1] - 0.75),
        (rnaseq[0], rnaseq[1] + 0.85),
        dashed=True, label="optional",
        label_offset=(0.0, 0.0),
    )

    # Legend strip at bottom
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

    return save(fig, "01_pipeline_overview.png")


# ---------------------------------------------------------------------------
# 02 -- Align stage internals (top-down)
# ---------------------------------------------------------------------------


def render_align_stage() -> Path:
    fig, ax = setup(figsize=(13, 14), xlim=(0, 13), ylim=(0, 14))

    ax.text(
        6.5, 13.5,
        "align stage  (mitoribopy.align)",
        ha="center", va="center",
        fontsize=18, fontweight="bold", color=TEXT_COL,
    )
    ax.text(
        6.5, 13.0,
        "Per-sample resolution + cutadapt + bowtie2 + dedup; UMI extracted into the read QNAME",
        ha="center", va="center",
        fontsize=11, style="italic", color="#555555",
    )

    # Steps top-to-bottom. Each step: (title, body, module label)
    steps = [
        ("FASTQ input",                    "(per sample)",                                                   COL_INPUT,  "user input"),
        ("per-sample resolution",          "kit + UMI + adapter + dedup strategy\nresolved independently per sample", COL_ALIGN,  "align.sample_resolve"),
        ("cutadapt trim",                  "kit-aware adapter; UMI extracted into\nQNAME (5' single-pass or 3' two-pass)", COL_ALIGN,  "align.trim"),
        ("contam subtract",                "bowtie2 vs rRNA / tRNA index;\nunaligned reads pass through",   COL_ALIGN,  "align.contam"),
        ("bowtie2 mt-align",               "bowtie2 vs mt-transcriptome\n--end-to-end --very-sensitive -L 18", COL_ALIGN,  "align.align"),
        ("MAPQ filter",                    "drop reads below --mapq;\nsuppresses NUMT cross-talk",          COL_ALIGN,  "align.bam_utils"),
        ("dedup",                          "umi-tools collapse on coord + UMI\nor SKIP for no-UMI samples", COL_ALIGN,  "align.dedup"),
        ("BAM -> BED6",                    "strand-aware export\nready for the rpf stage",                  COL_ALIGN,  "align.bam_utils"),
        ("align outputs",                  "bed/<sample>.bed\nread_counts.tsv\nkit_resolution.tsv\nrun_settings.json\n.sample_done/<sample>.json", COL_OUTPUT, "on disk"),
    ]

    box_w = 7.0
    box_h = 1.05
    centre_x = (13 - box_w) / 2
    spacing = 1.30  # vertical spacing between box centres

    centres = []
    for i, (title, body, fill, module) in enumerate(steps):
        y = 12.0 - i * spacing
        cx, cy = draw_box(
            ax, x=centre_x, y=y - box_h / 2, w=box_w, h=box_h,
            title=title, body=body, fill=fill,
            title_size=12, body_size=9,
        )
        # Module label in italics to the right
        ax.text(
            centre_x + box_w + 0.2, cy,
            module,
            ha="left", va="center",
            fontsize=9, style="italic", color="#555555",
        )
        centres.append((cx, cy))

    # Wire arrows top -> bottom
    for prev, nxt in zip(centres, centres[1:]):
        draw_arrow(ax, (prev[0], prev[1] - box_h / 2),
                       (nxt[0], nxt[1] + box_h / 2))

    return save(fig, "02_align_stage.png")


# ---------------------------------------------------------------------------
# 03 -- rpf stage internals (top-down)
# ---------------------------------------------------------------------------


def render_rpf_stage() -> Path:
    fig, ax = setup(figsize=(15, 14), xlim=(0, 15), ylim=(0, 14))

    ax.text(
        7.5, 13.5,
        "rpf stage  (mitoribopy.pipeline + analysis + plotting)",
        ha="center", va="center",
        fontsize=18, fontweight="bold", color=TEXT_COL,
    )
    ax.text(
        7.5, 13.0,
        "Offset selection + per-site translation profile + coverage plots",
        ha="center", va="center",
        fontsize=11, style="italic", color="#555555",
    )

    # Linear backbone (top -> bottom)
    backbone = [
        ("BED6 input",            "from align/bed/ or user-supplied",                               COL_INPUT,  ""),
        ("filter BED by RPF window", "drop reads outside [min, max] nt\n(short / monosome / disome)", COL_RPF,   "pipeline.steps"),
        ("offset enrichment",     "combined across samples +\nper-sample tables",                  COL_RPF,    "analysis.offset_enrichment"),
        ("offset selection",      "argmax + tie-break per read length\nper sample",                COL_RPF,    "analysis.offset_selection"),
    ]

    box_w = 6.5
    box_h = 1.05
    backbone_x = 1.0
    spacing = 1.30
    backbone_centres = []
    for i, (title, body, fill, module) in enumerate(backbone):
        y = 11.5 - i * spacing
        cx, cy = draw_box(
            ax, x=backbone_x, y=y - box_h / 2, w=box_w, h=box_h,
            title=title, body=body, fill=fill,
            title_size=12, body_size=9,
        )
        if module:
            ax.text(
                backbone_x + box_w + 0.2, cy, module,
                ha="left", va="center",
                fontsize=9, style="italic", color="#555555",
            )
        backbone_centres.append((cx, cy))

    for prev, nxt in zip(backbone_centres, backbone_centres[1:]):
        draw_arrow(ax, (prev[0], prev[1] - box_h / 2),
                       (nxt[0], nxt[1] + box_h / 2))

    # Branching outputs at the bottom of the backbone
    branch_top_y = backbone_centres[-1][1] - box_h / 2
    branch_y = branch_top_y - 1.0

    # Three branches: translation_profile, coverage_profile_plots, optional
    tp_x = 0.5
    cov_x = 5.5
    opt_x = 10.5

    tp = draw_box(
        ax, x=tp_x, y=branch_y - 1.6, w=4.0, h=1.6,
        title="translation_profile/<sample>/{p,a}/",
        body=(
            "footprint_density/    *_footprint_density.csv\n"
            "translating_frame/    frame_usage_*.csv\n"
            "codon_usage/          codon_usage_*.csv"
        ),
        fill=COL_RPF, title_size=10, body_size=8,
    )
    cov = draw_box(
        ax, x=cov_x, y=branch_y - 1.6, w=4.5, h=1.6,
        title="coverage_profile_plots/",
        body=(
            "read_coverage_*/         (site-independent)\n"
            "{p,a}/site_density_*/    (per-site density,\n"
            "                          incl. CDS-frame plots)"
        ),
        fill=COL_RPF, title_size=10, body_size=8,
    )
    opt = draw_box(
        ax, x=opt_x, y=branch_y - 1.6, w=4.0, h=1.6,
        title="optional modules",
        body=(
            "structure_density/   --structure_density\n"
            "codon_correlation/   --cor_plot --base_sample\n"
            "(SVG + 300 dpi PNG)"
        ),
        fill=COL_OPTIONAL, title_size=10, body_size=8,
        dashed=True,
    )

    # Connect branches
    fork = (backbone_centres[-1][0], branch_top_y - 0.4)
    draw_arrow(ax, (backbone_centres[-1][0], branch_top_y), fork, style="-")
    for branch_centre in (tp, cov, opt):
        draw_arrow(ax, fork, (branch_centre[0], branch_y - 0.0 + 0.0))

    # Output artefact at the very bottom
    out = draw_box(
        ax, x=4.5, y=0.3, w=6.0, h=1.0,
        title="rpf outputs",
        body="rpf_counts.tsv  +  run_settings.json (with reference_checksum)",
        fill=COL_OUTPUT, title_size=11, body_size=9,
    )
    draw_arrow(ax, (tp[0], branch_y - 1.6),  (out[0] - 1.5, out[1] + 0.5))
    draw_arrow(ax, (cov[0], branch_y - 1.6), (out[0],       out[1] + 0.5))
    draw_arrow(ax, (opt[0], branch_y - 1.6), (out[0] + 1.5, out[1] + 0.5),
               dashed=True)

    return save(fig, "03_rpf_stage.png")


# ---------------------------------------------------------------------------
# 04 -- rnaseq stage internals (top-down, tighter)
# ---------------------------------------------------------------------------


def render_rnaseq_stage() -> Path:
    fig, ax = setup(figsize=(13, 12), xlim=(0, 13), ylim=(0, 12))

    ax.text(
        6.5, 11.5,
        "rnaseq stage  (mitoribopy.rnaseq) -- optional",
        ha="center", va="center",
        fontsize=18, fontweight="bold", color=TEXT_COL,
    )
    ax.text(
        6.5, 11.0,
        "Translation efficiency + delta-TE from a paired DE table; "
        "SHA256 reference gate prevents mismatched runs.",
        ha="center", va="center",
        fontsize=11, style="italic", color="#555555",
    )

    # Two parallel inputs at the top
    de_in = draw_box(
        ax, x=0.5, y=8.6, w=5.0, h=1.6,
        title="DE table (external)",
        body=(
            "DESeq2 / Xtail / Anota2Seq output\n"
            "auto-detected by column headers"
        ),
        fill=COL_INPUT, title_size=11, body_size=9,
    )
    rpf_in = draw_box(
        ax, x=7.5, y=8.6, w=5.0, h=1.6,
        title="rpf outputs",
        body="rpf_counts.tsv  +  run_settings.json\n(reference_checksum recorded)",
        fill=COL_OUTPUT, title_size=11, body_size=9,
    )

    # Reference gate
    gate = draw_box(
        ax, x=4.0, y=6.6, w=5.0, h=1.4,
        title="SHA256 reference-consistency gate",
        body="hard-fail if Ribo-seq and DE references mismatch",
        fill=COL_RNASEQ, title_size=11, body_size=9,
    )
    draw_arrow(ax, (de_in[0],  de_in[1] - 0.8),  (gate[0] - 1.5, gate[1] + 0.7))
    draw_arrow(ax, (rpf_in[0], rpf_in[1] - 0.8), (gate[0] + 1.5, gate[1] + 0.7))

    # TE
    te = draw_box(
        ax, x=4.0, y=4.6, w=5.0, h=1.2,
        title="TE per (sample, gene)",
        body="te.tsv  --  rpf_count, mrna_abundance, te",
        fill=COL_RNASEQ, title_size=11, body_size=9,
    )
    draw_arrow(ax, (gate[0], gate[1] - 0.7), (te[0], te[1] + 0.6))

    # delta-TE
    dte = draw_box(
        ax, x=4.0, y=2.7, w=5.0, h=1.2,
        title="delta-TE per gene",
        body="delta_te.tsv  --  mrna_log2fc, rpf_log2fc, delta_te_log2, padj",
        fill=COL_RNASEQ, title_size=11, body_size=9,
    )
    draw_arrow(ax, (te[0], te[1] - 0.6), (dte[0], dte[1] + 0.6))

    # Plots
    plots = draw_box(
        ax, x=2.0, y=0.6, w=9.0, h=1.4,
        title="diagnostic plots",
        body=(
            "mrna_vs_rpf scatter (4-quadrant log2FC)   |   "
            "delta_te_volcano (delta-TE vs -log10 padj)"
        ),
        fill=COL_OUTPUT, title_size=11, body_size=9,
    )
    draw_arrow(ax, (dte[0], dte[1] - 0.6), (plots[0], plots[1] + 0.7))

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
