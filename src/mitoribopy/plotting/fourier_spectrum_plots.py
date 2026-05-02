"""Plot writers for the Wakigawa-style Fourier-spectrum analysis.

For each (sample, read_length) we render a stacked two-panel figure
(ORF on top, 3' UTR on bottom) — one trace per gene, x = period [nt],
y = amplitude. Two figures get written:

* ``<sample>_<read_length>nt_combined.{png,svg}`` — overlap-upstream
  genes (MT-ATP8, MT-ND4L) are EXCLUDED. This is the publication-style
  panel a reviewer reads first.
* ``<sample>_<read_length>nt_overlap_upstream.{png,svg}`` — only the
  overlap-upstream genes. Skipped when none are present in the input.
  Surfaced separately so reviewers can audit them deliberately
  ("analyze ATP8/ND4L, just don't combine them with the rest").

Each rendered PNG carries the canonical per-plot ``.metadata.json``
sidecar (Job 1 contract) so ``mitoribopy validate-figures`` can score
the file without re-running matplotlib.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .figure_validator import write_plot_metadata
from .style import apply_publication_style


apply_publication_style()


__all__ = [
    "render_fourier_spectrum_panels",
]


# Palette — Okabe-Ito plus a few extras so the typical 11-12 mt-mRNA
# overlay never has to recycle a colour. Keep the order stable so
# repeated runs colour the same gene the same way.
_GENE_PALETTE: tuple[str, ...] = (
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7", "#999999", "#000000", "#882255",
    "#117733", "#88CCEE", "#AA4499",
)


def _colour_for_gene(gene: str, ordered_genes: list[str]) -> str:
    try:
        idx = ordered_genes.index(gene)
    except ValueError:
        idx = abs(hash(gene)) % len(_GENE_PALETTE)
    return _GENE_PALETTE[idx % len(_GENE_PALETTE)]


def _plot_one_panel(
    ax,
    region_table: pd.DataFrame,
    *,
    title: str,
    ordered_genes: list[str],
) -> int:
    """Draw the per-gene amplitude curves on *ax*.

    Returns the count of trace points actually drawn (sum of
    ``len(period_nt)`` across all gene curves) so the caller can record
    it in the per-plot sidecar.
    """
    n_drawn = 0
    if region_table.empty:
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("Period [nt]")
        ax.set_ylabel("Amplitude")
        ax.text(
            0.5, 0.5, "no data", transform=ax.transAxes,
            ha="center", va="center", color="0.6", fontsize=10,
        )
        return n_drawn

    for gene, sub in region_table.groupby("gene", sort=False):
        sub = sub.sort_values("period_nt")
        ax.plot(
            sub["period_nt"].to_numpy(),
            sub["amplitude"].to_numpy(),
            color=_colour_for_gene(str(gene), ordered_genes),
            linewidth=1.2,
            label=str(gene),
        )
        n_drawn += int(len(sub))
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("Period [nt]")
    ax.set_ylabel("Amplitude")
    # Period axis: integer ticks across the published 2-10 range.
    period_lo = float(region_table["period_nt"].min())
    period_hi = float(region_table["period_nt"].max())
    integer_ticks = list(range(
        int(np.ceil(period_lo)), int(np.floor(period_hi)) + 1,
    ))
    if integer_ticks:
        ax.set_xticks(integer_ticks)
    ax.grid(True, alpha=0.25, linewidth=0.5)
    return n_drawn


def _render_two_panel(
    spectrum_subset: pd.DataFrame,
    *,
    sample: str,
    read_length: int,
    out_path: Path,
    title_suffix: str,
    source_data: str,
    plot_type: str,
    ordered_genes: list[str],
) -> int | None:
    """Write a stacked two-panel (ORF / 3' UTR) figure to *out_path*.

    Returns the per-plot ``n_points_drawn`` (integer) or ``None`` when
    the subset is empty (no figure written).
    """
    if spectrum_subset.empty:
        return None

    fig, axes = plt.subplots(2, 1, figsize=(6, 6), sharex=True)
    ax_orf, ax_utr = axes

    orf = spectrum_subset[spectrum_subset["region"] == "orf"]
    utr = spectrum_subset[spectrum_subset["region"] == "utr3"]

    n_orf = _plot_one_panel(
        ax_orf, orf,
        title=f"ORF (upstream of stop codon){title_suffix}",
        ordered_genes=ordered_genes,
    )
    n_utr = _plot_one_panel(
        ax_utr, utr,
        title="3' UTR (downstream of stop codon)",
        ordered_genes=ordered_genes,
    )

    # Shared legend on the right of the upper panel — keeps the data
    # area clean and matches Wakigawa's figure convention.
    handles, labels = ax_orf.get_legend_handles_labels()
    if handles:
        ax_orf.legend(
            handles, labels,
            loc="upper left", bbox_to_anchor=(1.02, 1.0),
            fontsize=8, frameon=False, borderaxespad=0.0,
        )

    fig.suptitle(
        f"{sample}  |  {read_length} nt", fontsize=12, fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 0.85, 0.96))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    svg_path = out_path.with_suffix(".svg")
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    n_drawn = n_orf + n_utr
    write_plot_metadata(
        out_path,
        plot_type=plot_type,
        stage="rpf",
        source_data=source_data,
        n_points_expected=n_drawn,
        n_points_drawn=n_drawn,
        formats=["png", "svg"],
        dpi=300,
        extra={
            "sample": str(sample),
            "read_length": int(read_length),
            "n_genes": int(spectrum_subset["gene"].nunique()),
        },
    )
    return n_drawn


def render_fourier_spectrum_panels(
    spectrum_table: pd.DataFrame,
    *,
    output_dir: Path | str,
    source_data_relpath: str = "fourier_spectrum.tsv",
) -> list[Path]:
    """Render every (sample, read_length) two-panel overlay figure.

    Writes per-sample subdirectories under *output_dir* with two
    figures per (sample, read_length): a ``*_combined.png`` (overlap-
    upstream genes excluded) and, when applicable, a
    ``*_overlap_upstream.png`` (only the excluded genes).

    Returns the list of PNG paths written.
    """
    output_dir = Path(output_dir)
    written: list[Path] = []
    if spectrum_table is None or spectrum_table.empty:
        return written

    # Stable gene-colour order across plots (alphabetical).
    ordered_genes = sorted(spectrum_table["gene"].astype(str).unique())

    group_cols = ["sample", "read_length"]
    for (sample, read_length), block in spectrum_table.groupby(group_cols):
        sample_dir = output_dir / str(sample)
        combined = block[~block["is_overlap_upstream_orf"]]
        overlap = block[block["is_overlap_upstream_orf"]]

        combined_png = sample_dir / f"{sample}_{int(read_length)}nt_combined.png"
        n = _render_two_panel(
            combined,
            sample=str(sample),
            read_length=int(read_length),
            out_path=combined_png,
            title_suffix="",
            source_data=source_data_relpath,
            plot_type="fourier_spectrum_combined",
            ordered_genes=ordered_genes,
        )
        if n is not None:
            written.append(combined_png)

        if not overlap.empty:
            overlap_png = (
                sample_dir / f"{sample}_{int(read_length)}nt_overlap_upstream.png"
            )
            n = _render_two_panel(
                overlap,
                sample=str(sample),
                read_length=int(read_length),
                out_path=overlap_png,
                title_suffix=" — overlap-upstream ORFs (analysed only, NOT combined)",
                source_data=source_data_relpath,
                plot_type="fourier_spectrum_overlap_upstream",
                ordered_genes=ordered_genes,
            )
            if n is not None:
                written.append(overlap_png)

    return written
