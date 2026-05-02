"""Plot writers for the metagene Fourier spectrum.

Three figures land per `(sample, read_length)` under
``<output>/qc/fourier_spectrum/<sample>/``:

* ``<sample>_<length>nt_combined.{png,svg}`` — combined-canonical
  metagene (every transcript that is NOT in the ATP8/ATP6 or ND4L/ND4
  overlap pair). Two stacked panels (orf_start top, orf_stop bottom),
  ONE aggregated trace per panel.
* ``<sample>_<length>nt_ATP86.{png,svg}`` — junction-bracketed
  bicistronic analysis. Top panel: ATP6 frame (orf_start window of
  ATP6, leaving the bicistronic junction). Bottom panel: ATP8 frame
  (orf_stop window of ATP8, leading INTO the junction). The two
  windows OVERLAP at the bicistronic junction so a period-3 peak in
  each panel says which reading frame is dominant in that region.
* ``<sample>_<length>nt_ND4L4.{png,svg}`` — same idea for the
  ND4L/ND4 overlap pair. The 4-nt overlap is too short to fall inside
  both windows; they flank the junction instead.

Each figure carries:
* The canonical per-plot ``.metadata.json`` sidecar (Job 1 contract)
  so ``mitoribopy validate-figures`` can score it without re-running
  matplotlib. The sidecar's ``panel_layout`` field records which
  transcript / region drives each panel.
* A vertical reference line at period = 3 nt.
* In-figure annotations with the spectral_ratio_3nt and snr_call for
  each panel — so a reviewer can read the headline number without
  opening the TSV.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..analysis.fourier_spectrum import (
    GENE_SETS,
    REGIONS,
    _ATP86_PANEL_TRANSCRIPT,
    _ND4L4_PANEL_TRANSCRIPT,
    snr_call_for_ratio,
)
from .figure_validator import write_plot_metadata
from .style import apply_publication_style


apply_publication_style()


__all__ = [
    "render_fourier_spectrum_panels",
]


# Trace colours per gene_set — Okabe-Ito-friendly, distinguishable.
_TRACE_COLOR: dict[str, str] = {
    "combined": "#0072B2",  # blue
    "ATP86":    "#D55E00",  # vermilion
    "ND4L4":    "#009E73",  # bluish green
}

# Per-region panel titles. For ATP86/ND4L4 we override with the
# transcript-specific title in the panel-rendering code.
_DEFAULT_PANEL_TITLES: dict[str, str] = {
    "orf_start": "Downstream of start codon (post-initiation)",
    "orf_stop":  "Upstream of stop codon (pre-termination)",
}

_SNR_TIER_DESCRIPTION: dict[str, str] = {
    "excellent":  "ratio >= 10x",
    "healthy":    "ratio >= 5x",
    "modest":     "ratio >= 2x",
    "broken":     "ratio < 2x",
    "no_signal":  "no data",
}


def _panel_title(*, gene_set: str, region: str) -> str:
    """Per-(gene_set, region) panel title for the plot."""
    if gene_set == "ATP86":
        tx = _ATP86_PANEL_TRANSCRIPT.get(region, "?")
        if region == "orf_start":
            return f"ATP6 frame — downstream of {tx} start codon"
        return f"ATP8 frame — upstream of {tx} stop codon"
    if gene_set == "ND4L4":
        tx = _ND4L4_PANEL_TRANSCRIPT.get(region, "?")
        if region == "orf_start":
            return f"ND4 frame — downstream of {tx} start codon"
        return f"ND4L frame — upstream of {tx} stop codon"
    return _DEFAULT_PANEL_TITLES.get(region, region)


def _plot_one_panel(
    ax,
    *,
    spectrum: pd.DataFrame,  # cols: period_nt, amplitude
    score_row: pd.Series | None,
    title: str,
    color: str,
) -> int:
    """Draw the metagene amplitude curve on *ax*. Returns n_points drawn."""
    if spectrum is None or spectrum.empty:
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("Period [nt]")
        ax.set_ylabel("Normalized amplitude")
        ax.text(
            0.5, 0.5, "no data", transform=ax.transAxes,
            ha="center", va="center", color="0.6", fontsize=10,
        )
        return 0

    sub = spectrum.sort_values("period_nt")
    periods = sub["period_nt"].to_numpy()
    amps = sub["amplitude"].to_numpy()
    ax.plot(periods, amps, color=color, linewidth=1.5)
    ax.fill_between(periods, 0, amps, color=color, alpha=0.10, linewidth=0)

    # Reference line at period = 3 nt (the headline frequency).
    ax.axvline(3.0, color="0.4", linestyle="--", linewidth=0.8, alpha=0.7)

    ax.set_title(title, fontsize=11)
    ax.set_xlabel("Period [nt]")
    ax.set_ylabel("Normalized amplitude")
    ax.set_xticks([2, 3, 4, 5, 6, 7, 8, 9, 10])
    ax.set_xlim(2.0, 10.0)
    if amps.size > 0:
        amp_max = float(np.nanmax(amps))
        if np.isfinite(amp_max) and amp_max > 0:
            ax.set_ylim(0, amp_max * 1.15)
    ax.grid(True, alpha=0.25, linewidth=0.5)

    if score_row is not None:
        ratio = float(score_row.get("spectral_ratio_3nt", float("nan")))
        snr_tier = str(score_row.get("snr_call", "no_signal"))
        ratio_local = float(score_row.get("spectral_ratio_3nt_local", float("nan")))
        snr_tier_local = str(score_row.get("snr_call_local", "no_signal"))
        n_genes = int(score_row.get("n_genes", 0))
        n_sites = int(score_row.get("n_sites_total", 0))
        text = (
            f"3-nt ratio (global, p=2-10): {ratio:.2f}x  ({snr_tier})\n"
            f"3-nt ratio (local, p=4-6):   {ratio_local:.2f}x  ({snr_tier_local})\n"
            f"n genes = {n_genes}  |  n sites = {n_sites:,}"
        )
        ax.text(
            0.97, 0.93, text,
            transform=ax.transAxes,
            ha="right", va="top", fontsize=8, family="monospace",
            bbox=dict(facecolor="white", edgecolor="0.7", alpha=0.85, pad=4),
        )
    return int(periods.size)


def _render_two_panel_figure(
    *,
    sample: str,
    read_length: int,
    gene_set: str,
    spectrum_block: pd.DataFrame,
    score_block: pd.DataFrame,
    out_path: Path,
    source_data: str,
) -> int | None:
    """Two-panel (orf_start top, orf_stop bottom) figure for one gene_set.

    Returns ``n_points_drawn`` total across both panels, or None if the
    block is empty (no figure written).
    """
    if spectrum_block.empty:
        return None

    fig, axes = plt.subplots(2, 1, figsize=(6, 6), sharex=True)
    color = _TRACE_COLOR.get(gene_set, "#444444")

    n_drawn = 0
    panel_layout: list[dict] = []
    for ax, region in zip(axes, REGIONS):
        sub_spec = spectrum_block[spectrum_block["region"] == region]
        sub_score = score_block[score_block["region"] == region]
        score_row = sub_score.iloc[0] if not sub_score.empty else None
        title = _panel_title(gene_set=gene_set, region=region)
        n_drawn += _plot_one_panel(
            ax,
            spectrum=sub_spec,
            score_row=score_row,
            title=title,
            color=color,
        )
        # Record what drove this panel — for the sidecar.
        if gene_set == "ATP86":
            tx = _ATP86_PANEL_TRANSCRIPT.get(region)
        elif gene_set == "ND4L4":
            tx = _ND4L4_PANEL_TRANSCRIPT.get(region)
        else:
            tx = "<combined>"
        panel_layout.append({
            "region": region,
            "transcript": tx,
            "n_genes": int(score_row["n_genes"]) if score_row is not None else 0,
            "n_sites_total": int(score_row["n_sites_total"]) if score_row is not None else 0,
            "spectral_ratio_3nt": (
                float(score_row["spectral_ratio_3nt"])
                if score_row is not None and pd.notna(score_row.get("spectral_ratio_3nt"))
                else None
            ),
            "spectral_ratio_3nt_local": (
                float(score_row["spectral_ratio_3nt_local"])
                if score_row is not None and pd.notna(score_row.get("spectral_ratio_3nt_local"))
                else None
            ),
        })

    fig.suptitle(
        f"{sample}  |  {read_length} nt  |  gene_set={gene_set}",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1.0, 0.96))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    svg_path = out_path.with_suffix(".svg")
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    write_plot_metadata(
        out_path,
        plot_type=f"fourier_spectrum_{gene_set}",
        stage="rpf",
        source_data=source_data,
        n_points_expected=n_drawn,
        n_points_drawn=n_drawn,
        formats=["png", "svg"],
        dpi=300,
        extra={
            "sample": str(sample),
            "read_length": int(read_length),
            "gene_set": str(gene_set),
            "panel_layout": panel_layout,
        },
    )
    return n_drawn


def render_fourier_spectrum_panels(
    spectrum_combined_table: pd.DataFrame,
    score_combined_table: pd.DataFrame,
    *,
    output_dir: Path | str,
    source_data_relpath: str = "fourier_spectrum_combined.tsv",
) -> list[Path]:
    """Render every (sample, read_length, gene_set) two-panel figure.

    Writes per-sample subdirectories under *output_dir* with up to three
    figures per (sample, read_length): one each for gene_set in
    {combined, ATP86, ND4L4}. Skips figures whose corresponding gene_set
    is empty.

    Returns the list of PNG paths written.
    """
    output_dir = Path(output_dir)
    written: list[Path] = []
    if spectrum_combined_table is None or spectrum_combined_table.empty:
        return written

    group_cols = ["sample", "read_length"]
    for (sample, read_length), block in spectrum_combined_table.groupby(group_cols):
        sample_dir = output_dir / str(sample)
        score_block_all = (
            score_combined_table[
                (score_combined_table["sample"] == sample)
                & (score_combined_table["read_length"] == int(read_length))
            ]
            if score_combined_table is not None and not score_combined_table.empty
            else pd.DataFrame()
        )
        for gene_set in GENE_SETS:
            sub_spec = block[block["gene_set"] == gene_set]
            sub_score = (
                score_block_all[score_block_all["gene_set"] == gene_set]
                if not score_block_all.empty else pd.DataFrame()
            )
            if sub_spec.empty:
                continue
            png = sample_dir / f"{sample}_{int(read_length)}nt_{gene_set}.png"
            n = _render_two_panel_figure(
                sample=str(sample),
                read_length=int(read_length),
                gene_set=str(gene_set),
                spectrum_block=sub_spec,
                score_block=sub_score,
                out_path=png,
                source_data=source_data_relpath,
            )
            if n is not None:
                written.append(png)

    return written
