"""Unified plotting entrypoints for MitoRiboPy outputs."""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import seaborn as sns

from ..console import log_info, log_warning
from .coverage_profile_plots import run_coverage_profile_plots
from .style import apply_publication_style
from .structure_density_export import run_structure_density_export

apply_publication_style()

FIVE_PRIME_COLOR = "#1f77b4"
THREE_PRIME_COLOR = "#ff7f0e"


def _validate_offset_mask_nt(offset_mask_nt: int) -> int:
    """Validate the near-anchor masking window."""
    offset_mask_nt = int(offset_mask_nt)
    if offset_mask_nt < 0:
        raise ValueError("offset_mask_nt must be zero or a positive integer.")
    return offset_mask_nt


def _is_masked_plot_offset(offset: int, offset_mask_nt: int) -> bool:
    """Return True when an offset lies in the masked near-anchor window."""
    return 0 < abs(int(offset)) <= offset_mask_nt


def plot_read_length_distribution(
    read_length_counts: pd.Series,
    sample_name: str,
    output_path: str | os.PathLike[str],
) -> None:
    """Render a per-sample read-length distribution bar chart."""
    plt.figure(figsize=(6, 4))
    plt.bar(read_length_counts.index, read_length_counts.values, width=0.6, color="skyblue")
    plt.xlabel("Read Length (nt)", fontsize=14, fontname="Arial", fontweight="bold")
    plt.ylabel("Read Count", fontsize=14, fontname="Arial")
    plt.title(
        f"Read Length Distribution ({sample_name})",
        fontsize=16,
        fontweight="bold",
        fontname="Arial",
    )
    plt.xticks(
        read_length_counts.index,
        rotation=45,
        fontsize=10,
        fontweight="bold",
        fontname="Arial",
    )
    plt.gca().yaxis.set_major_locator(MaxNLocator(nbins=4))
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    log_info("VIS", f"Read-length distribution plot saved => {output_path}")


def plot_unfiltered_read_length_heatmap(
    summary_csv_path: str,
    output_png_base: str,
    value_col: str = "count",
) -> None:
    """Plot unfiltered read-length heatmap from summary table."""
    if not os.path.exists(summary_csv_path):
        log_warning("VIS", f"{summary_csv_path} not found; skipping unfiltered heatmap.")
        return

    df = pd.read_csv(summary_csv_path)
    if df.empty:
        log_warning("VIS", "No unfiltered summary data found; skipping heatmap.")
        return

    pivot_df = df.pivot(index="sample_name", columns="read_length", values=value_col).fillna(0)

    plt.figure(figsize=(10, max(3, len(pivot_df) * 0.6)))
    ax = sns.heatmap(
        pivot_df,
        cmap=("Blues" if value_col.upper() == "RPM" else "Greens"),
        linewidths=0.5,
        cbar_kws={"label": value_col.upper()},
    )
    plt.title(
        f"Read-Length Heatmap ({value_col.upper()})",
        fontsize=16,
        fontweight="bold",
        fontname="Arial",
    )
    plt.xlabel("Read Length (nt)", fontsize=14, fontname="Arial")
    plt.ylabel("Sample Name", fontsize=14, fontname="Arial")

    cbar = ax.collections[0].colorbar
    cbar.ax.set_ylabel(value_col.upper(), fontsize=12, fontname="Arial")
    cbar.ax.tick_params(labelsize=10)

    plt.tight_layout()

    out_png = f"{output_png_base}_{value_col}.svg"
    plt.savefig(out_png)
    plt.close()
    log_info("VIS", f"Unfiltered read-length heatmap saved => {out_png}")


def _selected_offset_guides(
    selected_offset_row: pd.Series,
    plotted_offsets: set[int],
) -> list[dict[str, object]]:
    """Build vertical-guide definitions using the plotted axis convention."""
    guides: list[dict[str, object]] = []
    read_length = int(selected_offset_row["Read Length"])

    five_value = selected_offset_row.get("Most Enriched 5' Offset")
    if pd.notna(five_value):
        five_value = int(five_value)
        five_axis_offset = -five_value
        if five_axis_offset in plotted_offsets:
            guides.append(
                {
                    "axis_offset": five_axis_offset,
                    "color": FIVE_PRIME_COLOR,
                    "label": f"Selected 5' offset ({five_value} nt)",
                }
            )
        else:
            log_warning(
                "VIS",
                "Skipping 5' guide line for read length "
                f"{read_length}: plotted range does not include axis offset {five_axis_offset}.",
            )

    three_value = selected_offset_row.get("Most Enriched 3' Offset")
    if pd.notna(three_value):
        three_value = int(three_value)
        if three_value in plotted_offsets:
            guides.append(
                {
                    "axis_offset": three_value,
                    "color": THREE_PRIME_COLOR,
                    "label": f"Selected 3' offset ({three_value} nt)",
                }
            )
        else:
            log_warning(
                "VIS",
                "Skipping 3' guide line for read length "
                f"{read_length}: plotted range does not include axis offset {three_value}.",
            )

    return guides


def _add_selected_offset_guides(
    ax,
    data: pd.DataFrame,
    selected_offsets: pd.DataFrame | None,
    *,
    read_length: int,
) -> None:
    """Add selected-offset guide lines to a per-read-length line plot."""
    if selected_offsets is None:
        return

    row = selected_offsets[selected_offsets["Read Length"] == read_length]
    if row.empty:
        return

    guides = _selected_offset_guides(
        row.iloc[0],
        set(data["Offset"].astype(int)),
    )
    for guide in guides:
        ax.axvline(
            x=int(guide["axis_offset"]),
            color=str(guide["color"]),
            linestyle="--",
            linewidth=1.5,
            label=str(guide["label"]),
        )


def plot_offset_enrichment(
    summary_df: pd.DataFrame | None,
    align_to: str,
    plot_dir: str,
    plot_format: str = "svg",
    x_breaks: list[int] | None = None,
    line_plot_style: str = "combined",
    offset_limit: int = 20,
    offset_mask_nt: int = 5,
    selected_offsets: pd.DataFrame | None = None,
    offset_min: int = 11,
    offset_max: int = 20,
    five_offset_min: int | None = None,
    five_offset_max: int | None = None,
    three_offset_min: int | None = None,
    three_offset_max: int | None = None,
) -> None:
    """Plot offset heatmap and per-read-length line plots."""
    if summary_df is None or summary_df.empty:
        log_warning("VIS", "No offset-enrichment data available for plotting.")
        return

    offset_mask_nt = _validate_offset_mask_nt(offset_mask_nt)
    five_offset_min = offset_min if five_offset_min is None else int(five_offset_min)
    five_offset_max = offset_max if five_offset_max is None else int(five_offset_max)
    three_offset_min = offset_min if three_offset_min is None else int(three_offset_min)
    three_offset_max = offset_max if three_offset_max is None else int(three_offset_max)
    os.makedirs(plot_dir, exist_ok=True)
    log_info(
        "VIS",
        "Rendering offset enrichment plots "
        f"(align={align_to}, format={plot_format}, masked_nt={offset_mask_nt}).",
    )

    offset_positions = list(range(-offset_limit, offset_limit + 1))
    melted = summary_df.melt(
        id_vars=["Read Length"],
        value_vars=offset_positions,
        var_name="Offset",
        value_name="Count",
    )
    melted["Offset"] = melted["Offset"].astype(int)
    melted["Masked"] = melted["Offset"].apply(
        lambda offset: _is_masked_plot_offset(offset, offset_mask_nt)
    )

    plot_melted = melted.copy()
    plot_melted.loc[plot_melted["Masked"], "Count"] = np.nan

    heatmap_data = summary_df.set_index("Read Length")[offset_positions].astype(float)
    masked_columns = [
        offset for offset in offset_positions if _is_masked_plot_offset(offset, offset_mask_nt)
    ]
    if masked_columns:
        heatmap_data.loc[:, masked_columns] = np.nan
    heatmap_data = heatmap_data.reindex(sorted(heatmap_data.columns), axis=1)

    if x_breaks is None:
        x_breaks = list(range(-offset_limit, offset_limit + 1, 5))

    cmap = sns.color_palette("Reds", as_cmap=True)
    cmap.set_bad("#d9d9d9")

    plt.figure(figsize=(12, 3))
    sns.heatmap(
        heatmap_data,
        cmap=cmap,
        linewidths=0.5,
        cbar_kws={"label": "Read Count"},
        vmin=0,
        mask=heatmap_data.isna(),
    )
    plt.title(f"Offset Enrichment ({align_to.capitalize()} Codon)", fontsize=20, fontweight="bold")
    plt.xlabel("Offset Position (nt)", fontsize=16)
    plt.ylabel("Read Length (nt)", fontsize=16)
    if offset_mask_nt > 0:
        plt.figtext(
            0.99,
            0.01,
            f"Masked near-anchor bins: -{offset_mask_nt}..-1, +1..+{offset_mask_nt} nt",
            ha="right",
            fontsize=9,
        )
    plt.xticks(
        ticks=np.arange(len(heatmap_data.columns)) + 0.5,
        labels=heatmap_data.columns,
        rotation=45,
    )
    plt.tight_layout()
    heatmap_filename = f"offset_enrichment_heatmap_{align_to}.{plot_format}"
    heatmap_path = os.path.join(plot_dir, heatmap_filename)
    plt.savefig(heatmap_path)
    plt.close()
    log_info("VIS", f"Offset heatmap saved => {heatmap_path}")

    line_melted = plot_melted[plot_melted["Offset"] != 0].copy()
    line_melted["Offset Type"] = np.where(
        line_melted["Offset"] < 0, "5'-offset", "3'-offset"
    )
    if line_melted.empty:
        log_warning("VIS", "No non-zero offsets remain after masking; skipping line plots.")
        return

    line_palette = {"5'-offset": FIVE_PRIME_COLOR, "3'-offset": THREE_PRIME_COLOR}

    if line_plot_style == "combined":
        read_lengths = sorted(line_melted["Read Length"].unique())
        num_read_lengths = len(read_lengths)
        fig, axes = plt.subplots(
            nrows=num_read_lengths,
            ncols=1,
            figsize=(10, 3 * num_read_lengths),
            sharex=True,
        )
        if num_read_lengths == 1:
            axes = [axes]

        for ax, read_length in zip(axes, read_lengths):
            data = line_melted[line_melted["Read Length"] == read_length]
            sns.lineplot(
                x="Offset",
                y="Count",
                hue="Offset Type",
                data=data,
                marker="o",
                palette=line_palette,
                ax=ax,
            )
            ax.set_title(f"Read Length ({read_length} nt)", fontsize=16, fontweight="bold")
            ax.set_xlabel("Offset Position (nt)", fontsize=14)
            ax.set_ylabel("Read Count", fontsize=14)
            ax.grid(True, linestyle="--", alpha=0.5)
            ax.set_xticks(x_breaks)
            if offset_mask_nt > 0:
                ax.axvspan(-offset_mask_nt, offset_mask_nt, color="#d9d9d9", alpha=0.2)
            _add_selected_offset_guides(
                ax,
                data,
                selected_offsets,
                read_length=read_length,
            )
            ax.legend()

        plt.tight_layout()
        lineplot_filename = f"offset_enrichment_combined_lineplots_{align_to}.{plot_format}"
        lineplot_path = os.path.join(plot_dir, lineplot_filename)
        plt.savefig(lineplot_path)
        plt.close()
        log_info("VIS", f"Combined offset line plots saved => {lineplot_path}")
    elif line_plot_style == "separate":
        for read_length in sorted(line_melted["Read Length"].unique()):
            data = line_melted[line_melted["Read Length"] == read_length]
            fig, ax = plt.subplots(figsize=(10, 5))
            sns.lineplot(
                x="Offset",
                y="Count",
                hue="Offset Type",
                data=data,
                marker="o",
                palette=line_palette,
                ax=ax,
            )

            ax.set_title(
                f"Offset Enrichment (Read Len {read_length} nt) vs {align_to.capitalize()} Codon",
                fontsize=16,
            )
            ax.set_xlabel("Offset Position (nt)", fontsize=14)
            ax.set_ylabel("Read Count", fontsize=14)
            ax.set_xticks(ticks=x_breaks)
            ax.grid(True, linestyle="--", alpha=0.5)
            if offset_mask_nt > 0:
                ax.axvspan(-offset_mask_nt, offset_mask_nt, color="#d9d9d9", alpha=0.2)
            _add_selected_offset_guides(
                ax,
                data,
                selected_offsets,
                read_length=read_length,
            )
            ax.legend()

            lineplot_filename = (
                f"offset_enrichment_lineplot_read_{read_length}_{align_to}.{plot_format}"
            )
            lineplot_path = os.path.join(plot_dir, lineplot_filename)
            plt.savefig(lineplot_path)
            plt.close()
            log_info("VIS", f"Offset line plot for read length {read_length} saved => {lineplot_path}")
    else:
        raise ValueError("Invalid line_plot_style. Choose 'combined' or 'separate'.")

    log_info("VIS", f"All offset plots saved to {plot_dir}")


__all__ = [
    "plot_offset_enrichment",
    "plot_read_length_distribution",
    "plot_unfiltered_read_length_heatmap",
    "run_coverage_profile_plots",
    "run_structure_density_export",
]
