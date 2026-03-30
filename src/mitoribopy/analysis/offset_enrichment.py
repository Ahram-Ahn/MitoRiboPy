"""Offset enrichment analysis and plotting for mitochondrial Ribo-seq."""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ..plotting.style import apply_publication_style


apply_publication_style()


def _apply_site_shift(
    five_offset: int, three_offset: int, offset_site: str
) -> tuple[int, int]:
    """Convert P-site-based offsets to the requested site convention."""
    if offset_site == "a":
        return five_offset + 3, three_offset - 3
    return five_offset, three_offset


def compute_offsets(
    bed_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    align_to: str,
    manual_offset: int | None = None,
    offset_site: str = "p",
    codon_overlap_mode: str = "full",
) -> pd.DataFrame:
    """Compute read-level 5' and 3' offsets relative to start/stop codons."""
    results: list[tuple[int, int, int]] = []

    if offset_site not in {"p", "a"}:
        raise ValueError(f"Unsupported offset_site='{offset_site}'. Use 'p' or 'a'.")
    if codon_overlap_mode not in {"full", "any"}:
        raise ValueError(
            f"Unsupported codon_overlap_mode='{codon_overlap_mode}'. Use 'full' or 'any'."
        )

    if manual_offset is not None:
        for _, read in bed_df.iterrows():
            read_length = read["read_length"]
            results.append((read_length, int(manual_offset), 0))
        return pd.DataFrame(results, columns=["Read Length", "5' Offset", "3' Offset"])

    for _, annotation_row in annotation_df.iterrows():
        transcript_name = annotation_row["transcript"]
        codon_start = (
            int(annotation_row["start_codon"])
            if align_to == "start"
            else int(annotation_row["stop_codon"])
        )
        codon_end = codon_start + 2
        transcript_reads = bed_df[bed_df["chrom"] == transcript_name]

        for _, read in transcript_reads.iterrows():
            read_start = int(read["start"])
            # BED end is exclusive. Convert to inclusive coordinate.
            read_end = int(read["end"]) - 1
            read_length = read["read_length"]
            if read_end < read_start:
                continue

            has_any_overlap = (read_start <= codon_end) and (read_end >= codon_start)
            if not has_any_overlap:
                continue
            if codon_overlap_mode == "full" and not (
                read_start <= codon_start and read_end >= codon_end
            ):
                continue

            if codon_overlap_mode == "full":
                overlap_start = codon_start
                overlap_end = codon_end
            else:
                overlap_start = max(read_start, codon_start)
                overlap_end = min(read_end, codon_end)

            if align_to == "start":
                five_offset = overlap_start - read_start + 1
                three_offset = read_end - overlap_end + 1
            else:
                five_offset = overlap_start - read_start - 2
                three_offset = read_end - overlap_end + 4

            five_offset, three_offset = _apply_site_shift(
                five_offset=five_offset,
                three_offset=three_offset,
                offset_site=offset_site,
            )
            results.append((read_length, five_offset, three_offset))

    return pd.DataFrame(results, columns=["Read Length", "5' Offset", "3' Offset"])


def create_csv_for_offset_enrichment(
    bed_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    align_to: str,
    rpf_range: range,
    output_csv: str,
    offset_limit: int = 23,
    manual_offset: int | None = None,
    offset_site: str = "p",
    codon_overlap_mode: str = "full",
    strain: str = "y",
) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    """Create an offset enrichment summary table and write it as CSV."""
    _ = strain
    offsets_df = compute_offsets(
        bed_df=bed_df,
        annotation_df=annotation_df,
        align_to=align_to,
        manual_offset=manual_offset,
        offset_site=offset_site,
        codon_overlap_mode=codon_overlap_mode,
    )

    if offsets_df.empty:
        print("[create_csv_for_offset_enrichment] No reads overlapping codons.")
        return None, None

    offsets_df["File"] = "All"
    five_range = range(-offset_limit, 0)
    three_range = range(1, offset_limit + 1)
    summary_cols = list(range(-offset_limit, offset_limit + 1))
    summary_df = pd.DataFrame(0, index=rpf_range, columns=summary_cols)
    summary_df.index.name = "Read Length"

    for _, row in offsets_df.iterrows():
        read_length = row["Read Length"]
        five_offset = row["5' Offset"]
        three_offset = row["3' Offset"]
        five_axis_offset = -five_offset

        if five_axis_offset in five_range:
            summary_df.at[read_length, five_axis_offset] += 1
        if three_offset in three_range:
            summary_df.at[read_length, three_offset] += 1

    summary_df.reset_index(inplace=True)
    summary_df.to_csv(output_csv, index=False)
    print(f"[create_csv_for_offset_enrichment] Offset enrichment CSV => {output_csv}")
    return summary_df, offsets_df


def plot_offset_enrichment(
    summary_df: pd.DataFrame | None,
    align_to: str,
    plot_dir: str,
    plot_format: str = "svg",
    x_breaks: list[int] | None = None,
    line_plot_style: str = "combined",
    offset_limit: int = 20,
    p_site_offsets: pd.DataFrame | None = None,
    offset_min: int = 11,
    offset_max: int = 20,
) -> None:
    """Plot offset heatmap and per-read-length line plots."""
    if summary_df is None or summary_df.empty:
        print("[plot_offset_enrichment] No data available for plotting.")
        return

    print("[plot_offset_enrichment] Plotting offset enrichment...")
    os.makedirs(plot_dir, exist_ok=True)

    offset_positions = list(range(-offset_limit, offset_limit + 1))
    melted = summary_df.melt(
        id_vars=["Read Length"],
        value_vars=offset_positions,
        var_name="Offset",
        value_name="Count",
    )
    melted["Offset"] = melted["Offset"].astype(int)

    heatmap_data = melted.pivot_table(
        index="Read Length",
        columns="Offset",
        values="Count",
        aggfunc="sum",
        fill_value=0,
    )
    heatmap_data = heatmap_data.reindex(sorted(heatmap_data.columns), axis=1)

    if x_breaks is None:
        x_breaks = list(range(-offset_limit, offset_limit + 1, 5))

    plt.figure(figsize=(12, 3))
    sns.heatmap(
        heatmap_data,
        cmap="Reds",
        linewidths=0.5,
        cbar_kws={"label": "Read Count"},
        vmin=0,
    )
    plt.title(f"Offset Enrichment ({align_to.capitalize()} Codon)", fontsize=20, fontweight="bold")
    plt.xlabel("Offset Position (nt)", fontsize=16)
    plt.ylabel("Read Length (nt)", fontsize=16)
    plt.xticks(
        ticks=np.arange(len(heatmap_data.columns)) + 0.5,
        labels=heatmap_data.columns,
        rotation=45,
    )
    plt.tight_layout()
    heatmap_filename = f"offset_enrichment_heatmap_{align_to}.{plot_format}"
    plt.savefig(os.path.join(plot_dir, heatmap_filename))
    plt.close()
    print(f"[plot_offset_enrichment] Heatmap saved => {os.path.join(plot_dir, heatmap_filename)}")

    melted["Offset Type"] = np.where(melted["Offset"] < 0, "5'-offset", "3'-offset")

    if line_plot_style == "combined":
        read_lengths = sorted(melted["Read Length"].unique())
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
            data = melted[melted["Read Length"] == read_length]
            sns.lineplot(x="Offset", y="Count", hue="Offset Type", data=data, marker="o", ax=ax)
            ax.set_title(f"Read Length ({read_length} nt)", fontsize=16, fontweight="bold")
            ax.set_xlabel("Offset Position (nt)", fontsize=14)
            ax.set_ylabel("Read Count", fontsize=14)
            ax.grid(True, linestyle="--", alpha=0.5)
            ax.set_xticks(x_breaks)

            if p_site_offsets is not None:
                row = p_site_offsets[p_site_offsets["Read Length"] == read_length]
                if not row.empty:
                    if "Most Enriched 5' Offset" in row:
                        five_offset = -row["Most Enriched 5' Offset"].values[0]
                        if offset_min <= abs(five_offset) <= offset_max:
                            ax.axvline(
                                x=five_offset,
                                color="black",
                                linestyle="--",
                                label=f"5'-offset: {five_offset}",
                            )
                    if "Most Enriched 3' Offset" in row:
                        three_offset = row["Most Enriched 3' Offset"].values[0]
                        if offset_min <= abs(three_offset) <= offset_max:
                            ax.axvline(
                                x=three_offset,
                                color="black",
                                linestyle="--",
                                label=f"3'-offset: {three_offset}",
                            )
            ax.legend()

        plt.tight_layout()
        lineplot_filename = f"offset_enrichment_combined_lineplots_{align_to}.{plot_format}"
        plt.savefig(os.path.join(plot_dir, lineplot_filename))
        plt.close()
        print(
            "[plot_offset_enrichment] Combined line plots saved => "
            f"{os.path.join(plot_dir, lineplot_filename)}"
        )
    elif line_plot_style == "separate":
        for read_length in sorted(melted["Read Length"].unique()):
            data = melted[melted["Read Length"] == read_length]
            plt.figure(figsize=(10, 5))
            sns.lineplot(x="Offset", y="Count", hue="Offset Type", data=data, marker="o")

            plt.title(
                f"Offset Enrichment (Read Len {read_length} nt) vs {align_to.capitalize()} Codon",
                fontsize=16,
            )
            plt.xlabel("Offset Position (nt)", fontsize=14)
            plt.ylabel("Read Count", fontsize=14)
            plt.xticks(ticks=x_breaks)
            plt.grid(True, linestyle="--", alpha=0.5)

            if p_site_offsets is not None:
                row = p_site_offsets[p_site_offsets["Read Length"] == read_length]
                if not row.empty:
                    if "Most Enriched 5' Offset" in row:
                        five_offset = -row["Most Enriched 5' Offset"].values[0]
                        if offset_min <= abs(five_offset) <= offset_max:
                            plt.axvline(
                                x=five_offset,
                                color="black",
                                linestyle="--",
                                label=f"5'-offset: {five_offset}",
                            )
                    if "Most Enriched 3' Offset" in row:
                        three_offset = row["Most Enriched 3' Offset"].values[0]
                        if offset_min <= abs(three_offset) <= offset_max:
                            plt.axvline(
                                x=three_offset,
                                color="black",
                                linestyle="--",
                                label=f"3'-offset: {three_offset}",
                            )

            plt.legend()
            lineplot_filename = (
                f"offset_enrichment_lineplot_read_{read_length}_{align_to}.{plot_format}"
            )
            plt.savefig(os.path.join(plot_dir, lineplot_filename))
            plt.close()
            print(
                f"[plot_offset_enrichment] Line plot for read_len {read_length} => "
                f"{lineplot_filename}"
            )
    else:
        print("[plot_offset_enrichment] Invalid line_plot_style. Choose 'combined' or 'separate'.")

    print(f"[plot_offset_enrichment] All offset plots saved to {plot_dir}")

