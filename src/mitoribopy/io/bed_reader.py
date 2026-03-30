"""BED ingestion, read-length summaries, and heatmap plotting helpers."""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
import seaborn as sns

from ..plotting.style import apply_publication_style


apply_publication_style()


def process_bed_files(
    input_dir: str,
    output_dir: str,
    organism: str,
    annotation_df: pd.DataFrame,
    rpf_range: range | list[int],
) -> tuple[pd.DataFrame, list[str]]:
    """Filter BED files by read length, split by length, and write per-sample outputs."""
    _ = organism
    _ = annotation_df

    bed_files = [file_name for file_name in os.listdir(input_dir) if file_name.endswith(".bed")]
    if not bed_files:
        print("No .bed files found in the directory.")
        return pd.DataFrame(), []

    print(
        f"Found {len(bed_files)} BED files to process in '{input_dir}' "
        f"for filtering range {rpf_range}."
    )

    filtered_bed_dir = os.path.join(output_dir, "filtered_bed")
    os.makedirs(filtered_bed_dir, exist_ok=True)

    all_filtered_bed: list[pd.DataFrame] = []
    sample_dirs: list[str] = []
    read_length_summary: list[dict[str, object]] = []

    for bed_file in bed_files:
        if bed_file.startswith("filtered_"):
            print(f"Skipping already filtered file: {bed_file}")
            continue

        bed_file_path = os.path.join(input_dir, bed_file)
        try:
            bed_df = pd.read_csv(bed_file_path, sep="\t", header=None)
        except Exception as exc:
            print(f"Error reading {bed_file}: {exc}")
            continue

        if bed_df.shape[1] >= 3:
            bed_columns = ["chrom", "start", "end", "name", "score", "strand"][: bed_df.shape[1]]
        else:
            print(f"Unexpected BED file format in {bed_file}.")
            continue

        bed_df.columns = bed_columns
        bed_df["start"] = pd.to_numeric(bed_df["start"], errors="coerce")
        bed_df["end"] = pd.to_numeric(bed_df["end"], errors="coerce")
        bed_df.dropna(subset=["start", "end"], inplace=True)
        bed_df["start"] = bed_df["start"].astype(int)
        bed_df["end"] = bed_df["end"].astype(int)
        bed_df["read_length"] = bed_df["end"] - bed_df["start"]

        filtered_bed_df = bed_df[bed_df["read_length"].isin(rpf_range)].copy()
        if filtered_bed_df.empty:
            print(f"No reads of specified lengths in {bed_file} for rpf_range={rpf_range}.")
            continue

        sample_name = os.path.splitext(bed_file)[0]
        sample_dir = os.path.join(filtered_bed_dir, sample_name)
        os.makedirs(sample_dir, exist_ok=True)
        sample_dirs.append(sample_dir)

        for read_length in sorted(filtered_bed_df["read_length"].unique()):
            sub_df = filtered_bed_df[filtered_bed_df["read_length"] == read_length]
            out_bed_name = f"{sample_name}_{read_length}nt.bed"
            out_bed_path = os.path.join(sample_dir, out_bed_name)
            sub_df.to_csv(out_bed_path, sep="\t", header=False, index=False)

        read_length_counts = filtered_bed_df["read_length"].value_counts().sort_index()
        plot_path = os.path.join(output_dir, f"{sample_name}_read_length_distribution.svg")

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
        plt.savefig(plot_path)
        plt.close()
        print(f"[process_bed_files] Read length distribution plot saved => {plot_path}")

        read_length_summary.append(
            {
                "sample_name": sample_name,
                "read_length_counts": read_length_counts.to_dict(),
            }
        )
        all_filtered_bed.append(filtered_bed_df)

    if all_filtered_bed:
        concatenated_bed = pd.concat(all_filtered_bed, ignore_index=True)
        summary_csv_path = os.path.join(output_dir, "read_length_summary.csv")

        summary_rows: list[dict[str, int | str]] = []
        for entry in read_length_summary:
            sample = str(entry["sample_name"])
            for read_length, count in dict(entry["read_length_counts"]).items():
                summary_rows.append(
                    {
                        "sample_name": sample,
                        "read_length": int(read_length),
                        "count": int(count),
                    }
                )

        summary_df = pd.DataFrame(summary_rows)
        summary_df.to_csv(summary_csv_path, index=False)
        print(f"[process_bed_files] Filtered read length summary => {summary_csv_path}")
        return concatenated_bed, sample_dirs

    print("[process_bed_files] No reads found in the specified range => returning empty.")
    return pd.DataFrame(), []


def compute_unfiltered_read_length_summary(
    input_dir: str,
    output_csv: str,
    total_counts_map: dict[str, int],
    read_length_range: tuple[int, int] = (15, 40),
) -> None:
    """Build unfiltered per-sample read-length summary table with optional RPM."""
    bed_files = [file_name for file_name in os.listdir(input_dir) if file_name.endswith(".bed")]
    if not bed_files:
        print(f"[compute_unfiltered_read_length_summary] No .bed files in {input_dir}.")
        return

    all_rows: list[dict[str, int | float | str]] = []
    missing_norm_samples: set[str] = set()

    for bed_file in bed_files:
        bed_path = os.path.join(input_dir, bed_file)
        sample_name = os.path.splitext(bed_file)[0]

        try:
            df = pd.read_csv(bed_path, sep="\t", header=None)
        except Exception as exc:
            print(f"[compute_unfiltered_read_length_summary] Error reading {bed_file}: {exc}")
            continue

        if df.shape[1] >= 3:
            colnames = ["chrom", "start", "end", "name", "score", "strand"][: df.shape[1]]
        elif df.shape[1] == 2:
            colnames = ["chrom", "start"]
        else:
            print(
                f"[compute_unfiltered_read_length_summary] Unexpected BED format in {bed_file}, skipping."
            )
            continue

        df.columns = colnames
        df["start"] = pd.to_numeric(df["start"], errors="coerce")
        df["end"] = pd.to_numeric(df["end"], errors="coerce")
        df.dropna(subset=["start", "end"], inplace=True)
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)
        df["read_length"] = df["end"] - df["start"]

        mask = (df["read_length"] >= read_length_range[0]) & (
            df["read_length"] <= read_length_range[1]
        )
        filtered_df = df[mask].copy()
        if filtered_df.empty:
            continue

        length_counts = filtered_df["read_length"].value_counts()

        total_reads = total_counts_map.get(sample_name, None)
        if total_reads is None:
            total_reads = total_counts_map.get(sample_name.replace(".bed", ""), None)
        if total_reads is None:
            total_reads = total_counts_map.get(sample_name.lower(), None)
        if total_reads is None:
            total_reads = total_counts_map.get(sample_name.replace(".bed", "").lower(), None)
        if total_reads is None:
            missing_norm_samples.add(sample_name)

        for read_length, count in length_counts.items():
            if total_reads is not None and total_reads > 0:
                rpm = (count / total_reads) * 1e6
            else:
                rpm = 0
            all_rows.append(
                {
                    "sample_name": sample_name,
                    "read_length": int(read_length),
                    "count": int(count),
                    "RPM": float(rpm),
                }
            )

    if not all_rows:
        print("[compute_unfiltered_read_length_summary] No data after filtering by length range.")
        return

    result_df = pd.DataFrame(all_rows)
    result_df.to_csv(output_csv, index=False)
    print(f"[compute_unfiltered_read_length_summary] Unfiltered summary (25-50 nt) saved => {output_csv}")
    if missing_norm_samples:
        missing_list = ", ".join(sorted(missing_norm_samples))
        print(
            "[compute_unfiltered_read_length_summary] WARNING: No total read-count entry "
            f"for sample(s): {missing_list}. RPM was set to 0 for these sample(s)."
        )


def plot_unfiltered_read_length_heatmap(
    summary_csv_path: str,
    output_png_base: str,
    value_col: str = "count",
) -> None:
    """Plot unfiltered read-length heatmap from summary table."""
    if not os.path.exists(summary_csv_path):
        print(f"[plot_unfiltered_read_length_heatmap] {summary_csv_path} not found => skip heatmap.")
        return

    df = pd.read_csv(summary_csv_path)
    if df.empty:
        print("[plot_unfiltered_read_length_heatmap] No data in summary => skip.")
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
    print(f"[plot_unfiltered_read_length_heatmap] Heatmap saved => {out_png}")

