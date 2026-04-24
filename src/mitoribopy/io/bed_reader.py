"""BED ingestion, read-length summaries, and heatmap plotting helpers."""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
from ..console import iter_with_progress, log_info, log_warning


def process_bed_files(
    input_dir: str,
    output_dir: str,
    organism: str,
    annotation_df: pd.DataFrame,
    rpf_range: range | list[int],
    *,
    bam_mapq: int = 0,
) -> tuple[pd.DataFrame, list[str]]:
    """Filter BED / BAM inputs by read length and return in-memory per-sample data.

    Accepts two kinds of input files side-by-side in ``input_dir``:

    * ``*.bed`` - BED3 / BED4 / BED6 files (pre-existing path; unchanged).
    * ``*.bam`` - coordinate-sorted or unsorted BAM files. Each is
      converted to BED6 via
      :func:`mitoribopy.io.bam_reader.prepare_bam_inputs` into
      ``<output_dir>/bam_converted/``. The converted BEDs then flow
      through the same filtering + summary logic as native BEDs.

    ``bam_mapq`` (default ``0`` = off) applies a pysam MAPQ filter
    before the BAM -> BED6 conversion. Values around 10 are the
    recommended NUMT-suppression default used by ``mitoribopy align``.
    """
    _ = organism
    _ = annotation_df

    from ..plotting.visualization import plot_read_length_distribution
    from .bam_reader import prepare_bam_inputs

    converted_bed_paths = prepare_bam_inputs(
        input_dir=Path(input_dir),
        converted_dir=Path(output_dir) / "bam_converted",
        mapq_threshold=bam_mapq,
    )

    native_bed_files = sorted(
        file_name for file_name in os.listdir(input_dir) if file_name.endswith(".bed")
    )
    bed_sources: list[tuple[str, str]] = []  # (display_name, absolute_path)
    for bed_file in native_bed_files:
        bed_sources.append((bed_file, os.path.join(input_dir, bed_file)))
    for bed_path in converted_bed_paths:
        # Converted BEDs sit in a separate directory so their names do not
        # collide with the native-BED filenames above; their sample names
        # come from the BAM stem.
        bed_sources.append((bed_path.name, str(bed_path)))

    # Deterministic processing order by display name so sample ordering is
    # reproducible regardless of whether input was BAM, BED, or mixed.
    bed_sources.sort(key=lambda item: item[0])

    if not bed_sources:
        log_warning("BED", f"No .bed or .bam files found in {input_dir}.")
        return pd.DataFrame(), []

    log_info(
        "BED",
        f"Found {len(bed_sources)} input file(s) in {input_dir} "
        f"(native BED: {len(native_bed_files)}, BAM-converted: "
        f"{len(converted_bed_paths)}) for filtering range "
        f"{list(rpf_range)}.",
    )

    all_filtered_bed: list[pd.DataFrame] = []
    sample_names: list[str] = []
    read_length_summary: list[dict[str, object]] = []

    for bed_file, bed_file_path in iter_with_progress(
        bed_sources,
        component="BED",
        noun="BED file",
        labeler=lambda item: item[0],
    ):
        if bed_file.startswith("filtered_"):
            log_info("BED", f"Skipping previously filtered file: {bed_file}")
            continue
        try:
            bed_df = pd.read_csv(bed_file_path, sep="\t", header=None)
        except Exception as exc:
            log_warning("BED", f"Error reading {bed_file}: {exc}")
            continue

        if bed_df.shape[1] >= 3:
            bed_columns = ["chrom", "start", "end", "name", "score", "strand"][: bed_df.shape[1]]
        else:
            log_warning("BED", f"Unexpected BED format in {bed_file}; skipping.")
            continue

        bed_df.columns = bed_columns
        bed_df["start"] = pd.to_numeric(bed_df["start"], errors="coerce")
        bed_df["end"] = pd.to_numeric(bed_df["end"], errors="coerce")
        bed_df.dropna(subset=["start", "end"], inplace=True)
        bed_df["start"] = bed_df["start"].astype(int)
        bed_df["end"] = bed_df["end"].astype(int)
        # Drop malformed rows (end <= start) before computing read_length.
        malformed_mask = bed_df["end"] <= bed_df["start"]
        if malformed_mask.any():
            log_warning(
                "BED",
                f"{bed_file}: dropping {int(malformed_mask.sum())} BED row(s) with "
                "end <= start (malformed coordinates).",
            )
            bed_df = bed_df.loc[~malformed_mask].copy()
        bed_df["read_length"] = bed_df["end"] - bed_df["start"]

        filtered_bed_df = bed_df[bed_df["read_length"].isin(rpf_range)].copy()
        if filtered_bed_df.empty:
            log_warning("BED", f"No reads of requested lengths in {bed_file}; skipping.")
            continue

        sample_name = os.path.splitext(bed_file)[0]
        filtered_bed_df["sample_name"] = sample_name
        sample_names.append(sample_name)

        # Unfiltered histogram (full range found in the BED file) so the
        # plot can show the RPF window against the overall shape, not just
        # within the window itself.
        unfiltered_counts = bed_df["read_length"].value_counts().sort_index()
        plot_path = os.path.join(output_dir, f"{sample_name}_read_length_distribution.svg")
        plot_read_length_distribution(
            unfiltered_counts, sample_name, plot_path, rpf_range=rpf_range
        )
        # The filtered (RPF-window) histogram is still used for the summary CSV.
        read_length_counts = filtered_bed_df["read_length"].value_counts().sort_index()

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
        log_info("BED", f"Filtered read-length summary saved => {summary_csv_path}")
        return concatenated_bed, sample_names

    log_warning("BED", "No reads remained after length filtering; returning empty dataset.")
    return pd.DataFrame(), []


def compute_unfiltered_read_length_summary(
    input_dir: str,
    output_csv: str,
    total_counts_map: dict[str, int],
    read_length_range: tuple[int, int] = (15, 40),
) -> None:
    """Build unfiltered per-sample read-length summary table with optional RPM."""
    bed_files = sorted(
        file_name for file_name in os.listdir(input_dir) if file_name.endswith(".bed")
    )
    if not bed_files:
        log_warning("QC", f"No .bed files found in {input_dir}.")
        return

    all_rows: list[dict[str, int | float | str]] = []
    missing_norm_samples: set[str] = set()

    for bed_file in iter_with_progress(
        bed_files,
        component="QC",
        noun="BED file",
        labeler=str,
    ):
        bed_path = os.path.join(input_dir, bed_file)
        sample_name = os.path.splitext(bed_file)[0]

        try:
            df = pd.read_csv(bed_path, sep="\t", header=None)
        except Exception as exc:
            log_warning("QC", f"Error reading {bed_file}: {exc}")
            continue

        if df.shape[1] >= 3:
            colnames = ["chrom", "start", "end", "name", "score", "strand"][: df.shape[1]]
        elif df.shape[1] == 2:
            colnames = ["chrom", "start"]
        else:
            log_warning("QC", f"Unexpected BED format in {bed_file}; skipping.")
            continue

        df.columns = colnames
        df["start"] = pd.to_numeric(df["start"], errors="coerce")
        df["end"] = pd.to_numeric(df["end"], errors="coerce")
        df.dropna(subset=["start", "end"], inplace=True)
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)
        # Defensive: drop malformed rows where end <= start before computing
        # read_length. Such rows would produce zero or negative lengths that
        # could silently skew downstream QC.
        malformed_mask = df["end"] <= df["start"]
        if malformed_mask.any():
            log_warning(
                "QC",
                f"{bed_file}: dropping {int(malformed_mask.sum())} BED row(s) with "
                "end <= start (malformed coordinates).",
            )
            df = df.loc[~malformed_mask].copy()
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
        log_warning("QC", "No reads remained after unfiltered read-length QC filtering.")
        return

    result_df = pd.DataFrame(all_rows)
    result_df.to_csv(output_csv, index=False)
    log_info("QC", f"Unfiltered read-length summary saved => {output_csv}")
    if missing_norm_samples:
        missing_list = ", ".join(sorted(missing_norm_samples))
        log_warning(
            "QC",
            "No total read-count entry found for sample(s): "
            f"{missing_list}. RPM was set to 0 for these sample(s).",
        )


def plot_unfiltered_read_length_heatmap(
    summary_csv_path: str,
    output_png_base: str,
    value_col: str = "count",
) -> None:
    """Backward-compatible wrapper for the visualization module."""
    from ..plotting.visualization import plot_unfiltered_read_length_heatmap as render_heatmap

    return render_heatmap(
        summary_csv_path=summary_csv_path,
        output_png_base=output_png_base,
        value_col=value_col,
    )
