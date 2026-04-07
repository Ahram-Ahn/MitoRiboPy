"""Offset enrichment analysis and plotting for mitochondrial Ribo-seq."""

from __future__ import annotations

import pandas as pd
from ..console import log_info, log_warning
from ..data import resolve_sequence_name


def _validate_offset_mask_nt(offset_mask_nt: int) -> int:
    """Validate the near-anchor masking window."""
    offset_mask_nt = int(offset_mask_nt)
    if offset_mask_nt < 0:
        raise ValueError("offset_mask_nt must be zero or a positive integer.")
    return offset_mask_nt


def _is_masked_plot_offset(offset: int, offset_mask_nt: int) -> bool:
    """Return True when an offset lies in the masked near-anchor window."""
    return 0 < abs(int(offset)) <= offset_mask_nt


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

    available_sequence_names = bed_df["chrom"].dropna().astype(str).unique()
    for _, annotation_row in annotation_df.iterrows():
        transcript_name = resolve_sequence_name(annotation_row, available_sequence_names)
        if transcript_name is None:
            continue
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
    offset_mask_nt: int = 5,
    offset_site: str = "p",
    codon_overlap_mode: str = "full",
    strain: str = "y",
) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    """Create an offset enrichment summary table and write it as CSV."""
    _ = strain
    offset_mask_nt = _validate_offset_mask_nt(offset_mask_nt)
    offsets_df = compute_offsets(
        bed_df=bed_df,
        annotation_df=annotation_df,
        align_to=align_to,
        manual_offset=manual_offset,
        offset_site=offset_site,
        codon_overlap_mode=codon_overlap_mode,
    )

    if offsets_df.empty:
        log_warning("OFFSET", "No reads overlapping the anchor codons were found.")
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

        if five_axis_offset in five_range and not _is_masked_plot_offset(
            five_axis_offset, offset_mask_nt
        ):
            summary_df.at[read_length, five_axis_offset] += 1
        if three_offset in three_range and not _is_masked_plot_offset(
            three_offset, offset_mask_nt
        ):
            summary_df.at[read_length, three_offset] += 1

    summary_df.reset_index(inplace=True)
    summary_df.to_csv(output_csv, index=False)
    if offset_mask_nt > 0:
        log_info(
            "OFFSET",
            "Masked near-anchor offsets "
            f"-{offset_mask_nt}..-1 and +1..+{offset_mask_nt} nt from enrichment counts.",
        )
    log_info("OFFSET", f"Offset enrichment CSV saved => {output_csv}")
    return summary_df, offsets_df


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
    """Backward-compatible wrapper for the visualization module."""
    from ..plotting.visualization import plot_offset_enrichment as render_offset_enrichment

    return render_offset_enrichment(
        summary_df=summary_df,
        align_to=align_to,
        plot_dir=plot_dir,
        plot_format=plot_format,
        x_breaks=x_breaks,
        line_plot_style=line_plot_style,
        offset_limit=offset_limit,
        offset_mask_nt=offset_mask_nt,
        selected_offsets=selected_offsets,
        offset_min=offset_min,
        offset_max=offset_max,
        five_offset_min=five_offset_min,
        five_offset_max=five_offset_max,
        three_offset_min=three_offset_min,
        three_offset_max=three_offset_max,
    )
