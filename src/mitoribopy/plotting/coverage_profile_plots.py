"""
Generates coverage-profile plots (one figure per transcript) for each sample.

Two sets of figures are produced for every transcript in every run:

1. **RPM‑scaled** coverage (read‑depth normalised to 1 million mRNA‑aligned reads)
2. **Raw‑count** coverage (no scaling at all)

Each set contains a *read‑coverage* panel and a *P‑site* panel.
All samples share the same y‑axis per transcript so that visual
comparisons are not confounded by different scale bars.

The directory layout after a successful run looks like this::

    <output_dir>/
        read_coverage_rpm/   *.<plot_format>
        p_site_coverage_rpm/ *.<plot_format>
        read_coverage_raw/   *.<plot_format>
        p_site_coverage_raw/ *.<plot_format>

Usage
-----
Call the exported helper from *main.py*::

    run_coverage_profile_plots(
        sample_dirs,
        selected_offsets_df,
        offset_type,
        fasta_file,
        output_dir,
        args,
        annotation_df,
    )

`sample_dirs` must already be in the desired plotting order.
`selected_offsets_df` should have the columns:
    * Read Length
    * Most Enriched 5' Offset
    * Most Enriched 3' Offset

The bar width is fixed at 3 nt; coverage can be capped at the percentile
specified by `args.cap_percentile` (defaults to 0.999).
"""


from __future__ import annotations

import os
import re
import copy
from pathlib import Path
from typing import List, Tuple, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from ..console import iter_with_progress, log_info
from ..data import build_sequence_display_map
from .style import apply_publication_style

apply_publication_style()

# -----------------------------------------------------------------------------
# Public entry point
# -----------------------------------------------------------------------------

def run_coverage_profile_plots(
    sample_dirs: List[str],
    selected_offsets_df: pd.DataFrame,
    offset_type: str,
    fasta_file: str | Path,
    output_dir: str | Path,
    args,
    annotation_df: pd.DataFrame,
    filtered_bed_df: pd.DataFrame | None = None,
):
    """Create RPM and raw-count coverage-profile plots.

    Parameters
    ----------
    sample_dirs
        List of directories, each containing the sample's BED files.
        Their order defines the row order in the multi‑panel figures.
    selected_offsets_df
        DataFrame mapping read length to the selected 5' and 3' offsets.
    offset_type
        Either "5" or "3" (use 5′ or 3′ offsets).
    fasta_file
        FASTA file for sequence lengths.
    output_dir
        Parent directory in which the four sub‑directories will be created.
    args
        Parsed arguments namespace (expects .cap_percentile).
    annotation_df
        Gene annotation – *only* used to decide whether a BED record is
        CDS, hence not detailed here.
    """

    log_info("COVERAGE", "Generating coverage-profile plots (RPM + raw counts).")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_format = getattr(args, "plot_format", "svg")
    offset_site = str(getattr(args, "offset_site", "p")).lower()
    if offset_site not in {"p", "a"}:
        offset_site = "p"

    # ------------------------------------------------------------------
    # A) Parse the reference FASTA for transcript lengths
    # ------------------------------------------------------------------
    fasta_dict: Dict[str, SeqIO.SeqRecord] = SeqIO.to_dict(SeqIO.parse(str(fasta_file), "fasta"))
    sequence_display_map = build_sequence_display_map(annotation_df, fasta_dict.keys())

    # ------------------------------------------------------------------
    # B) Prepare containers keyed by sample → transcript → numpy array
    # ------------------------------------------------------------------
    read_cov_map: Dict[str, Dict[str, np.ndarray]] = {}
    psite_cov_map: Dict[str, Dict[str, np.ndarray]] = {}

    sample_list: List[Tuple[str, Path | None]] = [
        (Path(s).name, (None if filtered_bed_df is not None else Path(s))) for s in sample_dirs
    ]

    # ------------------------------------------------------------------
    # C) Gather per‑position coverage for every sample
    # ------------------------------------------------------------------
    for sample, sdir in iter_with_progress(
        sample_list,
        component="COVERAGE",
        noun="sample",
        labeler=lambda item: item[0],
    ):
        read_cov_map[sample] = {}
        psite_cov_map[sample] = {}

        if filtered_bed_df is not None:
            sample_bed_df = filtered_bed_df[filtered_bed_df["sample_name"] == sample].copy()
            read_length_blocks = [
                (int(read_length), sample_bed_df[sample_bed_df["read_length"] == read_length].copy())
                for read_length in sorted(sample_bed_df["read_length"].dropna().unique())
            ]
        else:
            bed_files = [f for f in os.listdir(sdir) if f.endswith(".bed")]
            read_length_blocks = []
            for bed in bed_files:
                m = re.search(r"_(\d+)nt\.bed$", bed)
                if not m:
                    continue
                rlen = int(m.group(1))
                bed_df = pd.read_csv(sdir / bed, sep="\t", header=None)
                if bed_df.shape[1] < 3:
                    continue
                bed_df.columns = ["chrom", "start", "end", "name", "score", "strand"][: bed_df.shape[1]]
                read_length_blocks.append((rlen, bed_df))

        for rlen, bed_df in read_length_blocks:

            # Fetch offset for this read length
            row = selected_offsets_df[selected_offsets_df["Read Length"] == rlen]
            if row.empty:
                continue
            if offset_type == "5":
                offset_val = row["Most Enriched 5' Offset"].values[0]
            else:
                offset_val = row["Most Enriched 3' Offset"].values[0]
            if pd.isna(offset_val):
                continue
            offset_val = int(offset_val)

            bed_df["start"] = pd.to_numeric(bed_df["start"], errors="coerce")
            bed_df["end"] = pd.to_numeric(bed_df["end"], errors="coerce")
            bed_df.dropna(subset=["start", "end"], inplace=True)
            bed_df["start"] = bed_df["start"].astype(int)
            bed_df["end"] = bed_df["end"].astype(int)
            bed_df = bed_df[bed_df["end"] > bed_df["start"]].copy()

            for chrom, local_df in bed_df.groupby("chrom", sort=False):
                if chrom not in fasta_dict:
                    continue

                arr_read = read_cov_map[sample].setdefault(
                    chrom, np.zeros(len(fasta_dict[chrom].seq), dtype=int)
                )
                for start, end in zip(local_df["start"], local_df["end"]):
                    arr_read[int(start): int(end)] += 1

                if offset_type == "5":
                    if offset_site == "p":
                        p_positions = local_df["start"] + offset_val
                    else:
                        p_positions = local_df["start"] + offset_val - 3
                else:  # 3′ offset
                    if offset_site == "p":
                        p_positions = local_df["end"] - offset_val - 1
                    else:
                        p_positions = local_df["end"] - offset_val - 4

                arr_psite = psite_cov_map[sample].setdefault(
                    chrom, np.zeros(len(fasta_dict[chrom].seq), dtype=int)
                )
                for p_pos in p_positions.astype(int):
                    if 0 <= p_pos < len(arr_read):
                        arr_psite[p_pos] += 1

    # ------------------------------------------------------------------
    # D) Preserve a deep copy **before** converting to RPM
    # ------------------------------------------------------------------
    raw_read_cov_map = copy.deepcopy(read_cov_map)
    raw_psite_cov_map = copy.deepcopy(psite_cov_map)

    # ------------------------------------------------------------------
    # E) Convert coverage → RPM (reads per million mRNA‑aligned)
    # ------------------------------------------------------------------
    # For simplicity we assume args.total_mrna_map exists (e.g. populated in main.py)
    if not hasattr(args, "total_mrna_map"):
        raise AttributeError("args.total_mrna_map missing – cannot compute RPM!")

    for sample in read_cov_map:
        scale = 1e6 / max(args.total_mrna_map.get(sample, 1), 1)
        for chrom in read_cov_map[sample]:
            read_cov_map[sample][chrom] = read_cov_map[sample][chrom] * scale
        for chrom in psite_cov_map[sample]:
            psite_cov_map[sample][chrom] = psite_cov_map[sample][chrom] * scale

    # Optionally cap extreme outliers to improve the visual dynamic range
    cap = args.cap_percentile
    if cap and 0 < cap < 1:
        for cov_map in (read_cov_map, psite_cov_map):
            for sample in cov_map:
                for chrom in cov_map[sample]:
                    arr = cov_map[sample][chrom]
                    lim = np.quantile(arr, cap)
                    arr[arr > lim] = lim

    # ------------------------------------------------------------------
    # F) Compute per‑transcript global maxima for y‑axis limits
    # ------------------------------------------------------------------
    def build_global_max(cov_map):
        gmax: Dict[str, float] = {}
        for sample in cov_map:
            for chrom, arr in cov_map[sample].items():
                gmax[chrom] = max(gmax.get(chrom, 0), float(arr.max()))
        return gmax

    read_global_max = build_global_max(read_cov_map)
    psite_global_max = build_global_max(psite_cov_map)
    raw_read_global_max = build_global_max(raw_read_cov_map)
    raw_psite_global_max = build_global_max(raw_psite_cov_map)

    # ------------------------------------------------------------------
    # G) Prepare output directories
    # ------------------------------------------------------------------
    read_dir_rpm = output_dir / "read_coverage_rpm"
    psite_dir_rpm = output_dir / "p_site_coverage_rpm"
    read_dir_raw = output_dir / "read_coverage_raw"
    psite_dir_raw = output_dir / "p_site_coverage_raw"
    for d in (read_dir_rpm, psite_dir_rpm, read_dir_raw, psite_dir_raw):
        d.mkdir(exist_ok=True)

    # ------------------------------------------------------------------
    # H) Helper to plot a multi‑panel coverage figure
    # ------------------------------------------------------------------
    def _plot_cov(tr: str, cov_map, ymax_map, out_dir: Path, title_suffix: str, bar_color: str):
        seqlen = len(fasta_dict[tr].seq)
        ymax = ymax_map.get(tr, 1)
        ymax = max(ymax, 1)

        fig, axes = plt.subplots(
            nrows=len(sample_list), ncols=1, figsize=(14, 3 * len(sample_list)), sharex=True
        )
        if len(sample_list) == 1:
            axes = [axes]

        for row, (sample, _) in enumerate(sample_list):
            ax = axes[row]
            arr = cov_map[sample].get(tr)
            if arr is None:
                ax.text(0.5, 0.5, f"No data for {sample}", ha="center", va="center", transform=ax.transAxes)
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            X = np.arange(seqlen)
            ax.bar(X, arr, width=3.0, color=bar_color)
            ax.set_ylim(0, ymax * 1.1)
            ax.set_ylabel(sample, fontsize=12, fontweight="bold")
            if row == len(sample_list) - 1:
                ax.set_xlabel("Nucleotide Position", fontsize=14, fontweight="bold")
            ax.set_xticks(np.arange(0, seqlen, 100))
            ax.tick_params(axis="x", rotation=45)

        fig.suptitle(
            f"{sequence_display_map.get(tr, tr)} - {title_suffix}",
            fontsize=16,
            fontweight="bold",
        )
        fig.tight_layout()
        out_path = out_dir / f"{tr}_{title_suffix.replace(' ', '_').lower()}.{plot_format}"
        fig.savefig(out_path)
        plt.close(fig)
        return out_path

    # ------------------------------------------------------------------
    # I) Iterate over all transcripts and create four figures each
    # ------------------------------------------------------------------
    transcripts = sorted({chrom for cov in read_cov_map.values() for chrom in cov})

    for tr in transcripts:
        # RPM
        out_r_rpm = _plot_cov(tr, read_cov_map, read_global_max, read_dir_rpm, "Read Coverage (RPM)", "steelblue")
        out_p_rpm = _plot_cov(tr, psite_cov_map, psite_global_max, psite_dir_rpm, "P-site Density (RPM)", "forestgreen")
        # RAW
        out_r_raw = _plot_cov(tr, raw_read_cov_map, raw_read_global_max, read_dir_raw, "Read Coverage (Raw Counts)", "steelblue")
        out_p_raw = _plot_cov(tr, raw_psite_cov_map, raw_psite_global_max, psite_dir_raw, "P-site Density (Raw Counts)", "forestgreen")
        log_info(
            "COVERAGE",
            f"{tr} -> rpm:({out_r_rpm.name},{out_p_rpm.name}) | raw:({out_r_raw.name},{out_p_raw.name})",
        )

    log_info("COVERAGE", "All coverage-profile plots generated.")
