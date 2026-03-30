# igv_style_plot.py
"""
Generates IGV‑style coverage plots (one figure per transcript) for each sample.

Two sets of figures are produced for every transcript in every run:

1. **RPM‑scaled** coverage (read‑depth normalised to 1 million mRNA‑aligned reads)
2. **Raw‑count** coverage (no scaling at all)

Each set contains a *read‑coverage* panel and a *P‑site* panel.
All samples share the same y‑axis per transcript so that visual
comparisons are not confounded by different scale bars.

The directory layout after a successful run looks like this::

    <output_dir>/
        read_coverage_rpm/   *.svg   # old behaviour (unchanged)
        p_site_coverage_rpm/ *.svg
        read_coverage_raw/   *.svg   # NEW
        p_site_coverage_raw/ *.svg   # NEW

Usage
-----
Call the exported helper from *main.py*::

    run_igv_style_plot(sample_dirs,
                       p_site_offsets_df,
                       offset_type,
                       fasta_file,
                       output_dir,
                       args,
                       annotation_df)

`sample_dirs` must already be in the desired plotting order.
`p_site_offsets_df` should have the columns:
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
from .style import apply_publication_style

apply_publication_style()

# -----------------------------------------------------------------------------
# Public entry point
# -----------------------------------------------------------------------------

def run_igv_style_plot(
    sample_dirs: List[str],
    p_site_offsets_df: pd.DataFrame,
    offset_type: str,
    fasta_file: str | Path,
    output_dir: str | Path,
    args,
    annotation_df: pd.DataFrame,
):
    """Create RPM **and** raw‑count IGV‑style plots.

    Parameters
    ----------
    sample_dirs
        List of directories, each containing the sample's BED files.
        Their order defines the row order in the multi‑panel figures.
    p_site_offsets_df
        DataFrame mapping read‑length → most‑enriched offset.
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

    print("[IGV] ▶︎ Generating IGV‑style plots (RPM + raw counts)…")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    offset_site = str(getattr(args, "offset_site", "p")).lower()
    if offset_site not in {"p", "a"}:
        offset_site = "p"

    # ------------------------------------------------------------------
    # A) Parse the reference FASTA for transcript lengths
    # ------------------------------------------------------------------
    fasta_dict: Dict[str, SeqIO.SeqRecord] = SeqIO.to_dict(SeqIO.parse(str(fasta_file), "fasta"))

    # ------------------------------------------------------------------
    # B) Prepare containers keyed by sample → transcript → numpy array
    # ------------------------------------------------------------------
    read_cov_map: Dict[str, Dict[str, np.ndarray]] = {}
    psite_cov_map: Dict[str, Dict[str, np.ndarray]] = {}

    sample_list: List[Tuple[str, Path]] = [
        (Path(s).name, Path(s)) for s in sample_dirs
    ]

    # ------------------------------------------------------------------
    # C) Gather per‑position coverage for every sample
    # ------------------------------------------------------------------
    for sample, sdir in sample_list:
        read_cov_map[sample] = {}
        psite_cov_map[sample] = {}

        bed_files = [f for f in os.listdir(sdir) if f.endswith(".bed")]
        for bed in bed_files:
            # Expect file names like  …_<READLEN>nt.bed
            m = re.search(r"_(\d+)nt\.bed$", bed)
            if not m:
                continue
            rlen = int(m.group(1))

            # Fetch offset for this read length
            row = p_site_offsets_df[p_site_offsets_df["Read Length"] == rlen]
            if row.empty:
                continue
            if offset_type == "5":
                offset_val = row["Most Enriched 5' Offset"].values[0]
            else:
                offset_val = row["Most Enriched 3' Offset"].values[0]
            if pd.isna(offset_val):
                continue
            offset_val = int(offset_val)

            # Parse BED
            with open(sdir / bed) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    chrom, start, end, *_ = line.rstrip().split("\t")
                    if chrom not in fasta_dict:
                        continue
                    start, end = int(start), int(end)
                    length = end - start

                    # Update read‑coverage
                    arr_read = read_cov_map[sample].setdefault(
                        chrom, np.zeros(len(fasta_dict[chrom].seq), dtype=int)
                    )
                    arr_read[start:end] += 1

                    # Update P‑site coverage
                    if offset_type == "5":
                        if offset_site == "p":
                            p_pos = start + offset_val
                        else:
                            # Offset points to A-site in this mode; convert to P-site.
                            p_pos = start + offset_val - 3
                    else:  # 3′ offset
                        if offset_site == "p":
                            p_pos = end - offset_val - 1
                        else:
                            # Offset points to A-site in this mode; convert to P-site.
                            p_pos = end - offset_val - 4
                    if 0 <= p_pos < len(arr_read):
                        arr_psite = psite_cov_map[sample].setdefault(
                            chrom, np.zeros(len(fasta_dict[chrom].seq), dtype=int)
                        )
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

        fig.suptitle(f"{tr} - {title_suffix}", fontsize=16, fontweight="bold")
        fig.tight_layout()
        out_path = out_dir / f"{tr}_{title_suffix.replace(' ', '_').lower()}.svg"
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
        print(f"[IGV] {tr} -> rpm:({out_r_rpm.name},{out_p_rpm.name}) | raw:({out_r_raw.name},{out_p_raw.name})")

    print("[IGV] ✔ All IGV‑style plots generated.")
