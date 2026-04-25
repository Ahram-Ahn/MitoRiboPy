"""
Generates coverage-profile plots (one figure per transcript) for each sample.

Two sets of figures are produced for every transcript in every run:

1. **RPM‑scaled** coverage (read‑depth normalised to 1 million mRNA‑aligned reads)
2. **Raw‑count** coverage (no scaling at all)

Each set contains a *read‑coverage* panel and a *P‑site* panel.
All samples share the same y‑axis per transcript so that visual
comparisons are not confounded by different scale bars.

Directory layout (v0.4.x, unified per-site)::

    <coverage_profile_plots>/
        read_coverage_rpm/         *.<plot_format>      (site-independent;
        read_coverage_raw/                               written ONCE even
        read_coverage_rpm_codon/                         when both sites are
        read_coverage_raw_codon/                         analysed)
        <site>/                                          # 'p' or 'a'
            site_density_rpm/      *.<plot_format>
            site_density_raw/
            site_density_rpm_codon/
            site_density_raw_codon/
            site_density_rpm_frame/                      (CDS only,
            site_density_raw_frame/                       frame-coloured)

Usage
-----
Call the exported helper from the package pipeline::

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
from ..console import iter_with_progress, log_info, log_warning
from ..data import build_sequence_display_map, resolve_sequence_name
from .style import apply_publication_style

apply_publication_style()


def _build_cds_lookup(
    annotation_df: pd.DataFrame,
    sequence_names: set[str],
) -> dict[str, tuple[int, int]]:
    """Map resolved transcript ids to CDS start/end (end-exclusive) bounds."""
    normalized = annotation_df.copy()
    if "start_codon" not in normalized.columns and {"l_utr5"}.issubset(normalized.columns):
        normalized["start_codon"] = normalized["l_utr5"]
    if "stop_codon" not in normalized.columns and {"l_tr", "l_utr3"}.issubset(normalized.columns):
        normalized["stop_codon"] = normalized["l_tr"] - normalized["l_utr3"] - 3

    lookup: dict[str, tuple[int, int]] = {}
    for _, row in normalized.iterrows():
        resolved_name = resolve_sequence_name(row, sequence_names)
        if resolved_name is None:
            continue
        cds_start = int(row["start_codon"])
        cds_end = int(row["stop_codon"]) + 3
        if cds_end <= cds_start:
            continue
        lookup[resolved_name] = (cds_start, cds_end)
    return lookup


def _build_codon_binned_map(
    cov_map: Dict[str, Dict[str, np.ndarray]],
    cds_lookup: dict[str, tuple[int, int]],
) -> Dict[str, Dict[str, np.ndarray]]:
    """Aggregate nucleotide-resolution coverage into CDS codon bins."""
    codon_map: Dict[str, Dict[str, np.ndarray]] = {}
    for sample, transcript_map in cov_map.items():
        codon_map[sample] = {}
        for transcript, arr in transcript_map.items():
            cds_bounds = cds_lookup.get(transcript)
            if cds_bounds is None:
                continue
            cds_start, cds_end = cds_bounds
            cds_end = min(cds_end, len(arr))
            cds_start = max(cds_start, 0)
            cds_len = cds_end - cds_start
            if cds_len < 3:
                continue
            usable_len = cds_len - (cds_len % 3)
            if usable_len < 3:
                continue
            cds_slice = np.asarray(arr[cds_start : cds_start + usable_len])
            codon_map[sample][transcript] = cds_slice.reshape(-1, 3).sum(axis=1)
    return codon_map


# Okabe-Ito colorblind-safe palette (frame 0 / 1 / 2).
# https://jfly.uni-koeln.de/color/
_FRAME_COLORS = ("#E69F00", "#56B4E9", "#009E73")
_FRAME_LABELS = ("Frame 0 (codon-locked)", "Frame 1", "Frame 2")


def _build_cds_nt_slice_map(
    cov_map: Dict[str, Dict[str, np.ndarray]],
    cds_lookup: dict[str, tuple[int, int]],
) -> Dict[str, Dict[str, tuple[np.ndarray, int]]]:
    """Extract the CDS nucleotide-resolution slice for each sample/transcript.

    Returns a map ``{sample: {transcript: (cds_slice_array, cds_start_nt)}}``.
    ``cds_start_nt`` is the 0-based nt coordinate that the slice starts at,
    so the caller can recover the original transcript coordinate if needed.
    UTR positions are excluded. This is the input to the frame-colored plot.
    """
    out: Dict[str, Dict[str, tuple[np.ndarray, int]]] = {}
    for sample, transcript_map in cov_map.items():
        out[sample] = {}
        for transcript, arr in transcript_map.items():
            cds_bounds = cds_lookup.get(transcript)
            if cds_bounds is None:
                continue
            cds_start, cds_end = cds_bounds
            cds_end = min(cds_end, len(arr))
            cds_start = max(cds_start, 0)
            cds_len = cds_end - cds_start
            if cds_len < 3:
                continue
            usable_len = cds_len - (cds_len % 3)
            if usable_len < 3:
                continue
            out[sample][transcript] = (
                np.asarray(arr[cds_start : cds_start + usable_len]),
                cds_start,
            )
    return out

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
    *,
    total_mrna_map: dict[str, int | float] | None = None,
    selected_offsets_by_sample: dict | None = None,
    site_override: str | None = None,
    read_coverage_dir: str | Path | None = None,
    write_read_coverage: bool = True,
):
    """Create RPM and raw-count coverage-profile plots.

    Parameters
    ----------
    sample_dirs
        List of directories, each containing the sample's BED files.
        Their order defines the row order in the multi-panel figures.
    selected_offsets_df
        DataFrame mapping read length to the selected 5' and 3' offsets.
    offset_type
        Either "5" or "3" (use 5' or 3' offsets).
    fasta_file
        FASTA file for sequence lengths.
    output_dir
        Parent directory in which the SITE-DEPENDENT subdirectories will
        be created (``site_density_rpm/`` etc.). For a per-site run this
        is typically ``<coverage_profile_plots>/p`` or ``.../a``.
    args
        Parsed arguments namespace (expects .cap_percentile).
    annotation_df
        Gene annotation - only used to decide whether a BED record is
        CDS, hence not detailed here.
    read_coverage_dir
        Optional separate parent for the SITE-INDEPENDENT
        ``read_coverage_*`` subdirectories. Defaults to ``output_dir``
        for backward-compatibility. The pipeline orchestrator passes
        ``<coverage_profile_plots>`` (one level above the per-site
        subdir) so the duplicate read_coverage figures are not written
        twice when both sites are analysed.
    write_read_coverage
        When False, the read_coverage_* directories and figures are
        skipped entirely. Used by the orchestrator on the second site
        to avoid recomputing site-independent plots.
    """

    log_info("COVERAGE", "Generating coverage-profile plots (RPM + raw counts).")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    read_coverage_dir = (
        Path(read_coverage_dir) if read_coverage_dir is not None else output_dir
    )
    if write_read_coverage:
        read_coverage_dir.mkdir(parents=True, exist_ok=True)
    plot_format = getattr(args, "plot_format", "svg")
    if site_override is not None:
        offset_site = str(site_override).lower()
    else:
        offset_site = str(getattr(args, "offset_site", "p")).lower()
    if offset_site not in {"p", "a"}:
        offset_site = "p"

    # ------------------------------------------------------------------
    # A) Parse the reference FASTA for transcript lengths
    # ------------------------------------------------------------------
    fasta_dict: Dict[str, SeqIO.SeqRecord] = SeqIO.to_dict(SeqIO.parse(str(fasta_file), "fasta"))
    sequence_display_map = build_sequence_display_map(annotation_df, fasta_dict.keys())
    cds_lookup = _build_cds_lookup(annotation_df, set(fasta_dict.keys()))

    # Warn loudly if any FASTA record has no matching annotation row -
    # coverage plots will silently skip that transcript's codon / frame
    # variants otherwise and the user might never notice.
    unmapped = [name for name in fasta_dict if name not in cds_lookup]
    if unmapped:
        log_warning(
            "COVERAGE",
            f"{len(unmapped)} FASTA record(s) have no matching annotation "
            "row; codon- and frame-resolution plots will be skipped for "
            "these. Missing: " + ", ".join(sorted(unmapped)[:10])
            + ("..." if len(unmapped) > 10 else ""),
        )

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

        sample_offsets_df = (
            selected_offsets_by_sample.get(sample, selected_offsets_df)
            if selected_offsets_by_sample
            else selected_offsets_df
        )
        if sample_offsets_df is None:
            continue
        for rlen, bed_df in read_length_blocks:

            # Fetch offset for this read length
            row = sample_offsets_df[sample_offsets_df["Read Length"] == rlen]
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
    # Prefer the explicit kwarg. Fall back to args.total_mrna_map with a
    # one-time deprecation warning, since it mutates the parsed argparse
    # namespace (Task 5c). When both are absent, fail loudly.
    resolved_totals: dict[str, int | float]
    if total_mrna_map is not None:
        resolved_totals = total_mrna_map
    elif hasattr(args, "total_mrna_map"):
        log_warning(
            "COVERAGE",
            "Reading total_mrna_map from args is deprecated; pass "
            "total_mrna_map=... as a keyword argument instead. "
            "The args fallback will be removed in v0.4.0.",
        )
        resolved_totals = args.total_mrna_map
    else:
        raise AttributeError(
            "total_mrna_map is required to compute RPM; pass it as a "
            "keyword argument to run_coverage_profile_plots."
        )

    for sample in read_cov_map:
        scale = 1e6 / max(resolved_totals.get(sample, 1), 1)
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
    # Site-INDEPENDENT (always under read_coverage_dir, which the
    # orchestrator points at <coverage_profile_plots>/ so the figures
    # are not duplicated for both sites).
    read_dir_rpm = read_coverage_dir / "read_coverage_rpm"
    read_dir_raw = read_coverage_dir / "read_coverage_raw"
    read_dir_rpm_codon = read_coverage_dir / "read_coverage_rpm_codon"
    read_dir_raw_codon = read_coverage_dir / "read_coverage_raw_codon"
    # Site-DEPENDENT (under output_dir, typically the per-site subdir
    # <coverage_profile_plots>/p or .../a). Renamed from
    # 'p_site_coverage_*' to 'site_density_*' because the parent dir
    # already names the site (P or A) and the same plot type is also
    # used in A-site mode.
    site_density_dir_rpm = output_dir / "site_density_rpm"
    site_density_dir_raw = output_dir / "site_density_raw"
    site_density_dir_rpm_codon = output_dir / "site_density_rpm_codon"
    site_density_dir_raw_codon = output_dir / "site_density_raw_codon"
    # Frame-coloured nucleotide-resolution plots scoped to the CDS only.
    # Frame-0 dominance is the canonical mt-Ribo-seq QC signature.
    site_density_dir_rpm_frame = output_dir / "site_density_rpm_frame"
    site_density_dir_raw_frame = output_dir / "site_density_raw_frame"

    site_density_dirs = (
        site_density_dir_rpm,
        site_density_dir_raw,
        site_density_dir_rpm_codon,
        site_density_dir_raw_codon,
        site_density_dir_rpm_frame,
        site_density_dir_raw_frame,
    )
    for d in site_density_dirs:
        d.mkdir(parents=True, exist_ok=True)
    if write_read_coverage:
        for d in (
            read_dir_rpm,
            read_dir_raw,
            read_dir_rpm_codon,
            read_dir_raw_codon,
        ):
            d.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # H) Helper to plot a multi‑panel coverage figure
    # ------------------------------------------------------------------
    def _plot_cov(
        tr: str,
        cov_map,
        ymax_map,
        out_dir: Path,
        title_suffix: str,
        bar_color: str,
        *,
        x_label: str = "Nucleotide Position",
        start_at_one: bool = False,
        rotation: int = 45,
    ):
        arr_ref = next(
            (cov_map[sample].get(tr) for sample, _ in sample_list if cov_map[sample].get(tr) is not None),
            None,
        )
        if arr_ref is None:
            return None
        seqlen = len(arr_ref)
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
            X = np.arange(1, seqlen + 1) if start_at_one else np.arange(seqlen)
            ax.bar(X, arr, width=0.95, color=bar_color, edgecolor="none")
            ax.set_ylim(0, ymax * 1.1)
            ax.set_ylabel(sample, fontsize=12, fontweight="bold")
            ax.grid(axis="y", alpha=0.2, linewidth=0.8)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            if row == len(sample_list) - 1:
                ax.set_xlabel(x_label, fontsize=14, fontweight="bold")
            if not start_at_one:
                ax.set_xticks(np.arange(0, seqlen, 100))
            else:
                tick_count = min(10, seqlen)
                tick_positions = np.linspace(1, seqlen, num=tick_count, dtype=int)
                ax.set_xticks(np.unique(tick_positions))
            ax.tick_params(axis="x", rotation=rotation)

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
    read_cov_map_codon = _build_codon_binned_map(read_cov_map, cds_lookup)
    psite_cov_map_codon = _build_codon_binned_map(psite_cov_map, cds_lookup)
    raw_read_cov_map_codon = _build_codon_binned_map(raw_read_cov_map, cds_lookup)
    raw_psite_cov_map_codon = _build_codon_binned_map(raw_psite_cov_map, cds_lookup)
    read_global_max_codon = build_global_max(read_cov_map_codon)
    psite_global_max_codon = build_global_max(psite_cov_map_codon)
    raw_read_global_max_codon = build_global_max(raw_read_cov_map_codon)
    raw_psite_global_max_codon = build_global_max(raw_psite_cov_map_codon)

    # CDS nt slices for the frame-colored P-site plots.
    psite_cds_slice_rpm = _build_cds_nt_slice_map(psite_cov_map, cds_lookup)
    psite_cds_slice_raw = _build_cds_nt_slice_map(raw_psite_cov_map, cds_lookup)

    def _frame_global_max(slice_map):
        gmax: Dict[str, float] = {}
        for sample in slice_map:
            for tr, (arr, _cds_start) in slice_map[sample].items():
                gmax[tr] = max(gmax.get(tr, 0), float(arr.max() if arr.size else 0))
        return gmax

    psite_frame_max_rpm = _frame_global_max(psite_cds_slice_rpm)
    psite_frame_max_raw = _frame_global_max(psite_cds_slice_raw)

    def _plot_frame_colored(
        tr: str,
        slice_map: Dict[str, Dict[str, tuple[np.ndarray, int]]],
        ymax_map: dict[str, float],
        out_dir: Path,
        title_suffix: str,
    ) -> Path | None:
        """One panel per sample; bars at nt resolution inside the CDS,
        colored by reading frame (0/1/2 relative to the CDS start).

        For ribosome-profiling, P-site density should strongly peak at
        frame 0 (codon-locked reads). The color partition makes that
        visible in a single look.
        """
        any_present = any(tr in slice_map[sample] for sample in slice_map)
        if not any_present:
            return None

        ymax = max(ymax_map.get(tr, 1), 1)
        fig, axes = plt.subplots(
            nrows=len(sample_list),
            ncols=1,
            figsize=(14, 3 * len(sample_list)),
            sharex=True,
        )
        if len(sample_list) == 1:
            axes = [axes]

        for row, (sample, _) in enumerate(sample_list):
            ax = axes[row]
            entry = slice_map.get(sample, {}).get(tr)
            if entry is None:
                ax.text(
                    0.5, 0.5, f"No data for {sample}", ha="center", va="center",
                    transform=ax.transAxes,
                )
                ax.set_xticks([]); ax.set_yticks([])
                continue
            arr, _cds_start = entry
            n = len(arr)
            x = np.arange(1, n + 1)  # codon 1 = start codon, nt 1 = first nt of start codon
            frames = np.arange(n) % 3  # 0, 1, 2, 0, 1, 2, ...
            for frame_idx in (0, 1, 2):
                mask = frames == frame_idx
                ax.bar(
                    x[mask],
                    arr[mask],
                    width=0.95,
                    color=_FRAME_COLORS[frame_idx],
                    edgecolor="none",
                    label=_FRAME_LABELS[frame_idx] if row == 0 else None,
                )
            ax.set_ylim(0, ymax * 1.1)
            ax.set_ylabel(sample, fontsize=12, fontweight="bold")
            ax.grid(axis="y", alpha=0.2, linewidth=0.8)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            if row == len(sample_list) - 1:
                ax.set_xlabel("CDS Nucleotide Position", fontsize=14, fontweight="bold")
            if row == 0:
                ax.legend(loc="upper right", framealpha=0.9, fontsize=10)

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

    primary_label = "P-site" if offset_site == "p" else "A-site"

    for tr in transcripts:
        # Site-DEPENDENT plots (always written under output_dir; that
        # parent encodes which site we are on, so the inner names use
        # the neutral 'site_density' prefix).
        out_p_rpm = _plot_cov(
            tr,
            psite_cov_map,
            psite_global_max,
            site_density_dir_rpm,
            f"{primary_label} Density (RPM)",
            "forestgreen",
        )
        out_p_raw = _plot_cov(
            tr,
            raw_psite_cov_map,
            raw_psite_global_max,
            site_density_dir_raw,
            f"{primary_label} Density (Raw Counts)",
            "forestgreen",
        )
        out_p_rpm_codon = _plot_cov(
            tr,
            psite_cov_map_codon,
            psite_global_max_codon,
            site_density_dir_rpm_codon,
            f"{primary_label} Density (RPM, Codon-binned CDS)",
            "forestgreen",
            x_label="CDS Codon Position",
            start_at_one=True,
            rotation=0,
        )
        out_p_raw_codon = _plot_cov(
            tr,
            raw_psite_cov_map_codon,
            raw_psite_global_max_codon,
            site_density_dir_raw_codon,
            f"{primary_label} Density (Raw Counts, Codon-binned CDS)",
            "forestgreen",
            x_label="CDS Codon Position",
            start_at_one=True,
            rotation=0,
        )
        out_p_rpm_frame = _plot_frame_colored(
            tr,
            psite_cds_slice_rpm,
            psite_frame_max_rpm,
            site_density_dir_rpm_frame,
            f"{primary_label} Density (RPM, CDS Frame Coloring)",
        )
        out_p_raw_frame = _plot_frame_colored(
            tr,
            psite_cds_slice_raw,
            psite_frame_max_raw,
            site_density_dir_raw_frame,
            f"{primary_label} Density (Raw Counts, CDS Frame Coloring)",
        )

        # Site-INDEPENDENT plots (skipped on the second site so the
        # orchestrator does not redo identical work).
        out_r_rpm = out_r_raw = out_r_rpm_codon = out_r_raw_codon = None
        if write_read_coverage:
            out_r_rpm = _plot_cov(
                tr, read_cov_map, read_global_max, read_dir_rpm,
                "Read Coverage (RPM)", "steelblue",
            )
            out_r_raw = _plot_cov(
                tr, raw_read_cov_map, raw_read_global_max, read_dir_raw,
                "Read Coverage (Raw Counts)", "steelblue",
            )
            out_r_rpm_codon = _plot_cov(
                tr,
                read_cov_map_codon,
                read_global_max_codon,
                read_dir_rpm_codon,
                "Read Coverage (RPM, Codon-binned CDS)",
                "steelblue",
                x_label="CDS Codon Position",
                start_at_one=True,
                rotation=0,
            )
            out_r_raw_codon = _plot_cov(
                tr,
                raw_read_cov_map_codon,
                raw_read_global_max_codon,
                read_dir_raw_codon,
                "Read Coverage (Raw Counts, Codon-binned CDS)",
                "steelblue",
                x_label="CDS Codon Position",
                start_at_one=True,
                rotation=0,
            )

        log_info(
            "COVERAGE",
            f"{tr} -> "
            + (
                f"rpm:({out_r_rpm.name if out_r_rpm else 'n/a'},{out_p_rpm.name}) | "
                f"raw:({out_r_raw.name if out_r_raw else 'n/a'},{out_p_raw.name})"
            )
            + (
                f" | codon:({out_r_rpm_codon.name if out_r_rpm_codon else 'n/a'},"
                f"{out_p_rpm_codon.name if out_p_rpm_codon else 'n/a'},"
                f"{out_r_raw_codon.name if out_r_raw_codon else 'n/a'},"
                f"{out_p_raw_codon.name if out_p_raw_codon else 'n/a'})"
            )
            + (
                f" | frame:({out_p_rpm_frame.name if out_p_rpm_frame else 'n/a'},"
                f"{out_p_raw_frame.name if out_p_raw_frame else 'n/a'})"
            ),
        )

    log_info("COVERAGE", "All coverage-profile plots generated.")
