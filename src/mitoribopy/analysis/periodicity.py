"""Periodicity + frame QC for mt-Ribo-seq.

Mt-Ribo-seq libraries are notoriously noisy on the periodicity axis:
mt-mRNAs are short (mean CDS ~700 nt for human), the 3-nt phasing is
muddier than nuclear cytosolic Ribo-seq, and visual inspection of a
single transcript is rarely conclusive. This module produces three
audit artefacts so reviewers can defend (or reject) a run's offset
choices:

1. **Frame summary.** For every sample, the fraction of P-site reads
   landing in frame 0, 1, 2 across every CDS position. A healthy
   library has frame-0 dominance; values are written verbatim instead
   of thresholded so reviewers see the actual numbers.

2. **Start-aligned metagene.** P-site density at codons 0..N from each
   transcript's start codon, summed across transcripts, per sample.
   Default window is 100 codons (300 nt). The community-standard
   "show me 3-nt periodicity" plot.

3. **Stop-aligned metagene.** Mirror of (2) aligned at the stop codon.
   For mt-mRNAs the start signal is sometimes muddy (unusual mt
   initiation context); the stop-aligned profile can be cleaner. Both
   are emitted so the user can pick whichever is more diagnostic.

Strand sanity is folded into the same QC bundle: it counts how many
reads were on the minus strand of the (already-oriented)
mt-transcriptome FASTA. A non-zero number suggests the wrong
``--library-strandedness`` setting upstream.

Public entry point: :func:`run_periodicity_qc`.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


__all__ = [
    "MetageneProfile",
    "StrandSanity",
    "compute_p_site_positions",
    "compute_metagene",
    "compute_strand_sanity",
    "run_periodicity_qc",
]


# Per-read-length inclusion thresholds. Defaults follow the spirit of
# riboWaltz's periodicity-threshold mode (Lauria et al. 2018) but are
# overridable from the CLI / config — different organisms and nucleases
# produce different mtRPF length distributions and asking the user to
# justify any deviation in the manifest is preferable to a silent
# drop.
_DEFAULT_MIN_CDS_READS_PER_LENGTH = 200
_DEFAULT_MIN_FRAME0_FRACTION = 0.50
_DEFAULT_MIN_FRAME0_DOMINANCE = 0.10
_DEFAULT_MAX_FRAME_ENTROPY = 1.45  # log2(3) ~= 1.585; 1.45 ~ "leaning periodic"


# Default windows. 100 codons (300 nt) is the figure most published
# mt-Ribo-seq papers show; long enough to see periodicity build up,
# short enough not to be smeared by transcript-to-transcript length
# heterogeneity (mt-mRNAs span ~200-1800 nt).
DEFAULT_WINDOW_NT = 300
DEFAULT_WINDOW_CODONS = DEFAULT_WINDOW_NT // 3


# ---------------------------------------------------------------------------
# Lightweight result containers (frozen for safe-by-default sharing)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MetageneProfile:
    """Per-sample, per-relative-position P-site density."""

    sample: str
    anchor: str  # 'start' | 'stop'
    positions: np.ndarray  # 1D int array (relative-to-anchor nt offset)
    density: np.ndarray  # 1D float array, same length as positions


@dataclass(frozen=True)
class StrandSanity:
    """Per-sample minus-strand fraction.

    The mt-transcriptome FASTA is already oriented in the sense
    direction of each transcript, so legitimate Ribo-seq reads must
    map to the plus strand. A non-zero ``minus_fraction`` is a strong
    signal that the wrong ``--library-strandedness`` was supplied
    upstream.
    """

    sample: str
    n_total: int
    n_minus: int
    minus_fraction: float


# ---------------------------------------------------------------------------
# Core computation — pure functions on DataFrames + offset tables
# ---------------------------------------------------------------------------


def _resolve_offsets_for_sample(
    selected_offsets_by_sample: dict[str, pd.DataFrame] | None,
    selected_offsets_combined: pd.DataFrame | None,
    sample: str,
) -> pd.DataFrame | None:
    """Return the per-length 5' offset table for ``sample``.

    Falls back to the combined table when no per-sample table is
    available (e.g. ``--offset_mode combined``). Returns ``None`` when
    neither is set.
    """
    if selected_offsets_by_sample and sample in selected_offsets_by_sample:
        return selected_offsets_by_sample[sample]
    return selected_offsets_combined


def compute_p_site_positions(
    bed_df: pd.DataFrame,
    *,
    sample: str,
    selected_offsets_by_sample: dict[str, pd.DataFrame] | None,
    selected_offsets_combined: pd.DataFrame | None,
    offset_type: str = "5",
    offset_site: str = "p",
) -> pd.DataFrame:
    """Annotate every BED row with its P-site coordinate.

    Returns a copy of the input subset to ``sample`` with an added
    ``P_site`` column (transcript-local 0-indexed nt position). Reads
    whose read length has no selected offset are dropped. Mirrors the
    arithmetic in :mod:`mitoribopy.analysis.translation_profile_analysis`
    so the QC and the plotted profiles agree by construction.
    """
    sub = bed_df[bed_df["sample_name"] == sample].copy()
    if sub.empty:
        return sub.assign(P_site=pd.Series(dtype=int))

    offsets = _resolve_offsets_for_sample(
        selected_offsets_by_sample, selected_offsets_combined, sample
    )
    if offsets is None or offsets.empty:
        return sub.iloc[0:0].assign(P_site=pd.Series(dtype=int))

    five_col = "Most Enriched 5' Offset"
    three_col = "Most Enriched 3' Offset"
    chosen_col = five_col if str(offset_type) == "5" else three_col
    if chosen_col not in offsets.columns:
        return sub.iloc[0:0].assign(P_site=pd.Series(dtype=int))

    offset_map = (
        offsets.set_index("Read Length")[chosen_col]
        .dropna()
        .astype(int)
        .to_dict()
    )
    sub = sub[sub["read_length"].isin(offset_map.keys())].copy()
    if sub.empty:
        return sub.assign(P_site=pd.Series(dtype=int))

    sub["start"] = pd.to_numeric(sub["start"], errors="coerce").astype("Int64")
    sub["end"] = pd.to_numeric(sub["end"], errors="coerce").astype("Int64")
    sub = sub.dropna(subset=["start", "end"])
    sub["start"] = sub["start"].astype(int)
    sub["end"] = sub["end"].astype(int)

    offset_series = sub["read_length"].map(offset_map).astype(int)
    site_normalised = str(offset_site).lower()
    if str(offset_type) == "5":
        # Forward-end placement; convert to canonical P-site space.
        if site_normalised == "p":
            sub["P_site"] = sub["start"] + offset_series
        else:
            sub["P_site"] = sub["start"] + offset_series - 3
    else:
        # 3'-end placement; BED is end-exclusive so use (end - 1).
        if site_normalised == "p":
            sub["P_site"] = sub["end"] - offset_series - 1
        else:
            sub["P_site"] = sub["end"] - offset_series - 4

    return sub


def compute_metagene(
    bed_with_psite: pd.DataFrame,
    annotation_df: pd.DataFrame,
    *,
    sample: str,
    anchor: str,
    window_nt: int = DEFAULT_WINDOW_NT,
) -> MetageneProfile:
    """Build a start- or stop-anchored P-site density metagene.

    For ``anchor='start'`` the relative-position axis is ``[0,
    window_nt)`` (P-site nt offset from the start codon). For
    ``anchor='stop'`` the axis is ``(-window_nt, 0]`` (P-site nt offset
    from the first nt of the stop codon).

    Density is summed across transcripts with no per-transcript
    normalisation (a longer / more-expressed transcript contributes
    more signal). This matches the convention used in published mt-
    Ribo-seq metagene plots and keeps the output interpretable as raw
    P-site read counts.
    """
    if anchor not in {"start", "stop"}:
        raise ValueError(f"anchor must be 'start' or 'stop', got {anchor!r}")

    if anchor == "start":
        positions = np.arange(0, window_nt, dtype=int)
    else:
        positions = np.arange(-window_nt + 1, 1, dtype=int)
    density = np.zeros(window_nt, dtype=float)

    if bed_with_psite.empty or "P_site" not in bed_with_psite.columns:
        return MetageneProfile(
            sample=sample, anchor=anchor, positions=positions, density=density,
        )

    ann = annotation_df.set_index("transcript")[
        ["start_codon", "stop_codon", "l_tr"]
    ]
    chrom_to_tx = annotation_df.set_index("sequence_name")["transcript"].to_dict()

    df = bed_with_psite.copy()
    df["transcript"] = df["chrom"].map(chrom_to_tx).fillna(df["chrom"])
    df = df[df["transcript"].isin(ann.index)]
    if df.empty:
        return MetageneProfile(
            sample=sample, anchor=anchor, positions=positions, density=density,
        )

    if anchor == "start":
        anchor_pos = df["transcript"].map(ann["start_codon"]).astype(int)
        rel = df["P_site"].astype(int) - anchor_pos
        mask = (rel >= 0) & (rel < window_nt)
        rel = rel[mask]
    else:
        anchor_pos = df["transcript"].map(ann["stop_codon"]).astype(int)
        rel = df["P_site"].astype(int) - anchor_pos
        # rel range we want: (-window_nt, 0]  (axis index = rel + window_nt - 1)
        mask = (rel > -window_nt) & (rel <= 0)
        rel = rel[mask]

    if rel.empty:
        return MetageneProfile(
            sample=sample, anchor=anchor, positions=positions, density=density,
        )

    if anchor == "start":
        idx = rel.to_numpy()
    else:
        idx = (rel + window_nt - 1).to_numpy()
    np.add.at(density, idx, 1.0)
    return MetageneProfile(
        sample=sample, anchor=anchor, positions=positions, density=density,
    )


def compute_strand_sanity(
    bed_df: pd.DataFrame, *, sample: str,
) -> StrandSanity:
    """Fraction of reads on the minus strand for one sample.

    Returns zeros when the BED has no ``strand`` column; the field is
    optional in standard BED3, and many upstream tools drop it.
    """
    sub = bed_df[bed_df["sample_name"] == sample]
    n_total = len(sub)
    if n_total == 0 or "strand" not in sub.columns:
        return StrandSanity(sample=sample, n_total=n_total,
                            n_minus=0, minus_fraction=0.0)
    minus_mask = sub["strand"].astype(str).str.strip() == "-"
    n_minus = int(minus_mask.sum())
    return StrandSanity(
        sample=sample, n_total=n_total, n_minus=n_minus,
        minus_fraction=n_minus / n_total if n_total else 0.0,
    )


# ---------------------------------------------------------------------------
# Orchestrator + writers
# ---------------------------------------------------------------------------


def _write_metagene_tsv(
    profiles: Iterable[MetageneProfile], path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as h:
        h.write("sample\tanchor\tposition\tdensity\n")
        for p in profiles:
            for pos, dens in zip(p.positions, p.density):
                h.write(f"{p.sample}\t{p.anchor}\t{int(pos)}\t{dens:.6g}\n")


def _write_strand_tsv(
    rows: Iterable[StrandSanity], path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as h:
        h.write("sample\tn_total\tn_minus\tminus_fraction\n")
        for r in rows:
            h.write(
                f"{r.sample}\t{r.n_total}\t{r.n_minus}\t{r.minus_fraction:.6g}\n"
            )


def _plot_metagene_panels(
    start_profiles: list[MetageneProfile],
    stop_profiles: list[MetageneProfile],
    *,
    out_path: Path,
) -> None:
    """Two-panel plot: start-anchored on top, stop-anchored on bottom.

    Bars are coloured by frame (0 / 1 / 2) so 3-nt periodicity is
    visible at a glance. Output is written as both PNG and SVG when the
    extension is one of those (the SVG sidecar pattern used elsewhere).
    """
    import matplotlib.pyplot as plt

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    samples = sorted({p.sample for p in start_profiles + stop_profiles})
    if not samples:
        return

    # Share Y per anchor column so panels are visually comparable
    # across samples. sharex stays False because start/stop windows
    # span different position ranges.
    fig, axes = plt.subplots(
        len(samples), 2, figsize=(12, 2.6 * len(samples)),
        sharex=False, sharey="col", squeeze=False,
    )
    frame_colors = {0: "#0072B2", 1: "#E69F00", 2: "#009E73"}
    start_by_sample = {p.sample: p for p in start_profiles}
    stop_by_sample = {p.sample: p for p in stop_profiles}

    legend_handles: list = []
    for row_i, sample in enumerate(samples):
        for col_i, (anchor, profiles_by_sample) in enumerate(
            (("start", start_by_sample), ("stop", stop_by_sample))
        ):
            ax = axes[row_i][col_i]
            profile = profiles_by_sample.get(sample)
            if profile is None:
                ax.set_visible(False)
                continue
            for f in (0, 1, 2):
                # Frame-by-frame colouring relative to the anchor
                # codon's first nt (anchor itself is in frame 0).
                # (pos % 3) handles negative positions on the stop
                # window correctly.
                frame_mask = (profile.positions % 3) == f
                bars = ax.bar(
                    profile.positions[frame_mask],
                    profile.density[frame_mask],
                    width=1.0,
                    color=frame_colors[f],
                    label=f"frame {f}",
                )
                if row_i == 0 and col_i == 0:
                    legend_handles.append(bars)
            ax.set_title(
                f"{sample} — {'start-aligned' if anchor == 'start' else 'stop-aligned'} P-site"
            )
            ax.set_xlabel(
                "nt from start codon" if anchor == "start"
                else "nt from stop codon"
            )
            if col_i == 0:
                ax.set_ylabel("P-site reads")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

    # Figure-level legend below all panels so it never overlaps data
    # in any subplot.
    if legend_handles:
        fig.legend(
            handles=legend_handles,
            labels=[h.get_label() for h in legend_handles],
            loc="lower center",
            bbox_to_anchor=(0.5, -0.02),
            ncol=3,
            frameon=False,
            fontsize=10,
        )

    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    if out_path.suffix.lower() == ".png":
        fig.savefig(out_path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def _plot_metagene_single_anchor(
    profiles: list[MetageneProfile],
    *,
    anchor: str,
    site_letter: str,
    out_path: Path,
) -> None:
    """Single-anchor metagene plot (one of start / stop only).

    Spec asks for ``metagene_start_p_site.svg`` and
    ``metagene_stop_p_site.svg`` as SEPARATE files (not the combined
    panel produced by :func:`_plot_metagene_panels`). One row per
    sample, bars frame-coloured.
    """
    import matplotlib.pyplot as plt

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    samples = sorted({p.sample for p in profiles})
    if not samples:
        return

    fig, axes = plt.subplots(
        len(samples), 1, figsize=(10, 2.6 * len(samples)),
        sharex=True, sharey=True, squeeze=False,
    )
    frame_colors = {0: "#0072B2", 1: "#E69F00", 2: "#009E73"}
    by_sample = {p.sample: p for p in profiles}
    legend_handles: list = []
    site_label = "P-site" if site_letter.lower() == "p" else "A-site"

    for row_i, sample in enumerate(samples):
        ax = axes[row_i][0]
        profile = by_sample.get(sample)
        if profile is None:
            ax.set_visible(False)
            continue
        for f in (0, 1, 2):
            frame_mask = (profile.positions % 3) == f
            bars = ax.bar(
                profile.positions[frame_mask],
                profile.density[frame_mask],
                width=1.0,
                color=frame_colors[f],
                label=f"frame {f}",
            )
            if row_i == 0:
                legend_handles.append(bars)
        ax.set_title(
            f"{sample} — {'start-aligned' if anchor == 'start' else 'stop-aligned'} {site_label}",
            fontsize=10,
        )
        ax.set_xlabel(
            "nt from start codon" if anchor == "start" else "nt from stop codon"
        )
        ax.set_ylabel(f"{site_label} reads")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    if legend_handles:
        fig.legend(
            handles=legend_handles,
            labels=[h.get_label() for h in legend_handles],
            loc="lower center",
            bbox_to_anchor=(0.5, -0.02),
            ncol=3,
            frameon=False,
            fontsize=10,
        )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    if out_path.suffix.lower() == ".png":
        fig.savefig(out_path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def run_periodicity_qc(
    *,
    bed_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    samples: Iterable[str],
    selected_offsets_by_sample: dict[str, pd.DataFrame] | None,
    selected_offsets_combined: pd.DataFrame | None,
    offset_type: str,
    offset_site: str,
    output_dir: Path,
    window_nt: int = DEFAULT_WINDOW_NT,
    plot: bool = True,
    fourier_window_nt: int | None = None,
    fourier_render_plots: bool = True,
    metagene_nt: int | None = None,
    # Accepted for backwards-compatible call sites; ignored under the
    # current Fourier-only QC contract.
    min_cds_reads_per_length: int = _DEFAULT_MIN_CDS_READS_PER_LENGTH,
    min_frame0_fraction: float = _DEFAULT_MIN_FRAME0_FRACTION,
    min_frame0_dominance: float = _DEFAULT_MIN_FRAME0_DOMINANCE,
    max_frame_entropy: float = _DEFAULT_MAX_FRAME_ENTROPY,
    exclude_start_codons: int = 6,
    exclude_stop_codons: int = 3,
    qc_thresholds: dict[str, float] | None = None,
) -> dict:
    """Compute start/stop metagenes + strand sanity + metagene Fourier QC.

    Writes:

    * ``periodicity_start.tsv`` / ``periodicity_stop.tsv`` — per-sample
      metagene profiles around start / stop codons (legacy filenames).
    * ``metagene_start.tsv`` / ``metagene_stop.tsv`` — same data, spec
      filenames.
    * ``strand_sanity.tsv`` — per-sample read-strand check.
    * ``periodicity_metagene.png/.svg`` — combined two-anchor metagene
      figure (when ``plot=True``).
    * ``metagene_{start,stop}_<site>_site.svg`` — single-anchor variants.

    Plus the metagene Fourier bundle written by
    :func:`mitoribopy.analysis.periodicity_qc.run_periodicity_qc_bundle`:
    ``fourier_spectrum_combined.tsv``, ``fourier_period3_score_combined.tsv``,
    ``fourier_spectrum/<sample>/*.{png,svg}``, ``periodicity.metadata.json``.

    The frame-fraction QC outputs (``frame_counts_*.tsv``, ``qc_summary.tsv``,
    ``gene_periodicity.tsv``, ``frame_fraction_heatmap.svg``,
    ``read_length_periodicity_barplot.svg``, ``gene_phase_score_dotplot.svg``)
    were removed in v0.8.0 — the Fourier ``spectral_ratio_3nt`` +
    ``snr_call`` columns now carry the headline QC story.
    """
    output_dir = Path(output_dir)
    effective_window_nt = int(metagene_nt) if metagene_nt is not None else int(window_nt)
    start_profiles: list[MetageneProfile] = []
    stop_profiles: list[MetageneProfile] = []
    strand_rows: list[StrandSanity] = []

    for sample in samples:
        psite = compute_p_site_positions(
            bed_df,
            sample=sample,
            selected_offsets_by_sample=selected_offsets_by_sample,
            selected_offsets_combined=selected_offsets_combined,
            offset_type=offset_type,
            offset_site=offset_site,
        )
        start_profiles.append(compute_metagene(
            psite, annotation_df,
            sample=sample, anchor="start", window_nt=effective_window_nt,
        ))
        stop_profiles.append(compute_metagene(
            psite, annotation_df,
            sample=sample, anchor="stop", window_nt=effective_window_nt,
        ))
        strand_rows.append(compute_strand_sanity(bed_df, sample=sample))

    _write_metagene_tsv(start_profiles, output_dir / "periodicity_start.tsv")
    _write_metagene_tsv(stop_profiles, output_dir / "periodicity_stop.tsv")
    _write_metagene_tsv(start_profiles, output_dir / "metagene_start.tsv")
    _write_metagene_tsv(stop_profiles, output_dir / "metagene_stop.tsv")
    _write_strand_tsv(strand_rows, output_dir / "strand_sanity.tsv")

    # Stack per-sample P-site coordinates so the Fourier bundle has all
    # the data in one BED-shaped DataFrame.
    psite_frames: list[pd.DataFrame] = []
    for sample in samples:
        psite = compute_p_site_positions(
            bed_df,
            sample=sample,
            selected_offsets_by_sample=selected_offsets_by_sample,
            selected_offsets_combined=selected_offsets_combined,
            offset_type=offset_type,
            offset_site=offset_site,
        )
        if not psite.empty:
            psite_frames.append(psite)
    bed_with_psite_concat = (
        pd.concat(psite_frames, ignore_index=True) if psite_frames else None
    )

    from .periodicity_qc import run_periodicity_qc_bundle as _qc_bundle

    fourier_kwargs = {}
    if fourier_window_nt is not None:
        fourier_kwargs["window_nt"] = int(fourier_window_nt)

    _qc_bundle(
        bed_with_psite=bed_with_psite_concat,
        annotation_df=annotation_df,
        samples=samples,
        output_dir=output_dir,
        site_type=str(offset_site).lower() or "p",  # type: ignore[arg-type]
        render_plots=fourier_render_plots,
        **fourier_kwargs,
    )

    if plot:
        _plot_metagene_panels(
            start_profiles, stop_profiles,
            out_path=output_dir / "periodicity_metagene.png",
        )
        site_letter = str(offset_site).lower() or "p"
        _plot_metagene_single_anchor(
            start_profiles,
            anchor="start",
            site_letter=site_letter,
            out_path=output_dir / f"metagene_start_{site_letter}_site.svg",
        )
        _plot_metagene_single_anchor(
            stop_profiles,
            anchor="stop",
            site_letter=site_letter,
            out_path=output_dir / f"metagene_stop_{site_letter}_site.svg",
        )

    return {
        "periodicity_start": start_profiles,
        "periodicity_stop": stop_profiles,
        "strand_sanity": strand_rows,
    }
