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

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


__all__ = [
    "FrameSummary",
    "MetageneProfile",
    "StrandSanity",
    "compute_p_site_positions",
    "compute_frame_summary",
    "compute_metagene",
    "compute_strand_sanity",
    "run_periodicity_qc",
]


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
class FrameSummary:
    """Per-sample frame fractions over every assigned P-site."""

    sample: str
    n_reads: int
    frame_0: float
    frame_1: float
    frame_2: float
    dominance: float  # frame_0 - max(frame_1, frame_2). >0.10 is healthy.


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


def compute_frame_summary(
    bed_with_psite: pd.DataFrame,
    annotation_df: pd.DataFrame,
    *,
    sample: str,
) -> FrameSummary:
    """Frame-0/1/2 fractions across every CDS-resident P-site.

    A read's frame is ``(P_site - start_codon) % 3`` evaluated only for
    reads whose P-site lies inside the annotated CDS (start_codon
    inclusive, stop_codon exclusive). Non-CDS reads are excluded so
    UTR-bound contamination cannot inflate one frame artificially.
    """
    if bed_with_psite.empty or "P_site" not in bed_with_psite.columns:
        return FrameSummary(sample=sample, n_reads=0,
                            frame_0=0.0, frame_1=0.0, frame_2=0.0,
                            dominance=0.0)

    ann = annotation_df.set_index("transcript")[["start_codon", "l_cds"]]
    chrom_to_tx = annotation_df.set_index("sequence_name")["transcript"].to_dict()

    df = bed_with_psite.copy()
    df["transcript"] = df["chrom"].map(chrom_to_tx).fillna(df["chrom"])
    df = df[df["transcript"].isin(ann.index)]
    if df.empty:
        return FrameSummary(sample=sample, n_reads=0,
                            frame_0=0.0, frame_1=0.0, frame_2=0.0,
                            dominance=0.0)

    starts = df["transcript"].map(ann["start_codon"]).astype(int)
    cds_lens = df["transcript"].map(ann["l_cds"]).astype(int)
    in_cds = (df["P_site"] >= starts) & (df["P_site"] < starts + cds_lens)
    df = df[in_cds]
    n_reads = len(df)
    if n_reads == 0:
        return FrameSummary(sample=sample, n_reads=0,
                            frame_0=0.0, frame_1=0.0, frame_2=0.0,
                            dominance=0.0)
    frames = ((df["P_site"] - starts.loc[df.index]) % 3).astype(int)
    counts = frames.value_counts().to_dict()
    f0 = counts.get(0, 0) / n_reads
    f1 = counts.get(1, 0) / n_reads
    f2 = counts.get(2, 0) / n_reads
    dominance = f0 - max(f1, f2)
    return FrameSummary(
        sample=sample, n_reads=n_reads,
        frame_0=f0, frame_1=f1, frame_2=f2, dominance=dominance,
    )


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


def _write_frame_summary_tsv(
    rows: Iterable[FrameSummary], path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as h:
        h.write("sample\tn_reads\tframe_0\tframe_1\tframe_2\tdominance\n")
        for r in rows:
            h.write(
                f"{r.sample}\t{r.n_reads}\t{r.frame_0:.6g}\t{r.frame_1:.6g}\t"
                f"{r.frame_2:.6g}\t{r.dominance:.6g}\n"
            )


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

    fig, axes = plt.subplots(
        len(samples), 2, figsize=(12, 2.6 * len(samples)),
        sharex=False, squeeze=False,
    )
    frame_colors = {0: "#0072B2", 1: "#E69F00", 2: "#009E73"}
    start_by_sample = {p.sample: p for p in start_profiles}
    stop_by_sample = {p.sample: p for p in stop_profiles}

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
                if anchor == "start":
                    frame_mask = (profile.positions % 3) == f
                else:
                    # Stop codon's first nt is also frame 0; positions
                    # run -window+1 .. 0, so (pos % 3) handles negatives
                    # as we want.
                    frame_mask = (profile.positions % 3) == f
                ax.bar(
                    profile.positions[frame_mask],
                    profile.density[frame_mask],
                    width=1.0,
                    color=frame_colors[f],
                    label=f"frame {f}" if (row_i == 0 and col_i == 0) else None,
                )
            ax.set_title(
                f"{sample} — {'start-aligned' if anchor == 'start' else 'stop-aligned'} P-site"
            )
            ax.set_xlabel(
                "nt from start codon" if anchor == "start"
                else "nt from stop codon"
            )
            ax.set_ylabel("P-site reads")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
    if axes[0][0].get_visible():
        axes[0][0].legend(loc="upper right", frameon=False)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    if out_path.suffix.lower() == ".png":
        fig.savefig(out_path.with_suffix(".svg"))
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
) -> dict:
    """Compute frame summary + start/stop metagenes + strand sanity per sample.

    Writes ``frame_summary.tsv``, ``periodicity_start.tsv``,
    ``periodicity_stop.tsv``, ``strand_sanity.tsv``, and (when
    ``plot=True``) ``periodicity_metagene.png/.svg`` under
    ``output_dir``. Returns a dict with the in-memory result objects so
    callers can inline assertions in tests without re-reading from disk.
    """
    output_dir = Path(output_dir)
    frame_rows: list[FrameSummary] = []
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
        frame_rows.append(compute_frame_summary(
            psite, annotation_df, sample=sample,
        ))
        start_profiles.append(compute_metagene(
            psite, annotation_df,
            sample=sample, anchor="start", window_nt=window_nt,
        ))
        stop_profiles.append(compute_metagene(
            psite, annotation_df,
            sample=sample, anchor="stop", window_nt=window_nt,
        ))
        strand_rows.append(compute_strand_sanity(bed_df, sample=sample))

    _write_frame_summary_tsv(frame_rows, output_dir / "frame_summary.tsv")
    _write_metagene_tsv(start_profiles, output_dir / "periodicity_start.tsv")
    _write_metagene_tsv(stop_profiles, output_dir / "periodicity_stop.tsv")
    _write_strand_tsv(strand_rows, output_dir / "strand_sanity.tsv")
    if plot:
        _plot_metagene_panels(
            start_profiles, stop_profiles,
            out_path=output_dir / "periodicity_metagene.png",
        )

    return {
        "frame_summary": frame_rows,
        "periodicity_start": start_profiles,
        "periodicity_stop": stop_profiles,
        "strand_sanity": strand_rows,
    }
