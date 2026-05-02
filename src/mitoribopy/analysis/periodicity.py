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
    "FrameSummary",
    "MetageneProfile",
    "StrandSanity",
    "compute_p_site_positions",
    "compute_frame_summary",
    "compute_frame_summary_by_length",
    "compute_metagene",
    "compute_strand_sanity",
    "select_read_lengths_by_periodicity",
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
    exclude_start_codons: int = 0,
    exclude_stop_codons: int = 0,
) -> FrameSummary:
    """Frame-0/1/2 fractions across every CDS-resident P-site.

    A read's frame is ``(P_site - start_codon) % 3`` evaluated only for
    reads whose P-site lies inside the annotated CDS (start_codon
    inclusive, stop_codon exclusive). Non-CDS reads are excluded so
    UTR-bound contamination cannot inflate one frame artificially.

    ``exclude_start_codons`` and ``exclude_stop_codons`` mask the
    initiation- and termination-proximal codons (in nt, 3 * codons)
    to keep initiation pause and stop-codon stacking out of the frame
    fraction. Defaults are 0 to preserve historical pooled numbers;
    callers (``run_periodicity_qc``, the standalone CLI) can pass the
    spec defaults of 6 / 3.
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
    edge_lo = starts + 3 * int(exclude_start_codons)
    edge_hi = starts + cds_lens - 3 * int(exclude_stop_codons)
    in_cds = (df["P_site"] >= edge_lo) & (df["P_site"] < edge_hi)
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


def compute_frame_summary_by_length(
    bed_with_psite: pd.DataFrame,
    annotation_df: pd.DataFrame,
    *,
    sample: str,
    min_cds_reads_per_length: int = _DEFAULT_MIN_CDS_READS_PER_LENGTH,
    min_frame0_fraction: float = _DEFAULT_MIN_FRAME0_FRACTION,
    min_frame0_dominance: float = _DEFAULT_MIN_FRAME0_DOMINANCE,
    max_frame_entropy: float = _DEFAULT_MAX_FRAME_ENTROPY,
    exclude_start_codons: int = 0,
    exclude_stop_codons: int = 0,
) -> pd.DataFrame:
    """Per-read-length frame fractions, dominance, entropy and inclusion call.

    Stratifies frame phasing by ``read_length`` instead of pooling. This
    is the QC layer the publication review (and riboWaltz) demand: a
    pooled frame-0 number can look healthy when one bad length class is
    contaminating an otherwise clean library, and a length that looks
    pooled-bad can be salvaged by reporting its own offset.

    Returns a dataframe with one row per length actually present in
    ``bed_with_psite``. Length classes with no CDS reads are still
    reported (``include_for_downstream=False``, ``exclusion_reason=
    "no_cds_reads"``) so the manifest distinguishes "absent" from "tested
    and rejected".
    """
    cols = [
        "sample_id",
        "read_length",
        "n_reads_total",
        "n_reads_cds",
        "frame0_fraction",
        "frame1_fraction",
        "frame2_fraction",
        "dominant_frame",
        "frame0_dominance",
        "periodicity_score",
        "frame_entropy",
        "include_for_downstream",
        "exclusion_reason",
    ]
    if (
        bed_with_psite.empty
        or "P_site" not in bed_with_psite.columns
        or "read_length" not in bed_with_psite.columns
    ):
        return pd.DataFrame(columns=cols)

    ann = annotation_df.set_index("transcript")[["start_codon", "l_cds"]]
    chrom_to_tx = annotation_df.set_index("sequence_name")["transcript"].to_dict()

    df = bed_with_psite.copy()
    df["transcript"] = df["chrom"].map(chrom_to_tx).fillna(df["chrom"])
    df = df[df["transcript"].isin(ann.index)]

    rows: list[dict] = []
    if df.empty:
        return pd.DataFrame(columns=cols)

    starts = df["transcript"].map(ann["start_codon"]).astype(int)
    cds_lens = df["transcript"].map(ann["l_cds"]).astype(int)
    edge_lo = starts + 3 * int(exclude_start_codons)
    edge_hi = starts + cds_lens - 3 * int(exclude_stop_codons)
    in_cds_mask = (df["P_site"] >= edge_lo) & (df["P_site"] < edge_hi)
    df_cds = df.loc[in_cds_mask].copy()
    df_cds["frame"] = (
        (df_cds["P_site"].astype(int) - starts.loc[df_cds.index]) % 3
    ).astype(int)

    all_lengths = sorted({int(x) for x in bed_with_psite["read_length"].unique()})
    cds_groups = dict(tuple(df_cds.groupby("read_length")))
    total_groups = (
        bed_with_psite.groupby("read_length").size().to_dict()
    )

    for read_length in all_lengths:
        n_total = int(total_groups.get(read_length, 0))
        sub = cds_groups.get(read_length)
        if sub is None or sub.empty:
            rows.append({
                "sample_id": sample,
                "read_length": int(read_length),
                "n_reads_total": n_total,
                "n_reads_cds": 0,
                "frame0_fraction": 0.0,
                "frame1_fraction": 0.0,
                "frame2_fraction": 0.0,
                "dominant_frame": -1,
                "frame0_dominance": 0.0,
                "periodicity_score": 0.0,
                "frame_entropy": float("nan"),
                "include_for_downstream": False,
                "exclusion_reason": "no_cds_reads",
            })
            continue

        n_cds = int(len(sub))
        counts = sub["frame"].value_counts().reindex([0, 1, 2], fill_value=0)
        f0, f1, f2 = (counts.loc[i] / n_cds for i in (0, 1, 2))
        f0, f1, f2 = float(f0), float(f1), float(f2)
        dominance = f0 - max(f1, f2)
        dominant_frame = int(counts.idxmax())
        entropy = -sum(p * np.log2(p) for p in (f0, f1, f2) if p > 0)

        if n_cds < min_cds_reads_per_length:
            include = False
            reason = "low_count"
        elif f0 < min_frame0_fraction:
            include = False
            reason = "weak_periodicity"
        elif dominance < min_frame0_dominance:
            include = False
            reason = "ambiguous_dominance"
        elif entropy > max_frame_entropy:
            include = False
            reason = "high_entropy"
        else:
            include = True
            reason = "none"

        rows.append({
            "sample_id": sample,
            "read_length": int(read_length),
            "n_reads_total": n_total,
            "n_reads_cds": n_cds,
            "frame0_fraction": f0,
            "frame1_fraction": f1,
            "frame2_fraction": f2,
            "dominant_frame": dominant_frame,
            "frame0_dominance": float(dominance),
            "periodicity_score": float(dominance),
            "frame_entropy": float(entropy),
            "include_for_downstream": bool(include),
            "exclusion_reason": reason,
        })

    return pd.DataFrame(rows, columns=cols)


def select_read_lengths_by_periodicity(
    frame_by_length: pd.DataFrame,
) -> dict[str, list[int]]:
    """Group included read lengths by sample.

    Convenience wrapper that returns ``{sample_id: [read_length, ...]}``
    for rows where ``include_for_downstream`` is true. Returns an empty
    dict when the input is empty (caller decides whether that is a
    fall-back-to-pooled situation).
    """
    if frame_by_length.empty:
        return {}
    out: dict[str, list[int]] = {}
    keep = frame_by_length[frame_by_length["include_for_downstream"].astype(bool)]
    for sample_id, group in keep.groupby("sample_id"):
        out[str(sample_id)] = sorted(int(x) for x in group["read_length"])
    return out


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


def _write_frame_by_length_tsv(
    frame_tables: list[pd.DataFrame], path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not frame_tables:
        path.write_text(
            "sample_id\tread_length\tn_reads_total\tn_reads_cds\t"
            "frame0_fraction\tframe1_fraction\tframe2_fraction\t"
            "dominant_frame\tframe0_dominance\tperiodicity_score\t"
            "frame_entropy\tinclude_for_downstream\texclusion_reason\n",
            encoding="utf-8",
        )
        return
    combined = pd.concat(frame_tables, ignore_index=True)
    combined.to_csv(path, sep="\t", index=False, na_rep="")


def _write_length_inclusion_decisions_tsv(
    frame_tables: list[pd.DataFrame], path: Path,
) -> None:
    """Distil ``frame_by_length.tsv`` to the include/exclude decision.

    A reviewer-friendly subset that pairs every (sample, read_length)
    with its ``include`` boolean and the human-readable reason. Mirrors
    riboWaltz's per-length QC table.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    if not frame_tables:
        path.write_text(
            "sample_id\tread_length\tn_reads_cds\tframe0_fraction\t"
            "frame0_dominance\tinclude_for_downstream\texclusion_reason\n",
            encoding="utf-8",
        )
        return
    combined = pd.concat(frame_tables, ignore_index=True)
    keep = [
        "sample_id",
        "read_length",
        "n_reads_cds",
        "frame0_fraction",
        "frame0_dominance",
        "include_for_downstream",
        "exclusion_reason",
    ]
    combined[keep].to_csv(path, sep="\t", index=False, na_rep="")


def _plot_frame_by_length_heatmap(
    frame_tables: list[pd.DataFrame], out_path: Path,
) -> None:
    """Frame-fraction heatmap with read-length on Y, frame on X.

    For each sample produces one column showing frame 0/1/2 fractions
    across every read length encountered. A length that fails the
    inclusion filter is annotated with a marker so the reader can tell
    a "low frame-0 because rejected" cell from a "low frame-0, kept"
    cell.
    """
    import matplotlib.pyplot as plt

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if not frame_tables:
        return
    combined = pd.concat(frame_tables, ignore_index=True)
    if combined.empty:
        return
    samples = sorted(combined["sample_id"].astype(str).unique())
    n_samples = len(samples)
    # Build the layout with an explicit colorbar column so the cbar
    # cannot steal width from the last data panel.
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(3.4 * n_samples + 1.2, 4.2))
    width_ratios = [1.0] * n_samples + [0.06]
    gs = GridSpec(
        nrows=1,
        ncols=n_samples + 1,
        width_ratios=width_ratios,
        wspace=0.18,
        figure=fig,
    )
    data_axes = []
    im = None
    for col_i, sample in enumerate(samples):
        ax = fig.add_subplot(
            gs[0, col_i],
            sharey=data_axes[0] if data_axes else None,
        )
        data_axes.append(ax)
        sub = combined[combined["sample_id"].astype(str) == sample].sort_values(
            "read_length"
        )
        if sub.empty:
            ax.set_visible(False)
            continue
        matrix = sub[["frame0_fraction", "frame1_fraction", "frame2_fraction"]].to_numpy()
        im = ax.imshow(
            matrix, aspect="auto", cmap="viridis",
            vmin=0.0, vmax=1.0,
            extent=(-0.5, 2.5, sub["read_length"].max() + 0.5,
                    sub["read_length"].min() - 0.5),
        )
        # Mark excluded rows with a hatch / outline so the reader can
        # distinguish "low signal because excluded" from "low signal
        # included".
        for _, row in sub.iterrows():
            if not bool(row["include_for_downstream"]):
                ax.add_patch(
                    plt.Rectangle(
                        (-0.5, float(row["read_length"]) - 0.5), 3, 1,
                        fill=False, edgecolor="red", linewidth=0.8,
                    )
                )
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["F0", "F+1", "F+2"])
        ax.set_title(sample, fontsize=10)
        if col_i == 0:
            ax.set_ylabel("read length (nt)")
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
    if im is not None:
        cax = fig.add_subplot(gs[0, -1])
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label("frame fraction")
    fig.suptitle("Per-read-length frame distribution (red border = excluded)")
    fig.subplots_adjust(top=0.88, left=0.08, right=0.94, bottom=0.12)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    if out_path.suffix.lower() == ".png":
        fig.savefig(out_path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


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
    min_cds_reads_per_length: int = _DEFAULT_MIN_CDS_READS_PER_LENGTH,
    min_frame0_fraction: float = _DEFAULT_MIN_FRAME0_FRACTION,
    min_frame0_dominance: float = _DEFAULT_MIN_FRAME0_DOMINANCE,
    max_frame_entropy: float = _DEFAULT_MAX_FRAME_ENTROPY,
    exclude_start_codons: int = 0,
    exclude_stop_codons: int = 0,
    compute_phase_score: bool = False,
    qc_thresholds: dict[str, float] | None = None,
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
    frame_by_length_tables: list[pd.DataFrame] = []

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
            exclude_start_codons=exclude_start_codons,
            exclude_stop_codons=exclude_stop_codons,
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
        frame_by_length_tables.append(
            compute_frame_summary_by_length(
                psite,
                annotation_df,
                sample=sample,
                min_cds_reads_per_length=min_cds_reads_per_length,
                min_frame0_fraction=min_frame0_fraction,
                min_frame0_dominance=min_frame0_dominance,
                max_frame_entropy=max_frame_entropy,
                exclude_start_codons=exclude_start_codons,
                exclude_stop_codons=exclude_stop_codons,
            )
        )

    _write_frame_summary_tsv(frame_rows, output_dir / "frame_summary.tsv")
    _write_metagene_tsv(start_profiles, output_dir / "periodicity_start.tsv")
    _write_metagene_tsv(stop_profiles, output_dir / "periodicity_stop.tsv")
    _write_strand_tsv(strand_rows, output_dir / "strand_sanity.tsv")

    by_length_dir = output_dir / "by_length"
    _write_frame_by_length_tsv(
        frame_by_length_tables, by_length_dir / "frame_by_length.tsv",
    )
    _write_length_inclusion_decisions_tsv(
        frame_by_length_tables,
        by_length_dir / "length_inclusion_decisions.tsv",
    )

    # v0.6.2: emit the publication-grade QC bundle
    # (frame_counts_by_sample_length.tsv + qc_summary.tsv +
    # qc_summary.md + gene_periodicity.tsv) alongside the
    # by_length/* legacy outputs.
    from .periodicity_qc import run_periodicity_qc_bundle as _qc_bundle

    bed_with_psite_concat: pd.DataFrame | None = None
    if frame_by_length_tables:
        # Re-run compute_p_site_positions per sample and stack so the
        # gene-level table has the per-read coordinates.
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
        if psite_frames:
            bed_with_psite_concat = pd.concat(psite_frames, ignore_index=True)

    combined_frame_by_length = (
        pd.concat(frame_by_length_tables, ignore_index=True)
        if frame_by_length_tables
        else pd.DataFrame()
    )
    _qc_bundle(
        frame_by_length=combined_frame_by_length,
        bed_with_psite=bed_with_psite_concat,
        annotation_df=annotation_df,
        samples=samples,
        output_dir=output_dir,
        site_type=str(offset_site).lower() or "p",
        thresholds=qc_thresholds,
        compute_phase_score=compute_phase_score,
        exclude_start_codons=exclude_start_codons,
        exclude_stop_codons=exclude_stop_codons,
    )

    # Persist the thresholds the inclusion calls were made under so a
    # reviewer can verify them from disk without re-running the CLI.
    by_length_dir.mkdir(parents=True, exist_ok=True)
    (by_length_dir / "periodicity.metadata.json").write_text(
        json.dumps(
            {
                "thresholds": {
                    "min_cds_reads_per_length": int(min_cds_reads_per_length),
                    "min_frame0_fraction": float(min_frame0_fraction),
                    "min_frame0_dominance": float(min_frame0_dominance),
                    "max_frame_entropy": float(max_frame_entropy),
                },
                "exclude_start_codons": int(exclude_start_codons),
                "exclude_stop_codons": int(exclude_stop_codons),
                "phase_score_enabled": bool(compute_phase_score),
                "frame_formula": "(P_site_nt - CDS_start_nt) % 3",
                "frame_0_definition": (
                    "assigned P-site lies in the annotated coding frame"
                ),
                "site": "P-site",
                "offset_type": str(offset_type),
                "offset_site": str(offset_site),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    if plot:
        _plot_metagene_panels(
            start_profiles, stop_profiles,
            out_path=output_dir / "periodicity_metagene.png",
        )
        _plot_frame_by_length_heatmap(
            frame_by_length_tables,
            out_path=by_length_dir / "frame_by_length_heatmap.png",
        )

    return {
        "frame_summary": frame_rows,
        "periodicity_start": start_profiles,
        "periodicity_stop": stop_profiles,
        "strand_sanity": strand_rows,
        "frame_by_length": (
            pd.concat(frame_by_length_tables, ignore_index=True)
            if frame_by_length_tables
            else pd.DataFrame()
        ),
    }
