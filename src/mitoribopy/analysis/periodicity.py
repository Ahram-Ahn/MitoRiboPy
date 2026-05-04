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
    "PerTranscriptStrandSanity",
    "compute_p_site_positions",
    "compute_metagene",
    "compute_strand_sanity",
    "compute_per_transcript_strand_sanity",
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
    """Per-sample, per-relative-position P-site density.

    The ``normalize`` field records how per-transcript signals were
    aggregated into ``density``:

    * ``"per_gene_unit_mean"`` (the default) — each transcript's
      per-position density is divided by its own mean before being
      averaged with the others. Removes the depth-weighting bias where
      one high-expression transcript dominates the metagene shape.
    * ``"none"`` — raw per-position counts summed across transcripts;
      proportional to depth × expression. Kept as an opt-in for
      preserving older depth-weighted metagene semantics.

    ``n_transcripts`` is the number of transcripts that contributed
    (the qualifying-tracks count after any zero-mean filter under
    ``per_gene_unit_mean``). Zero means the profile is the empty
    placeholder.
    """

    sample: str
    anchor: str  # 'start' | 'stop'
    positions: np.ndarray  # 1D int array (relative-to-anchor nt offset)
    density: np.ndarray  # 1D float array, same length as positions
    normalize: str = "per_gene_unit_mean"
    n_transcripts: int = 0


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


@dataclass(frozen=True)
class PerTranscriptStrandSanity:
    """Per-(sample, transcript) minus-strand fraction.

    The per-sample summary (:class:`StrandSanity`) can hide a localized
    antisense bleed-through on one mt-mRNA: a single misannotated
    transcript or an unmasked NUMT can produce many minus-strand reads
    on one chromosome while every other chromosome is clean. Splitting
    the audit by transcript surfaces that pattern at a glance.
    """

    sample: str
    transcript: str
    n_total: int
    n_minus: int
    minus_fraction: float


# ---------------------------------------------------------------------------
# Core computation — pure functions on DataFrames + offset tables
# ---------------------------------------------------------------------------


_OFFSET_FALLBACK_REPORTED: set[tuple[str, int]] = set()


def _reset_offset_fallback_reporter() -> None:
    """Test hook: clear the per-process dedup set for the fallback warning."""
    _OFFSET_FALLBACK_REPORTED.clear()


def _resolve_offsets_for_sample(
    selected_offsets_by_sample: dict[str, pd.DataFrame] | None,
    selected_offsets_combined: pd.DataFrame | None,
    sample: str,
    *,
    record_warning: bool = True,
    stage: str = "PERIODICITY",
) -> pd.DataFrame | None:
    """Return the per-length 5' offset table for ``sample``.

    Falls back to the combined table when no per-sample table is
    available (e.g. ``--offset_mode combined``, or a per-sample table
    that simply does not include this sample). Returns ``None`` when
    neither source is populated.

    When the per-sample dictionary IS supplied but the requested sample
    is missing from it, we record a structured warning
    (``E_OFFSET_FALLBACK_USED``) so a reviewer reading
    ``warnings.tsv`` / ``run_manifest.json`` can see at a glance that
    the run silently degraded from per-sample to combined offsets for
    that sample. The opposite case — no per-sample dict at all,
    intentional combined-mode runs — is silent. Per (stage, sample)
    deduplication keeps the log free of repeats when the resolver is
    called many times for the same sample (start metagene + stop
    metagene + Fourier path in one run).
    """
    if selected_offsets_by_sample is None:
        return selected_offsets_combined
    if sample in selected_offsets_by_sample:
        return selected_offsets_by_sample[sample]
    if record_warning and selected_offsets_combined is not None:
        key = (str(stage), id(selected_offsets_by_sample))
        # Dedup per (stage, dict-instance, sample) so multi-call sites
        # do not flood warnings.tsv.
        sample_key = (str(stage), hash((id(selected_offsets_by_sample), sample)))
        if sample_key not in _OFFSET_FALLBACK_REPORTED:
            _OFFSET_FALLBACK_REPORTED.add(sample_key)
            from ..errors import E_OFFSET_FALLBACK_USED
            from ..io import warnings_log

            warnings_log.record(
                stage=stage,
                sample_id=str(sample),
                code=E_OFFSET_FALLBACK_USED,
                severity="warn",
                message=(
                    f"Per-sample offsets requested but sample {sample!r} is "
                    f"missing from the per-sample table; falling back to the "
                    f"combined-across-samples offsets. Downstream codon-usage "
                    f"and periodicity outputs may be biased toward the "
                    f"average-of-other-samples offset choice."
                ),
                suggested_action=(
                    "Inspect the offset_diagnostics/ outputs for this "
                    "sample, and either pool more reads (e.g. concatenate "
                    "technical replicates before alignment) or rerun with "
                    "--offset_mode combined and document that choice in the "
                    "manuscript methods."
                ),
            )
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


_VALID_METAGENE_NORMALIZATIONS = ("per_gene_unit_mean", "none")


def compute_metagene(
    bed_with_psite: pd.DataFrame,
    annotation_df: pd.DataFrame,
    *,
    sample: str,
    anchor: str,
    window_nt: int = DEFAULT_WINDOW_NT,
    normalize: str = "per_gene_unit_mean",
) -> MetageneProfile:
    """Build a start- or stop-anchored P-site density metagene.

    For ``anchor='start'`` the relative-position axis is ``[0,
    window_nt)`` (P-site nt offset from the start codon). For
    ``anchor='stop'`` the axis is ``(-window_nt, 0]`` (P-site nt offset
    from the first nt of the stop codon).

    Aggregation modes (``normalize``):

    * ``"per_gene_unit_mean"`` (default) — each contributing
      transcript's per-position density is divided by its own mean
      across the window before being averaged with the others. The
      resulting metagene reflects the *shape* shared across transcripts,
      not the depth × expression product, so a single high-expression
      transcript no longer dominates. This matches the per-gene
      normalisation already used in the Fourier-spectrum path
      (:mod:`mitoribopy.analysis.fourier_spectrum`); using it here too
      keeps the headline TSV / plot consistent with the QC story.
    * ``"none"`` — raw per-position counts summed across transcripts;
      proportional to depth × expression. Kept as an opt-in for
      preserving the older depth-weighted interpretation verbatim.

    Migration note: scripts that parsed ``metagene_*.tsv`` and assumed
    integer-count semantics now see fractional means. Pass
    ``normalize='none'`` to recover the old behaviour.
    """
    if anchor not in {"start", "stop"}:
        raise ValueError(f"anchor must be 'start' or 'stop', got {anchor!r}")
    norm = str(normalize)
    if norm not in _VALID_METAGENE_NORMALIZATIONS:
        raise ValueError(
            f"normalize must be one of {_VALID_METAGENE_NORMALIZATIONS}, "
            f"got {normalize!r}"
        )

    if anchor == "start":
        positions = np.arange(0, window_nt, dtype=int)
    else:
        positions = np.arange(-window_nt + 1, 1, dtype=int)
    density = np.zeros(window_nt, dtype=float)

    if bed_with_psite.empty or "P_site" not in bed_with_psite.columns:
        return MetageneProfile(
            sample=sample, anchor=anchor, positions=positions, density=density,
            normalize=norm, n_transcripts=0,
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
            normalize=norm, n_transcripts=0,
        )

    if anchor == "start":
        anchor_pos_all = df["transcript"].map(ann["start_codon"]).astype(int)
        rel_all = df["P_site"].astype(int) - anchor_pos_all
        mask_window = (rel_all >= 0) & (rel_all < window_nt)
    else:
        anchor_pos_all = df["transcript"].map(ann["stop_codon"]).astype(int)
        rel_all = df["P_site"].astype(int) - anchor_pos_all
        # rel range we want: (-window_nt, 0]  (axis index = rel + window_nt - 1)
        mask_window = (rel_all > -window_nt) & (rel_all <= 0)

    df_window = df.loc[mask_window].copy()
    df_window["_rel"] = rel_all.loc[mask_window]
    if df_window.empty:
        return MetageneProfile(
            sample=sample, anchor=anchor, positions=positions, density=density,
            normalize=norm, n_transcripts=0,
        )

    if anchor == "start":
        df_window["_idx"] = df_window["_rel"].astype(int)
    else:
        df_window["_idx"] = (df_window["_rel"] + window_nt - 1).astype(int)

    if norm == "none":
        np.add.at(density, df_window["_idx"].to_numpy(), 1.0)
        return MetageneProfile(
            sample=sample, anchor=anchor, positions=positions, density=density,
            normalize=norm,
            n_transcripts=int(df_window["transcript"].nunique()),
        )

    # per_gene_unit_mean: build per-transcript density vectors, divide
    # each by its own mean, then average across qualifying transcripts.
    n_qualified = 0
    accumulator = np.zeros(window_nt, dtype=float)
    for tx, sub in df_window.groupby("transcript", sort=False):
        per_tx = np.zeros(window_nt, dtype=float)
        np.add.at(per_tx, sub["_idx"].to_numpy(), 1.0)
        mean_per_tx = float(per_tx.mean())
        if mean_per_tx <= 0:
            continue
        accumulator += per_tx / mean_per_tx
        n_qualified += 1

    if n_qualified > 0:
        density = accumulator / n_qualified
    return MetageneProfile(
        sample=sample, anchor=anchor, positions=positions, density=density,
        normalize=norm, n_transcripts=n_qualified,
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


def compute_per_transcript_strand_sanity(
    bed_df: pd.DataFrame, *, sample: str,
) -> list[PerTranscriptStrandSanity]:
    """Per-transcript minus-strand fraction for one sample.

    The per-sample summary collapses an antisense bleed-through
    localised to one mt-mRNA into the all-transcripts mean. Splitting
    by transcript surfaces a single misannotated reference contig or
    an unmasked NUMT pulling reads onto the wrong strand. Returns an
    empty list when the BED has no rows for ``sample`` or no
    ``strand`` column.
    """
    sub = bed_df[bed_df["sample_name"] == sample]
    if sub.empty or "strand" not in sub.columns or "chrom" not in sub.columns:
        return []
    minus = sub["strand"].astype(str).str.strip() == "-"
    rows: list[PerTranscriptStrandSanity] = []
    for chrom, group in sub.groupby("chrom", sort=True):
        n_total = int(len(group))
        n_minus = int(minus.loc[group.index].sum())
        rows.append(PerTranscriptStrandSanity(
            sample=str(sample),
            transcript=str(chrom),
            n_total=n_total,
            n_minus=n_minus,
            minus_fraction=(n_minus / n_total) if n_total else 0.0,
        ))
    return rows


# ---------------------------------------------------------------------------
# Orchestrator + writers
# ---------------------------------------------------------------------------


_DEPRECATED_PERIODICITY_TSVS: frozenset[str] = frozenset({
    "periodicity_start.tsv",
    "periodicity_stop.tsv",
})


def _write_metagene_tsv(
    profiles: Iterable[MetageneProfile], path: Path,
) -> None:
    """Write a metagene TSV with the normalisation mode in a header comment.

    The header comment carries the normalisation policy so downstream
    consumers can detect at a glance whether ``density`` is a
    per-transcript-mean-then-mean (``per_gene_unit_mean``) or a
    raw-position-count sum (``none``). pandas / csv readers that pass
    ``comment="#"`` ignore the line; readers that don't see a leading
    comment line that's still safe to skip.

    The legacy ``periodicity_{start,stop}.tsv`` filenames are still
    written verbatim for backwards compatibility and carry identical
    data to the spec-compliant ``metagene_{start,stop}.tsv`` filenames.
    A ``# DEPRECATED`` comment is added to the legacy paths so a human
    reading the file sees the warning.
    """
    profiles_list = list(profiles)
    path.parent.mkdir(parents=True, exist_ok=True)
    norm_modes = sorted({p.normalize for p in profiles_list}) or ["per_gene_unit_mean"]
    is_deprecated = path.name in _DEPRECATED_PERIODICITY_TSVS
    with path.open("w", encoding="utf-8") as h:
        if is_deprecated:
            h.write(
                "# DEPRECATED: periodicity_{start,stop}.tsv will be removed "
                "in v1.0.0. Read metagene_{start,stop}.tsv instead "
                "(identical contents).\n"
            )
        h.write(f"# normalize: {','.join(norm_modes)}\n")
        h.write("sample\tanchor\tposition\tdensity\tnormalize\tn_transcripts\n")
        for p in profiles_list:
            for pos, dens in zip(p.positions, p.density):
                h.write(
                    f"{p.sample}\t{p.anchor}\t{int(pos)}\t{dens:.6g}"
                    f"\t{p.normalize}\t{int(p.n_transcripts)}\n"
                )


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


def _write_per_transcript_strand_tsv(
    rows: Iterable[PerTranscriptStrandSanity], path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as h:
        h.write("sample\ttranscript\tn_total\tn_minus\tminus_fraction\n")
        for r in rows:
            h.write(
                f"{r.sample}\t{r.transcript}\t{r.n_total}\t{r.n_minus}"
                f"\t{r.minus_fraction:.6g}\n"
            )


def _metagene_y_label(profile: MetageneProfile | None, *, site_letter: str) -> str:
    """Pick the right y-axis label given the profile's normalisation mode.

    With per-gene unit-mean normalisation, the y-axis carries
    fractional density values that are NOT raw counts; mislabelling
    them as "reads" would mislead a reviewer about the scale.
    """
    site_label = "P-site" if str(site_letter).lower() == "p" else "A-site"
    if profile is None:
        return f"{site_label} density"
    if profile.normalize == "per_gene_unit_mean":
        return f"{site_label} density (per-gene unit-mean)"
    return f"{site_label} reads"


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
                ax.set_ylabel(_metagene_y_label(profile, site_letter="p"))
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
        ax.set_ylabel(_metagene_y_label(profile, site_letter=site_letter))
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
    # How the per-transcript signals are aggregated into the
    # start- / stop-anchored metagene density. See compute_metagene().
    metagene_normalize: str = "per_gene_unit_mean",
    # Statistical hardening. All have library defaults; pass
    # `compute_stats=False` to skip bootstrap + permutation entirely.
    fourier_n_bootstrap: int | None = None,
    fourier_n_permutations: int | None = None,
    fourier_ci_alpha: float | None = None,
    fourier_random_seed: int | None = None,
    fourier_compute_stats: bool = True,
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
    are not emitted by the current package; the Fourier ``spectral_ratio_3nt``
    + ``snr_call`` columns now carry the headline QC story.
    """
    output_dir = Path(output_dir)
    effective_window_nt = int(metagene_nt) if metagene_nt is not None else int(window_nt)
    start_profiles: list[MetageneProfile] = []
    stop_profiles: list[MetageneProfile] = []
    strand_rows: list[StrandSanity] = []
    per_tx_strand_rows: list[PerTranscriptStrandSanity] = []

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
            normalize=str(metagene_normalize),
        ))
        stop_profiles.append(compute_metagene(
            psite, annotation_df,
            sample=sample, anchor="stop", window_nt=effective_window_nt,
            normalize=str(metagene_normalize),
        ))
        strand_rows.append(compute_strand_sanity(bed_df, sample=sample))
        per_tx_strand_rows.extend(
            compute_per_transcript_strand_sanity(bed_df, sample=sample)
        )

    _write_metagene_tsv(start_profiles, output_dir / "periodicity_start.tsv")
    _write_metagene_tsv(stop_profiles, output_dir / "periodicity_stop.tsv")
    _write_metagene_tsv(start_profiles, output_dir / "metagene_start.tsv")
    _write_metagene_tsv(stop_profiles, output_dir / "metagene_stop.tsv")
    _write_strand_tsv(strand_rows, output_dir / "strand_sanity.tsv")
    _write_per_transcript_strand_tsv(
        per_tx_strand_rows, output_dir / "strand_sanity_per_transcript.tsv",
    )

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
    if fourier_n_bootstrap is not None:
        fourier_kwargs["n_bootstrap"] = int(fourier_n_bootstrap)
    if fourier_n_permutations is not None:
        fourier_kwargs["n_permutations"] = int(fourier_n_permutations)
    if fourier_ci_alpha is not None:
        fourier_kwargs["ci_alpha"] = float(fourier_ci_alpha)
    if fourier_random_seed is not None:
        fourier_kwargs["random_seed"] = int(fourier_random_seed)

    _qc_bundle(
        bed_with_psite=bed_with_psite_concat,
        annotation_df=annotation_df,
        samples=samples,
        output_dir=output_dir,
        site_type=str(offset_site).lower() or "p",  # type: ignore[arg-type]
        render_plots=fourier_render_plots,
        compute_stats=bool(fourier_compute_stats),
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
        "strand_sanity_per_transcript": per_tx_strand_rows,
    }
