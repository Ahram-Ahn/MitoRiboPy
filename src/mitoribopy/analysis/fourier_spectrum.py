"""Wakigawa-style discrete-Fourier-transform periodicity analysis.

Replaces the legacy single-scalar ``fft_period3_power`` path. The
contract here mirrors the published mt-Ribo-seq analysis from Wakigawa
et al. (Mol Cell, 2025) — the only published DFT recipe that
distinguishes translating mitoribosomes from background by comparing
the period-3 amplitude of the **ORF window upstream of the canonical
stop codon** with the matched **3′ UTR window downstream of it**.

Why this replaces the old path
------------------------------

The legacy ``calculate_fft_period3_ratio`` collapsed each (sample, gene)
into ``power[k=1/3] / median(background)`` over the whole CDS. Three
hidden weaknesses motivated the rewrite:

1. **Single number.** A reviewer cannot see whether the peak is
   actually at period 3 or just leaked in from neighbours. The full
   amplitude curve is the auditable artefact.
2. **No 3′ UTR negative control.** The strongest defensible claim is
   "ORF periodic AND 3′ UTR flat". Without the matched downstream DFT
   the pipeline cannot verify it.
3. **Whole-CDS window.** Long mt-mRNAs get dominated by transcript-
   length-dependent low-frequency noise. Wakigawa fixes this with a
   fixed 100-nt window anchored at the stop codon — comparable across
   genes regardless of CDS length.

Window convention
-----------------

For a transcript with annotated ``stop_codon`` (the 0-based first nt
of the stop codon trinucleotide):

* ORF window  = ``[stop_codon - W, stop_codon)``     (W nt strictly upstream)
* 3′ UTR      = ``[stop_codon + 3, stop_codon + 3 + W)``  (W nt strictly downstream)

The 3′ UTR window starts at ``stop_codon + 3`` so the stop-codon
trinucleotide itself is never included in either panel. The default W
is 100 nt — Wakigawa's published choice.

Read-length stratification
--------------------------

The DFT is computed per (sample, read_length) — Wakigawa publishes
"32 nt" panels because pooling lengths can mask a single bad-length
class that's contaminating a healthy library. Every read length present
in the input BED gets its own row in the long-format table; the plot
layer renders one figure per (sample, read_length).

Overlap-upstream gene policy
----------------------------

For human mt-mRNAs, the stop codons of MT-ATP8 and MT-ND4L sit
*inside* the downstream MT-ATP6 and MT-ND4 ORFs respectively. Their
"3′ UTR" 100-nt window is therefore not really 3′ UTR — it's another
ORF being translated, so it can show period-3 signal that has nothing
to do with stop-codon termination behaviour. These two genes are
flagged via :func:`is_overlap_upstream_orf` and EXCLUDED from the
combined per-(sample, length, region) summary. They are still scored
per-gene so a reviewer can inspect them deliberately.

Note that fused-FASTA spellings (``ATP86``, ``ND4L4``) are NOT
excluded — those represent the whole combined region whose canonical
stop codon is at the natural transcript end, so the 3′ UTR window is
legitimate.

Amplitude normalization
-----------------------

The published figure plots a one-sided amplitude spectrum on the
``|FFT[k]| * 2 / N`` scale (k > 0; the standard convention for
estimating sinusoidal amplitude from a finite signal). This is the
same scale matplotlib / scipy users see by default and is depth-
sensitive on raw counts; per-gene scaling (e.g. divide by mean window
density) is intentionally NOT applied here — Wakigawa overlays absolute
amplitudes per gene so high-coverage genes are visually dominant, and
flattening by depth would erase that signal. Downstream consumers that
want depth-normalised scores can divide by ``n_sites``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal

import numpy as np
import pandas as pd


__all__ = [
    "DEFAULT_WINDOW_NT",
    "DEFAULT_PERIOD_RANGE",
    "OVERLAP_UPSTREAM_GENES",
    "Region",
    "Site",
    "FourierWindow",
    "compute_amplitude_spectrum",
    "build_window_coverage",
    "build_fourier_spectrum_table",
    "build_period3_score_table",
    "build_period3_summary_table",
    "is_overlap_upstream_orf",
]


DEFAULT_WINDOW_NT: int = 100
"""Width in nucleotides of the ORF / 3' UTR window (Wakigawa default)."""

DEFAULT_PERIOD_RANGE: tuple[float, float] = (2.0, 10.0)
"""Period range (nt) the published figures span on the x-axis."""

# Genes whose annotated stop codon sits inside another mt-mRNA ORF.
# Their 3' UTR window is contaminated by translation of the downstream
# ORF (MT-ATP6 for MT-ATP8; MT-ND4 for MT-ND4L). They are scored
# per-gene but excluded from combined-genes summaries and combined-
# overlay plots. Fused-FASTA names (e.g. "ATP86") are deliberately NOT
# in this set: those entries cover the whole combined region and have a
# legitimate 3' UTR downstream of the natural transcript-end stop codon.
OVERLAP_UPSTREAM_GENES: frozenset[str] = frozenset({
    "MT-ATP8", "MTATP8", "ATP8",
    "MT-ND4L", "MTND4L", "ND4L",
})


Region = Literal["orf", "utr3"]
Site = Literal["a", "p", "5p"]


def is_overlap_upstream_orf(name: str | None) -> bool:
    """Return True only for genes whose stop codon is INSIDE another ORF.

    Matching is case-insensitive and treats ``_`` as ``-``. See
    :data:`OVERLAP_UPSTREAM_GENES` for the recognised spellings.
    """
    if name is None:
        return False
    s = str(name).strip().upper().replace("_", "-")
    return s in OVERLAP_UPSTREAM_GENES


# ---------------------------------------------------------------------------
# Window definition
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FourierWindow:
    """Resolved nt-coordinate window for one (transcript, region) pair.

    ``lo`` is inclusive, ``hi`` is exclusive — a half-open interval over
    transcript-local nt coordinates. ``valid`` is False when the window
    falls off either end of the transcript (the caller drops these).
    """

    transcript: str
    region: Region
    lo: int
    hi: int
    valid: bool

    @property
    def width(self) -> int:
        return max(0, self.hi - self.lo)


def _resolve_window(
    *,
    transcript: str,
    region: Region,
    stop_codon: int,
    transcript_length: int,
    window_nt: int,
) -> FourierWindow:
    """Compute the nt-coordinate window for one (transcript, region)."""
    if region == "orf":
        lo = stop_codon - window_nt
        hi = stop_codon
    elif region == "utr3":
        lo = stop_codon + 3
        hi = stop_codon + 3 + window_nt
    else:  # pragma: no cover - Literal guards this at the boundary
        raise ValueError(f"region must be 'orf' or 'utr3', got {region!r}")
    valid = lo >= 0 and hi <= transcript_length and (hi - lo) == window_nt
    return FourierWindow(
        transcript=transcript, region=region, lo=lo, hi=hi, valid=valid,
    )


# ---------------------------------------------------------------------------
# Per-window coverage + spectrum
# ---------------------------------------------------------------------------


def build_window_coverage(
    site_positions: np.ndarray,
    *,
    window: FourierWindow,
) -> np.ndarray:
    """Return a length-W coverage vector for the given window.

    *site_positions* is a 1-D integer array of P-site or A-site nt
    positions (transcript-local, same coordinate system as
    ``window.lo``/``window.hi``). Positions outside the window are
    silently dropped.
    """
    width = window.width
    out = np.zeros(width, dtype=float)
    if width == 0 or site_positions.size == 0:
        return out
    in_window = (site_positions >= window.lo) & (site_positions < window.hi)
    if not np.any(in_window):
        return out
    idx = site_positions[in_window].astype(np.int64) - window.lo
    np.add.at(out, idx, 1.0)
    return out


def compute_amplitude_spectrum(
    coverage: np.ndarray,
    *,
    period_range: tuple[float, float] = DEFAULT_PERIOD_RANGE,
) -> pd.DataFrame:
    """Return ``DataFrame[period_nt, amplitude]`` from a 1-D coverage signal.

    Subtracts the mean (so the DC component is dropped), runs
    :func:`numpy.fft.rfft`, converts frequency bins (cycles per nt) to
    period bins (nt per cycle), and restricts the result to the
    requested period range (default 2-10 nt).

    Amplitude is the single-sided spectrum ``|FFT[k]| * 2 / N`` for
    k > 0 — the convention that gives back the peak amplitude of a
    pure sinusoid in the input. See module docstring for why we do
    NOT depth-normalise.

    Returns an empty DataFrame for inputs with fewer than 2 *
    ``ceil(period_range[1])`` samples or zero variance — there isn't
    enough signal to estimate a spectrum at that resolution.
    """
    x = np.asarray(coverage, dtype=float)
    n = x.size
    period_lo, period_hi = float(period_range[0]), float(period_range[1])
    if n < int(np.ceil(period_hi)) * 2:
        return pd.DataFrame(columns=["period_nt", "amplitude"])
    x_centered = x - float(np.nanmean(x))
    if not np.isfinite(x_centered).any() or float(np.nanvar(x_centered)) <= 1e-12:
        return pd.DataFrame(columns=["period_nt", "amplitude"])
    spectrum = np.abs(np.fft.rfft(np.nan_to_num(x_centered)))
    if n > 0:
        spectrum = spectrum * (2.0 / n)
    freqs = np.fft.rfftfreq(n, d=1.0)
    if freqs.size < 2:
        return pd.DataFrame(columns=["period_nt", "amplitude"])
    nonzero = freqs > 0
    periods = 1.0 / freqs[nonzero]
    amplitudes = spectrum[nonzero]
    keep = (periods >= period_lo) & (periods <= period_hi)
    return pd.DataFrame({
        "period_nt": periods[keep],
        "amplitude": amplitudes[keep],
    })


# ---------------------------------------------------------------------------
# Per (sample, length, gene, region) long-format table
# ---------------------------------------------------------------------------


_ANNOTATION_REQUIRED_COLS: tuple[str, ...] = (
    "transcript", "sequence_name", "stop_codon", "l_tr",
)
_BED_REQUIRED_COLS: tuple[str, ...] = (
    "sample_name", "chrom", "read_length", "P_site",
)


def _resolve_site_offset(site: Site) -> int:
    """Translate a site-name into a nt offset added to ``P_site``.

    Wakigawa uses A-site assignment, so the default is +3.
    """
    site = str(site).lower()  # type: ignore[assignment]
    if site == "p":
        return 0
    if site == "a":
        return 3
    if site == "5p":
        # Caller must supply ``5p_site`` instead of ``P_site``; offset 0
        # is recorded in metadata for provenance only.
        return 0
    raise ValueError(f"site must be 'a', 'p', or '5p', got {site!r}")


def build_fourier_spectrum_table(
    bed_with_psite_and_gene: pd.DataFrame,
    annotation_df: pd.DataFrame,
    *,
    samples: Iterable[str],
    window_nt: int = DEFAULT_WINDOW_NT,
    period_range: tuple[float, float] = DEFAULT_PERIOD_RANGE,
    site: Site = "a",
) -> pd.DataFrame:
    """Long-format ``[sample, read_length, gene, region, period_nt, amplitude]`` table.

    Iterates every (sample, read_length, transcript, region) cell,
    builds the windowed coverage, and emits the amplitude spectrum.
    Genes whose window falls off the transcript ends are dropped.

    Adds an ``is_overlap_upstream_orf`` boolean column so downstream
    consumers can split the combined-data view from the diagnostic
    view without re-deriving the predicate.
    """
    cols = [
        "sample", "read_length", "gene", "transcript_id",
        "region", "n_sites", "n_nt", "period_nt", "amplitude",
        "is_overlap_upstream_orf",
    ]
    if bed_with_psite_and_gene is None or bed_with_psite_and_gene.empty:
        return pd.DataFrame(columns=cols)
    missing = [c for c in _BED_REQUIRED_COLS if c not in bed_with_psite_and_gene.columns]
    if missing:
        raise ValueError(
            f"bed_with_psite_and_gene is missing required column(s): {missing}"
        )
    if annotation_df is None or annotation_df.empty:
        return pd.DataFrame(columns=cols)
    ann_missing = [c for c in _ANNOTATION_REQUIRED_COLS if c not in annotation_df.columns]
    if ann_missing:
        raise ValueError(
            f"annotation_df is missing required column(s): {ann_missing}"
        )

    samples_set = {str(s) for s in samples}
    if not samples_set:
        return pd.DataFrame(columns=cols)

    site_offset = _resolve_site_offset(site)

    ann = annotation_df.set_index("transcript")[["stop_codon", "l_tr"]]
    chrom_to_tx = annotation_df.set_index("sequence_name")["transcript"].to_dict()

    df = bed_with_psite_and_gene.copy()
    df["transcript"] = df["chrom"].map(chrom_to_tx).fillna(df["chrom"])
    df = df[df["transcript"].isin(ann.index)]
    df = df[df["sample_name"].astype(str).isin(samples_set)]
    if df.empty:
        return pd.DataFrame(columns=cols)
    df["site_pos"] = df["P_site"].astype(int) + site_offset
    df["read_length"] = df["read_length"].astype(int)

    rows: list[dict] = []
    group_cols = ["sample_name", "transcript", "read_length"]
    for (sample, transcript, read_length), group in df.groupby(group_cols):
        try:
            stop_codon = int(ann.loc[transcript, "stop_codon"])
            l_tr = int(ann.loc[transcript, "l_tr"])
        except (KeyError, TypeError, ValueError):
            continue

        site_positions = group["site_pos"].astype(int).to_numpy()
        for region in ("orf", "utr3"):
            window = _resolve_window(
                transcript=transcript,
                region=region,
                stop_codon=stop_codon,
                transcript_length=l_tr,
                window_nt=window_nt,
            )
            if not window.valid:
                continue
            coverage = build_window_coverage(site_positions, window=window)
            n_sites = float(coverage.sum())
            if n_sites <= 0:
                continue
            spectrum = compute_amplitude_spectrum(
                coverage, period_range=period_range,
            )
            if spectrum.empty:
                continue
            is_overlap = bool(is_overlap_upstream_orf(transcript))
            for _, srow in spectrum.iterrows():
                rows.append({
                    "sample": str(sample),
                    "read_length": int(read_length),
                    "gene": str(transcript),
                    "transcript_id": str(transcript),
                    "region": str(region),
                    "n_sites": int(n_sites),
                    "n_nt": int(window.width),
                    "period_nt": float(srow["period_nt"]),
                    "amplitude": float(srow["amplitude"]),
                    "is_overlap_upstream_orf": is_overlap,
                })
    return pd.DataFrame(rows, columns=cols)


# ---------------------------------------------------------------------------
# Derived per-gene scalars
# ---------------------------------------------------------------------------


def _nearest_period_bin(periods: np.ndarray, target: float) -> int:
    """Index of the period bin closest to *target* (defaults to nt=3)."""
    return int(np.argmin(np.abs(periods - target)))


def build_period3_score_table(
    spectrum_table: pd.DataFrame,
    *,
    target_period: float = 3.0,
    reference_periods: tuple[float, ...] = (4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0),
) -> pd.DataFrame:
    """Per (sample, read_length, gene, region) period-3 amplitude + score.

    For each spectrum cell:

    * ``period3_amplitude`` is the amplitude at the period-bin closest
      to *target_period* (default 3 nt).
    * ``period3_score`` is ``period3_amplitude / mean(amplitude in
      period bins closest to each value of *reference_periods*)``.
      This is the new auditable scalar that replaces the legacy
      ``fft_period3_power`` column. NaN when the reference mean is 0.

    The ``is_overlap_upstream_orf`` column is propagated unchanged so
    callers can filter out the contaminated genes when aggregating.
    """
    cols = [
        "sample", "read_length", "gene", "transcript_id", "region",
        "n_sites", "n_nt", "period3_amplitude", "period3_score",
        "is_overlap_upstream_orf",
    ]
    if spectrum_table is None or spectrum_table.empty:
        return pd.DataFrame(columns=cols)

    rows: list[dict] = []
    group_cols = ["sample", "read_length", "gene", "region"]
    for keys, group in spectrum_table.groupby(group_cols, sort=False):
        sample, read_length, gene, region = keys
        periods = group["period_nt"].to_numpy()
        amps = group["amplitude"].to_numpy()
        if periods.size == 0:
            continue
        target_bin = _nearest_period_bin(periods, target_period)
        period3_amp = float(amps[target_bin])
        ref_indices = sorted({
            _nearest_period_bin(periods, p) for p in reference_periods
        })
        # Drop the target bin itself in case a reference period coincides.
        ref_indices = [i for i in ref_indices if i != target_bin]
        if not ref_indices:
            score = float("nan")
        else:
            ref_mean = float(np.mean(amps[ref_indices]))
            score = period3_amp / ref_mean if ref_mean > 0 else float("nan")
        first_row = group.iloc[0]
        rows.append({
            "sample": str(sample),
            "read_length": int(read_length),
            "gene": str(gene),
            "transcript_id": str(first_row["transcript_id"]),
            "region": str(region),
            "n_sites": int(first_row["n_sites"]),
            "n_nt": int(first_row["n_nt"]),
            "period3_amplitude": period3_amp,
            "period3_score": score,
            "is_overlap_upstream_orf": bool(first_row["is_overlap_upstream_orf"]),
        })
    return pd.DataFrame(rows, columns=cols)


def build_period3_summary_table(
    score_table: pd.DataFrame,
) -> pd.DataFrame:
    """Per (sample, read_length, region) combined summary.

    Two parallel views per cell:

    * ``*_combined`` columns aggregate over genes EXCLUDING any with
      ``is_overlap_upstream_orf=True``. This is the headline number
      that should track translation phasing without contamination from
      the downstream ORF.
    * ``*_overlap_upstream`` columns aggregate over the excluded genes
      only — surfaced separately so reviewers can check that the
      contaminated 3′ UTR window in those genes really does show
      translation-like signal (which is the biological expectation).
    """
    cols = [
        "sample", "read_length", "region",
        "n_genes_combined", "median_period3_score_combined",
        "max_period3_score_combined", "n_genes_overlap_upstream",
        "median_period3_score_overlap_upstream",
        "max_period3_score_overlap_upstream",
    ]
    if score_table is None or score_table.empty:
        return pd.DataFrame(columns=cols)

    rows: list[dict] = []
    group_cols = ["sample", "read_length", "region"]
    for keys, group in score_table.groupby(group_cols, sort=False):
        sample, read_length, region = keys
        combined = group[~group["is_overlap_upstream_orf"]]
        overlap = group[group["is_overlap_upstream_orf"]]
        rows.append({
            "sample": str(sample),
            "read_length": int(read_length),
            "region": str(region),
            "n_genes_combined": int(len(combined)),
            "median_period3_score_combined": (
                float(combined["period3_score"].median())
                if not combined.empty else float("nan")
            ),
            "max_period3_score_combined": (
                float(combined["period3_score"].max())
                if not combined.empty else float("nan")
            ),
            "n_genes_overlap_upstream": int(len(overlap)),
            "median_period3_score_overlap_upstream": (
                float(overlap["period3_score"].median())
                if not overlap.empty else float("nan")
            ),
            "max_period3_score_overlap_upstream": (
                float(overlap["period3_score"].max())
                if not overlap.empty else float("nan")
            ),
        })
    return pd.DataFrame(rows, columns=cols)
