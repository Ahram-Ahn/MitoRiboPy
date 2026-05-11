"""Metagene Fourier periodicity analysis.

Replaces the legacy per-gene FFT overlay (`fft_period3_power.tsv` and the
v0.7.x sum-coverage path). The contract is **aggregate first, then DFT**:
per-gene normalised tracks are averaged into a single metagene before a
single direct DFT is evaluated on the result.

Why aggregate-then-DFT instead of per-gene overlay
--------------------------------------------------

* **Initiation pile-ups create flat broadband spectra.** The first
  ~5 codons after every AUG carry massive ribosome occupancy
  (initiation stalling). The Fourier transform of an impulse is a flat
  spectrum across all frequencies — which smears period-3 energy
  uniformly across periods 2-10 in the ORF figure. Skipping those
  first 5 codons (15 nt) of every transcript removes the contaminant.

* **Stop-codon stalling does the same on the 3' end.** Especially severe
  for noncanonical-stop ORFs (MT-CO1: AGA; MT-ND6: AGG). Including the
  stop trinucleotide adds another impulse, so the last codon (3 nt)
  before the stop is dropped.

* **Per-gene 100 nt is structurally underpowered.** 33 codons of signal
  per gene is too short for clean DFT on noisy data. Aggregating ~9
  ORFs into one metagene before the DFT averages out gene-specific
  noise while reinforcing codon-locked 3-nt phasing.

* **Coverage magnitude dominates raw overlays.** The highest-expression
  gene (MT-ND6 in many libraries) buries every other trace. Per-gene
  unit-mean normalisation equalises contributions before aggregation.

* **Hann window + mean-centering kill spectral leakage.** Without them,
  the DC component leaks into nearby frequencies and creates an upward
  amplitude ramp toward periods 8-10 (the "rising tail" you see in the
  raw figure).

Window convention
-----------------

For a transcript with annotated ``start_codon`` and ``stop_codon`` (each
the 0-based first nt of its trinucleotide):

* ``orf_start`` = ``[start_codon + 15, start_codon + 15 + W)``
  (skip 5 codons of initiation peak, then W nt of elongation signal)
* ``orf_stop``  = ``[stop_codon - 3 - W, stop_codon - 3)``
  (skip 1 codon before stop, then W nt of elongation signal upstream)

Default ``W = 99 nt`` = 33 codons exactly. Avoids periods 3-leakage that
non-multiple-of-3 windows produce.

Per-gene processing pipeline
----------------------------

For each (sample, read_length, transcript, region):

  1. Build P-/A-site coverage vector of length W in window-local coords.
  2. Filter: skip if ``mean(coverage) < 0.1`` or ``sum(coverage) < 30``
     (low-coverage genes contribute noise without signal).
  3. ``x = coverage / mean(coverage)``  -> unit-mean track.
  4. ``x = x - mean(x)``                -> mean-centred.
  5. ``x = x * np.hanning(W)``          -> Hann-windowed (reduces leakage).

Aggregation
-----------

For each (sample, read_length, region, gene_set):

  metagene = element-wise mean across qualifying per-gene tracks.

Three gene_sets are scored:

* ``combined`` — every transcript in the dataset that is NOT part of
  the ATP8/ATP6 or ND4L/ND4 overlap pair. Typically 9 mt-mRNAs in
  human (MT-CO1, MT-CO2, MT-CO3, MT-CYB, MT-ND1, MT-ND2, MT-ND3,
  MT-ND5, MT-ND6).
* ``ATP86`` — junction-bracketed analysis of the bicistronic
  ATP8/ATP6 region. The ``orf_stop`` panel uses the ATP8 transcript
  (window leads up to ATP8's stop at the start of the overlap); the
  ``orf_start`` panel uses the ATP6 transcript (window leaves ATP6's
  start, just past the overlap). The two windows OVERLAP at the
  bicistronic junction, so a period-3 peak in each panel says which
  reading frame is dominant in that region.
* ``ND4L4`` — same idea for the ND4L/ND4 bicistronic pair. The four-
  nt overlap is too short to bracket inside both windows; the windows
  flank the junction instead.

Direct DFT (period-grid evaluation)
-----------------------------------

Standard ``np.fft.rfft`` evaluates at frequency bins
``k / N`` for integer ``k``, which means the period-3 bin lands at
``N / round(N / 3)`` — never exactly 3 except when ``N`` is a multiple
of 3. We sidestep that by computing the DFT directly at arbitrary
periods via

    z          = sum_n x[n] * exp(-2j * pi * n / period)
    amplitude  = |z| / sqrt(sum(x ** 2) + 1e-12)

This evaluates exactly at ``period = 3.0`` and any other period the
user requests. The denominator is the L2 norm of x, so amplitude is
in ``[0, 1]`` and is comparable across windows of different length.

The headline scalar is the **3-nt spectral ratio**:

    spectral_ratio_3nt = amp(3.0) / median( amp(p) for p in 2..10
                                            excluding 2.8 <= p <= 3.2 )

Healthy library: ``ratio >= 5x``. Excellent: ``ratio >= 10x``.
Frame-shifted / unrecycled: ``ratio < 2x``. The thresholds map to the
``snr_call`` column ({"excellent", "healthy", "modest", "broken"}).

Annotation contract
-------------------

This module reads ``transcript``, ``sequence_name``, ``sequence_aliases``,
``start_codon``, ``stop_codon``, and ``l_tr`` from the annotation table.
The chrom-to-transcript lookup checks both ``sequence_name`` and the
semicolon-separated ``sequence_aliases``, so BED reads from a fused-
FASTA chromosome (``ATP86``, ``ND4L4``) resolve to BOTH constituent
transcripts (ATP8 + ATP6, or ND4L + ND4) for analysis. This is what
makes the ATP86 / ND4L4 figures work without code-side aliasing.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal

import numpy as np
import pandas as pd
from scipy import signal as _sp_signal


__all__ = [
    "DEFAULT_WINDOW_NT",
    "DEFAULT_DROP_CODONS_AFTER_START",
    "DEFAULT_DROP_CODONS_BEFORE_STOP",
    "DEFAULT_LOCAL_BACKGROUND_RANGE",
    "DEFAULT_PERIOD_GRID",
    "DEFAULT_MIN_MEAN_COVERAGE",
    "DEFAULT_MIN_TOTAL_COUNTS",
    "REGIONS",
    "GENE_SETS",
    "OVERLAP_PAIR_TRANSCRIPTS",
    "ATP86_GENE_SET",
    "ND4L4_GENE_SET",
    "Region",
    "GeneSet",
    "Site",
    "FourierWindow",
    "build_chrom_to_transcripts",
    "extract_per_gene_normalized_tracks",
    "build_metagene",
    "direct_dft_amplitude",
    "compute_spectrum_grid",
    "compute_spectral_ratio_3nt",
    "compute_spectral_ratio_3nt_local",
    "snr_call_for_ratio",
    "DEFAULT_BOOTSTRAP_N",
    "DEFAULT_PERMUTATIONS_N",
    "DEFAULT_CI_ALPHA",
    "DEFAULT_RANDOM_SEED",
    "MIN_GENES_FOR_BOOTSTRAP_CI",
    "build_basis_matrix",
    "bootstrap_period3_ci",
    "circular_shift_permutation_p",
    "build_fourier_spectrum_combined_table",
    "build_fourier_period3_score_combined_table",
]


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_WINDOW_NT: int = 99
"""33 codons of signal (multiple of 3) — clean integer-bin period-3."""

DEFAULT_DROP_CODONS_AFTER_START: int = 5
"""Drop 5 codons (15 nt) after AUG to skip the initiation peak."""

DEFAULT_DROP_CODONS_BEFORE_STOP: int = 1
"""Drop 1 codon (3 nt) before stop to skip the termination peak."""

DEFAULT_MIN_MEAN_COVERAGE: float = 0.1
"""Skip per-gene windows whose mean coverage is below this threshold."""

DEFAULT_MIN_TOTAL_COUNTS: int = 30
"""Skip per-gene windows whose total site count is below this threshold."""

DEFAULT_PERIOD_GRID: np.ndarray = np.round(np.arange(2.0, 10.05, 0.05), 2)
"""Period grid (nt) for the publication-style spectrum plot."""

DEFAULT_LOCAL_BACKGROUND_RANGE: tuple[float, float] = (4.0, 6.0)
"""Codon-scale period window used for the local-background ratio.

The default global-background ratio (``spectral_ratio_3nt``) divides
amp(3) by the median of amplitudes over the full 2-10 nt range
(excluding 2.8-3.2). When the metagene carries strong long-period
structure (5'->3' density gradients, broad pause clusters), that long-
period signal inflates the global median and can mask a genuine
period-3 peak. The local ratio (``spectral_ratio_3nt_local``) divides
amp(3) by the median over a NARROW codon-scale window (default
periods 4-6 nt) so the score reflects "peak vs. nearby noise" rather
than "peak vs. everything else." Report both columns and let the
reviewer compare.
"""

DEFAULT_BOOTSTRAP_N: int = 200
"""Bootstrap iterations for percentile CI over per-gene tracks."""

DEFAULT_PERMUTATIONS_N: int = 200
"""Circular-shift permutations for the empirical null on spectral ratios."""

DEFAULT_CI_ALPHA: float = 0.10
"""Two-sided alpha for the percentile bootstrap CI (default 90% CI)."""

DEFAULT_RANDOM_SEED: int = 42
"""Default RNG seed so bootstrap / permutation outputs are reproducible.

The seed is recorded in ``periodicity.metadata.json`` so a reviewer can
reproduce the exact CI / p-value bounds. Pass a different value via
``mitoribopy periodicity --random-seed`` to spot-check that conclusions
do not depend on this single draw.
"""

MIN_GENES_FOR_BOOTSTRAP_CI: int = 3
"""Below this many qualifying per-gene tracks we refuse to emit a CI.

A bootstrap CI from < 3 genes is misleadingly narrow (the bootstrap
cannot resample variability it has not seen). Below this threshold the
``ci_method`` column is set to ``"skipped_too_few_genes"`` and the CI
columns are NaN. The point estimate (``spectral_ratio_3nt``) is still
emitted; the user retains the ability to look at the metagene figure
and judge it visually.
"""

# Bicistronic-overlap transcripts (and their fused-FASTA spellings).
# Membership in this set excludes a transcript from the ``combined``
# gene_set and routes it to either the ``ATP86`` or ``ND4L4`` gene_set
# for separate analysis. Keys are case-insensitive; ``_`` is treated as
# ``-`` at lookup.
ATP86_GENE_SET: frozenset[str] = frozenset({
    "ATP8", "MTATP8", "MT-ATP8",
    "ATP6", "MTATP6", "MT-ATP6",
    "ATP86", "MTATP86", "MT-ATP86",
})

ND4L4_GENE_SET: frozenset[str] = frozenset({
    "ND4L", "MTND4L", "MT-ND4L",
    "ND4", "MTND4", "MT-ND4",
    "ND4L4", "MTND4L4", "MT-ND4L4",
})

OVERLAP_PAIR_TRANSCRIPTS: frozenset[str] = ATP86_GENE_SET | ND4L4_GENE_SET


Region = Literal["orf_start", "orf_stop"]
GeneSet = Literal["combined", "ATP86", "ND4L4"]
Site = Literal["a", "p", "5p"]


REGIONS: tuple[Region, Region] = ("orf_start", "orf_stop")
GENE_SETS: tuple[GeneSet, GeneSet, GeneSet] = ("combined", "ATP86", "ND4L4")


# ATP86 panel mapping: which transcript drives each panel of the ATP86
# figure. The orf_stop panel uses ATP8's window (upstream of the ATP8
# stop codon, leading INTO the bicistronic junction). The orf_start
# panel uses ATP6's window (downstream of the ATP6 start codon, leaving
# the junction). With the standard W=99 / drop_5_codons_after_start /
# drop_1_codon_before_stop defaults, the two windows overlap at the
# junction: ATP8 orf_stop = [103, 202), ATP6 orf_start = [177, 276),
# overlap = [177, 202). The user's interest in "the overlapping junction"
# is satisfied here.
_ATP86_PANEL_TRANSCRIPT: dict[str, str] = {
    "orf_stop": "ATP8",
    "orf_start": "ATP6",
}

# ND4L4 panel mapping. The 4-nt ND4L/ND4 overlap (nt 290-294) is too
# short to fall inside both windows; with W=99 the orf_stop of ND4L
# ends at 291 and orf_start of ND4 begins at 305, flanking the junction.
_ND4L4_PANEL_TRANSCRIPT: dict[str, str] = {
    "orf_stop": "ND4L",
    "orf_start": "ND4",
}


def _normalize_gene_name(name: str | None) -> str:
    if name is None:
        return ""
    return str(name).strip().upper().replace("_", "-")


def _is_in_set(name: str | None, gene_set: frozenset[str]) -> bool:
    return _normalize_gene_name(name) in gene_set


# ---------------------------------------------------------------------------
# Window definition
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FourierWindow:
    """Resolved nt-coordinate window for one (transcript, region)."""

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
    start_codon: int,
    stop_codon: int,
    transcript_length: int,
    window_nt: int,
    drop_codons_after_start: int,
    drop_codons_before_stop: int,
) -> FourierWindow:
    """Compute the half-open ``[lo, hi)`` nt window for one (tx, region)."""
    if region == "orf_start":
        skip = drop_codons_after_start * 3
        lo = start_codon + skip
        hi = lo + window_nt
        valid = (
            lo >= 0
            and hi <= transcript_length
            and hi <= stop_codon
            and (hi - lo) == window_nt
        )
    elif region == "orf_stop":
        skip = drop_codons_before_stop * 3
        hi = stop_codon - skip
        lo = hi - window_nt
        valid = (
            lo >= 0
            and lo >= start_codon + drop_codons_after_start * 3
            and hi <= transcript_length
            and (hi - lo) == window_nt
        )
    else:  # pragma: no cover - Literal guards this at the boundary
        raise ValueError(
            "region must be 'orf_start' or 'orf_stop', got " + repr(region)
        )
    return FourierWindow(
        transcript=transcript, region=region, lo=lo, hi=hi, valid=valid,
    )


# ---------------------------------------------------------------------------
# Chrom -> transcripts mapping (handles fused-FASTA aliases)
# ---------------------------------------------------------------------------


def _split_aliases(value) -> list[str]:
    if value is None:
        return []
    s = str(value).strip()
    if not s or s.lower() == "nan":
        return []
    return [token.strip() for token in s.split(";") if token.strip()]


def build_chrom_to_transcripts(annotation_df: pd.DataFrame) -> dict[str, list[str]]:
    """Map every BED chromosome name to the list of transcripts on it.

    The mapping checks ``sequence_name`` first, then each entry in the
    semicolon-separated ``sequence_aliases`` field. This is what allows
    a fused-FASTA chromosome (``ATP86``) to resolve to BOTH constituent
    transcripts (ATP8 + ATP6) for analysis.
    """
    mapping: dict[str, list[str]] = {}
    for _, row in annotation_df.iterrows():
        transcript = str(row["transcript"]).strip()
        if not transcript:
            continue
        chroms: list[str] = [str(row.get("sequence_name", "")).strip()]
        chroms.extend(_split_aliases(row.get("sequence_aliases")))
        # Always allow the transcript name itself as a chrom (some BEDs
        # name the chromosome after the transcript directly).
        chroms.append(transcript)
        for chrom in chroms:
            chrom = chrom.strip()
            if not chrom:
                continue
            mapping.setdefault(chrom, [])
            if transcript not in mapping[chrom]:
                mapping[chrom].append(transcript)
    return mapping


# ---------------------------------------------------------------------------
# Per-gene coverage extraction + normalisation
# ---------------------------------------------------------------------------


_BED_REQUIRED_COLS: tuple[str, ...] = (
    "sample_name", "chrom", "read_length", "P_site",
)
_ANNOTATION_REQUIRED_COLS: tuple[str, ...] = (
    "transcript", "sequence_name", "start_codon", "stop_codon", "l_tr",
)


def _resolve_site_offset(site: Site) -> int:
    site = str(site).lower()  # type: ignore[assignment]
    if site == "p":
        return 0
    if site == "a":
        return 3
    if site == "5p":
        return 0
    raise ValueError("site must be 'a', 'p', or '5p', got " + repr(site))


def _build_window_coverage(
    site_positions: np.ndarray, *, window: FourierWindow,
) -> np.ndarray:
    """Return a length-W coverage vector for the given window."""
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


def _normalize_and_window(
    coverage: np.ndarray, *, hann: np.ndarray,
) -> np.ndarray:
    """Per-gene preprocessing: unit-mean normalise, linear detrend, Hann window.

    Why linear detrend (instead of plain mean-centre): a 5'->3' density
    gradient or a single broad pause cluster shows up as a slow linear
    ramp inside the W-nt window. After Hann tapering, that ramp is
    concentrated in the long-period FFT bins (the Hann main lobe is
    ~4/W wide in normalised frequency; for W=99 nt that's a frequency
    band corresponding to periods >~ 25 nt). The result is the "rising
    tail" toward periods 7-10 in the raw figure. Subtracting a
    least-squares-fitted linear trend (which also implicitly removes
    the DC offset) before windowing kills that contamination at the
    source.
    """
    mean_cov = float(np.mean(coverage))
    if mean_cov <= 0:
        return np.zeros_like(coverage)
    x = coverage / mean_cov
    if x.size >= 2:
        x = _sp_signal.detrend(x, type="linear")
    else:
        x = x - float(np.mean(x))
    return x * hann


@dataclass(frozen=True)
class _PerGeneTrack:
    """One qualifying per-(sample, read_length, gene_set, region, transcript) track."""
    sample: str
    read_length: int
    gene_set: GeneSet
    region: Region
    transcript: str
    n_sites: int
    mean_coverage: float
    raw_coverage: np.ndarray
    normalized: np.ndarray  # unit-mean, mean-centred, Hann-windowed


def extract_per_gene_normalized_tracks(
    bed_with_psite_and_gene: pd.DataFrame,
    annotation_df: pd.DataFrame,
    *,
    samples: Iterable[str],
    window_nt: int = DEFAULT_WINDOW_NT,
    drop_codons_after_start: int = DEFAULT_DROP_CODONS_AFTER_START,
    drop_codons_before_stop: int = DEFAULT_DROP_CODONS_BEFORE_STOP,
    min_mean_coverage: float = DEFAULT_MIN_MEAN_COVERAGE,
    min_total_counts: int = DEFAULT_MIN_TOTAL_COUNTS,
    site: Site = "a",
) -> list[_PerGeneTrack]:
    """Extract per-(sample, length, gene_set, region, transcript) tracks.

    Each track is the preprocessed coverage vector (unit-mean normalised,
    mean-centred, Hann-windowed) ready to be aggregated into a metagene.
    Tracks failing the min-coverage / min-count filters are dropped
    silently and not returned.
    """
    missing = [c for c in _BED_REQUIRED_COLS if c not in bed_with_psite_and_gene.columns]
    if missing:
        raise ValueError(
            "bed_with_psite_and_gene is missing required column(s): " + repr(missing)
        )
    ann_missing = [c for c in _ANNOTATION_REQUIRED_COLS if c not in annotation_df.columns]
    if ann_missing:
        raise ValueError(
            "annotation_df is missing required column(s): " + repr(ann_missing)
        )

    samples_set = {str(s) for s in samples}
    if not samples_set:
        return []

    site_offset = _resolve_site_offset(site)
    chrom_to_txs = build_chrom_to_transcripts(annotation_df)
    ann = annotation_df.set_index("transcript")[
        ["start_codon", "stop_codon", "l_tr"]
    ]

    df = bed_with_psite_and_gene.copy()
    df = df[df["sample_name"].astype(str).isin(samples_set)]
    if df.empty:
        return []
    df["site_pos"] = df["P_site"].astype(int) + site_offset
    df["read_length"] = df["read_length"].astype(int)

    hann = np.hanning(window_nt)

    tracks: list[_PerGeneTrack] = []
    group_cols = ["sample_name", "chrom", "read_length"]
    for (sample, chrom, read_length), group in df.groupby(group_cols):
        chrom_str = str(chrom)
        transcripts = chrom_to_txs.get(chrom_str, [])
        if not transcripts:
            continue
        site_positions = group["site_pos"].astype(int).to_numpy()
        for transcript in transcripts:
            try:
                start_codon = int(ann.loc[transcript, "start_codon"])
                stop_codon = int(ann.loc[transcript, "stop_codon"])
                l_tr = int(ann.loc[transcript, "l_tr"])
            except (KeyError, TypeError, ValueError):
                continue

            gene_set = _classify_transcript(transcript)

            for region in REGIONS:
                # Skip regions that don't apply to this gene_set/transcript
                # combination (ATP86/ND4L4 use only one transcript per panel).
                if not _transcript_drives_region(gene_set, region, transcript):
                    continue

                window = _resolve_window(
                    transcript=transcript,
                    region=region,
                    start_codon=start_codon,
                    stop_codon=stop_codon,
                    transcript_length=l_tr,
                    window_nt=window_nt,
                    drop_codons_after_start=drop_codons_after_start,
                    drop_codons_before_stop=drop_codons_before_stop,
                )
                if not window.valid:
                    continue

                raw = _build_window_coverage(site_positions, window=window)
                n_sites = int(raw.sum())
                mean_cov = float(np.mean(raw))
                if mean_cov < float(min_mean_coverage):
                    continue
                if n_sites < int(min_total_counts):
                    continue

                normalized = _normalize_and_window(raw, hann=hann)
                tracks.append(_PerGeneTrack(
                    sample=str(sample),
                    read_length=int(read_length),
                    gene_set=gene_set,
                    region=region,
                    transcript=transcript,
                    n_sites=n_sites,
                    mean_coverage=mean_cov,
                    raw_coverage=raw,
                    normalized=normalized,
                ))
    return tracks


def _classify_transcript(transcript: str) -> GeneSet:
    """Return the gene_set this transcript belongs to."""
    if _is_in_set(transcript, ATP86_GENE_SET):
        return "ATP86"
    if _is_in_set(transcript, ND4L4_GENE_SET):
        return "ND4L4"
    return "combined"


def _transcript_drives_region(
    gene_set: GeneSet, region: Region, transcript: str,
) -> bool:
    """Whether *transcript* contributes to (gene_set, region) in the figure.

    ATP86 figure: orf_stop panel uses ATP8 only; orf_start panel uses
    ATP6 only. ND4L4 figure: orf_stop uses ND4L; orf_start uses ND4.
    Combined figure: every non-overlap-pair transcript drives both
    panels.
    """
    if gene_set == "combined":
        return True
    panel_map = (
        _ATP86_PANEL_TRANSCRIPT if gene_set == "ATP86" else _ND4L4_PANEL_TRANSCRIPT
    )
    expected = panel_map.get(region)
    if expected is None:
        return False
    return _normalize_gene_name(transcript) == _normalize_gene_name(expected)


# ---------------------------------------------------------------------------
# Metagene aggregation + DFT
# ---------------------------------------------------------------------------


def build_metagene(tracks: list[_PerGeneTrack]) -> np.ndarray:
    """Element-wise mean across qualifying per-gene normalised tracks.

    Returns a length-W vector. Returns an empty array if *tracks* is
    empty.
    """
    if not tracks:
        return np.array([], dtype=float)
    stacked = np.stack([t.normalized for t in tracks], axis=0)
    return np.mean(stacked, axis=0)


def direct_dft_amplitude(x: np.ndarray, *, period: float) -> float:
    """Direct DFT amplitude at an arbitrary period (not bin-snapped).

    Formula: ``|sum_n x[n] * exp(-2j*pi*n / period)| / sqrt(sum(x**2))``.
    Amplitude is in [0, 1]: 1.0 means x is a pure sinusoid at that
    period; 0 means no projection onto that frequency.

    Returns NaN for empty / zero-norm inputs.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    if n == 0:
        return float("nan")
    norm = float(np.sqrt(np.sum(x ** 2)))
    if not np.isfinite(norm) or norm <= 0:
        return float("nan")
    indices = np.arange(n)
    z = np.sum(x * np.exp(-2j * np.pi * indices / float(period)))
    return float(np.abs(z) / norm)


def compute_spectrum_grid(
    metagene: np.ndarray,
    *,
    periods: np.ndarray = DEFAULT_PERIOD_GRID,
) -> pd.DataFrame:
    """Direct DFT amplitudes evaluated on a dense period grid.

    Returns ``DataFrame[period_nt, amplitude]`` with one row per period
    in *periods*. Returns an empty DataFrame for empty input.
    """
    if metagene.size == 0:
        return pd.DataFrame(columns=["period_nt", "amplitude"])
    amps = np.array([direct_dft_amplitude(metagene, period=p) for p in periods])
    return pd.DataFrame({
        "period_nt": np.asarray(periods, dtype=float),
        "amplitude": amps,
    })


def compute_spectral_ratio_3nt(
    metagene: np.ndarray,
    *,
    periods: np.ndarray = DEFAULT_PERIOD_GRID,
    target_period: float = 3.0,
    exclude_window: tuple[float, float] = (2.8, 3.2),
) -> tuple[float, float]:
    """Return ``(amp_at_3nt, spectral_ratio_3nt)``.

    * ``amp_at_3nt`` is the direct DFT amplitude exactly at period 3.0.
    * ``spectral_ratio_3nt`` is ``amp_at_3nt / median(background)``
      where background is the amplitudes at all *periods* outside the
      exclusion window around 3 nt.

    Returns ``(nan, nan)`` for empty input or zero background.
    """
    if metagene.size == 0:
        return float("nan"), float("nan")
    amp3 = direct_dft_amplitude(metagene, period=float(target_period))
    grid = compute_spectrum_grid(metagene, periods=periods)
    if grid.empty:
        return amp3, float("nan")
    lo, hi = float(exclude_window[0]), float(exclude_window[1])
    background = grid[
        (grid["period_nt"] < lo) | (grid["period_nt"] > hi)
    ]["amplitude"].to_numpy()
    background = background[np.isfinite(background)]
    if background.size == 0:
        return amp3, float("nan")
    median_bg = float(np.median(background))
    if median_bg <= 0:
        return amp3, float("nan")
    return amp3, amp3 / median_bg


def compute_spectral_ratio_3nt_local(
    metagene: np.ndarray,
    *,
    periods: np.ndarray = DEFAULT_PERIOD_GRID,
    target_period: float = 3.0,
    local_range: tuple[float, float] = DEFAULT_LOCAL_BACKGROUND_RANGE,
) -> tuple[float, float]:
    """Return ``(amp_at_3nt, local_spectral_ratio_3nt)``.

    Like :func:`compute_spectral_ratio_3nt` but the background is the
    median amplitude in a NARROW codon-scale window (default periods
    4-6 nt) instead of the full 2-10 nt range. This isolates the
    period-3 peak from competing low-frequency biology (5'->3' density
    gradients, broad pause clusters) that lives at longer periods and
    can artificially inflate the global-background ratio.

    Returns ``(nan, nan)`` for empty input or zero background.
    """
    if metagene.size == 0:
        return float("nan"), float("nan")
    amp3 = direct_dft_amplitude(metagene, period=float(target_period))
    grid = compute_spectrum_grid(metagene, periods=periods)
    if grid.empty:
        return amp3, float("nan")
    lo, hi = float(local_range[0]), float(local_range[1])
    local_bg = grid[
        (grid["period_nt"] >= lo) & (grid["period_nt"] <= hi)
    ]["amplitude"].to_numpy()
    local_bg = local_bg[np.isfinite(local_bg)]
    if local_bg.size == 0:
        return amp3, float("nan")
    median_local = float(np.median(local_bg))
    if median_local <= 0:
        return amp3, float("nan")
    return amp3, amp3 / median_local


def snr_call_for_ratio(ratio: float) -> str:
    """Map a 3-nt spectral ratio to a publication-grade QC tier.

    Thresholds:
      * >= 10x  : excellent
      * >=  5x  : healthy
      * >=  2x  : modest
      * <   2x  : broken (likely frame-shifted or wrong offset)
    """
    if not np.isfinite(ratio):
        return "no_signal"
    r = float(ratio)
    if r >= 10.0:
        return "excellent"
    if r >= 5.0:
        return "healthy"
    if r >= 2.0:
        return "modest"
    return "broken"


# ---------------------------------------------------------------------------
# Vectorised DFT basis + bootstrap / permutation helpers
# ---------------------------------------------------------------------------


def build_basis_matrix(window_nt: int, periods: np.ndarray) -> np.ndarray:
    """Pre-compute the direct-DFT basis matrix ``B[k, n] = exp(-2j*pi*n/p_k)``.

    Returns a ``(P, N)`` complex array. Multiplying a length-N signal by
    ``B.T`` produces the per-period complex amplitudes in one matmul,
    avoiding the per-period Python loop in :func:`direct_dft_amplitude`.
    Sharing the matrix across many bootstrap / permutation iterations is
    what makes the statistical hardening cheap enough to run on by
    default.
    """
    n = np.arange(int(window_nt), dtype=float)
    periods_arr = np.asarray(periods, dtype=float)
    exponent = -2j * np.pi * np.outer(1.0 / periods_arr, n)
    return np.exp(exponent)


def _amplitudes_for_signals(
    signals: np.ndarray, basis: np.ndarray,
) -> np.ndarray:
    """Direct DFT amplitudes for a batch of length-N signals.

    ``signals`` shape ``(M, N)``; ``basis`` shape ``(P, N)``. Returns a
    ``(M, P)`` array of normalised amplitudes (``|sum| / ||x||_2``).
    Rows whose L2 norm is 0 emit NaN amplitudes.
    """
    z = signals @ basis.T
    norms = np.sqrt(np.sum(signals ** 2, axis=1))
    norms = np.where(norms > 0, norms, np.nan)
    return np.abs(z) / norms[:, None]


def _ratios_from_amplitudes(
    amps: np.ndarray, *,
    idx_target: int,
    idx_global_bg: np.ndarray,
    idx_local_bg: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """For each signal (row) return ``(global_ratio, local_ratio)``."""
    target = amps[:, idx_target]
    global_med = np.nanmedian(amps[:, idx_global_bg], axis=1) if idx_global_bg.size else np.full(amps.shape[0], np.nan)
    local_med = np.nanmedian(amps[:, idx_local_bg], axis=1) if idx_local_bg.size else np.full(amps.shape[0], np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        global_ratio = np.where(global_med > 0, target / global_med, np.nan)
        local_ratio = np.where(local_med > 0, target / local_med, np.nan)
    return global_ratio, local_ratio


def _percentile_ci(values: np.ndarray, alpha: float) -> tuple[float, float]:
    """Two-sided percentile CI from a finite-value sample."""
    finite = values[np.isfinite(values)]
    if finite.size < 2:
        return float("nan"), float("nan")
    half = 100.0 * float(alpha) / 2.0
    lo = float(np.percentile(finite, half))
    hi = float(np.percentile(finite, 100.0 - half))
    return lo, hi


def bootstrap_period3_ci(
    tracks: list["_PerGeneTrack"],
    *,
    periods: np.ndarray = DEFAULT_PERIOD_GRID,
    target_period: float = 3.0,
    exclude_window: tuple[float, float] = (2.8, 3.2),
    local_range: tuple[float, float] = DEFAULT_LOCAL_BACKGROUND_RANGE,
    n_bootstrap: int = DEFAULT_BOOTSTRAP_N,
    ci_alpha: float = DEFAULT_CI_ALPHA,
    rng: np.random.Generator | None = None,
    basis: np.ndarray | None = None,
) -> dict:
    """Nonparametric bootstrap-over-genes CI for amp(3) and the spectral ratios.

    Resamples the per-gene normalised tracks with replacement, recomputes
    the metagene + amp(3) + global / local ratios, and returns the
    two-sided percentile CI for each. Returns NaN CI columns when fewer
    than :data:`MIN_GENES_FOR_BOOTSTRAP_CI` genes are available — a CI
    drawn from < 3 genes is misleadingly tight.

    Notes
    -----
    The bootstrap operates on the already-preprocessed
    (unit-mean → linear-detrend → Hann-windowed) per-gene tracks, so the
    resulting CI reflects gene-to-gene variability of the metagene
    spectrum. It does NOT capture per-read sampling noise; the
    :func:`circular_shift_permutation_p` companion captures the
    significance of the *peak location* itself.
    """
    n_genes = len(tracks)
    out = {
        "amp_3nt_ci_low": float("nan"),
        "amp_3nt_ci_high": float("nan"),
        "spectral_ratio_3nt_ci_low": float("nan"),
        "spectral_ratio_3nt_ci_high": float("nan"),
        "spectral_ratio_3nt_local_ci_low": float("nan"),
        "spectral_ratio_3nt_local_ci_high": float("nan"),
        "n_bootstrap": int(n_bootstrap),
        "ci_alpha": float(ci_alpha),
        "ci_method": "percentile_over_genes",
    }
    if n_genes < MIN_GENES_FOR_BOOTSTRAP_CI or n_bootstrap < 1:
        out["ci_method"] = "skipped_too_few_genes"
        return out

    rng = rng if rng is not None else np.random.default_rng(DEFAULT_RANDOM_SEED)
    window_nt = tracks[0].normalized.size
    if window_nt == 0:
        out["ci_method"] = "skipped_empty_window"
        return out

    periods_arr = np.asarray(periods, dtype=float)
    if basis is None:
        basis = build_basis_matrix(window_nt, periods_arr)

    target_idx = int(np.argmin(np.abs(periods_arr - float(target_period))))
    lo, hi = float(exclude_window[0]), float(exclude_window[1])
    idx_global_bg = np.where((periods_arr < lo) | (periods_arr > hi))[0]
    local_lo, local_hi = float(local_range[0]), float(local_range[1])
    idx_local_bg = np.where(
        (periods_arr >= local_lo) & (periods_arr <= local_hi)
    )[0]

    stacked = np.stack([t.normalized for t in tracks], axis=0)  # (G, W)
    sample_idx = rng.integers(0, n_genes, size=(int(n_bootstrap), n_genes))
    boot_metagenes = np.mean(stacked[sample_idx], axis=1)  # (B, W)
    amps = _amplitudes_for_signals(boot_metagenes, basis)  # (B, P)
    target_amps = amps[:, target_idx]
    global_ratios, local_ratios = _ratios_from_amplitudes(
        amps, idx_target=target_idx,
        idx_global_bg=idx_global_bg, idx_local_bg=idx_local_bg,
    )

    out["amp_3nt_ci_low"], out["amp_3nt_ci_high"] = _percentile_ci(
        target_amps, ci_alpha,
    )
    out["spectral_ratio_3nt_ci_low"], out["spectral_ratio_3nt_ci_high"] = (
        _percentile_ci(global_ratios, ci_alpha)
    )
    (
        out["spectral_ratio_3nt_local_ci_low"],
        out["spectral_ratio_3nt_local_ci_high"],
    ) = _percentile_ci(local_ratios, ci_alpha)
    return out


def circular_shift_permutation_p(
    tracks: list["_PerGeneTrack"],
    *,
    observed_ratio: float,
    observed_ratio_local: float,
    periods: np.ndarray = DEFAULT_PERIOD_GRID,
    target_period: float = 3.0,
    exclude_window: tuple[float, float] = (2.8, 3.2),
    local_range: tuple[float, float] = DEFAULT_LOCAL_BACKGROUND_RANGE,
    n_permutations: int = DEFAULT_PERMUTATIONS_N,
    rng: np.random.Generator | None = None,
    basis: np.ndarray | None = None,
) -> dict:
    """Circular-shift permutation null for the spectral ratios.

    For each permutation iteration: every per-gene RAW coverage track is
    independently shifted by a uniformly random integer offset in
    ``[0, W-1]``; the per-gene unit-mean → linear-detrend → Hann-window
    pipeline is re-applied; the resulting metagene spectral ratios are
    computed. The Laplace-smoothed empirical p-value is the fraction of
    permutations whose ratio is >= the observed ratio,
    ``(k + 1) / (n_permutations + 1)``.

    The null hypothesis is "no codon-phase coherence across genes" —
    circular shifts preserve each gene's marginal coverage shape but
    randomise its position relative to the codon grid. A p-value near
    0 means the observed peak at period 3 is unlikely under that null.
    """
    n_genes = len(tracks)
    out = {
        "permutation_p": float("nan"),
        "permutation_p_local": float("nan"),
        "n_permutations": int(n_permutations),
        "null_method": "circular_shift_per_gene",
    }
    if n_genes < 1 or n_permutations < 1:
        out["null_method"] = "skipped_no_tracks"
        return out
    rng = rng if rng is not None else np.random.default_rng(
        DEFAULT_RANDOM_SEED + 1
    )
    window_nt = tracks[0].raw_coverage.size
    if window_nt == 0:
        out["null_method"] = "skipped_empty_window"
        return out

    periods_arr = np.asarray(periods, dtype=float)
    if basis is None:
        basis = build_basis_matrix(window_nt, periods_arr)
    target_idx = int(np.argmin(np.abs(periods_arr - float(target_period))))
    lo, hi = float(exclude_window[0]), float(exclude_window[1])
    idx_global_bg = np.where((periods_arr < lo) | (periods_arr > hi))[0]
    local_lo, local_hi = float(local_range[0]), float(local_range[1])
    idx_local_bg = np.where(
        (periods_arr >= local_lo) & (periods_arr <= local_hi)
    )[0]
    hann = np.hanning(window_nt)

    raw = np.stack([t.raw_coverage for t in tracks], axis=0)  # (G, W)
    # Vectorise the per-gene circular shift for one permutation iter:
    # gather indices `(arange(W) - shift) mod W` give the shifted vector.
    base_idx = np.arange(window_nt)
    n_perm = int(n_permutations)
    metagenes = np.empty((n_perm, window_nt), dtype=float)
    for perm_i in range(n_perm):
        shifts = rng.integers(0, window_nt, size=n_genes)
        gather = (base_idx[None, :] - shifts[:, None]) % window_nt
        shifted = np.take_along_axis(raw, gather, axis=1)
        # Per-gene normalize + linear detrend + Hann (matches the live
        # path in _normalize_and_window).
        means = shifted.mean(axis=1)
        keep = means > 0
        normalized = np.zeros_like(shifted)
        if np.any(keep):
            unit = shifted[keep] / means[keep, None]
            unit_dt = _sp_signal.detrend(unit, axis=1, type="linear")
            normalized[keep] = unit_dt * hann[None, :]
        metagenes[perm_i] = np.mean(normalized, axis=0)

    amps = _amplitudes_for_signals(metagenes, basis)
    null_global, null_local = _ratios_from_amplitudes(
        amps, idx_target=target_idx,
        idx_global_bg=idx_global_bg, idx_local_bg=idx_local_bg,
    )

    if np.isfinite(observed_ratio):
        n_ge = int(np.sum(np.isfinite(null_global) & (null_global >= float(observed_ratio))))
        out["permutation_p"] = (n_ge + 1) / (n_perm + 1)
    if np.isfinite(observed_ratio_local):
        n_ge = int(np.sum(np.isfinite(null_local) & (null_local >= float(observed_ratio_local))))
        out["permutation_p_local"] = (n_ge + 1) / (n_perm + 1)
    return out


# ---------------------------------------------------------------------------
# Aggregate-table builders
# ---------------------------------------------------------------------------


_SPECTRUM_COMBINED_COLS: tuple[str, ...] = (
    "sample", "read_length", "gene_set", "region",
    "n_genes", "n_sites_total", "n_nt",
    "period_nt", "amplitude",
)

_SCORE_COMBINED_COLS: tuple[str, ...] = (
    "sample", "read_length", "gene_set", "region",
    "n_genes", "n_sites_total", "n_nt",
    "amp_at_3nt", "background_amp_median", "spectral_ratio_3nt",
    "local_background_amp_median", "spectral_ratio_3nt_local",
    "snr_call", "snr_call_local",
    # Statistical hardening — bootstrap CI + permutation null.
    "amp_3nt_ci_low", "amp_3nt_ci_high",
    "spectral_ratio_3nt_ci_low", "spectral_ratio_3nt_ci_high",
    "spectral_ratio_3nt_local_ci_low", "spectral_ratio_3nt_local_ci_high",
    "permutation_p", "permutation_p_local",
    "n_bootstrap", "n_permutations", "ci_alpha",
    "ci_method", "null_method",
    "transcripts",
)


def _group_tracks(
    tracks: list[_PerGeneTrack],
) -> dict[tuple[str, int, GeneSet, Region], list[_PerGeneTrack]]:
    grouped: dict[tuple[str, int, GeneSet, Region], list[_PerGeneTrack]] = {}
    for t in tracks:
        key = (t.sample, t.read_length, t.gene_set, t.region)
        grouped.setdefault(key, []).append(t)
    return grouped


def build_fourier_spectrum_combined_table(
    tracks: list[_PerGeneTrack],
    *,
    periods: np.ndarray = DEFAULT_PERIOD_GRID,
) -> pd.DataFrame:
    """Build the per-(sample, length, gene_set, region) spectrum table.

    Columns: sample, read_length, gene_set, region, n_genes,
    n_sites_total, n_nt, period_nt, amplitude.

    One row per (sample, read_length, gene_set, region, period_nt). The
    metagene is the element-wise mean of qualifying per-gene normalised
    tracks; amplitude is the direct DFT amplitude at the period.
    """
    if not tracks:
        return pd.DataFrame(columns=_SPECTRUM_COMBINED_COLS)

    rows: list[dict] = []
    grouped = _group_tracks(tracks)
    for (sample, read_length, gene_set, region), group in grouped.items():
        metagene = build_metagene(group)
        if metagene.size == 0:
            continue
        spectrum = compute_spectrum_grid(metagene, periods=periods)
        n_genes = len(group)
        n_sites_total = int(sum(t.n_sites for t in group))
        n_nt = int(metagene.size)
        for _, srow in spectrum.iterrows():
            amp = float(srow["amplitude"])
            if not np.isfinite(amp):
                continue
            rows.append({
                "sample": sample,
                "read_length": int(read_length),
                "gene_set": str(gene_set),
                "region": str(region),
                "n_genes": int(n_genes),
                "n_sites_total": n_sites_total,
                "n_nt": n_nt,
                "period_nt": float(srow["period_nt"]),
                "amplitude": amp,
            })
    return pd.DataFrame(rows, columns=list(_SPECTRUM_COMBINED_COLS))


def build_fourier_period3_score_combined_table(
    tracks: list[_PerGeneTrack],
    *,
    periods: np.ndarray = DEFAULT_PERIOD_GRID,
    target_period: float = 3.0,
    exclude_window: tuple[float, float] = (2.8, 3.2),
    local_range: tuple[float, float] = DEFAULT_LOCAL_BACKGROUND_RANGE,
    n_bootstrap: int = DEFAULT_BOOTSTRAP_N,
    n_permutations: int = DEFAULT_PERMUTATIONS_N,
    ci_alpha: float = DEFAULT_CI_ALPHA,
    random_seed: int = DEFAULT_RANDOM_SEED,
    compute_stats: bool = True,
) -> pd.DataFrame:
    """Per-(sample, length, gene_set, region) period-3 amplitude + ratios + CI + p.

    Two complementary scalars are reported per row:

    * ``spectral_ratio_3nt`` (global) — ``amp(3) / median(amp at 2..10
      excluding 2.8..3.2)``. The classical metric. Sensitive to long-
      period structure that inflates the background median.
    * ``spectral_ratio_3nt_local`` (codon-scale neighbourhood) —
      ``amp(3) / median(amp at 4..6)``. Robust to long-period biology
      (5'->3' density gradients, broad pause clusters); reflects "peak
      vs. nearby noise" rather than "peak vs. everything else."

    Each ratio gets its own ``snr_call`` / ``snr_call_local`` tier
    ({"excellent", "healthy", "modest", "broken", "no_signal"}) per
    the thresholds in :func:`snr_call_for_ratio`.

    Statistical hardening: when ``compute_stats=True`` (the
    default) every row also carries

    * ``amp_3nt_ci_{low,high}``, ``spectral_ratio_3nt_ci_{low,high}``,
      ``spectral_ratio_3nt_local_ci_{low,high}`` — two-sided percentile
      CI from a nonparametric bootstrap-over-genes (default 200 iters,
      90 % CI). Skipped (NaN columns, ``ci_method="skipped_too_few_genes"``)
      when fewer than :data:`MIN_GENES_FOR_BOOTSTRAP_CI` genes are
      available.
    * ``permutation_p`` / ``permutation_p_local`` — Laplace-smoothed
      empirical p-values from a per-gene circular-shift permutation null
      (default 200 iters). The null hypothesis is "no codon-phase
      coherence across genes"; small p means the period-3 peak is
      unlikely under that null.
    * ``n_bootstrap`` / ``n_permutations`` / ``ci_alpha`` /
      ``ci_method`` / ``null_method`` — audit fields recording the
      exact procedure that produced the CI / p columns.

    The bootstrap and permutation steps share a single deterministic
    RNG seeded with ``random_seed`` (default 42), so re-running the same
    inputs produces byte-identical CI / p columns.

    The ``transcripts`` column is a semicolon-joined list of the
    transcripts that contributed to the metagene — auditability.
    """
    if not tracks:
        return pd.DataFrame(columns=_SCORE_COMBINED_COLS)

    rows: list[dict] = []
    grouped = _group_tracks(tracks)
    boot_rng = np.random.default_rng(int(random_seed))
    perm_rng = np.random.default_rng(int(random_seed) + 1)
    # Pre-compute a basis matrix once per window-size — every track in
    # one run shares the same window_nt, so this is a single allocation
    # that's reused across every (sample, length, gene_set, region).
    basis_cache: dict[int, np.ndarray] = {}
    periods_arr = np.asarray(periods, dtype=float)
    for (sample, read_length, gene_set, region), group in grouped.items():
        metagene = build_metagene(group)
        if metagene.size == 0:
            continue
        amp3, ratio = compute_spectral_ratio_3nt(
            metagene,
            periods=periods,
            target_period=target_period,
            exclude_window=exclude_window,
        )
        _, ratio_local = compute_spectral_ratio_3nt_local(
            metagene,
            periods=periods,
            target_period=target_period,
            local_range=local_range,
        )
        # Compute background medians exposed in output for audit.
        grid = compute_spectrum_grid(metagene, periods=periods)
        lo, hi = float(exclude_window[0]), float(exclude_window[1])
        bg = grid[(grid["period_nt"] < lo) | (grid["period_nt"] > hi)]["amplitude"]
        bg_med = float(np.median(bg.dropna())) if not bg.dropna().empty else float("nan")
        local_lo, local_hi = float(local_range[0]), float(local_range[1])
        local_bg = grid[
            (grid["period_nt"] >= local_lo) & (grid["period_nt"] <= local_hi)
        ]["amplitude"]
        local_bg_med = (
            float(np.median(local_bg.dropna()))
            if not local_bg.dropna().empty else float("nan")
        )

        # Statistical hardening — bootstrap CI + permutation null.
        if compute_stats:
            window_nt = group[0].normalized.size
            basis = basis_cache.get(window_nt)
            if basis is None:
                basis = build_basis_matrix(window_nt, periods_arr)
                basis_cache[window_nt] = basis
            stats_boot = bootstrap_period3_ci(
                group,
                periods=periods_arr,
                target_period=target_period,
                exclude_window=exclude_window,
                local_range=local_range,
                n_bootstrap=int(n_bootstrap),
                ci_alpha=float(ci_alpha),
                rng=boot_rng,
                basis=basis,
            )
            stats_perm = circular_shift_permutation_p(
                group,
                observed_ratio=float(ratio),
                observed_ratio_local=float(ratio_local),
                periods=periods_arr,
                target_period=target_period,
                exclude_window=exclude_window,
                local_range=local_range,
                n_permutations=int(n_permutations),
                rng=perm_rng,
                basis=basis,
            )
        else:
            stats_boot = {
                "amp_3nt_ci_low": float("nan"),
                "amp_3nt_ci_high": float("nan"),
                "spectral_ratio_3nt_ci_low": float("nan"),
                "spectral_ratio_3nt_ci_high": float("nan"),
                "spectral_ratio_3nt_local_ci_low": float("nan"),
                "spectral_ratio_3nt_local_ci_high": float("nan"),
                "n_bootstrap": 0,
                "ci_alpha": float(ci_alpha),
                "ci_method": "disabled",
            }
            stats_perm = {
                "permutation_p": float("nan"),
                "permutation_p_local": float("nan"),
                "n_permutations": 0,
                "null_method": "disabled",
            }

        rows.append({
            "sample": sample,
            "read_length": int(read_length),
            "gene_set": str(gene_set),
            "region": str(region),
            "n_genes": int(len(group)),
            "n_sites_total": int(sum(t.n_sites for t in group)),
            "n_nt": int(metagene.size),
            "amp_at_3nt": float(amp3),
            "background_amp_median": bg_med,
            "spectral_ratio_3nt": float(ratio),
            "local_background_amp_median": local_bg_med,
            "spectral_ratio_3nt_local": float(ratio_local),
            "snr_call": snr_call_for_ratio(ratio),
            "snr_call_local": snr_call_for_ratio(ratio_local),
            "amp_3nt_ci_low": stats_boot["amp_3nt_ci_low"],
            "amp_3nt_ci_high": stats_boot["amp_3nt_ci_high"],
            "spectral_ratio_3nt_ci_low": stats_boot["spectral_ratio_3nt_ci_low"],
            "spectral_ratio_3nt_ci_high": stats_boot["spectral_ratio_3nt_ci_high"],
            "spectral_ratio_3nt_local_ci_low": stats_boot["spectral_ratio_3nt_local_ci_low"],
            "spectral_ratio_3nt_local_ci_high": stats_boot["spectral_ratio_3nt_local_ci_high"],
            "permutation_p": stats_perm["permutation_p"],
            "permutation_p_local": stats_perm["permutation_p_local"],
            "n_bootstrap": int(stats_boot["n_bootstrap"]),
            "n_permutations": int(stats_perm["n_permutations"]),
            "ci_alpha": float(stats_boot["ci_alpha"]),
            "ci_method": str(stats_boot["ci_method"]),
            "null_method": str(stats_perm["null_method"]),
            "transcripts": ";".join(sorted(t.transcript for t in group)),
        })
    return pd.DataFrame(rows, columns=list(_SCORE_COMBINED_COLS))
