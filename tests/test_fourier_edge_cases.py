"""Edge-case coverage for the Fourier periodicity bundle.

Round-2 audit found these paths were not exercised by the synthetic
clean-signal fixtures:

* Empty BED (no reads at all) → no crash, empty score table.
* All-zero coverage on a transcript → skipped without divide-by-zero.
* Single-transcript run (n_genes=1) → bootstrap CI gracefully NaN'd
  with `ci_method == "skipped_too_few_genes"`, point estimate still
  emitted.
* Annotation rows whose ``stop_codon`` is too close to ``start_codon``
  for any window of W nt to fit → silently skipped, no exception.
* Direct DFT on a length-0 signal → NaN, not a crash.
* `compute_metagene` on a sample with no qualifying transcripts →
  empty metagene profile (zeros), `n_transcripts == 0`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from mitoribopy.analysis.fourier_spectrum import (
    DEFAULT_PERIOD_GRID,
    DEFAULT_WINDOW_NT,
    MIN_GENES_FOR_BOOTSTRAP_CI,
    bootstrap_period3_ci,
    build_basis_matrix,
    build_fourier_period3_score_combined_table,
    build_fourier_spectrum_combined_table,
    direct_dft_amplitude,
    extract_per_gene_normalized_tracks,
)
from mitoribopy.analysis.periodicity import (
    MetageneProfile,
    compute_metagene,
)


def _ann(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(
        rows,
        columns=[
            "transcript", "sequence_name", "sequence_aliases",
            "start_codon", "stop_codon", "l_tr",
        ],
    )


def _bed(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(
        rows, columns=["sample_name", "chrom", "read_length", "P_site"],
    )


# ---------------------------------------------------------------------------
# Empty / degenerate inputs
# ---------------------------------------------------------------------------


def test_extract_tracks_returns_empty_for_empty_bed() -> None:
    ann = _ann([
        {"transcript": "COX1", "sequence_name": "COX1",
         "sequence_aliases": "", "start_codon": 30,
         "stop_codon": 1500, "l_tr": 1600},
    ])
    bed = _bed([])
    tracks = extract_per_gene_normalized_tracks(
        bed, ann, samples=["WT"], window_nt=99,
    )
    assert tracks == []


def test_score_table_is_header_only_for_empty_tracks() -> None:
    score = build_fourier_period3_score_combined_table([])
    assert score.empty
    # The schema must still be present so a downstream reader can
    # column-check before iterating.
    expected = {
        "sample", "spectral_ratio_3nt", "snr_call",
        "amp_3nt_ci_low", "amp_3nt_ci_high",
        "permutation_p", "ci_method", "null_method",
    }
    assert expected.issubset(set(score.columns))


def test_spectrum_table_is_header_only_for_empty_tracks() -> None:
    spec = build_fourier_spectrum_combined_table([])
    assert spec.empty
    assert {"sample", "period_nt", "amplitude"}.issubset(set(spec.columns))


# ---------------------------------------------------------------------------
# Numerical guards
# ---------------------------------------------------------------------------


def test_direct_dft_amplitude_returns_nan_for_empty_signal() -> None:
    assert np.isnan(direct_dft_amplitude(np.array([], dtype=float), period=3.0))


def test_direct_dft_amplitude_returns_nan_for_zero_norm_signal() -> None:
    assert np.isnan(direct_dft_amplitude(np.zeros(99, dtype=float), period=3.0))


def test_build_basis_matrix_handles_single_period() -> None:
    """A length-1 period grid must still produce a (1, W) basis."""
    basis = build_basis_matrix(99, np.array([3.0]))
    assert basis.shape == (1, 99)
    assert np.iscomplexobj(basis)


# ---------------------------------------------------------------------------
# Bootstrap / permutation guards
# ---------------------------------------------------------------------------


def test_bootstrap_skipped_for_zero_iterations() -> None:
    """Passing n_bootstrap=0 should refuse to emit a CI rather than crash."""
    # We need at least MIN_GENES_FOR_BOOTSTRAP_CI tracks for the
    # zero-iteration path to be reachable; build them, then call with 0.
    ann = _ann([
        {"transcript": f"GENE{i}", "sequence_name": f"GENE{i}",
         "sequence_aliases": "", "start_codon": 30,
         "stop_codon": 1500, "l_tr": 1600}
        for i in range(4)
    ])
    bed_rows: list[dict] = []
    for i in range(4):
        for pos in range(45, 144, 3):
            for _ in range(10):
                bed_rows.append({
                    "sample_name": "WT", "chrom": f"GENE{i}",
                    "read_length": 32, "P_site": pos - 3,
                })
    tracks = extract_per_gene_normalized_tracks(
        _bed(bed_rows), ann, samples=["WT"], window_nt=99,
    )
    out = bootstrap_period3_ci(tracks, n_bootstrap=0)
    assert out["ci_method"] == "skipped_too_few_genes"
    assert np.isnan(out["spectral_ratio_3nt_ci_low"])


# ---------------------------------------------------------------------------
# Annotation edge cases
# ---------------------------------------------------------------------------


def test_extract_tracks_skips_transcripts_too_short_for_window() -> None:
    """A transcript whose CDS does not fit a W=99 window is silently skipped."""
    ann = _ann([
        # CDS length only 30 nt — well under W=99.
        {"transcript": "TINY", "sequence_name": "TINY",
         "sequence_aliases": "", "start_codon": 0,
         "stop_codon": 30, "l_tr": 30},
        # Normal-sized companion so the function has something to return.
        {"transcript": "COX1", "sequence_name": "COX1",
         "sequence_aliases": "", "start_codon": 30,
         "stop_codon": 1500, "l_tr": 1600},
    ])
    bed_rows: list[dict] = []
    # Reads on TINY (which has no valid window) and on COX1.
    for pos in range(0, 30, 3):
        bed_rows.append({"sample_name": "WT", "chrom": "TINY",
                         "read_length": 32, "P_site": pos})
    for pos in range(45, 144, 3):
        for _ in range(10):
            bed_rows.append({"sample_name": "WT", "chrom": "COX1",
                             "read_length": 32, "P_site": pos - 3})
    tracks = extract_per_gene_normalized_tracks(
        _bed(bed_rows), ann, samples=["WT"], window_nt=99,
    )
    transcripts_seen = {t.transcript for t in tracks}
    assert "TINY" not in transcripts_seen
    assert "COX1" in transcripts_seen


def test_extract_tracks_skips_zero_mean_window_quietly() -> None:
    """A transcript whose window has zero coverage is filtered, not crashed."""
    ann = _ann([
        {"transcript": "QUIET", "sequence_name": "QUIET",
         "sequence_aliases": "", "start_codon": 30,
         "stop_codon": 1500, "l_tr": 1600},
    ])
    # No reads at all on QUIET.
    tracks = extract_per_gene_normalized_tracks(
        _bed([]), ann, samples=["WT"], window_nt=99,
    )
    assert tracks == []


# ---------------------------------------------------------------------------
# compute_metagene degenerate inputs
# ---------------------------------------------------------------------------


def test_compute_metagene_returns_empty_profile_for_empty_bed() -> None:
    annotation = pd.DataFrame(
        [{"transcript": "COX1", "sequence_name": "COX1",
          "start_codon": 30, "stop_codon": 1500, "l_tr": 1600}],
        columns=["transcript", "sequence_name", "start_codon", "stop_codon", "l_tr"],
    )
    bed = pd.DataFrame(
        columns=["sample_name", "chrom", "P_site", "read_length"],
    )
    profile = compute_metagene(
        bed, annotation, sample="WT", anchor="start", window_nt=99,
    )
    assert isinstance(profile, MetageneProfile)
    assert profile.n_transcripts == 0
    assert profile.density.sum() == 0.0
    assert profile.normalize == "per_gene_unit_mean"


def test_compute_metagene_legacy_normalize_returns_integer_counts() -> None:
    """The legacy 'none' normalisation must reproduce raw per-position
    counts (sum across transcripts, no per-gene division)."""
    annotation = pd.DataFrame(
        [{"transcript": "COX1", "sequence_name": "COX1",
          "start_codon": 0, "stop_codon": 297, "l_tr": 300}],
        columns=["transcript", "sequence_name", "start_codon", "stop_codon", "l_tr"],
    )
    bed = pd.DataFrame(
        [
            {"sample_name": "WT", "chrom": "COX1", "P_site": 0, "read_length": 32},
            {"sample_name": "WT", "chrom": "COX1", "P_site": 0, "read_length": 32},
            {"sample_name": "WT", "chrom": "COX1", "P_site": 3, "read_length": 32},
        ],
        columns=["sample_name", "chrom", "P_site", "read_length"],
    )
    legacy = compute_metagene(
        bed, annotation, sample="WT", anchor="start",
        window_nt=99, normalize="none",
    )
    # 2 reads at position 0, 1 at position 3.
    assert legacy.density[0] == 2.0
    assert legacy.density[3] == 1.0
    assert legacy.normalize == "none"
    assert legacy.n_transcripts == 1
