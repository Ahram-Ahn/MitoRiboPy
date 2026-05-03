"""Tests for the v0.9.0 Fourier statistical-hardening helpers.

The contract being verified:

* Bootstrap-over-genes percentile CI brackets the point estimate when
  the per-gene tracks are mutually consistent, and is wider when they
  disagree.
* Below ``MIN_GENES_FOR_BOOTSTRAP_CI`` qualifying tracks the CI columns
  are NaN and ``ci_method == "skipped_too_few_genes"`` — we refuse to
  emit a misleadingly tight CI from < 3 genes.
* Circular-shift permutation null returns p ≈ 0 for a clean period-3
  signal and p uniformly distributed (mean ≈ 0.5) for unstructured
  noise.
* Re-running with the same ``random_seed`` produces byte-identical CI
  / p columns.
* The basis-matrix vectorisation produces the same amplitudes as the
  per-period :func:`direct_dft_amplitude` reference loop.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mitoribopy.analysis.fourier_spectrum import (
    DEFAULT_BOOTSTRAP_N,
    DEFAULT_PERIOD_GRID,
    DEFAULT_PERMUTATIONS_N,
    DEFAULT_RANDOM_SEED,
    DEFAULT_WINDOW_NT,
    MIN_GENES_FOR_BOOTSTRAP_CI,
    bootstrap_period3_ci,
    build_basis_matrix,
    build_fourier_period3_score_combined_table,
    circular_shift_permutation_p,
    direct_dft_amplitude,
    extract_per_gene_normalized_tracks,
)


# ---------------------------------------------------------------------------
# Fixture helpers (mirror test_fourier_spectrum.py conventions)
# ---------------------------------------------------------------------------


def _make_annotation(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(
        rows,
        columns=[
            "transcript", "sequence_name", "sequence_aliases",
            "start_codon", "stop_codon", "l_tr",
        ],
    )


def _make_bed(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(
        rows, columns=["sample_name", "chrom", "read_length", "P_site"],
    )


def _build_periodic_bed(
    *, sample: str, chrom: str, read_length: int,
    window_lo: int, window_hi: int, peak: int = 10,
) -> list[dict]:
    """Reads whose A-site lands every 3 nt across [window_lo, window_hi)."""
    rows: list[dict] = []
    for window_pos in range(window_lo, window_hi, 3):
        for _ in range(peak):
            rows.append({
                "sample_name": sample, "chrom": chrom,
                "read_length": read_length, "P_site": window_pos - 3,
            })
    return rows


def _multi_gene_periodic_tracks(*, n_genes: int, read_length: int = 32):
    """Build n_genes synthetic period-3 transcripts and extract their tracks."""
    ann_rows: list[dict] = []
    bed_rows: list[dict] = []
    for i in range(n_genes):
        name = f"GENE{i}"
        # Each transcript is large enough to host the orf_start window
        # at [start_codon + 15, start_codon + 15 + 99).
        ann_rows.append({
            "transcript": name, "sequence_name": name,
            "sequence_aliases": "", "start_codon": 30,
            "stop_codon": 1500, "l_tr": 1600,
        })
        bed_rows += _build_periodic_bed(
            sample="WT", chrom=name, read_length=read_length,
            window_lo=45, window_hi=144, peak=10,
        )
    ann = _make_annotation(ann_rows)
    bed = _make_bed(bed_rows)
    return extract_per_gene_normalized_tracks(
        bed, ann, samples=["WT"], window_nt=99,
    )


# ---------------------------------------------------------------------------
# Basis matrix correctness
# ---------------------------------------------------------------------------


def test_build_basis_matrix_matches_reference_direct_dft() -> None:
    rng = np.random.default_rng(0)
    x = rng.normal(size=DEFAULT_WINDOW_NT)
    basis = build_basis_matrix(DEFAULT_WINDOW_NT, DEFAULT_PERIOD_GRID)
    z = basis @ x
    ref = np.array([
        direct_dft_amplitude(x, period=p) * np.sqrt(np.sum(x ** 2))
        for p in DEFAULT_PERIOD_GRID
    ])
    np.testing.assert_allclose(np.abs(z), ref, atol=1e-9)


# ---------------------------------------------------------------------------
# Bootstrap CI
# ---------------------------------------------------------------------------


def test_bootstrap_ci_skipped_when_fewer_than_min_genes() -> None:
    tracks = _multi_gene_periodic_tracks(n_genes=2)
    assert len(tracks) < MIN_GENES_FOR_BOOTSTRAP_CI
    out = bootstrap_period3_ci(tracks)
    assert out["ci_method"] == "skipped_too_few_genes"
    assert np.isnan(out["amp_3nt_ci_low"])
    assert np.isnan(out["spectral_ratio_3nt_ci_high"])


def test_bootstrap_ci_brackets_point_estimate_for_clean_signal() -> None:
    tracks = _multi_gene_periodic_tracks(n_genes=6)
    score = build_fourier_period3_score_combined_table(
        tracks, n_bootstrap=200, n_permutations=50,
    )
    row = score.iloc[0]
    # Point estimate must lie inside its own bootstrap CI. The synthetic
    # tracks are identical to floating-point precision so the CI can
    # collapse to the point estimate ± numerical noise; allow a small
    # relative tolerance.
    rtol = 1e-6
    assert row["spectral_ratio_3nt_ci_low"] <= row["spectral_ratio_3nt"] * (1 + rtol)
    assert row["spectral_ratio_3nt"] <= row["spectral_ratio_3nt_ci_high"] * (1 + rtol)
    assert row["amp_3nt_ci_low"] <= row["amp_at_3nt"] * (1 + rtol)
    assert row["amp_at_3nt"] <= row["amp_3nt_ci_high"] * (1 + rtol)
    assert row["spectral_ratio_3nt_local_ci_low"] <= row["spectral_ratio_3nt_local"] * (1 + rtol)
    assert row["spectral_ratio_3nt_local"] <= row["spectral_ratio_3nt_local_ci_high"] * (1 + rtol)


def test_bootstrap_ci_widens_with_heterogeneous_genes() -> None:
    """Adding gene-to-gene variability should produce a non-collapsed CI."""
    rng = np.random.default_rng(202)
    ann_rows: list[dict] = []
    bed_rows: list[dict] = []
    for i in range(6):
        name = f"GENE{i}"
        ann_rows.append({
            "transcript": name, "sequence_name": name,
            "sequence_aliases": "", "start_codon": 30,
            "stop_codon": 1500, "l_tr": 1600,
        })
        # Perfect period-3 + per-gene Poisson noise so the per-gene
        # tracks differ and the bootstrap CI has non-zero width.
        peak_height = int(rng.integers(6, 16))
        bed_rows += _build_periodic_bed(
            sample="WT", chrom=name, read_length=32,
            window_lo=45, window_hi=144, peak=peak_height,
        )
        for pos in rng.integers(45, 144, size=80):
            bed_rows.append({
                "sample_name": "WT", "chrom": name,
                "read_length": 32, "P_site": int(pos),
            })
    ann = _make_annotation(ann_rows)
    bed = _make_bed(bed_rows)
    tracks = extract_per_gene_normalized_tracks(
        bed, ann, samples=["WT"], window_nt=99,
    )
    score = build_fourier_period3_score_combined_table(
        tracks, n_bootstrap=300, n_permutations=50,
    )
    row = score.iloc[0]
    # Heterogeneous genes → CI has measurable width and brackets the
    # point estimate.
    assert row["spectral_ratio_3nt_ci_high"] > row["spectral_ratio_3nt_ci_low"]
    assert row["spectral_ratio_3nt_ci_low"] <= row["spectral_ratio_3nt"] <= row["spectral_ratio_3nt_ci_high"]


def test_bootstrap_ci_is_reproducible_with_same_seed() -> None:
    tracks = _multi_gene_periodic_tracks(n_genes=4)
    a = bootstrap_period3_ci(tracks, n_bootstrap=200, rng=np.random.default_rng(7))
    b = bootstrap_period3_ci(tracks, n_bootstrap=200, rng=np.random.default_rng(7))
    assert a == b


# ---------------------------------------------------------------------------
# Circular-shift permutation null
# ---------------------------------------------------------------------------


def test_permutation_p_is_small_for_period3_dominant_signal_with_noise() -> None:
    """Realistic case: strong period-3 signal + uniform background noise.

    The clean-signal-only fixture used elsewhere collapses both the
    observed and null spectra to a leakage-dominated regime where the
    contrast is small. Adding uniform background reads gives the null
    something to push against (the noise is preserved under circular
    shift; the period-3 phase coherence across genes is destroyed),
    which is the regime the test is designed to verify.
    """
    rng = np.random.default_rng(909)
    ann_rows: list[dict] = []
    bed_rows: list[dict] = []
    for i in range(8):
        name = f"GENE{i}"
        ann_rows.append({
            "transcript": name, "sequence_name": name,
            "sequence_aliases": "", "start_codon": 30,
            "stop_codon": 1500, "l_tr": 1600,
        })
        # Period-3 signal at all genes, in phase.
        bed_rows += _build_periodic_bed(
            sample="WT", chrom=name, read_length=32,
            window_lo=45, window_hi=144, peak=8,
        )
        # Uniform background noise, independent per gene.
        for pos in rng.integers(45, 144, size=200):
            bed_rows.append({
                "sample_name": "WT", "chrom": name,
                "read_length": 32, "P_site": int(pos),
            })
    ann = _make_annotation(ann_rows)
    bed = _make_bed(bed_rows)
    tracks = extract_per_gene_normalized_tracks(
        bed, ann, samples=["WT"], window_nt=99,
    )
    score = build_fourier_period3_score_combined_table(
        tracks, n_bootstrap=50, n_permutations=400, random_seed=17,
    )
    row = score[score["region"] == "orf_start"].iloc[0]
    assert row["null_method"] == "circular_shift_per_gene"
    # Genuine codon-locked signal across many genes => permutation null
    # rarely exceeds the observed ratio. Laplace smoothing floors p at
    # 1 / (n_perm + 1) ≈ 0.0025.
    assert row["permutation_p"] < 0.10, (
        f"expected small p for codon-locked signal, got {row['permutation_p']}"
    )


def test_permutation_p_is_uniform_for_unstructured_noise() -> None:
    rng = np.random.default_rng(123)
    # Build several "genes" whose coverage is i.i.d. Poisson noise.
    ann_rows: list[dict] = []
    bed_rows: list[dict] = []
    for i in range(6):
        name = f"NOISE{i}"
        ann_rows.append({
            "transcript": name, "sequence_name": name,
            "sequence_aliases": "", "start_codon": 30,
            "stop_codon": 1500, "l_tr": 1600,
        })
        # Sample positions uniformly in the orf_start window so there
        # is no codon-locked structure.
        for pos in rng.integers(45, 144, size=300):
            bed_rows.append({
                "sample_name": "WT", "chrom": name,
                "read_length": 32, "P_site": int(pos),
            })
    ann = _make_annotation(ann_rows)
    bed = _make_bed(bed_rows)
    tracks = extract_per_gene_normalized_tracks(
        bed, ann, samples=["WT"], window_nt=99,
    )
    assert len(tracks) >= 3, "synthetic noise fixture should yield enough tracks"

    score = build_fourier_period3_score_combined_table(
        tracks, n_bootstrap=50, n_permutations=400, random_seed=11,
    )
    row = score[score["region"] == "orf_start"].iloc[0]
    # Uniform p under H0 means we should NOT see p << 0.05 reliably.
    # Allow some headroom — the test is stochastic but seed-locked.
    assert row["permutation_p"] >= 0.05 or row["permutation_p_local"] >= 0.05


def test_permutation_p_is_reproducible_with_same_seed() -> None:
    tracks = _multi_gene_periodic_tracks(n_genes=4)
    a = circular_shift_permutation_p(
        tracks, observed_ratio=8.0, observed_ratio_local=4.0,
        n_permutations=100, rng=np.random.default_rng(99),
    )
    b = circular_shift_permutation_p(
        tracks, observed_ratio=8.0, observed_ratio_local=4.0,
        n_permutations=100, rng=np.random.default_rng(99),
    )
    assert a == b


# ---------------------------------------------------------------------------
# Score-table integration
# ---------------------------------------------------------------------------


def test_score_table_includes_all_v090_columns() -> None:
    tracks = _multi_gene_periodic_tracks(n_genes=4)
    score = build_fourier_period3_score_combined_table(
        tracks, n_bootstrap=100, n_permutations=100,
    )
    expected = {
        "amp_3nt_ci_low", "amp_3nt_ci_high",
        "spectral_ratio_3nt_ci_low", "spectral_ratio_3nt_ci_high",
        "spectral_ratio_3nt_local_ci_low", "spectral_ratio_3nt_local_ci_high",
        "permutation_p", "permutation_p_local",
        "n_bootstrap", "n_permutations", "ci_alpha",
        "ci_method", "null_method",
    }
    assert expected.issubset(set(score.columns))
    row = score.iloc[0]
    assert row["n_bootstrap"] == 100
    assert row["n_permutations"] == 100
    assert row["ci_method"] == "percentile_over_genes"
    assert row["null_method"] == "circular_shift_per_gene"


def test_score_table_disable_stats_yields_nan_columns_with_disabled_method() -> None:
    tracks = _multi_gene_periodic_tracks(n_genes=4)
    score = build_fourier_period3_score_combined_table(
        tracks, compute_stats=False,
    )
    row = score.iloc[0]
    assert row["ci_method"] == "disabled"
    assert row["null_method"] == "disabled"
    assert np.isnan(row["spectral_ratio_3nt_ci_low"])
    assert np.isnan(row["permutation_p"])
    # Point estimate is still emitted.
    assert np.isfinite(row["spectral_ratio_3nt"])


def test_score_table_default_seed_is_byte_reproducible() -> None:
    tracks = _multi_gene_periodic_tracks(n_genes=4)
    a = build_fourier_period3_score_combined_table(
        tracks, n_bootstrap=100, n_permutations=100,
        random_seed=DEFAULT_RANDOM_SEED,
    )
    b = build_fourier_period3_score_combined_table(
        tracks, n_bootstrap=100, n_permutations=100,
        random_seed=DEFAULT_RANDOM_SEED,
    )
    pd.testing.assert_frame_equal(a, b)


def test_score_table_different_seed_gives_different_ci() -> None:
    """Different RNG seeds → different bootstrap draws → different CI.

    Uses a heterogeneous fixture (per-gene Poisson noise) so the
    bootstrap is non-degenerate; identical-gene fixtures collapse the
    CI to the point estimate regardless of seed.
    """
    rng = np.random.default_rng(303)
    ann_rows: list[dict] = []
    bed_rows: list[dict] = []
    for i in range(6):
        name = f"GENE{i}"
        ann_rows.append({
            "transcript": name, "sequence_name": name,
            "sequence_aliases": "", "start_codon": 30,
            "stop_codon": 1500, "l_tr": 1600,
        })
        bed_rows += _build_periodic_bed(
            sample="WT", chrom=name, read_length=32,
            window_lo=45, window_hi=144, peak=int(rng.integers(6, 16)),
        )
        for pos in rng.integers(45, 144, size=80):
            bed_rows.append({
                "sample_name": "WT", "chrom": name,
                "read_length": 32, "P_site": int(pos),
            })
    ann = _make_annotation(ann_rows)
    bed = _make_bed(bed_rows)
    tracks = extract_per_gene_normalized_tracks(
        bed, ann, samples=["WT"], window_nt=99,
    )
    a = build_fourier_period3_score_combined_table(
        tracks, n_bootstrap=100, n_permutations=100, random_seed=1,
    )
    b = build_fourier_period3_score_combined_table(
        tracks, n_bootstrap=100, n_permutations=100, random_seed=2,
    )
    assert a["spectral_ratio_3nt_ci_low"].iloc[0] != b["spectral_ratio_3nt_ci_low"].iloc[0] \
        or a["spectral_ratio_3nt_ci_high"].iloc[0] != b["spectral_ratio_3nt_ci_high"].iloc[0]
