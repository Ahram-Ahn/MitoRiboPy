"""Tests for the v0.9.0 per-gene metagene normalization.

The contract being verified:

* The default ``normalize='per_gene_unit_mean'`` removes the depth-
  weighting bias: a single high-expression transcript no longer
  dominates the metagene shape.
* The legacy ``normalize='none'`` still reproduces the < v0.9.0
  raw-position-count semantics for users that pinned to those numbers.
* The TSV writer emits a ``# normalize: ...`` header comment + the new
  ``normalize`` and ``n_transcripts`` columns so consumers can
  unambiguously distinguish the two flavours.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mitoribopy.analysis.periodicity import (
    MetageneProfile,
    _write_metagene_tsv,
    compute_metagene,
)


def _ann(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(
        rows,
        columns=["transcript", "sequence_name", "start_codon", "stop_codon", "l_tr"],
    )


def _bed(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(
        rows, columns=["sample_name", "chrom", "read_length", "P_site"],
    )


# ---------------------------------------------------------------------------
# Per-gene unit-mean removes depth-weighting bias
# ---------------------------------------------------------------------------


def test_per_gene_unit_mean_does_not_let_one_high_expression_gene_dominate() -> None:
    """A 100x-expressed transcript dominates the legacy sum but not the
    per-gene unit-mean metagene."""
    annotation = _ann([
        {"transcript": "TX_HIGH", "sequence_name": "TX_HIGH",
         "start_codon": 0, "stop_codon": 297, "l_tr": 300},
        {"transcript": "TX_LOW", "sequence_name": "TX_LOW",
         "start_codon": 0, "stop_codon": 297, "l_tr": 300},
    ])

    # TX_HIGH: 100 reads at every codon position 0, 3, 6 ... 99 (frame 0).
    # TX_LOW: 1 read per codon position, but on frame 2 (positions 2, 5, 8 ...).
    rows: list[dict] = []
    for codon in range(0, 30):
        for _ in range(100):
            rows.append({
                "sample_name": "WT", "chrom": "TX_HIGH",
                "read_length": 32, "P_site": 3 * codon,
            })
        rows.append({
            "sample_name": "WT", "chrom": "TX_LOW",
            "read_length": 32, "P_site": 3 * codon + 2,
        })
    bed = _bed(rows)

    legacy = compute_metagene(
        bed, annotation, sample="WT", anchor="start",
        window_nt=99, normalize="none",
    )
    per_gene = compute_metagene(
        bed, annotation, sample="WT", anchor="start",
        window_nt=99, normalize="per_gene_unit_mean",
    )

    # In the legacy aggregation, frame-0 mass (TX_HIGH) is 100x larger
    # than frame-2 mass (TX_LOW), so frame_2 / total << 0.01.
    legacy_f0 = legacy.density[(legacy.positions % 3) == 0].sum()
    legacy_f2 = legacy.density[(legacy.positions % 3) == 2].sum()
    legacy_total = legacy.density.sum()
    assert legacy_f0 / legacy_total > 0.95
    assert legacy_f2 / legacy_total < 0.05

    # In the per-gene-unit-mean aggregation, both transcripts
    # contribute equal *normalized* mass, so frame-0 and frame-2 mass
    # should be roughly equal (each ≈ 50% of the total signal).
    pg_f0 = per_gene.density[(per_gene.positions % 3) == 0].sum()
    pg_f2 = per_gene.density[(per_gene.positions % 3) == 2].sum()
    pg_total = per_gene.density.sum()
    # Healthy de-biasing: each frame contributes between 30-70% of total.
    assert 0.30 < pg_f0 / pg_total < 0.70, f"frame-0 fraction = {pg_f0/pg_total:.3f}"
    assert 0.30 < pg_f2 / pg_total < 0.70, f"frame-2 fraction = {pg_f2/pg_total:.3f}"


def test_compute_metagene_records_normalize_and_n_transcripts() -> None:
    annotation = _ann([
        {"transcript": "TX1", "sequence_name": "TX1",
         "start_codon": 0, "stop_codon": 297, "l_tr": 300},
        {"transcript": "TX2", "sequence_name": "TX2",
         "start_codon": 0, "stop_codon": 297, "l_tr": 300},
    ])
    bed = _bed([
        {"sample_name": "WT", "chrom": "TX1", "read_length": 32, "P_site": 0},
        {"sample_name": "WT", "chrom": "TX2", "read_length": 32, "P_site": 3},
    ])
    p = compute_metagene(
        bed, annotation, sample="WT", anchor="start", window_nt=99,
    )
    assert p.normalize == "per_gene_unit_mean"
    assert p.n_transcripts == 2


def test_compute_metagene_rejects_unknown_normalize() -> None:
    bed = _bed([])
    annotation = _ann([])
    with pytest.raises(ValueError, match="normalize must be one of"):
        compute_metagene(
            bed, annotation, sample="WT", anchor="start",
            window_nt=99, normalize="unknown_mode",
        )


# ---------------------------------------------------------------------------
# TSV writer carries the normalization mode
# ---------------------------------------------------------------------------


def test_write_metagene_tsv_carries_normalize_header_and_columns(tmp_path: Path) -> None:
    profiles = [
        MetageneProfile(
            sample="WT",
            anchor="start",
            positions=np.array([0, 1, 2], dtype=int),
            density=np.array([0.5, 0.0, 0.5], dtype=float),
            normalize="per_gene_unit_mean",
            n_transcripts=2,
        ),
    ]
    path = tmp_path / "metagene_start.tsv"
    _write_metagene_tsv(profiles, path)

    content = path.read_text(encoding="utf-8")
    lines = content.splitlines()
    assert lines[0] == "# normalize: per_gene_unit_mean"
    assert lines[1] == "sample\tanchor\tposition\tdensity\tnormalize\tn_transcripts"
    # Three data rows.
    assert len(lines) == 5

    # pandas can read the file back when comment="#" is passed.
    df = pd.read_csv(path, sep="\t", comment="#")
    assert list(df.columns) == [
        "sample", "anchor", "position", "density", "normalize", "n_transcripts",
    ]
    assert (df["normalize"] == "per_gene_unit_mean").all()
    assert (df["n_transcripts"] == 2).all()


def test_legacy_periodicity_tsvs_carry_deprecation_header(tmp_path: Path) -> None:
    """The legacy filename writes a DEPRECATED comment so consumers see it."""
    profiles = [
        MetageneProfile(
            sample="WT", anchor="start",
            positions=np.array([0, 1, 2], dtype=int),
            density=np.array([0.1, 0.0, 0.1], dtype=float),
            normalize="per_gene_unit_mean", n_transcripts=2,
        ),
    ]
    legacy_path = tmp_path / "periodicity_start.tsv"
    _write_metagene_tsv(profiles, legacy_path)
    text = legacy_path.read_text(encoding="utf-8")
    assert text.startswith("# DEPRECATED:")
    assert "metagene_{start,stop}.tsv" in text

    # Spec-compliant filename does NOT carry the deprecation warning.
    spec_path = tmp_path / "metagene_start.tsv"
    _write_metagene_tsv(profiles, spec_path)
    spec_text = spec_path.read_text(encoding="utf-8")
    assert not spec_text.startswith("# DEPRECATED")

    # Both files round-trip through pandas with comment="#".
    df_legacy = pd.read_csv(legacy_path, sep="\t", comment="#")
    df_spec = pd.read_csv(spec_path, sep="\t", comment="#")
    pd.testing.assert_frame_equal(df_legacy, df_spec)
