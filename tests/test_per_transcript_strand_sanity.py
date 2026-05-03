"""Tests for the v0.9.0 per-(sample, transcript) strand-sanity audit."""

from __future__ import annotations

import pandas as pd

from mitoribopy.analysis.periodicity import (
    PerTranscriptStrandSanity,
    _write_per_transcript_strand_tsv,
    compute_per_transcript_strand_sanity,
)


def _bed(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(
        rows, columns=["sample_name", "chrom", "P_site", "read_length", "strand"],
    )


def test_per_transcript_strand_sanity_localizes_antisense_bleed() -> None:
    """A single misannotated transcript with all minus reads should
    surface as a per-transcript anomaly, even when the global
    minus-strand fraction is small."""
    rows: list[dict] = []
    # 1000 plus-strand reads on COX1 — clean.
    for pos in range(1000):
        rows.append({
            "sample_name": "WT", "chrom": "COX1",
            "P_site": pos, "read_length": 32, "strand": "+",
        })
    # 1000 plus-strand reads on ND1 — clean.
    for pos in range(1000):
        rows.append({
            "sample_name": "WT", "chrom": "ND1",
            "P_site": pos, "read_length": 32, "strand": "+",
        })
    # 50 minus-strand reads on ND6 — local antisense issue.
    for pos in range(50):
        rows.append({
            "sample_name": "WT", "chrom": "ND6",
            "P_site": pos, "read_length": 32, "strand": "-",
        })
    bed = _bed(rows)

    rows_out = compute_per_transcript_strand_sanity(bed, sample="WT")
    by_chrom = {r.transcript: r for r in rows_out}

    assert by_chrom["COX1"].minus_fraction == 0.0
    assert by_chrom["ND1"].minus_fraction == 0.0
    assert by_chrom["ND6"].minus_fraction == 1.0

    # Counts add up to total reads on each chrom.
    assert by_chrom["COX1"].n_total == 1000
    assert by_chrom["ND6"].n_minus == 50


def test_per_transcript_strand_sanity_returns_empty_for_missing_strand_column() -> None:
    bed = pd.DataFrame({
        "sample_name": ["WT"], "chrom": ["COX1"],
        "P_site": [10], "read_length": [32],
    })
    assert compute_per_transcript_strand_sanity(bed, sample="WT") == []


def test_per_transcript_strand_sanity_skips_other_samples() -> None:
    bed = _bed([
        {"sample_name": "WT", "chrom": "COX1", "P_site": 0, "read_length": 32, "strand": "+"},
        {"sample_name": "WT", "chrom": "COX1", "P_site": 1, "read_length": 32, "strand": "-"},
        {"sample_name": "KO", "chrom": "ND1", "P_site": 0, "read_length": 32, "strand": "+"},
    ])
    out = compute_per_transcript_strand_sanity(bed, sample="WT")
    assert {r.sample for r in out} == {"WT"}
    assert {r.transcript for r in out} == {"COX1"}


def test_per_transcript_strand_tsv_has_expected_schema(tmp_path) -> None:
    rows = [
        PerTranscriptStrandSanity(
            sample="WT", transcript="COX1",
            n_total=100, n_minus=2, minus_fraction=0.02,
        ),
        PerTranscriptStrandSanity(
            sample="WT", transcript="ND6",
            n_total=20, n_minus=20, minus_fraction=1.0,
        ),
    ]
    path = tmp_path / "strand_sanity_per_transcript.tsv"
    _write_per_transcript_strand_tsv(rows, path)

    df = pd.read_csv(path, sep="\t")
    assert list(df.columns) == [
        "sample", "transcript", "n_total", "n_minus", "minus_fraction",
    ]
    assert df.loc[df["transcript"] == "ND6", "minus_fraction"].iloc[0] == 1.0
    assert df.loc[df["transcript"] == "COX1", "minus_fraction"].iloc[0] == 0.02
