"""Tests for the UMI QC writer added in the publication-readiness branch."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from mitoribopy.align.dedup import (
    UmiQCRow,
    build_umi_qc_row,
    recommend_umi_method,
    write_umi_qc_tsv,
)


def test_recommend_unique_for_short_umi() -> None:
    method, warn = recommend_umi_method(umi_length=6, duplicate_fraction=0.5)
    assert method == "unique"
    assert warn == "short_umi_collision_risk"


def test_recommend_directional_for_long_umi_high_dup() -> None:
    method, warn = recommend_umi_method(umi_length=12, duplicate_fraction=0.4)
    assert method == "directional"
    assert warn == "none"


def test_recommend_unique_for_long_umi_low_dup() -> None:
    method, warn = recommend_umi_method(umi_length=10, duplicate_fraction=0.05)
    assert method == "unique"
    assert warn == "none"


def test_build_umi_qc_row_marks_skipped_umi() -> None:
    row = build_umi_qc_row(
        sample_id="S1",
        umi_length=8,
        umi_position="5p",
        dedup_strategy="skip",
        dedup_method="skip",
        pre_count=1000,
        post_count=1000,
    )
    assert row.umi_present is True
    assert row.warning_code == "umi_present_but_skipped"
    assert row.duplicate_fraction == 0.0


def test_build_umi_qc_row_no_umi_no_warning() -> None:
    row = build_umi_qc_row(
        sample_id="S1",
        umi_length=0,
        umi_position=None,
        dedup_strategy="skip",
        dedup_method="skip",
        pre_count=1000,
        post_count=1000,
    )
    assert row.umi_present is False
    assert row.warning_code == "no_umi"


def test_duplicate_fraction_is_complement_of_kept_ratio() -> None:
    row = build_umi_qc_row(
        sample_id="S1",
        umi_length=8,
        umi_position="5p",
        dedup_strategy="umi-tools",
        dedup_method="unique",
        pre_count=1000,
        post_count=750,
    )
    assert abs(row.duplicate_fraction - 0.25) < 1e-9


def test_write_umi_qc_tsv_round_trips(tmp_path: Path) -> None:
    rows = [
        build_umi_qc_row(
            sample_id="S1",
            umi_length=8,
            umi_position="5p",
            dedup_strategy="umi-tools",
            dedup_method="unique",
            pre_count=1000,
            post_count=900,
        ),
        build_umi_qc_row(
            sample_id="S2",
            umi_length=0,
            umi_position=None,
            dedup_strategy="skip",
            dedup_method="skip",
            pre_count=2000,
            post_count=2000,
        ),
    ]
    out = tmp_path / "umi_qc.tsv"
    write_umi_qc_tsv(rows, out)
    df = pd.read_csv(out, sep="\t", dtype={"umi_present": str})
    assert list(df.columns) == [
        "sample_id",
        "umi_present",
        "umi_length",
        "umi_position",
        "n_reads_pre_dedup",
        "n_reads_post_dedup",
        "duplicate_fraction",
        "dedup_strategy",
        "dedup_method",
        "warning_code",
    ]
    assert len(df) == 2
    s1 = df[df["sample_id"] == "S1"].iloc[0]
    assert s1["umi_present"] == "true"
    assert int(s1["n_reads_post_dedup"]) == 900
