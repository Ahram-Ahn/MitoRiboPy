"""Tests for the publication-readiness codon-correlation refactor."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mitoribopy.analysis import codon_correlation as cc


def _make_codon_csv(path: Path, codons: list[tuple[str, str, str, float]]) -> None:
    """Write a minimal codon-usage CSV in the format the function expects."""
    df = pd.DataFrame(codons, columns=["Codon", "AA", "Category", "CoverageDivFreq"])
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def _build_two_sample_run(tmp_path: Path) -> Path:
    """Two samples with deliberately different codon densities so the
    log2_density_rpm metric and the MA plot have something to show."""
    base = tmp_path / "translation_profile"
    codons = [
        ("AAA", "K", "basic", 100.0),
        ("AAC", "N", "polar", 80.0),
        ("AGG", "R", "basic", 60.0),
        ("CTG", "L", "non-polar", 40.0),
        ("CTC", "L", "non-polar", 30.0),
        ("ATG", "M", "start", 20.0),
        ("TGA", "*", "stop", 5.0),
    ]
    perturbed = [
        # AAA codon dropped 4-fold (clear outlier in MA)
        ("AAA", "K", "basic", 25.0),
        ("AAC", "N", "polar", 80.0),
        ("AGG", "R", "basic", 60.0),
        ("CTG", "L", "non-polar", 40.0),
        ("CTC", "L", "non-polar", 30.0),
        ("ATG", "M", "start", 20.0),
        ("TGA", "*", "stop", 5.0),
    ]
    for sample, vals in (("WT", codons), ("KO", perturbed)):
        _make_codon_csv(
            base / sample / "codon_usage" / "p_site_codon_usage_total.csv",
            vals,
        )
        _make_codon_csv(
            base / sample / "codon_usage" / "a_site_codon_usage_total.csv",
            vals,
        )
    return base


def test_default_writes_metrics_tsv_and_metadata(tmp_path: Path) -> None:
    base = _build_two_sample_run(tmp_path)
    out_dir = tmp_path / "cor_out"
    cc.run_codon_correlation(
        translation_profile_dir=str(base),
        samples=["WT", "KO"],
        base_sample="WT",
        output_dir=str(out_dir),
        site="p",
        mask_method="none",
    )
    metrics_path = out_dir / "codon_correlation_metrics.tsv"
    assert metrics_path.is_file()
    df = pd.read_csv(metrics_path, sep="\t")
    expected_cols = {
        "site",
        "base_sample",
        "compare_sample",
        "version",
        "Codon",
        "base_metric",
        "sample_metric",
        "log2_fold_change",
        "mean_log2_density",
        "support_min_raw",
        "label_score",
        "include_primary",
        "robust_residual",
    }
    assert expected_cols.issubset(df.columns)
    # AAA was perturbed 4× downward — its log2 fold change should be
    # noticeably negative.
    aaa = df[(df["Codon"] == "AAA") & (df["compare_sample"] == "KO")]
    assert not aaa.empty
    assert (aaa["log2_fold_change"] < -0.5).all()

    metadata_path = out_dir / "codon_correlation.metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    assert metadata["metric"] == "log2_density_rpm"
    assert metadata["regression"] == "theil_sen"
    assert metadata["warnings"] == []


def test_raw_count_metric_emits_warning_and_qc_dir(tmp_path: Path) -> None:
    base = _build_two_sample_run(tmp_path)
    out_dir = tmp_path / "cor_raw"
    cc.run_codon_correlation(
        translation_profile_dir=str(base),
        samples=["WT", "KO"],
        base_sample="WT",
        output_dir=str(out_dir),
        site="p",
        mask_method="none",
        metric="raw_count",
        regression="ols",
    )
    metadata = json.loads(
        (out_dir / "codon_correlation.metadata.json").read_text(encoding="utf-8")
    )
    assert "W_CODON_RAW_COUNT_PRIMARY" in metadata["warnings"]
    # Plot should land under raw_count_qc/ when metric=raw_count.
    qc_dir = out_dir / "raw_count_qc"
    assert any(qc_dir.glob("*_vs_*_*.svg"))


def test_label_score_prefers_well_supported_outliers() -> None:
    log2_fc = pd.Series([3.0, 3.0, 0.5])
    support_min = pd.Series([100, 1, 100])
    scores = cc._label_score(log2_fc, support_min)
    # Same |log2 FC| but the high-support codon outranks the low-support one.
    assert scores.iloc[0] > scores.iloc[1]
    # The well-supported codon with smaller |FC| can still be higher than
    # a low-support 3.0; that's what we want — label the supported one.
    assert scores.iloc[2] > scores.iloc[1]


def test_invalid_metric_raises() -> None:
    with pytest.raises(ValueError, match="metric="):
        cc.run_codon_correlation(
            translation_profile_dir="/dev/null",
            samples=["a", "b"],
            base_sample="a",
            metric="not_a_metric",
        )


def test_invalid_regression_raises() -> None:
    with pytest.raises(ValueError, match="regression="):
        cc.run_codon_correlation(
            translation_profile_dir="/dev/null",
            samples=["a", "b"],
            base_sample="a",
            regression="rolls_royce",
        )
