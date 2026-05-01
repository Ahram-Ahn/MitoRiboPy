"""Lock-in tests for the v0.6.2 fourth-edit cleanup release."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import mitoribopy
from mitoribopy.align.dedup import (
    canonicalize_dedup_strategy,
    resolve_dedup_strategy,
)
from mitoribopy.analysis.periodicity_qc import (
    QC_THRESHOLDS_DEFAULT,
    build_frame_counts_by_sample_length,
    build_qc_summary,
    calculate_entropy_bias,
    calculate_fft_period3_ratio,
    calculate_frame_enrichment,
    calculate_phase_score,
)
from mitoribopy.config.migrate import migrate
from mitoribopy.io.warning_codes import WARNING_CODES, lookup
from mitoribopy.pipeline.resource_plan import plan_parallelism


# ---------------------------------------------------------------------------
# Version
# ---------------------------------------------------------------------------


def test_package_version_is_0_6_2():
    assert mitoribopy.__version__ == "0.6.2"


# ---------------------------------------------------------------------------
# Dedup strategy canonicalisation (P0-08)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw,canonical",
    [
        ("umi-tools", "umi_coordinate"),
        ("umi_tools", "umi_coordinate"),
        ("UMI-Tools", "umi_coordinate"),
        ("umi_coordinate", "umi_coordinate"),
        ("auto", "auto"),
        ("skip", "skip"),
    ],
)
def test_canonicalize_dedup_strategy_aliases(raw, canonical):
    assert canonicalize_dedup_strategy(raw) == canonical


def test_resolve_dedup_strategy_accepts_canonical_name():
    # canonical umi_coordinate with UMI -> dispatch returns "umi-tools"
    # at the implementation layer.
    assert resolve_dedup_strategy("umi_coordinate", umi_length=8) == "umi-tools"
    # auto resolves the same way via the canonicaliser.
    assert resolve_dedup_strategy("umi-tools", umi_length=8) == "umi-tools"


def test_resolve_dedup_strategy_rejects_canonical_without_umis():
    with pytest.raises(ValueError, match="umi_coordinate"):
        resolve_dedup_strategy("umi_coordinate", umi_length=0)


# ---------------------------------------------------------------------------
# Migration: legacy umi-tools + allow_pseudo_replicates rewrites
# ---------------------------------------------------------------------------


def test_migrate_rewrites_dedup_strategy():
    raw = {"align": {"dedup_strategy": "umi-tools"}}
    out, log = migrate(raw)
    assert out["align"]["dedup_strategy"] == "umi_coordinate"
    assert any("umi-tools" in line and "umi_coordinate" in line for line in log)


def test_migrate_rewrites_dedup_strategy_in_per_sample_overrides():
    raw = {
        "align": {
            "samples": [
                {"name": "s1", "dedup_strategy": "umi_tools"},
                {"name": "s2", "dedup_strategy": "skip"},
            ],
        }
    }
    out, log = migrate(raw)
    samples = out["align"]["samples"]
    assert samples[0]["dedup_strategy"] == "umi_coordinate"
    assert samples[1]["dedup_strategy"] == "skip"
    assert any("samples[0]" in line for line in log)


def test_migrate_renames_short_pseudo_replicate_key():
    raw = {"rnaseq": {"allow_pseudo_replicates": True}}
    out, log = migrate(raw)
    assert "allow_pseudo_replicates_for_demo_not_publication" in out["rnaseq"]
    assert out["rnaseq"]["allow_pseudo_replicates_for_demo_not_publication"] is True
    assert "allow_pseudo_replicates" not in out["rnaseq"]


# ---------------------------------------------------------------------------
# Top-level execution: cascade + resource_plan.json (P0-05/06)
# ---------------------------------------------------------------------------


def test_apply_execution_block_cascades_to_align():
    from mitoribopy.cli.all_ import _apply_execution_block

    cfg = {
        "execution": {"threads": 16, "parallel_samples": "auto", "memory_gb": 8.0},
        "align": {},
    }
    _apply_execution_block(cfg, cli_threads=None)
    assert cfg["align"]["threads"] == 16
    assert cfg["align"]["max_parallel_samples"] == "auto"
    assert cfg["align"]["memory_gb"] == 8.0


def test_apply_execution_block_cli_threads_seeds_execution():
    from mitoribopy.cli.all_ import _apply_execution_block

    cfg = {"align": {}}
    _apply_execution_block(cfg, cli_threads=20)
    assert cfg["execution"]["threads"] == 20
    assert cfg["align"]["threads"] == 20


def test_apply_execution_block_explicit_stage_value_wins():
    from mitoribopy.cli.all_ import _apply_execution_block

    cfg = {
        "execution": {"threads": 16},
        "align": {"threads": 4},
    }
    _apply_execution_block(cfg, cli_threads=None)
    assert cfg["align"]["threads"] == 4


def test_resource_plan_round_trip(tmp_path: Path):
    from mitoribopy.pipeline.resource_plan import write_resource_plan

    plan = plan_parallelism(
        n_samples=8, requested_threads=20, requested_parallel="auto",
        memory_gb=None,
    )
    path = write_resource_plan(plan, tmp_path)
    assert path.exists()
    import json
    data = json.loads(path.read_text())
    assert data["n_samples"] == 8
    assert data["total_threads"] == 20
    assert data["parallel_samples"] >= 1


# ---------------------------------------------------------------------------
# Periodicity QC (P spec)
# ---------------------------------------------------------------------------


def test_entropy_bias_is_zero_for_uniform_and_one_for_dominated():
    assert calculate_entropy_bias((1.0, 0.0, 0.0)) == pytest.approx(1.0)
    assert calculate_entropy_bias((1 / 3, 1 / 3, 1 / 3)) == pytest.approx(0.0, abs=1e-6)


def test_frame_enrichment_picks_expected_frame():
    # All in frame 0 -> infinite enrichment (other frames are zero).
    assert calculate_frame_enrichment((1.0, 0.0, 0.0)) == float("inf")
    # Uniform -> enrichment of 1.
    assert calculate_frame_enrichment((1 / 3, 1 / 3, 1 / 3)) == pytest.approx(1.0)
    # Frame 1 dominant when expected_frame=1.
    val = calculate_frame_enrichment((0.1, 0.8, 0.1), expected_frame=1)
    assert val > 1.0


def test_phase_score_consistent_dominance():
    # Every codon dominates frame 0 -> score == 1.
    codons = [(10, 1, 1)] * 20
    score = calculate_phase_score(codons)
    assert score == pytest.approx(1.0, abs=1e-6)


def test_phase_score_uniform_low():
    codons = [(1, 1, 1)] * 20
    score = calculate_phase_score(codons)
    # Uniform per codon yields a zero vector that we skip -> NaN.
    assert score != score or score < 0.5  # NaN fall-through allowed


def test_fft_period3_ratio_high_for_periodic_signal():
    import numpy as np

    coverage = np.tile([10.0, 1.0, 1.0], 50)
    ratio = calculate_fft_period3_ratio(coverage)
    assert ratio > 5.0


def test_fft_period3_ratio_nan_for_short_signal():
    import math

    val = calculate_fft_period3_ratio([1.0, 2.0, 3.0])
    assert math.isnan(val)


def test_build_frame_counts_assigns_qc_call():
    df = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "read_length": 30,
                "n_reads_total": 5000,
                "n_reads_cds": 5000,
                "frame0_fraction": 0.7, "frame1_fraction": 0.2, "frame2_fraction": 0.1,
                "dominant_frame": 0, "frame0_dominance": 0.5,
                "periodicity_score": 0.5, "frame_entropy": 0.9,
                "include_for_downstream": True, "exclusion_reason": "none",
            },
            {
                "sample_id": "S1",
                "read_length": 31,
                "n_reads_total": 100,
                "n_reads_cds": 100,
                "frame0_fraction": 0.5, "frame1_fraction": 0.3, "frame2_fraction": 0.2,
                "dominant_frame": 0, "frame0_dominance": 0.2,
                "periodicity_score": 0.2, "frame_entropy": 1.4,
                "include_for_downstream": False, "exclusion_reason": "low_count",
            },
        ]
    )
    annotated = build_frame_counts_by_sample_length(df)
    calls = sorted(annotated["qc_call"].tolist())
    assert "good" in calls
    assert "low_depth" in calls


def test_build_qc_summary_picks_best_length():
    df = pd.DataFrame(
        [
            {
                "sample_id": "S1", "read_length": 30,
                "n_reads_total": 5000, "n_reads_cds": 5000,
                "frame0_fraction": 0.75, "frame1_fraction": 0.15, "frame2_fraction": 0.10,
                "dominant_frame": 0, "frame0_dominance": 0.6,
                "periodicity_score": 0.6, "frame_entropy": 0.8,
                "include_for_downstream": True, "exclusion_reason": "none",
                "expected_frame": 0, "expected_frame_fraction": 0.75,
                "expected_frame_enrichment": 6.0,
                "entropy_bias": 0.4, "qc_call": "good",
            },
        ]
    )
    summary = build_qc_summary(df)
    assert len(summary) == 1
    assert summary.iloc[0]["overall_qc_call"] == "good"
    assert summary.iloc[0]["best_read_length"] == 30
    assert summary.iloc[0]["n_total_sites"] == 5000


def test_qc_thresholds_default_present():
    for k in (
        "good_frame_fraction", "warn_frame_fraction",
        "min_reads_per_length", "min_reads_per_gene", "min_phase_score_good",
    ):
        assert k in QC_THRESHOLDS_DEFAULT


# ---------------------------------------------------------------------------
# Warning code registry (P1-06)
# ---------------------------------------------------------------------------


def test_warning_code_registry_has_required_codes():
    expected = {
        "E_CONFIG_UNKNOWN_KEY",
        "W_UMI_SHORT_COLLISION_RISK",
        "W_UMI_HIGH_DUPLICATE_FRACTION",
        "W_UMI_DEDUP_SKIPPED_WITH_UMI_PRESENT",
        "W_CODON_RAW_COUNT_PRIMARY",
        "W_RNASEQ_FROM_FASTQ_EXPLORATORY",
        "W_RNASEQ_PSEUDO_REPLICATE",
        "E_REFERENCE_CHECKSUM_MISMATCH",
    }
    have = {entry.code for entry in WARNING_CODES}
    missing = expected - have
    assert not missing, f"warning code registry is missing: {missing}"


def test_warning_code_lookup_returns_entry():
    entry = lookup("W_UMI_SHORT_COLLISION_RISK")
    assert entry is not None
    assert entry.severity == "warn"
    assert entry.stage == "align"
    assert entry.sample_scoped is True


def test_warning_code_lookup_unknown_returns_none():
    assert lookup("Z_NOT_REAL") is None
