"""Per-sample offset selection regression tests."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from mitoribopy.analysis.offset_enrichment import (
    build_per_sample_summaries,
    compute_offsets,
    create_csv_for_offset_enrichment,
)
from mitoribopy.analysis.offset_selection import determine_p_site_offsets
from mitoribopy.data.transcript_annotations import human_annotation_df


def _human_annotation_with_codon_positions() -> pd.DataFrame:
    annotation = human_annotation_df.copy()
    annotation["start_codon"] = annotation["l_utr5"]
    annotation["stop_codon"] = annotation["l_tr"] - annotation["l_utr3"] - 3
    return annotation


def _bed_with_two_samples_different_offsets() -> pd.DataFrame:
    """Build a BED-like frame where two samples bracket the ND1 stop codon
    at slightly different read-start positions, producing distinct
    per-sample 5' offsets at RL30."""
    rows: list[dict] = []
    # Sample S1: offset peak at one position
    rows.extend(
        {"chrom": "ND1", "start": 940, "end": 970, "read_length": 30, "sample_name": "S1"}
        for _ in range(40)
    )
    # Sample S2: offset peak shifted by 2 nt (different starts)
    rows.extend(
        {"chrom": "ND1", "start": 942, "end": 972, "read_length": 30, "sample_name": "S2"}
        for _ in range(40)
    )
    return pd.DataFrame(rows)


def test_compute_offsets_carries_sample_name(tmp_path: Path) -> None:
    bed_df = _bed_with_two_samples_different_offsets()
    annotation = _human_annotation_with_codon_positions()
    offsets = compute_offsets(
        bed_df=bed_df,
        annotation_df=annotation,
        align_to="stop",
        offset_site="p",
    )
    assert "sample_name" in offsets.columns
    samples = set(offsets["sample_name"].unique())
    assert samples == {"S1", "S2"}


def test_per_sample_summaries_are_independent(tmp_path: Path) -> None:
    bed_df = _bed_with_two_samples_different_offsets()
    annotation = _human_annotation_with_codon_positions()
    _, offsets = create_csv_for_offset_enrichment(
        bed_df=bed_df,
        annotation_df=annotation,
        align_to="stop",
        rpf_range=range(30, 31),
        output_csv=str(tmp_path / "combined.csv"),
        offset_limit=25,
        offset_site="p",
        codon_overlap_mode="full",
        strain="h",
    )
    assert offsets is not None
    summaries = build_per_sample_summaries(
        offsets, range(30, 31), offset_limit=25, offset_mask_nt=5
    )
    assert set(summaries) == {"S1", "S2"}
    # The two summaries should not be element-wise identical because the
    # samples bracket the stop codon at distinct read-start positions.
    s1 = summaries["S1"].drop(columns=["Read Length"]).iloc[0].to_numpy()
    s2 = summaries["S2"].drop(columns=["Read Length"]).iloc[0].to_numpy()
    assert (s1 != s2).any(), "expected per-sample summaries to differ"


def test_determine_p_site_offsets_per_sample_picks_independent_values(tmp_path: Path) -> None:
    bed_df = _bed_with_two_samples_different_offsets()
    annotation = _human_annotation_with_codon_positions()
    _, offsets = create_csv_for_offset_enrichment(
        bed_df=bed_df,
        annotation_df=annotation,
        align_to="stop",
        rpf_range=range(30, 31),
        output_csv=str(tmp_path / "combined.csv"),
        offset_limit=25,
        offset_site="p",
        codon_overlap_mode="full",
        strain="h",
    )
    assert offsets is not None

    per_sample_picks: dict[str, pd.DataFrame] = {}
    for sample, sub in offsets.groupby("sample_name"):
        out = tmp_path / f"sel_{sample}.csv"
        picked = determine_p_site_offsets(
            offsets_df=sub,
            align_to="stop",
            out_file=str(out),
            offset_min=10,
            offset_max=22,
            offset_site="p",
            selection_reference="p_site",
        )
        assert picked is not None
        per_sample_picks[sample] = picked

    s1_5 = per_sample_picks["S1"]["Most Enriched 5' Offset"].iloc[0]
    s2_5 = per_sample_picks["S2"]["Most Enriched 5' Offset"].iloc[0]
    # Sample S2's reads start 2 nt later, so its 5' P-site offset should
    # land 2 nt lower than S1's.
    assert s1_5 != s2_5
    assert s1_5 - s2_5 == 2


# ---------- Phase 2.1: offset confidence scores ----------------------------


def test_offset_confidence_columns_present_in_selected_table(tmp_path: Path) -> None:
    """The selected-offsets CSV gains diagnostic columns (n_reads_5,
    top_count_5, second_best_*, delta_score_5, enrichment_score_5,
    confidence_5; same for the 3' side) so reviewers can audit each
    pick instead of trusting it on faith."""
    bed_df = _bed_with_two_samples_different_offsets()
    annotation = _human_annotation_with_codon_positions()
    _, offsets = create_csv_for_offset_enrichment(
        bed_df=bed_df,
        annotation_df=annotation,
        align_to="stop",
        rpf_range=range(30, 31),
        output_csv=str(tmp_path / "combined.csv"),
        offset_limit=25,
        offset_site="p",
        codon_overlap_mode="full",
        strain="h",
    )
    out = tmp_path / "selected.csv"
    picked = determine_p_site_offsets(
        offsets_df=offsets,
        align_to="stop",
        out_file=str(out),
        offset_min=10,
        offset_max=22,
        offset_site="p",
        selection_reference="p_site",
    )
    assert picked is not None

    expected_cols = {
        "Read Length", "Most Enriched 5' Offset", "Most Enriched 3' Offset",
        "n_reads_5", "top_count_5", "second_best_count_5",
        "second_best_offset_5", "delta_score_5", "enrichment_score_5",
        "confidence_5",
        "n_reads_3", "top_count_3", "second_best_count_3",
        "second_best_offset_3", "delta_score_3", "enrichment_score_3",
        "confidence_3",
    }
    missing = expected_cols - set(picked.columns)
    assert not missing, f"expected confidence columns missing: {missing}"
    # The CSV on disk carries the same columns.
    on_disk = pd.read_csv(out)
    assert expected_cols.issubset(set(on_disk.columns))


def test_confidence_label_thresholds_classify_as_documented() -> None:
    """The threshold table in offset_selection._confidence_label is the
    package's defensibility argument; this test pins the documented
    rules so silent threshold drift is caught in CI."""
    from mitoribopy.analysis.offset_selection import _confidence_label

    # Strong peak with plenty of reads -> high.
    assert _confidence_label(
        n_reads=1000, top_count=600, second_best_count=200,
        enrichment_score=0.6,
    ) == "high"
    # Decent enrichment but only borderline read count -> medium.
    assert _confidence_label(
        n_reads=80, top_count=30, second_best_count=20,
        enrichment_score=0.375,
    ) == "medium"
    # Strong fraction but barely any reads -> low (n_reads gate).
    assert _confidence_label(
        n_reads=10, top_count=8, second_best_count=1,
        enrichment_score=0.8,
    ) == "low"
    # Empty length window -> insufficient.
    assert _confidence_label(
        n_reads=0, top_count=0, second_best_count=0, enrichment_score=0.0,
    ) == "insufficient"
    # High-tier delta_ratio fails because runner-up is too close.
    assert _confidence_label(
        n_reads=1000, top_count=400, second_best_count=350,
        enrichment_score=0.4,
    ) == "medium"  # enrichment + reads enough for medium, but delta_ratio kills high


def test_high_confidence_emitted_for_clean_signal(tmp_path: Path) -> None:
    """A synthetic offsets DataFrame with one strong, well-separated
    peak should yield confidence_5 == 'high'."""
    # Build a fake offsets table for a single read length: 1000 reads
    # cluster at offset 12, with thin tails at 11 and 13. That's a
    # textbook "high confidence" situation for the human monosome.
    rows = []
    for _ in range(1000):
        rows.append({"Read Length": 30, "5' Offset": 12, "3' Offset": -18})
    for _ in range(50):
        rows.append({"Read Length": 30, "5' Offset": 11, "3' Offset": -19})
    for _ in range(50):
        rows.append({"Read Length": 30, "5' Offset": 13, "3' Offset": -17})
    offsets = pd.DataFrame(rows)
    out = tmp_path / "clean.csv"
    picked = determine_p_site_offsets(
        offsets_df=offsets,
        align_to="stop",
        out_file=str(out),
        offset_min=10,
        offset_max=22,
        offset_site="p",
        selection_reference="p_site",
    )
    assert picked is not None
    row = picked.iloc[0]
    assert row["Most Enriched 5' Offset"] == 12
    assert row["confidence_5"] == "high"
    assert row["n_reads_5"] == 1100
    assert row["top_count_5"] == 1000
    assert row["second_best_offset_5"] in (11, 13)
    assert row["delta_score_5"] == 950
