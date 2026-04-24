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
