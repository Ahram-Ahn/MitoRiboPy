"""Tests for per-read-length periodicity QC."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from mitoribopy.analysis.periodicity import (
    compute_frame_summary_by_length,
    compute_p_site_positions,
    run_periodicity_qc,
    select_read_lengths_by_periodicity,
)


def _annotation(transcript: str, l_utr5: int, l_cds: int, l_utr3: int) -> pd.DataFrame:
    return pd.DataFrame([{
        "transcript": transcript,
        "sequence_name": transcript,
        "display_name": transcript,
        "l_tr": l_utr5 + l_cds + l_utr3,
        "l_utr5": l_utr5,
        "l_cds": l_cds,
        "l_utr3": l_utr3,
        "start_codon": l_utr5,
        "stop_codon": l_utr5 + l_cds,
    }])


def _selected_offsets(read_length_to_offset: dict[int, int]) -> pd.DataFrame:
    rows = []
    for rl, off in read_length_to_offset.items():
        rows.append({
            "Read Length": rl,
            "Most Enriched 5' Offset": off,
            "Most Enriched 3' Offset": rl - off - 3,
        })
    return pd.DataFrame(rows)


def _make_psite_bed_with_lengths(
    sample: str,
    transcript: str,
    l_utr5: int,
    *,
    good_length_reads: list[int],
    bad_length_reads: list[int],
) -> pd.DataFrame:
    """Build a BED + offsets such that one read length is frame-0
    dominant and another is uniformly distributed across frames."""
    # Generate raw bed reads at frame-0 starts for the "good" length and
    # at random frames for the "bad" length. Offset is fixed at 12 nt
    # so P_site = start + 12 — pick start = (utr5 + frame_target * 1 - 12)
    # so the resulting P_site lands in the requested frame.
    rows = []
    rng = np.random.default_rng(0)
    for idx, target_frame in enumerate(good_length_reads):
        # P_site == l_utr5 + target_frame  => start + 12 == that
        start = l_utr5 + 3 * idx + target_frame - 12
        rows.append({
            "sample_name": sample,
            "chrom": transcript,
            "start": start,
            "end": start + 30,
            "read_length": 30,
        })
    for idx, target_frame in enumerate(bad_length_reads):
        start = l_utr5 + 3 * idx + target_frame - 12
        rows.append({
            "sample_name": sample,
            "chrom": transcript,
            "start": start,
            "end": start + 32,
            "read_length": 32,
        })
    return pd.DataFrame(rows)


def test_frame_by_length_distinguishes_good_and_bad_lengths() -> None:
    ann = _annotation("ND1", l_utr5=10, l_cds=900, l_utr3=10)
    # 30-mer reads: 100% in frame 0
    good = [0] * 240
    # 32-mer reads: uniformly spread across frames 0/1/2
    bad = ([0, 1, 2]) * 80
    bed = _make_psite_bed_with_lengths(
        "S1", "ND1", 10, good_length_reads=good, bad_length_reads=bad,
    )
    offsets = _selected_offsets({30: 12, 32: 12})
    psite = compute_p_site_positions(
        bed,
        sample="S1",
        selected_offsets_by_sample={"S1": offsets},
        selected_offsets_combined=None,
    )

    by_length = compute_frame_summary_by_length(
        psite, ann, sample="S1",
    )
    assert set(by_length["read_length"].astype(int)) >= {30, 32}
    row30 = by_length[by_length["read_length"] == 30].iloc[0]
    row32 = by_length[by_length["read_length"] == 32].iloc[0]

    assert row30["frame0_fraction"] > 0.9
    assert bool(row30["include_for_downstream"]) is True
    assert row32["frame0_fraction"] < 0.5
    assert bool(row32["include_for_downstream"]) is False
    assert row32["exclusion_reason"] in {"weak_periodicity", "ambiguous_dominance"}


def test_select_read_lengths_returns_only_included() -> None:
    df = pd.DataFrame([
        {
            "sample_id": "S1", "read_length": 28,
            "include_for_downstream": True, "exclusion_reason": "none",
        },
        {
            "sample_id": "S1", "read_length": 32,
            "include_for_downstream": False, "exclusion_reason": "low_count",
        },
        {
            "sample_id": "S2", "read_length": 30,
            "include_for_downstream": True, "exclusion_reason": "none",
        },
    ])
    chosen = select_read_lengths_by_periodicity(df)
    assert chosen == {"S1": [28], "S2": [30]}


def test_run_periodicity_qc_writes_by_length_outputs(tmp_path: Path) -> None:
    ann = _annotation("ND1", l_utr5=10, l_cds=900, l_utr3=10)
    good = [0] * 240
    bed = _make_psite_bed_with_lengths(
        "S1", "ND1", 10, good_length_reads=good, bad_length_reads=[],
    )
    offsets = _selected_offsets({30: 12})
    out_dir = tmp_path / "qc"
    run_periodicity_qc(
        bed_df=bed,
        annotation_df=ann,
        samples=["S1"],
        selected_offsets_by_sample={"S1": offsets},
        selected_offsets_combined=None,
        offset_type="5",
        offset_site="p",
        output_dir=out_dir,
        plot=False,
    )
    by_length_dir = out_dir / "by_length"
    assert (by_length_dir / "frame_by_length.tsv").is_file()
    assert (by_length_dir / "length_inclusion_decisions.tsv").is_file()
    metadata = json.loads(
        (by_length_dir / "periodicity.metadata.json").read_text(encoding="utf-8")
    )
    assert metadata["frame_formula"] == "(P_site_nt - CDS_start_nt) % 3"
    # Sanity: the table has at least one length with the included flag.
    df = pd.read_csv(by_length_dir / "frame_by_length.tsv", sep="\t")
    assert (df["include_for_downstream"].astype(bool)).any()
