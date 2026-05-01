"""Phase 2.2 — periodicity + frame + strand QC."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mitoribopy.analysis.periodicity import (
    DEFAULT_WINDOW_NT,
    compute_frame_summary,
    compute_metagene,
    compute_p_site_positions,
    compute_strand_sanity,
    run_periodicity_qc,
)


def _annotation(transcript: str, l_utr5: int, l_cds: int, l_utr3: int) -> pd.DataFrame:
    """Single-row annotation matching the package's normalised schema."""
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


def _selected_offsets(offset_5: int, read_length: int = 30) -> pd.DataFrame:
    """One row of the selected-offsets table for a single read length."""
    return pd.DataFrame([{
        "Read Length": read_length,
        "Most Enriched 5' Offset": offset_5,
        "Most Enriched 3' Offset": read_length - offset_5 - 3,
    }])


# ---------- compute_p_site_positions ----------------------------------------


def test_compute_p_site_uses_per_sample_offset_for_5prime_p_site() -> None:
    """5' offset of 12 → P_site = start + 12 in P-site placement mode."""
    bed = pd.DataFrame([
        {"sample_name": "S1", "chrom": "ND1", "start": 100, "end": 130, "read_length": 30},
        {"sample_name": "S1", "chrom": "ND1", "start": 200, "end": 230, "read_length": 30},
        # Different read length — must be dropped (no offset for it).
        {"sample_name": "S1", "chrom": "ND1", "start": 300, "end": 332, "read_length": 32},
        # Other sample — must be excluded entirely from this sample's view.
        {"sample_name": "S2", "chrom": "ND1", "start": 100, "end": 130, "read_length": 30},
    ])
    psite = compute_p_site_positions(
        bed,
        sample="S1",
        selected_offsets_by_sample={"S1": _selected_offsets(12)},
        selected_offsets_combined=None,
        offset_type="5",
        offset_site="p",
    )
    assert list(psite["P_site"]) == [112, 212]


def test_compute_p_site_falls_back_to_combined_when_no_per_sample() -> None:
    bed = pd.DataFrame([
        {"sample_name": "S1", "chrom": "ND1", "start": 50, "end": 80, "read_length": 30},
    ])
    psite = compute_p_site_positions(
        bed, sample="S1",
        selected_offsets_by_sample=None,
        selected_offsets_combined=_selected_offsets(11),
        offset_type="5",
        offset_site="p",
    )
    assert list(psite["P_site"]) == [61]


def test_compute_p_site_returns_empty_when_no_offset_table() -> None:
    bed = pd.DataFrame([
        {"sample_name": "S1", "chrom": "ND1", "start": 10, "end": 40, "read_length": 30},
    ])
    psite = compute_p_site_positions(
        bed, sample="S1",
        selected_offsets_by_sample=None,
        selected_offsets_combined=None,
        offset_type="5",
        offset_site="p",
    )
    assert psite.empty


# ---------- compute_frame_summary ------------------------------------------


def test_frame_summary_perfect_frame0_signal() -> None:
    """Reads whose P-sites all sit on codon-1 starts of a CDS that
    begins at position 0 must yield frame_0 == 1.0."""
    ann = _annotation("ND1", l_utr5=0, l_cds=300, l_utr3=10)
    psite = pd.DataFrame([
        {"sample_name": "S1", "chrom": "ND1", "P_site": 0, "read_length": 30, "start": 0, "end": 0},
        {"sample_name": "S1", "chrom": "ND1", "P_site": 3, "read_length": 30, "start": 0, "end": 0},
        {"sample_name": "S1", "chrom": "ND1", "P_site": 6, "read_length": 30, "start": 0, "end": 0},
        {"sample_name": "S1", "chrom": "ND1", "P_site": 9, "read_length": 30, "start": 0, "end": 0},
    ])
    fs = compute_frame_summary(psite, ann, sample="S1")
    assert fs.n_reads == 4
    assert fs.frame_0 == pytest.approx(1.0)
    assert fs.dominance == pytest.approx(1.0)


def test_frame_summary_balanced_signal_low_dominance() -> None:
    """3 reads in each frame -> dominance == 0."""
    ann = _annotation("ND1", l_utr5=0, l_cds=300, l_utr3=10)
    psite = pd.DataFrame([
        {"sample_name": "S1", "chrom": "ND1", "P_site": p, "read_length": 30, "start": 0, "end": 0}
        for p in (0, 1, 2, 3, 4, 5, 6, 7, 8)
    ])
    fs = compute_frame_summary(psite, ann, sample="S1")
    assert fs.n_reads == 9
    assert fs.frame_0 == pytest.approx(1 / 3)
    assert fs.dominance == pytest.approx(0.0)


def test_frame_summary_excludes_utr_reads() -> None:
    """UTR-bound reads must not count toward frame fractions — that's
    why frame summary uses CDS positions only."""
    ann = _annotation("ND1", l_utr5=10, l_cds=99, l_utr3=10)
    psite = pd.DataFrame([
        # 5'-UTR (P_site < l_utr5): should be ignored.
        {"sample_name": "S1", "chrom": "ND1", "P_site": 0, "read_length": 30, "start": 0, "end": 0},
        {"sample_name": "S1", "chrom": "ND1", "P_site": 5, "read_length": 30, "start": 0, "end": 0},
        # CDS frame 0:
        {"sample_name": "S1", "chrom": "ND1", "P_site": 10, "read_length": 30, "start": 0, "end": 0},
        {"sample_name": "S1", "chrom": "ND1", "P_site": 13, "read_length": 30, "start": 0, "end": 0},
    ])
    fs = compute_frame_summary(psite, ann, sample="S1")
    assert fs.n_reads == 2  # only the 2 CDS reads counted
    assert fs.frame_0 == pytest.approx(1.0)


# ---------- compute_metagene -----------------------------------------------


def test_metagene_start_window_and_3nt_periodicity() -> None:
    """Reads at codon starts should produce density at positions 0, 3,
    6, … in the start-aligned metagene."""
    ann = _annotation("ND1", l_utr5=0, l_cds=300, l_utr3=10)
    psite = pd.DataFrame([
        {"sample_name": "S1", "chrom": "ND1", "P_site": p, "read_length": 30, "start": 0, "end": 0}
        for p in (0, 3, 6, 9, 12, 15)
    ])
    profile = compute_metagene(psite, ann, sample="S1", anchor="start", window_nt=18)
    # All density at positions divisible by 3.
    nonzero_positions = profile.positions[profile.density > 0]
    assert list(nonzero_positions) == [0, 3, 6, 9, 12, 15]
    assert profile.density.sum() == 6.0
    assert len(profile.positions) == 18


def test_metagene_stop_window_uses_negative_relative_axis() -> None:
    """A stop-aligned metagene's positions are -window+1 .. 0 (the
    first nt of the stop codon is position 0)."""
    ann = _annotation("ND1", l_utr5=0, l_cds=300, l_utr3=10)
    stop = ann["stop_codon"].iloc[0]  # 300
    psite = pd.DataFrame([
        # P_site exactly on the stop codon -> rel = 0.
        {"sample_name": "S1", "chrom": "ND1", "P_site": stop, "read_length": 30, "start": 0, "end": 0},
        # 3 nt upstream of stop -> rel = -3.
        {"sample_name": "S1", "chrom": "ND1", "P_site": stop - 3, "read_length": 30, "start": 0, "end": 0},
        # 6 nt upstream -> rel = -6.
        {"sample_name": "S1", "chrom": "ND1", "P_site": stop - 6, "read_length": 30, "start": 0, "end": 0},
        # 1 nt downstream of stop -> rel = +1, OUTSIDE the (-w, 0] window.
        {"sample_name": "S1", "chrom": "ND1", "P_site": stop + 1, "read_length": 30, "start": 0, "end": 0},
    ])
    profile = compute_metagene(psite, ann, sample="S1", anchor="stop", window_nt=12)
    # Axis runs from -11 .. 0 inclusive (12 entries).
    assert profile.positions[0] == -11
    assert profile.positions[-1] == 0
    nonzero_positions = profile.positions[profile.density > 0]
    assert list(nonzero_positions) == [-6, -3, 0]
    assert profile.density.sum() == 3.0  # the +1 read was excluded


def test_metagene_empty_input_returns_zero_profile() -> None:
    ann = _annotation("ND1", l_utr5=0, l_cds=300, l_utr3=10)
    profile = compute_metagene(
        pd.DataFrame(), ann, sample="S1", anchor="start", window_nt=DEFAULT_WINDOW_NT,
    )
    assert profile.density.sum() == 0.0
    assert len(profile.positions) == DEFAULT_WINDOW_NT


def test_metagene_rejects_invalid_anchor() -> None:
    ann = _annotation("ND1", l_utr5=0, l_cds=300, l_utr3=10)
    with pytest.raises(ValueError, match="anchor"):
        compute_metagene(pd.DataFrame(), ann, sample="S1", anchor="middle")


# ---------- compute_strand_sanity ------------------------------------------


def test_strand_sanity_counts_minus_strand_reads() -> None:
    bed = pd.DataFrame([
        {"sample_name": "S1", "chrom": "ND1", "start": 0, "end": 30, "read_length": 30, "strand": "+"},
        {"sample_name": "S1", "chrom": "ND1", "start": 50, "end": 80, "read_length": 30, "strand": "-"},
        {"sample_name": "S1", "chrom": "ND1", "start": 90, "end": 120, "read_length": 30, "strand": "+"},
        {"sample_name": "S2", "chrom": "ND1", "start": 0, "end": 30, "read_length": 30, "strand": "-"},
    ])
    s1 = compute_strand_sanity(bed, sample="S1")
    assert s1.n_total == 3
    assert s1.n_minus == 1
    assert s1.minus_fraction == pytest.approx(1 / 3)
    s2 = compute_strand_sanity(bed, sample="S2")
    assert s2.minus_fraction == 1.0


def test_strand_sanity_returns_zero_when_no_strand_column() -> None:
    bed = pd.DataFrame([
        {"sample_name": "S1", "chrom": "ND1", "start": 0, "end": 30, "read_length": 30},
    ])
    s1 = compute_strand_sanity(bed, sample="S1")
    assert s1.n_total == 1
    assert s1.n_minus == 0


# ---------- run_periodicity_qc (end-to-end) --------------------------------


def test_run_periodicity_qc_writes_all_artefacts(tmp_path: Path) -> None:
    """The orchestrator should drop frame_summary.tsv,
    periodicity_start.tsv, periodicity_stop.tsv, strand_sanity.tsv,
    and a metagene PNG/SVG into output_dir."""
    ann = _annotation("ND1", l_utr5=0, l_cds=300, l_utr3=10)
    bed = pd.DataFrame([
        {"sample_name": "S1", "chrom": "ND1", "start": p - 12, "end": p - 12 + 30,
         "read_length": 30, "strand": "+"}
        for p in (0, 3, 6, 9, 12)  # P-sites all in frame 0
    ])
    out = tmp_path / "qc"
    result = run_periodicity_qc(
        bed_df=bed,
        annotation_df=ann,
        samples=["S1"],
        selected_offsets_by_sample={"S1": _selected_offsets(12)},
        selected_offsets_combined=None,
        offset_type="5",
        offset_site="p",
        output_dir=out,
        window_nt=30,
        plot=True,
    )

    for name in (
        "frame_summary.tsv", "periodicity_start.tsv",
        "periodicity_stop.tsv", "strand_sanity.tsv",
        "periodicity_metagene.png", "periodicity_metagene.svg",
    ):
        assert (out / name).exists(), f"missing {name} in QC output"

    # Spot-check the on-disk frame summary matches the in-memory result.
    body = (out / "frame_summary.tsv").read_text()
    assert "S1\t5\t1\t0\t0\t1" in body  # frame_0=1, dominance=1, 5 reads
    fs = result["frame_summary"][0]
    assert fs.frame_0 == pytest.approx(1.0)
    # Both metagenes are populated.
    start = result["periodicity_start"][0]
    stop = result["periodicity_stop"][0]
    assert start.density.sum() == 5.0
    # P_site positions 0/3/6/9/12 are all in the first 30 nt of CDS but
    # NOT within 30 nt of the stop codon (stop is at position 300), so
    # the stop-aligned metagene picks up zero reads with this fixture.
    assert stop.density.sum() == 0.0
