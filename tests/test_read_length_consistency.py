"""Regression test: the read-length CSV summary and the value that feeds
the distribution plot must always agree, row-for-row.

Background: in earlier builds the reported (CSV) and plotted distribution
appeared to diverge by ~3 nt. Current code uses ``end - start`` in both
paths with no offset in between, and this test locks that invariant in.

See reports/read_length_discrepancy_report.md for the full diagnosis.
"""

from __future__ import annotations

import pandas as pd
import pytest

from mitoribopy.io.bed_reader import (
    compute_unfiltered_read_length_summary,
    process_bed_files,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


# Known-by-construction read-length histogram. Keeping this tight so a
# future maintainer can verify the expected values by hand.
EXPECTED_HISTOGRAM: dict[int, int] = {
    28: 1,   # exactly one read of length 28
    29: 2,   # two reads of length 29
    30: 5,   # the peak
    31: 3,
    32: 1,
}


def _write_bed(path, histogram: dict[int, int]) -> None:
    """Write a BED6 file whose per-length row count matches ``histogram``.

    All rows are on transcript 'TX1', forward strand. Starts are walked
    forward so starts never collide — the coordinates themselves do not
    matter for the length summary.
    """
    rows = []
    cursor = 0
    for length, count in histogram.items():
        for i in range(count):
            rows.append(
                f"TX1\t{cursor}\t{cursor + length}\tr{cursor}\t42\t+\n"
            )
            cursor += length + 1  # separate rows so starts stay unique
    path.write_text("".join(rows), encoding="utf-8")


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_filtered_summary_matches_bed_arithmetic(tmp_path):
    """process_bed_files must emit ``end - start`` counts verbatim."""
    bed_dir = tmp_path / "bed"
    out_dir = tmp_path / "out"
    bed_dir.mkdir()
    out_dir.mkdir()
    _write_bed(bed_dir / "sampleA.bed", EXPECTED_HISTOGRAM)

    rpf_range = list(range(28, 33))  # keep every length in EXPECTED_HISTOGRAM
    annotation_df = pd.DataFrame()   # process_bed_files does not consume this
    concatenated, samples = process_bed_files(
        input_dir=str(bed_dir),
        output_dir=str(out_dir),
        organism="h",
        annotation_df=annotation_df,
        rpf_range=rpf_range,
    )

    assert samples == ["sampleA"]
    summary_df = pd.read_csv(out_dir / "read_length_summary.csv")
    got = dict(zip(summary_df["read_length"], summary_df["count"]))
    assert got == EXPECTED_HISTOGRAM, (
        f"Filtered summary disagrees with end - start histogram.\n"
        f"  expected: {EXPECTED_HISTOGRAM}\n"
        f"  got:      {got}"
    )

    # The in-memory DataFrame the plot sees must encode the same lengths.
    plot_histogram = (
        concatenated["read_length"].value_counts().sort_index().to_dict()
    )
    assert plot_histogram == EXPECTED_HISTOGRAM


def test_unfiltered_summary_matches_bed_arithmetic(tmp_path):
    """compute_unfiltered_read_length_summary must also agree verbatim."""
    bed_dir = tmp_path / "bed"
    bed_dir.mkdir()
    out_csv = tmp_path / "unfiltered.csv"
    _write_bed(bed_dir / "sampleA.bed", EXPECTED_HISTOGRAM)

    compute_unfiltered_read_length_summary(
        input_dir=str(bed_dir),
        output_csv=str(out_csv),
        total_counts_map={"sampleA": 1_000_000},
        read_length_range=(15, 50),
    )
    assert out_csv.exists()

    summary_df = pd.read_csv(out_csv)
    got = dict(zip(summary_df["read_length"], summary_df["count"]))
    assert got == EXPECTED_HISTOGRAM


def test_filtered_and_unfiltered_paths_agree_on_shared_window(tmp_path):
    """Both paths, restricted to the same window, must return identical counts.

    This is the regression guard for the originally-reported ~3 nt shift.
    If either path ever re-introduces an offset, this test fails.
    """
    bed_dir = tmp_path / "bed"
    out_dir = tmp_path / "out"
    bed_dir.mkdir()
    out_dir.mkdir()
    _write_bed(bed_dir / "sampleA.bed", EXPECTED_HISTOGRAM)

    # Filtered path
    process_bed_files(
        input_dir=str(bed_dir),
        output_dir=str(out_dir),
        organism="h",
        annotation_df=pd.DataFrame(),
        rpf_range=list(EXPECTED_HISTOGRAM.keys()),
    )
    filtered_df = pd.read_csv(out_dir / "read_length_summary.csv")

    # Unfiltered path (same numeric window)
    unfilt_csv = tmp_path / "unfilt.csv"
    compute_unfiltered_read_length_summary(
        input_dir=str(bed_dir),
        output_csv=str(unfilt_csv),
        total_counts_map={},
        read_length_range=(
            min(EXPECTED_HISTOGRAM),
            max(EXPECTED_HISTOGRAM),
        ),
    )
    unfilt_df = pd.read_csv(unfilt_csv)

    filtered_hist = dict(
        zip(filtered_df["read_length"], filtered_df["count"])
    )
    unfilt_hist = dict(zip(unfilt_df["read_length"], unfilt_df["count"]))
    assert filtered_hist == unfilt_hist == EXPECTED_HISTOGRAM
