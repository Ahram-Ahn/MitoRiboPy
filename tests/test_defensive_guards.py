"""Targeted defensive-guard tests for Task 9.

Each test drives a predictable failure mode and asserts the pipeline
either recovers with a clear warning or fails with an actionable error
(no silent corruption, no bare stack trace).

Note: ``mitoribopy``'s console logger sets ``propagate=False`` and attaches
its StreamHandler at module-load time (before pytest replaces sys.stdout),
so neither ``capsys``, ``capfd``, nor ``caplog`` sees the records cleanly.
We attach our own in-memory handler to the LOGGER for the test's scope.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import pytest

from mitoribopy.cli import common as cli_common
from mitoribopy.console import LOGGER
from mitoribopy.io.bed_reader import (
    compute_unfiltered_read_length_summary,
    process_bed_files,
)


@pytest.fixture
def captured_log():
    """Yield a list of log messages captured from the ``mitoribopy`` logger
    during the test's scope."""
    buffer: list[str] = []

    class ListHandler(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            buffer.append(self.format(record))

    handler = ListHandler()
    handler.setFormatter(logging.Formatter("%(message)s"))
    LOGGER.addHandler(handler)
    try:
        yield buffer
    finally:
        LOGGER.removeHandler(handler)


# ---------- --threads <= 0 (CLI guard) ---------------------------------------


def test_threads_zero_raises_systemexit_with_message() -> None:
    args = type("ns", (), {"threads": 0})()
    with pytest.raises(SystemExit) as exc:
        cli_common.apply_common_arguments(args)
    assert "--threads must be a positive integer" in str(exc.value)


def test_threads_negative_raises_systemexit() -> None:
    args = type("ns", (), {"threads": -1})()
    with pytest.raises(SystemExit):
        cli_common.apply_common_arguments(args)


# ---------- malformed BED rows (end <= start) -------------------------------


def _write_bed(path: Path, rows: list[tuple[str, int, int]]) -> None:
    path.write_text(
        "\n".join(
            f"{chrom}\t{start}\t{end}\tr{idx}\t42\t+"
            for idx, (chrom, start, end) in enumerate(rows)
        )
        + "\n",
        encoding="utf-8",
    )


def test_process_bed_files_drops_malformed_rows_with_warning(tmp_path, captured_log) -> None:
    """end <= start must be dropped with a clear WARNING, not propagate as
    a zero/negative read_length into downstream filters."""
    bed = tmp_path / "s.bed"
    _write_bed(
        bed,
        [
            ("TX1", 0, 30),    # ok
            ("TX1", 10, 10),   # malformed: length 0
            ("TX1", 50, 40),   # malformed: end < start
            ("TX1", 100, 130),  # ok
        ],
    )
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    _, samples = process_bed_files(
        input_dir=str(tmp_path),
        output_dir=str(out_dir),
        organism="h",
        annotation_df=pd.DataFrame(),
        rpf_range=[30],
    )

    assert samples == ["s"]
    summary = pd.read_csv(out_dir / "read_length_summary.csv")
    # Only the two well-formed rows survive.
    assert summary["count"].sum() == 2
    assert all(summary["read_length"] == 30)

    warned = [msg for msg in captured_log if "malformed" in msg]
    assert warned, "Expected WARNING about malformed BED rows"
    assert "[BED] WARNING:" in warned[0]


def test_compute_unfiltered_summary_drops_malformed_rows(tmp_path, captured_log) -> None:
    bed = tmp_path / "s.bed"
    _write_bed(
        bed,
        [
            ("TX1", 0, 30),
            ("TX1", 10, 10),  # malformed
            ("TX1", 100, 130),
        ],
    )
    out_csv = tmp_path / "unfiltered.csv"
    compute_unfiltered_read_length_summary(
        input_dir=str(tmp_path),
        output_csv=str(out_csv),
        total_counts_map={},
        read_length_range=(15, 50),
    )
    assert out_csv.is_file()
    summary = pd.read_csv(out_csv)
    assert summary["count"].sum() == 2
    warned = [msg for msg in captured_log if "malformed" in msg]
    assert warned, "Expected WARNING about malformed BED rows"
    assert "[QC] WARNING:" in warned[0]


# ---------- FASTA sequences with no annotation match ------------------------


def test_coverage_plots_warn_on_unmapped_fasta_record(tmp_path, captured_log) -> None:
    """Coverage plots should WARN (not crash, not silently continue) when a
    FASTA record has no corresponding annotation row - the codon/frame
    variants will have been skipped for that transcript."""
    from types import SimpleNamespace
    from mitoribopy.plotting.coverage_profile_plots import run_coverage_profile_plots

    fasta_path = tmp_path / "tiny.fa"
    # Two records: ND1 (annotated) and GHOST (not annotated).
    fasta_path.write_text(
        ">ND1\n" + ("A" * 120) + "\n>GHOST\n" + ("C" * 90) + "\n",
        encoding="utf-8",
    )
    annotation_df = pd.DataFrame(
        [
            {
                "transcript": "ND1",
                "sequence_name": "ND1",
                "l_tr": 120,
                "l_utr5": 9,
                "l_utr3": 9,
            }
        ]
    )
    filtered_bed_df = pd.DataFrame(
        [
            {
                "chrom": "ND1",
                "start": 9, "end": 39,
                "read_length": 30, "sample_name": "s1",
            }
        ]
    )
    selected_offsets_df = pd.DataFrame(
        [{"Read Length": 30, "Most Enriched 5' Offset": 12, "Most Enriched 3' Offset": 18}]
    )
    args = SimpleNamespace(
        plot_format="svg", offset_site="p", total_mrna_map={"s1": 100}, cap_percentile=0.999,
    )

    run_coverage_profile_plots(
        sample_dirs=[str(tmp_path / "s1")],
        selected_offsets_df=selected_offsets_df,
        offset_type="5",
        fasta_file=fasta_path,
        output_dir=tmp_path / "cov",
        args=args,
        annotation_df=annotation_df,
        filtered_bed_df=filtered_bed_df,
    )
    warned = [msg for msg in captured_log if "no matching annotation" in msg]
    assert warned, "Expected WARNING about GHOST having no annotation"
    assert "GHOST" in warned[0]
