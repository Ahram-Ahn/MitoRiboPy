"""Tests for the standalone ``mitoribopy periodicity`` subcommand."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mitoribopy.cli import periodicity as periodicity_cli


def _make_site_table(
    *,
    sample: str = "WT",
    gene: str = "MT-CO1",
    transcript: str = "MT-CO1",
    cds_start: int = 0,
    cds_end: int = 600,
    read_length: int = 31,
    n_codons: int = 100,
    frame: int = 0,
    extra_uniform: int = 0,
) -> pd.DataFrame:
    """Build a site table with strong frame-0 dominance plus optional noise."""
    rows: list[dict] = []
    # Place sites starting at codon 10 so they survive the default
    # exclude_start_codons=6 mask.
    for i in range(n_codons):
        site_pos = cds_start + 30 + 3 * i + frame
        if site_pos >= cds_end - 9:  # respect exclude_stop_codons=3
            break
        rows.append({
            "sample": sample,
            "gene": gene,
            "transcript_id": transcript,
            "read_length": read_length,
            "site_type": "p",
            "site_pos": site_pos,
            "cds_start": cds_start,
            "cds_end": cds_end,
        })
    for i in range(extra_uniform):
        rows.append({
            "sample": sample,
            "gene": gene,
            "transcript_id": transcript,
            "read_length": read_length,
            "site_type": "p",
            "site_pos": cds_start + 30 + i,  # uniform across frames
            "cds_start": cds_start,
            "cds_end": cds_end,
        })
    return pd.DataFrame(rows)


def test_strong_periodicity_yields_good_call(tmp_path: Path) -> None:
    site_table = tmp_path / "sites.tsv"
    df = _make_site_table(cds_end=10_000, n_codons=2000)
    df.to_csv(site_table, sep="\t", index=False)
    out = tmp_path / "periodicity_out"
    rc = periodicity_cli.run([
        "--site-table", str(site_table),
        "--output", str(out),
        "--site", "p",
    ])
    assert rc == 0
    qc = pd.read_csv(out / "qc_summary.tsv", sep="\t")
    assert (qc["overall_qc_call"] == "good").any()
    assert (out / "frame_counts_by_sample_length.tsv").exists()
    assert (out / "qc_summary.md").exists()
    meta = json.loads((out / "periodicity.metadata.json").read_text())
    assert meta["frame_formula"] == "(site_pos - cds_start) % 3"


def test_uniform_input_yields_poor_or_low_depth(tmp_path: Path) -> None:
    # Uniform across all three frames -> no dominance.
    site_table = tmp_path / "sites.tsv"
    df = _make_site_table(n_codons=0, extra_uniform=900)
    df.to_csv(site_table, sep="\t", index=False)
    out = tmp_path / "out"
    rc = periodicity_cli.run([
        "--site-table", str(site_table),
        "--output", str(out),
        "--min-reads-per-length", "100",
    ])
    assert rc == 0
    qc = pd.read_csv(out / "qc_summary.tsv", sep="\t")
    assert qc["overall_qc_call"].iloc[0] in ("poor", "warn", "low_depth")


def test_missing_required_column_fails(tmp_path: Path) -> None:
    bad = tmp_path / "bad.tsv"
    pd.DataFrame({"sample": ["x"], "gene": ["MT-CO1"]}).to_csv(
        bad, sep="\t", index=False,
    )
    rc = periodicity_cli.run([
        "--site-table", str(bad),
        "--output", str(tmp_path / "out"),
    ])
    assert rc == 2


def test_count_column_is_honored(tmp_path: Path) -> None:
    site_table = tmp_path / "sites.tsv"
    df = _make_site_table(n_codons=10)
    df["count"] = 100
    df.to_csv(site_table, sep="\t", index=False)
    out = tmp_path / "out"
    rc = periodicity_cli.run([
        "--site-table", str(site_table),
        "--output", str(out),
        "--min-reads-per-length", "500",  # only passes if counts are summed
    ])
    assert rc == 0
    by_length = pd.read_csv(out / "frame_counts_by_sample_length.tsv", sep="\t")
    # 10 rows * count=100 = 1000 ≥ 500, so qc_call should be 'good'.
    assert (by_length["qc_call"] == "good").all()


def test_help_lists_periodicity_subcommand() -> None:
    from mitoribopy.cli import _SUBCOMMANDS, _SUBCOMMAND_SUMMARIES

    assert "periodicity" in _SUBCOMMANDS
    names = [name for name, _ in _SUBCOMMAND_SUMMARIES]
    assert "periodicity" in names
