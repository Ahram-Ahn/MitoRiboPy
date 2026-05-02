"""Tests for the standalone ``mitoribopy periodicity`` subcommand.

The subcommand now runs the metagene Fourier bundle on a pre-
assigned site table. The contract being verified:

* A strong period-3 signal in the input yields a high spectral_ratio_3nt
  and an ``snr_call`` of ``healthy`` or ``excellent``.
* The expected TSV outputs are written and the metadata sidecar records
  the Fourier knobs.
* Missing required columns produce a non-zero exit.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mitoribopy.cli import periodicity as periodicity_cli


def _make_periodic_site_table(
    *,
    sample: str = "WT",
    gene: str = "MT-CO1",
    transcript: str = "MT-CO1",
    cds_start: int = 0,
    cds_end: int = 1500,
    read_length: int = 31,
    peak_per_codon: int = 5,
) -> pd.DataFrame:
    """Synthetic site table with one read every 3 nt (perfect frame-0)."""
    rows: list[dict] = []
    for site_pos in range(cds_start + 18, cds_end - 6, 3):
        for _ in range(peak_per_codon):
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
    return pd.DataFrame(rows)


def test_strong_periodicity_yields_healthy_or_better_snr_call(tmp_path: Path) -> None:
    site_table = tmp_path / "sites.tsv"
    df = _make_periodic_site_table(cds_end=2000, peak_per_codon=10)
    df.to_csv(site_table, sep="\t", index=False)
    out = tmp_path / "periodicity_out"
    rc = periodicity_cli.run([
        "--site-table", str(site_table),
        "--output", str(out),
        "--site", "p",
    ])
    assert rc == 0
    score = pd.read_csv(out / "fourier_period3_score_combined.tsv", sep="\t")
    assert not score.empty
    # The combined panel for orf_start should hit at least 'healthy' on
    # a perfect period-3 input.
    combined = score[
        (score["gene_set"] == "combined") & (score["region"] == "orf_start")
    ]
    assert not combined.empty
    assert combined["snr_call"].iloc[0] in {"healthy", "excellent"}


def test_outputs_and_metadata_are_written(tmp_path: Path) -> None:
    site_table = tmp_path / "sites.tsv"
    df = _make_periodic_site_table(cds_end=2000, peak_per_codon=10)
    df.to_csv(site_table, sep="\t", index=False)
    out = tmp_path / "out"
    rc = periodicity_cli.run([
        "--site-table", str(site_table),
        "--output", str(out),
        "--site", "p",
    ])
    assert rc == 0
    assert (out / "fourier_spectrum_combined.tsv").is_file()
    assert (out / "fourier_period3_score_combined.tsv").is_file()
    meta = json.loads((out / "periodicity.metadata.json").read_text())
    assert meta["method"] == "metagene_dft"
    assert meta["fourier_window_nt"] == 99
    assert meta["regions"] == ["orf_start", "orf_stop"]
    assert meta["gene_sets"] == ["combined", "ATP86", "ND4L4"]


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


def test_window_override_is_recorded_in_metadata(tmp_path: Path) -> None:
    site_table = tmp_path / "sites.tsv"
    df = _make_periodic_site_table(cds_end=2000, peak_per_codon=10)
    df.to_csv(site_table, sep="\t", index=False)
    out = tmp_path / "out"
    rc = periodicity_cli.run([
        "--site-table", str(site_table),
        "--output", str(out),
        "--fourier-window-nt", "150",
    ])
    assert rc == 0
    meta = json.loads((out / "periodicity.metadata.json").read_text())
    assert meta["fourier_window_nt"] == 150


def test_help_lists_periodicity_subcommand() -> None:
    from mitoribopy.cli import _SUBCOMMANDS, _SUBCOMMAND_SUMMARIES

    assert "periodicity" in _SUBCOMMANDS
    names = [name for name, _ in _SUBCOMMAND_SUMMARIES]
    assert "periodicity" in names
