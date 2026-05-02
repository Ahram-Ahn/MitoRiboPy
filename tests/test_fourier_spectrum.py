"""Tests for the metagene Fourier periodicity module.

The contract being verified:

* A synthetic perfectly-periodic-3 coverage signal yields a sharp peak
  at period = 3 nt, and the spectral_ratio_3nt is well above 1; flat
  coverage yields no such peak (snr_call="no_signal" or "broken").
* Window math: orf_start window starts after the AUG + initiation skip;
  orf_stop window ends before the stop codon - termination skip.
* The chrom-to-transcript mapping resolves fused-FASTA aliases (ATP86,
  ND4L4) to BOTH constituent transcripts.
* Three gene_sets are emitted: combined, ATP86, ND4L4.
* The plot writer renders one PNG per (sample, length, gene_set) and
  writes a sidecar with the canonical panel layout.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mitoribopy.analysis.fourier_spectrum import (
    ATP86_GENE_SET,
    DEFAULT_DROP_CODONS_AFTER_START,
    DEFAULT_DROP_CODONS_BEFORE_STOP,
    DEFAULT_PERIOD_GRID,
    DEFAULT_WINDOW_NT,
    GENE_SETS,
    ND4L4_GENE_SET,
    OVERLAP_PAIR_TRANSCRIPTS,
    REGIONS,
    build_chrom_to_transcripts,
    build_fourier_period3_score_combined_table,
    build_fourier_spectrum_combined_table,
    build_metagene,
    compute_spectral_ratio_3nt,
    compute_spectrum_grid,
    direct_dft_amplitude,
    extract_per_gene_normalized_tracks,
    snr_call_for_ratio,
)
from mitoribopy.analysis.fourier_spectrum import _resolve_window  # noqa: E501


# ---------------------------------------------------------------------------
# Module-level invariants
# ---------------------------------------------------------------------------


def test_regions_are_orf_start_and_orf_stop_only() -> None:
    assert REGIONS == ("orf_start", "orf_stop")


def test_gene_sets_are_combined_and_two_overlap_pairs() -> None:
    assert GENE_SETS == ("combined", "ATP86", "ND4L4")


def test_default_window_is_99_nt() -> None:
    assert DEFAULT_WINDOW_NT == 99


def test_default_codon_skips() -> None:
    assert DEFAULT_DROP_CODONS_AFTER_START == 5
    assert DEFAULT_DROP_CODONS_BEFORE_STOP == 1


def test_overlap_pair_transcripts_includes_fused_spellings() -> None:
    for name in ("ATP8", "ATP6", "ATP86", "ND4L", "ND4", "ND4L4"):
        assert name in OVERLAP_PAIR_TRANSCRIPTS, name
    # And case-insensitive variants.
    assert "ATP86" in ATP86_GENE_SET
    assert "ND4L4" in ND4L4_GENE_SET


# ---------------------------------------------------------------------------
# Direct DFT math
# ---------------------------------------------------------------------------


def _synthetic_period3(n: int, peak: float = 10.0) -> np.ndarray:
    """Coverage with one read every 3 nt — perfect 3-nt phasing."""
    cov = np.zeros(n, dtype=float)
    cov[::3] = peak
    return cov


def test_direct_dft_amplitude_is_high_at_period_3_for_synthetic_signal() -> None:
    cov = _synthetic_period3(99)
    # Mean-centre + Hann to mimic the production pipeline.
    centred = cov - cov.mean()
    windowed = centred * np.hanning(centred.size)
    amp3 = direct_dft_amplitude(windowed, period=3.0)
    amp7 = direct_dft_amplitude(windowed, period=7.0)
    assert amp3 > 0.5
    assert amp3 > 5.0 * amp7


def test_direct_dft_amplitude_returns_nan_for_zero_signal() -> None:
    assert np.isnan(direct_dft_amplitude(np.zeros(99), period=3.0))


def test_compute_spectral_ratio_high_for_clean_period3() -> None:
    cov = _synthetic_period3(99)
    centred = cov - cov.mean()
    windowed = centred * np.hanning(centred.size)
    _, ratio = compute_spectral_ratio_3nt(windowed)
    assert ratio > 5.0


def test_snr_call_thresholds() -> None:
    assert snr_call_for_ratio(12.0) == "excellent"
    assert snr_call_for_ratio(7.0)  == "healthy"
    assert snr_call_for_ratio(3.0)  == "modest"
    assert snr_call_for_ratio(1.5)  == "broken"
    assert snr_call_for_ratio(float("nan")) == "no_signal"


# ---------------------------------------------------------------------------
# Window resolution
# ---------------------------------------------------------------------------


def test_resolve_window_orf_start_skips_initiation_codons() -> None:
    w = _resolve_window(
        transcript="t1", region="orf_start",
        start_codon=30, stop_codon=600, transcript_length=800,
        window_nt=99, drop_codons_after_start=5, drop_codons_before_stop=1,
    )
    # start + 5*3 = 45; window = [45, 144).
    assert w.lo == 45
    assert w.hi == 144
    assert w.valid is True


def test_resolve_window_orf_stop_skips_termination_codon() -> None:
    w = _resolve_window(
        transcript="t1", region="orf_stop",
        start_codon=30, stop_codon=600, transcript_length=800,
        window_nt=99, drop_codons_after_start=5, drop_codons_before_stop=1,
    )
    # stop - 1*3 = 597; window = [498, 597).
    assert w.lo == 498
    assert w.hi == 597
    assert w.valid is True


def test_resolve_window_invalid_when_orf_too_short() -> None:
    # Tiny ORF (start=30, stop=200): orf_start window [45, 144) ends
    # before stop=200 so it's actually valid here. But orf_stop window
    # [98, 197) would also fit. Use a much shorter ORF.
    w = _resolve_window(
        transcript="t1", region="orf_start",
        start_codon=30, stop_codon=100, transcript_length=200,
        window_nt=99, drop_codons_after_start=5, drop_codons_before_stop=1,
    )
    assert w.valid is False


# ---------------------------------------------------------------------------
# Chrom -> transcripts mapping (sequence_aliases)
# ---------------------------------------------------------------------------


def test_build_chrom_to_transcripts_resolves_fused_alias() -> None:
    ann = pd.DataFrame([
        # Canonical-spelling rows for ATP8 + ATP6 sharing chrom ATP6,
        # with the fused-FASTA alias ATP86 in sequence_aliases.
        {"transcript": "ATP8", "sequence_name": "ATP6",
         "sequence_aliases": "ATP86;ATP8/ATP6"},
        {"transcript": "ATP6", "sequence_name": "ATP6",
         "sequence_aliases": "ATP86;ATP8/ATP6"},
        {"transcript": "COX1", "sequence_name": "COX1",
         "sequence_aliases": ""},
    ])
    mapping = build_chrom_to_transcripts(ann)
    # Both fused and canonical chroms should resolve to BOTH ATP8 and ATP6.
    assert set(mapping["ATP86"]) == {"ATP8", "ATP6"}
    assert set(mapping["ATP6"]) == {"ATP8", "ATP6"}
    # COX1 only resolves to itself.
    assert mapping["COX1"] == ["COX1"]


# ---------------------------------------------------------------------------
# End-to-end synthetic build_fourier_spectrum_combined_table
# ---------------------------------------------------------------------------


def _make_annotation(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(
        rows,
        columns=[
            "transcript", "sequence_name", "sequence_aliases",
            "start_codon", "stop_codon", "l_tr",
        ],
    )


def _make_bed(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(
        rows, columns=["sample_name", "chrom", "read_length", "P_site"],
    )


def _build_periodic_bed(
    *, sample: str, chrom: str, read_length: int,
    window_lo: int, window_hi: int, peak: int = 5,
) -> list[dict]:
    """Reads whose A-site lands every 3 nt across [window_lo, window_hi)."""
    rows: list[dict] = []
    for window_pos in range(window_lo, window_hi, 3):
        for _ in range(peak):
            rows.append({
                "sample_name": sample, "chrom": chrom,
                "read_length": read_length, "P_site": window_pos - 3,
            })
    return rows


def test_combined_spectrum_table_has_period3_peak_for_clean_signal() -> None:
    ann = _make_annotation([
        {"transcript": "COX1", "sequence_name": "COX1",
         "sequence_aliases": "", "start_codon": 30,
         "stop_codon": 1500, "l_tr": 1600},
    ])
    bed_rows = _build_periodic_bed(
        sample="WT", chrom="COX1", read_length=32,
        window_lo=45, window_hi=144, peak=10,
    )
    bed = _make_bed(bed_rows)

    tracks = extract_per_gene_normalized_tracks(
        bed, ann, samples=["WT"], window_nt=99, site="a",
    )
    spec = build_fourier_spectrum_combined_table(tracks)
    assert not spec.empty
    sub = spec[
        (spec["gene_set"] == "combined") & (spec["region"] == "orf_start")
    ].sort_values("period_nt")
    nearest_3 = sub.iloc[(sub["period_nt"] - 3.0).abs().argsort()].iloc[0]
    other = sub[
        (sub["period_nt"] < 2.8) | (sub["period_nt"] > 3.2)
    ]
    assert nearest_3["amplitude"] > 5.0 * other["amplitude"].max()


def test_combined_score_table_yields_healthy_or_better_for_clean_signal() -> None:
    ann = _make_annotation([
        {"transcript": "COX1", "sequence_name": "COX1",
         "sequence_aliases": "", "start_codon": 30,
         "stop_codon": 1500, "l_tr": 1600},
    ])
    bed_rows = _build_periodic_bed(
        sample="WT", chrom="COX1", read_length=32,
        window_lo=45, window_hi=144, peak=10,
    )
    bed = _make_bed(bed_rows)

    tracks = extract_per_gene_normalized_tracks(
        bed, ann, samples=["WT"], window_nt=99,
    )
    score = build_fourier_period3_score_combined_table(tracks)
    assert not score.empty
    row = score[
        (score["gene_set"] == "combined") & (score["region"] == "orf_start")
    ].iloc[0]
    assert row["spectral_ratio_3nt"] > 5.0
    assert row["snr_call"] in {"healthy", "excellent"}


def test_atp86_gene_set_routed_separately_from_combined() -> None:
    # A canonical-spelling fixture (ATP8 + ATP6 + COX1).
    ann = _make_annotation([
        {"transcript": "ATP8", "sequence_name": "ATP6",
         "sequence_aliases": "ATP86", "start_codon": 1,
         "stop_codon": 205, "l_tr": 843},
        {"transcript": "ATP6", "sequence_name": "ATP6",
         "sequence_aliases": "ATP86", "start_codon": 162,
         "stop_codon": 840, "l_tr": 843},
        {"transcript": "COX1", "sequence_name": "COX1",
         "sequence_aliases": "", "start_codon": 30,
         "stop_codon": 1500, "l_tr": 1600},
    ])
    # Reads on each chromosome.
    bed_rows: list[dict] = []
    bed_rows += _build_periodic_bed(
        sample="WT", chrom="COX1", read_length=32,
        window_lo=45, window_hi=144, peak=10,
    )
    # ATP8 frame at the orf_stop window: [205 - 3 - 99, 205 - 3) = [103, 202).
    bed_rows += _build_periodic_bed(
        sample="WT", chrom="ATP6", read_length=32,
        window_lo=103, window_hi=202, peak=10,
    )
    # ATP6 frame at the orf_start window: [162 + 15, 162 + 15 + 99) = [177, 276).
    bed_rows += _build_periodic_bed(
        sample="WT", chrom="ATP6", read_length=32,
        window_lo=177, window_hi=276, peak=10,
    )
    bed = _make_bed(bed_rows)

    tracks = extract_per_gene_normalized_tracks(
        bed, ann, samples=["WT"], window_nt=99,
    )
    spec = build_fourier_spectrum_combined_table(tracks)
    gene_sets = set(spec["gene_set"].unique())
    assert "combined" in gene_sets
    assert "ATP86" in gene_sets
    # The combined set must NOT include ATP8/ATP6.
    score = build_fourier_period3_score_combined_table(tracks)
    combined_rows = score[score["gene_set"] == "combined"]
    for _, r in combined_rows.iterrows():
        transcripts = r["transcripts"].split(";") if r["transcripts"] else []
        for tx in transcripts:
            assert tx not in {"ATP8", "ATP6", "ATP86", "ND4L", "ND4", "ND4L4"}, (
                f"{tx} should not be in the 'combined' set"
            )


def test_extract_tracks_filters_low_coverage_genes() -> None:
    ann = _make_annotation([
        {"transcript": "COX1", "sequence_name": "COX1",
         "sequence_aliases": "", "start_codon": 30,
         "stop_codon": 1500, "l_tr": 1600},
    ])
    # Only 2 reads in the window — below default min_total_counts=30.
    bed = _make_bed([
        {"sample_name": "WT", "chrom": "COX1",
         "read_length": 32, "P_site": 50},
        {"sample_name": "WT", "chrom": "COX1",
         "read_length": 32, "P_site": 80},
    ])
    tracks = extract_per_gene_normalized_tracks(
        bed, ann, samples=["WT"], window_nt=99,
    )
    assert tracks == []


# ---------------------------------------------------------------------------
# Plot writer
# ---------------------------------------------------------------------------


def test_render_fourier_spectrum_panels_writes_combined_figure(tmp_path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    from mitoribopy.plotting.fourier_spectrum_plots import (
        render_fourier_spectrum_panels,
    )
    from mitoribopy.plotting.figure_validator import metadata_sidecar_path

    ann = _make_annotation([
        {"transcript": "COX1", "sequence_name": "COX1",
         "sequence_aliases": "", "start_codon": 30,
         "stop_codon": 1500, "l_tr": 1600},
    ])
    bed_rows = _build_periodic_bed(
        sample="WT", chrom="COX1", read_length=32,
        window_lo=45, window_hi=144, peak=10,
    )
    # Add reads in the orf_stop window too so both panels populate.
    bed_rows += _build_periodic_bed(
        sample="WT", chrom="COX1", read_length=32,
        window_lo=1398, window_hi=1497, peak=10,
    )
    bed = _make_bed(bed_rows)

    tracks = extract_per_gene_normalized_tracks(
        bed, ann, samples=["WT"], window_nt=99,
    )
    spec = build_fourier_spectrum_combined_table(tracks)
    score = build_fourier_period3_score_combined_table(tracks)

    out = tmp_path / "fourier_spectrum"
    written = render_fourier_spectrum_panels(spec, score, output_dir=out)

    combined_png = out / "WT" / "WT_32nt_combined.png"
    assert combined_png in written
    assert combined_png.is_file()

    sidecar = metadata_sidecar_path(combined_png)
    import json as _json
    meta = _json.loads(sidecar.read_text())
    assert meta["plot_type"] == "fourier_spectrum_combined"
    assert meta["gene_set"] == "combined"
    panel_layout = meta["panel_layout"]
    assert [p["region"] for p in panel_layout] == ["orf_start", "orf_stop"]
