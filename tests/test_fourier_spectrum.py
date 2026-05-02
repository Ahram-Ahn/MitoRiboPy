"""Tests for the Wakigawa-style Fourier-spectrum periodicity module.

The contract being verified:

* A synthetic perfectly-periodic-3 coverage signal yields a sharp peak
  at period = 3 nt; flat coverage yields no such peak.
* Window math: ORF window is strictly upstream of the stop codon, 3'
  UTR window is strictly downstream — no off-by-three at the boundary.
* Read-length stratification produces independent spectra per length.
* The overlap-upstream predicate flags ATP8 / ND4L (not their fused
  spellings) and they are excluded from the combined summary while
  remaining present in the per-gene table.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mitoribopy.analysis.fourier_spectrum import (
    DEFAULT_PERIOD_RANGE,
    DEFAULT_WINDOW_NT,
    OVERLAP_UPSTREAM_GENES,
    build_fourier_spectrum_table,
    build_period3_score_table,
    build_period3_summary_table,
    build_window_coverage,
    compute_amplitude_spectrum,
    is_overlap_upstream_orf,
)
from mitoribopy.analysis.fourier_spectrum import _resolve_window  # noqa: E501 - private helper, tested intentionally


# ---------------------------------------------------------------------------
# Overlap-upstream predicate
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,expected",
    [
        ("MT-ATP8", True),
        ("MTATP8", True),
        ("ATP8", True),
        ("mt-atp8", True),
        ("MT_ATP8", True),
        ("MT-ND4L", True),
        ("ND4L", True),
        # NOT excluded: the downstream genes are fine on their own.
        ("MT-ATP6", False),
        ("MT-ND4", False),
        # NOT excluded: fused-FASTA entries cover the whole region.
        ("ATP86", False),
        ("ND4L4", False),
        ("MT-ATP86", False),
        ("MT-ND4L4", False),
        # Negative controls.
        ("MT-CO1", False),
        ("MT-CYB", False),
        (None, False),
        ("", False),
    ],
)
def test_is_overlap_upstream_orf_recognises_the_right_spellings(
    name, expected,
) -> None:
    assert is_overlap_upstream_orf(name) is expected


def test_overlap_upstream_set_excludes_fused_names() -> None:
    # Guard against future maintainers re-adding "ATP86" by mistake.
    assert "ATP86" not in OVERLAP_UPSTREAM_GENES
    assert "ND4L4" not in OVERLAP_UPSTREAM_GENES


# ---------------------------------------------------------------------------
# Spectrum math
# ---------------------------------------------------------------------------


def _synthetic_period3(window_nt: int, peak: float = 10.0) -> np.ndarray:
    """Coverage with one read every 3 nt — perfect 3-nt phasing."""
    cov = np.zeros(window_nt, dtype=float)
    cov[::3] = peak
    return cov


def test_amplitude_spectrum_peaks_at_period_three_for_synthetic_signal() -> None:
    cov = _synthetic_period3(99, peak=10.0)  # multiple of 3 -> exact bin
    spec = compute_amplitude_spectrum(cov)
    assert not spec.empty
    nearest_3 = spec.iloc[(spec["period_nt"] - 3.0).abs().argsort()].iloc[0]
    # The period-3 bin should hold most of the amplitude mass.
    other = spec[spec["period_nt"] != nearest_3["period_nt"]]
    assert nearest_3["amplitude"] > 5.0 * other["amplitude"].max()


def test_amplitude_spectrum_returns_empty_for_flat_signal() -> None:
    cov = np.full(100, 5.0)
    spec = compute_amplitude_spectrum(cov)
    # Mean-subtracting a constant signal leaves zero variance — empty
    # spectrum is the documented behaviour.
    assert spec.empty


def test_amplitude_spectrum_returns_empty_for_short_signal() -> None:
    # Need at least 2*ceil(period_hi) samples; default period_hi=10 -> 20.
    cov = _synthetic_period3(15, peak=5.0)
    spec = compute_amplitude_spectrum(cov)
    assert spec.empty


def test_amplitude_spectrum_period_range_is_respected() -> None:
    cov = _synthetic_period3(99, peak=10.0)
    spec = compute_amplitude_spectrum(cov, period_range=(2.5, 4.5))
    assert not spec.empty
    assert spec["period_nt"].min() >= 2.5
    assert spec["period_nt"].max() <= 4.5


# ---------------------------------------------------------------------------
# Window resolution
# ---------------------------------------------------------------------------


def test_resolve_window_orf_is_strictly_upstream_of_stop_codon() -> None:
    w = _resolve_window(
        transcript="t1", region="orf",
        stop_codon=300, transcript_length=400, window_nt=100,
    )
    assert w.lo == 200
    assert w.hi == 300  # stop_codon itself NOT included
    assert w.valid is True


def test_resolve_window_utr3_is_strictly_downstream_of_stop_codon() -> None:
    w = _resolve_window(
        transcript="t1", region="utr3",
        stop_codon=300, transcript_length=500, window_nt=100,
    )
    # Wakigawa convention: window starts AFTER the stop codon trinucleotide.
    assert w.lo == 303
    assert w.hi == 403
    assert w.valid is True


def test_resolve_window_marks_invalid_when_off_transcript_end() -> None:
    # 3' UTR runs off the end.
    w = _resolve_window(
        transcript="t1", region="utr3",
        stop_codon=400, transcript_length=420, window_nt=100,
    )
    assert w.valid is False


def test_resolve_window_marks_invalid_when_off_transcript_start() -> None:
    # ORF runs before transcript start (implausible but defensive).
    w = _resolve_window(
        transcript="t1", region="orf",
        stop_codon=50, transcript_length=200, window_nt=100,
    )
    assert w.valid is False


def test_build_window_coverage_drops_positions_outside_window() -> None:
    w = _resolve_window(
        transcript="t1", region="orf",
        stop_codon=300, transcript_length=500, window_nt=100,
    )
    positions = np.array([199, 200, 250, 299, 300, 305])
    coverage = build_window_coverage(positions, window=w)
    # Window is [200, 300). Positions 199 / 300 / 305 are dropped.
    assert int(coverage.sum()) == 3  # 200, 250, 299
    assert coverage[0] == 1.0  # nt 200 -> index 0
    assert coverage[50] == 1.0  # nt 250 -> index 50
    assert coverage[99] == 1.0  # nt 299 -> index 99


# ---------------------------------------------------------------------------
# End-to-end synthetic build_fourier_spectrum_table
# ---------------------------------------------------------------------------


def _make_annotation(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["transcript", "sequence_name", "stop_codon", "l_tr"])


def _make_bed(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(
        rows, columns=["sample_name", "chrom", "read_length", "P_site"],
    )


def _build_periodic_bed(
    *, sample: str, chrom: str, read_length: int,
    window_lo: int, window_hi: int, peak: int = 5,
) -> list[dict]:
    """Return BED rows whose A-site lands every 3 nt across the window.

    A-site = P_site + 3, so we set P_site = window_pos - 3 to land the
    A-site exactly at window_pos.
    """
    rows = []
    for window_pos in range(window_lo, window_hi, 3):
        for _ in range(peak):
            rows.append({
                "sample_name": sample, "chrom": chrom,
                "read_length": read_length, "P_site": window_pos - 3,
            })
    return rows


def test_build_fourier_spectrum_table_returns_period3_peak_in_orf() -> None:
    ann = _make_annotation([
        {"transcript": "GENE1", "sequence_name": "GENE1",
         "stop_codon": 300, "l_tr": 500},
    ])
    # Use window_nt=99 (multiple of 3) so the period-3 energy lands in
    # one FFT bin instead of leaking across two neighbours. This is a
    # test-fixture choice; the default window remains 100 nt
    # (Wakigawa's published value), where some leakage is unavoidable
    # but the period-3 score still distinguishes ORF from 3' UTR.
    # ORF window = [201, 300); 3' UTR window = [303, 402).
    bed_rows = _build_periodic_bed(
        sample="WT", chrom="GENE1", read_length=32,
        window_lo=201, window_hi=300, peak=5,
    )
    bed = _make_bed(bed_rows)

    table = build_fourier_spectrum_table(
        bed, ann, samples=["WT"], window_nt=99, site="a",
    )

    orf = table[table["region"] == "orf"]
    assert not orf.empty
    nearest_3 = orf.iloc[(orf["period_nt"] - 3.0).abs().argsort()].iloc[0]
    other = orf[orf["period_nt"] != nearest_3["period_nt"]]
    assert nearest_3["amplitude"] > 10.0 * other["amplitude"].max()

    # 3' UTR has no reads -> no rows for that region (n_sites == 0 drop).
    utr3 = table[table["region"] == "utr3"]
    assert utr3.empty


def test_build_fourier_spectrum_table_stratifies_by_read_length() -> None:
    ann = _make_annotation([
        {"transcript": "GENE1", "sequence_name": "GENE1",
         "stop_codon": 300, "l_tr": 500},
    ])
    bed_rows = _build_periodic_bed(
        sample="WT", chrom="GENE1", read_length=32,
        window_lo=200, window_hi=300,
    )
    bed_rows += _build_periodic_bed(
        sample="WT", chrom="GENE1", read_length=29,
        window_lo=200, window_hi=300,
    )
    bed = _make_bed(bed_rows)

    table = build_fourier_spectrum_table(bed, ann, samples=["WT"])
    lengths = sorted(table["read_length"].unique())
    assert lengths == [29, 32]


def test_build_fourier_spectrum_table_flags_overlap_upstream_genes() -> None:
    ann = _make_annotation([
        {"transcript": "MT-ATP8", "sequence_name": "MT-ATP8",
         "stop_codon": 300, "l_tr": 500},
        {"transcript": "MT-CO1", "sequence_name": "MT-CO1",
         "stop_codon": 300, "l_tr": 500},
    ])
    bed_rows = _build_periodic_bed(
        sample="WT", chrom="MT-ATP8", read_length=32,
        window_lo=200, window_hi=300,
    )
    bed_rows += _build_periodic_bed(
        sample="WT", chrom="MT-CO1", read_length=32,
        window_lo=200, window_hi=300,
    )
    bed = _make_bed(bed_rows)
    table = build_fourier_spectrum_table(bed, ann, samples=["WT"])

    atp8_rows = table[table["gene"] == "MT-ATP8"]
    co1_rows = table[table["gene"] == "MT-CO1"]
    assert atp8_rows["is_overlap_upstream_orf"].all()
    assert not co1_rows["is_overlap_upstream_orf"].any()


# ---------------------------------------------------------------------------
# Per-gene scalar + combined summary
# ---------------------------------------------------------------------------


def test_period3_score_table_returns_high_score_for_clean_signal() -> None:
    ann = _make_annotation([
        {"transcript": "GENE1", "sequence_name": "GENE1",
         "stop_codon": 300, "l_tr": 500},
    ])
    bed_rows = _build_periodic_bed(
        sample="WT", chrom="GENE1", read_length=32,
        window_lo=200, window_hi=300, peak=5,
    )
    bed = _make_bed(bed_rows)
    spec = build_fourier_spectrum_table(bed, ann, samples=["WT"])
    scores = build_period3_score_table(spec)

    assert not scores.empty
    orf = scores[(scores["gene"] == "GENE1") & (scores["region"] == "orf")]
    assert len(orf) == 1
    # period3_score = period3_amplitude / mean(reference periods).
    # For a pure period-3 signal the reference bins are nearly zero, so
    # the score should be very large.
    assert float(orf["period3_score"].iloc[0]) > 5.0


def test_period3_summary_excludes_overlap_upstream_genes_from_combined() -> None:
    ann = _make_annotation([
        {"transcript": "MT-ATP8", "sequence_name": "MT-ATP8",
         "stop_codon": 300, "l_tr": 500},
        {"transcript": "MT-CO1", "sequence_name": "MT-CO1",
         "stop_codon": 300, "l_tr": 500},
        {"transcript": "MT-CO2", "sequence_name": "MT-CO2",
         "stop_codon": 300, "l_tr": 500},
    ])
    bed_rows: list[dict] = []
    for chrom in ("MT-ATP8", "MT-CO1", "MT-CO2"):
        bed_rows += _build_periodic_bed(
            sample="WT", chrom=chrom, read_length=32,
            window_lo=200, window_hi=300, peak=5,
        )
    bed = _make_bed(bed_rows)
    spec = build_fourier_spectrum_table(bed, ann, samples=["WT"])
    scores = build_period3_score_table(spec)
    summary = build_period3_summary_table(scores)

    orf_summary = summary[
        (summary["region"] == "orf") & (summary["sample"] == "WT")
    ]
    assert len(orf_summary) == 1
    row = orf_summary.iloc[0]
    # Combined: 2 genes (MT-CO1 + MT-CO2). Excluded: 1 (MT-ATP8).
    assert int(row["n_genes_combined"]) == 2
    assert int(row["n_genes_overlap_upstream"]) == 1
    assert float(row["median_period3_score_combined"]) > 5.0
    assert float(row["median_period3_score_overlap_upstream"]) > 5.0


def test_build_fourier_spectrum_raises_on_bad_bed_columns() -> None:
    bed = pd.DataFrame({"foo": [1]})
    ann = _make_annotation([
        {"transcript": "GENE1", "sequence_name": "GENE1",
         "stop_codon": 300, "l_tr": 500},
    ])
    with pytest.raises(ValueError, match="bed_with_psite_and_gene"):
        build_fourier_spectrum_table(bed, ann, samples=["WT"])


def test_build_fourier_spectrum_returns_empty_for_no_samples() -> None:
    ann = _make_annotation([
        {"transcript": "GENE1", "sequence_name": "GENE1",
         "stop_codon": 300, "l_tr": 500},
    ])
    bed = _make_bed([
        {"sample_name": "WT", "chrom": "GENE1",
         "read_length": 32, "P_site": 200},
    ])
    out = build_fourier_spectrum_table(bed, ann, samples=[])
    assert out.empty


# ---------------------------------------------------------------------------
# Plot writer (Job 1 sidecar contract)
# ---------------------------------------------------------------------------


def test_render_fourier_spectrum_panels_writes_combined_and_overlap(tmp_path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    from mitoribopy.plotting.fourier_spectrum_plots import (
        render_fourier_spectrum_panels,
    )
    from mitoribopy.plotting.figure_validator import metadata_sidecar_path

    ann = _make_annotation([
        {"transcript": "MT-ATP8", "sequence_name": "MT-ATP8",
         "stop_codon": 300, "l_tr": 500},
        {"transcript": "MT-CO1", "sequence_name": "MT-CO1",
         "stop_codon": 300, "l_tr": 500},
    ])
    bed_rows: list[dict] = []
    for chrom in ("MT-ATP8", "MT-CO1"):
        bed_rows += _build_periodic_bed(
            sample="WT", chrom=chrom, read_length=32,
            window_lo=201, window_hi=300, peak=5,
        )
    bed = _make_bed(bed_rows)
    spec = build_fourier_spectrum_table(
        bed, ann, samples=["WT"], window_nt=99,
    )
    out_dir = tmp_path / "fourier_spectrum"
    written = render_fourier_spectrum_panels(spec, output_dir=out_dir)

    combined_png = out_dir / "WT" / "WT_32nt_combined.png"
    overlap_png = out_dir / "WT" / "WT_32nt_overlap_upstream.png"
    assert combined_png in written
    assert overlap_png in written
    assert combined_png.is_file()
    assert overlap_png.is_file()

    # Per-plot sidecars from Job 1 contract.
    for png in (combined_png, overlap_png):
        sidecar = metadata_sidecar_path(png)
        assert sidecar.is_file()
        import json as _json
        meta = _json.loads(sidecar.read_text())
        for key in ("plot_type", "stage", "source_data",
                    "n_points_expected", "n_points_drawn"):
            assert key in meta, f"sidecar missing {key} in {sidecar}"
        assert meta["stage"] == "rpf"
        assert meta["source_data"] == "fourier_spectrum.tsv"

    # The combined plot must NOT carry MT-ATP8 traces; the overlap-
    # upstream plot must carry exactly MT-ATP8 (and nothing else).
    import json as _json
    combined_meta = _json.loads(metadata_sidecar_path(combined_png).read_text())
    overlap_meta = _json.loads(metadata_sidecar_path(overlap_png).read_text())
    assert combined_meta["n_genes"] == 1  # only MT-CO1
    assert overlap_meta["n_genes"] == 1  # only MT-ATP8
