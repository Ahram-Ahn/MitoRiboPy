"""Refinements to the 3-nt periodicity QC bundle.

Covers:
* Depth-aware ``best_read_length`` selection in :func:`build_qc_summary`.
* Correct ``best_read_length_dominant_fraction`` (max of frame fractions).
* ``exclude_start_codons`` / ``exclude_stop_codons`` plumbed through
  :func:`compute_frame_summary` and :func:`compute_frame_summary_by_length`.
* Human mt-mRNA overlap-pair annotation in :func:`build_gene_periodicity`.
"""

from __future__ import annotations

import pandas as pd
import pytest

from mitoribopy.analysis.periodicity import (
    compute_frame_summary,
    compute_frame_summary_by_length,
)
from mitoribopy.analysis.periodicity_qc import (
    HUMAN_MT_OVERLAPPING_GENES,
    build_gene_periodicity,
    build_qc_summary,
    is_known_overlap_gene,
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


# ---------- build_qc_summary refinements -----------------------------------


def _by_length_row(
    *,
    sample: str,
    read_length: int,
    n: int,
    f0: float,
    f1: float,
    f2: float,
    qc_call: str,
) -> dict:
    return {
        "sample_id": sample,
        "read_length": read_length,
        "n_reads_total": n,
        "n_reads_cds": n,
        "frame0_fraction": f0,
        "frame1_fraction": f1,
        "frame2_fraction": f2,
        "dominant_frame": int(max(range(3), key=lambda i: (f0, f1, f2)[i])),
        "frame0_dominance": f0 - max(f1, f2),
        "periodicity_score": f0 - max(f1, f2),
        "frame_entropy": 0.0,
        "include_for_downstream": True,
        "exclusion_reason": "none",
        "expected_frame": 0,
        "expected_frame_fraction": f0,
        "expected_frame_enrichment": 0.0,
        "entropy_bias": 0.0,
        "qc_call": qc_call,
    }


def test_qc_summary_prefers_deep_length_over_shallow_lucky_one() -> None:
    """A 5-read length with frame0=1.0 must NOT win over a 5000-read length
    with frame0=0.75 once the depth threshold is honored."""
    df = pd.DataFrame([
        _by_length_row(
            sample="S1", read_length=31, n=5_000,
            f0=0.75, f1=0.15, f2=0.10, qc_call="good",
        ),
        _by_length_row(
            sample="S1", read_length=22, n=5,
            f0=1.00, f1=0.00, f2=0.00, qc_call="low_depth",
        ),
    ])
    summary = build_qc_summary(df)
    assert int(summary.iloc[0]["best_read_length"]) == 31
    assert summary.iloc[0]["best_read_length_expected_frame_fraction"] == pytest.approx(0.75)


def test_qc_summary_dominant_fraction_is_max_of_frame_fractions() -> None:
    """When dominant_frame != expected_frame, dominant_fraction must be the
    actual max(frame fraction), not zero (the previous and/or-ternary bug
    could collapse to the wrong value in this case)."""
    df = pd.DataFrame([
        _by_length_row(
            sample="S1", read_length=31, n=5_000,
            f0=0.20, f1=0.70, f2=0.10, qc_call="poor",
        ),
    ])
    summary = build_qc_summary(df)
    row = summary.iloc[0]
    assert int(row["best_read_length_dominant_frame"]) == 1
    assert float(row["best_read_length_dominant_fraction"]) == pytest.approx(0.70)
    # And expected-frame fraction is still f0 (not the max).
    assert float(row["best_read_length_expected_frame_fraction"]) == pytest.approx(0.20)


def test_qc_summary_falls_back_to_all_lengths_when_none_clear_depth() -> None:
    """If every length is below min_reads_per_length, the best-length
    pick must still produce a row instead of crashing on an empty pool."""
    df = pd.DataFrame([
        _by_length_row(
            sample="S1", read_length=22, n=10,
            f0=0.50, f1=0.30, f2=0.20, qc_call="low_depth",
        ),
        _by_length_row(
            sample="S1", read_length=31, n=20,
            f0=0.55, f1=0.30, f2=0.15, qc_call="low_depth",
        ),
    ])
    summary = build_qc_summary(df)
    assert len(summary) == 1
    assert int(summary.iloc[0]["best_read_length"]) == 31  # higher f0


# ---------- exclude_start_codons / exclude_stop_codons ---------------------


def _row(p: int, sample: str = "S1", chrom: str = "ND1", read_length: int = 30) -> dict:
    return {"sample_name": sample, "chrom": chrom, "P_site": p,
            "read_length": read_length, "start": 0, "end": 0}


def test_compute_frame_summary_honors_exclude_start_codons() -> None:
    """A burst of P-sites in the first 6 codons should be masked when
    exclude_start_codons=6 and visible otherwise."""
    ann = _annotation("ND1", l_utr5=0, l_cds=300, l_utr3=10)
    initiation_burst = [_row(p) for p in (1, 4, 7)]  # first 6 codons, frame 1
    clean_codons = [_row(p) for p in (18, 21, 24, 27, 30, 33)]  # past mask, frame 0
    psite = pd.DataFrame(initiation_burst + clean_codons)
    masked = compute_frame_summary(
        psite, ann, sample="S1",
        exclude_start_codons=6, exclude_stop_codons=0,
    )
    assert masked.n_reads == 6
    assert masked.frame_0 == pytest.approx(1.0)

    unmasked = compute_frame_summary(psite, ann, sample="S1")
    assert unmasked.n_reads == 9  # both bursts counted
    assert unmasked.frame_0 < 1.0


def test_compute_frame_summary_by_length_honors_codon_edges() -> None:
    ann = _annotation("ND1", l_utr5=0, l_cds=300, l_utr3=10)
    # 6 reads at frame-1 in the masked initiation zone, 6 reads at frame-0
    # past it (within the same length class).
    masked_zone = [_row(p) for p in (1, 4, 7, 10, 13, 16)]
    clean_zone = [_row(p) for p in (18, 21, 24, 27, 30, 33)]
    psite = pd.DataFrame(masked_zone + clean_zone)
    masked = compute_frame_summary_by_length(
        psite, ann, sample="S1",
        exclude_start_codons=6, exclude_stop_codons=0,
        min_cds_reads_per_length=1,  # don't reject for low count in this test
    )
    row = masked[masked["read_length"] == 30].iloc[0]
    assert int(row["n_reads_cds"]) == 6  # only the 6 post-mask reads
    assert float(row["frame0_fraction"]) == pytest.approx(1.0)


# ---------- overlap-pair helper --------------------------------------------


def test_overlap_set_covers_known_human_mt_mrna_pairs() -> None:
    for s in ("MT-ATP8", "MT-ATP6", "MT-ND4L", "MT-ND4"):
        assert is_known_overlap_gene(s)
    # Common spelling variants
    assert is_known_overlap_gene("mtatp8")
    assert is_known_overlap_gene("MT_ATP6")
    # Fused-overlap transcript names used by some FASTAs
    # (e.g. input_data/human-mt-mRNA.fasta in this repo).
    assert is_known_overlap_gene("ATP86")
    assert is_known_overlap_gene("ND4L4")
    assert is_known_overlap_gene("atp86")
    # Non-overlap gene
    assert not is_known_overlap_gene("MT-CO1")
    assert not is_known_overlap_gene("COX1")
    assert not is_known_overlap_gene(None)
    # The constant is a frozenset and contains the canonical names.
    assert "MT-ATP8" in HUMAN_MT_OVERLAPPING_GENES
    assert "ATP86" in HUMAN_MT_OVERLAPPING_GENES
    assert isinstance(HUMAN_MT_OVERLAPPING_GENES, frozenset)


def test_build_gene_periodicity_marks_overlap_pairs() -> None:
    """When annotate_overlap=True, MT-ATP8 should be flagged is_overlap_pair=True
    while MT-CO1 should not."""
    ann = pd.concat([
        _annotation("MT-ATP8", l_utr5=0, l_cds=210, l_utr3=0),
        _annotation("MT-CO1", l_utr5=0, l_cds=300, l_utr3=0),
    ], ignore_index=True)
    bed = pd.DataFrame([
        {"sample_name": "S1", "chrom": "MT-ATP8", "P_site": p, "read_length": 30}
        for p in range(0, 60, 3)
    ] + [
        {"sample_name": "S1", "chrom": "MT-CO1", "P_site": p, "read_length": 30}
        for p in range(0, 60, 3)
    ])
    table = build_gene_periodicity(
        bed, ann, samples=["S1"], annotate_overlap=True,
    )
    assert "is_overlap_pair" in table.columns
    atp8 = table[table["gene"] == "MT-ATP8"].iloc[0]
    co1 = table[table["gene"] == "MT-CO1"].iloc[0]
    assert bool(atp8["is_overlap_pair"]) is True
    assert bool(co1["is_overlap_pair"]) is False


# ---------- spec outputs added in 2026-05-02 refactor ----------------------


def _bed_for(sample: str, transcript: str, p_sites: list[int],
             read_length: int = 30) -> pd.DataFrame:
    return pd.DataFrame([
        {"sample_name": sample, "chrom": transcript, "P_site": int(p),
         "read_length": int(read_length)}
        for p in p_sites
    ])


def test_build_frame_counts_by_gene_schema_and_values() -> None:
    """Spec output frame_counts_by_gene.tsv: per-(sample, gene, read_length)
    frame counts. Distinct from gene_periodicity.tsv which pools across
    read lengths."""
    from mitoribopy.analysis.periodicity_qc import build_frame_counts_by_gene

    ann = _annotation("ND1", l_utr5=0, l_cds=300, l_utr3=10)
    # Two read lengths on the same gene; length 30 is frame-0 dominant,
    # length 32 is uniform.
    bed = pd.concat([
        _bed_for("S1", "ND1", [21, 24, 27, 30, 33, 36], read_length=30),
        _bed_for("S1", "ND1", [22, 23, 25, 26, 28, 29], read_length=32),
    ], ignore_index=True)
    table = build_frame_counts_by_gene(
        bed, ann, samples=["S1"], annotate_overlap=True, site_type="p",
    )
    # Required spec columns.
    must_have = {
        "sample", "gene", "transcript_id", "site_type", "read_length",
        "n_frame0", "n_frame1", "n_frame2", "n_total",
        "frame0_fraction", "frame1_fraction", "frame2_fraction",
        "expected_frame", "expected_frame_fraction",
        "expected_frame_enrichment", "dominant_frame",
        "dominant_frame_fraction", "entropy_bias", "qc_call",
        "is_overlap_pair",
    }
    assert must_have.issubset(set(table.columns)), \
        f"missing: {must_have - set(table.columns)}"
    # One row per (sample, gene, read_length).
    assert len(table) == 2
    rl30 = table[table["read_length"] == 30].iloc[0]
    rl32 = table[table["read_length"] == 32].iloc[0]
    assert int(rl30["n_total"]) == 6
    assert float(rl30["frame0_fraction"]) == pytest.approx(1.0)
    assert int(rl32["n_total"]) == 6
    assert float(rl32["frame0_fraction"]) < 0.5  # uniform spread
    # site_type stamped through.
    assert all(table["site_type"] == "p")


def test_fourier_spectrum_replaces_legacy_fft_period3_table() -> None:
    """The Wakigawa-style spectrum must score a strong period-3 signal
    cleanly. Detailed coverage (window math, overlap-upstream policy,
    plot writers) is in tests/test_fourier_spectrum.py — this is the
    smoke test that the legacy build_fft_period3_table replacement is
    importable and returns a populated table for a perfect signal.
    """
    from mitoribopy.analysis.fourier_spectrum import (
        build_fourier_spectrum_table,
        build_period3_score_table,
    )

    # Use an upstream-of-stop window (the Wakigawa convention).
    ann = pd.DataFrame([
        {"transcript": "LONG", "sequence_name": "LONG",
         "stop_codon": 300, "l_tr": 500},
    ])
    # A-site = P_site + 3, so we want A-site every 3 nt in [201, 300).
    psite_positions = list(range(198, 300, 3))
    bed = pd.DataFrame([
        {"sample_name": "S1", "chrom": "LONG",
         "read_length": 30, "P_site": p}
        for p in psite_positions
    ])
    spec = build_fourier_spectrum_table(
        bed, ann, samples=["S1"], window_nt=99, site="a",
    )
    assert not spec.empty
    score = build_period3_score_table(spec)
    long_row = score[(score["gene"] == "LONG") & (score["region"] == "orf")]
    assert len(long_row) == 1
    assert float(long_row["period3_score"].iloc[0]) > 5.0


def test_run_periodicity_qc_writes_every_spec_output(tmp_path) -> None:
    """End-to-end: run_periodicity_qc on a small synthetic fixture must
    produce every TSV + plot the spec asks for, alongside the legacy
    artefacts."""
    from mitoribopy.analysis.periodicity import (
        compute_p_site_positions, run_periodicity_qc,
    )

    ann = _annotation("ND1", l_utr5=0, l_cds=300, l_utr3=10)
    # 60 reads on length 30, P-sites at 30, 33, ... (frame 0) past the
    # 6-codon mask zone, plus a few in the masked zone that should be
    # dropped by the spec defaults.
    starts = list(range(30, 270, 3))
    bed = pd.DataFrame([
        {"sample_name": "S1", "chrom": "ND1",
         "start": p - 12, "end": p - 12 + 30,
         "read_length": 30, "strand": "+"}
        for p in starts
    ])
    offsets = pd.DataFrame([{
        "Read Length": 30,
        "Most Enriched 5' Offset": 12,
        "Most Enriched 3' Offset": 30 - 12 - 3,
    }])
    out = tmp_path / "qc"
    run_periodicity_qc(
        bed_df=bed,
        annotation_df=ann,
        samples=["S1"],
        selected_offsets_by_sample={"S1": offsets},
        selected_offsets_combined=None,
        offset_type="5",
        offset_site="p",
        output_dir=out,
        plot=True,
        min_cds_reads_per_length=1,  # tiny fixture
        compute_fourier_spectrum=True,
        fourier_window_nt=99,  # multiple-of-3 -> exact period-3 bin
        fourier_render_plots=False,  # fixture is too sparse for an overlay plot
        metagene_nt=60,
    )
    spec_outputs = [
        # TSVs
        "frame_counts_by_sample_length.tsv",
        "frame_counts_by_gene.tsv",
        "gene_periodicity.tsv",
        "metagene_start.tsv",
        "metagene_stop.tsv",
        "fourier_spectrum.tsv",
        "fourier_period3_score.tsv",
        "fourier_period3_summary.tsv",
        # Plots
        "frame_fraction_heatmap.svg",
        "read_length_periodicity_barplot.svg",
        "metagene_start_p_site.svg",
        "metagene_stop_p_site.svg",
        "gene_phase_score_dotplot.svg",
        # QC summaries
        "qc_summary.tsv",
        "qc_summary.md",
    ]
    missing = [f for f in spec_outputs if not (out / f).is_file()]
    assert not missing, f"spec outputs missing: {missing}"


def test_periodicity_metadata_records_new_fields(tmp_path) -> None:
    """periodicity.metadata.json must record exclude_start_codons,
    exclude_stop_codons, phase_score_enabled, fourier_spectrum_enabled,
    fourier_window_nt, metagene_nt — so a reviewer can re-derive every
    number from disk."""
    import json
    from mitoribopy.analysis.periodicity import run_periodicity_qc

    ann = _annotation("ND1", l_utr5=0, l_cds=300, l_utr3=10)
    bed = pd.DataFrame([
        {"sample_name": "S1", "chrom": "ND1", "start": 18, "end": 48,
         "read_length": 30, "strand": "+"}
    ])
    offsets = pd.DataFrame([{
        "Read Length": 30,
        "Most Enriched 5' Offset": 12,
        "Most Enriched 3' Offset": 30 - 12 - 3,
    }])
    out = tmp_path / "qc"
    run_periodicity_qc(
        bed_df=bed, annotation_df=ann, samples=["S1"],
        selected_offsets_by_sample={"S1": offsets},
        selected_offsets_combined=None,
        offset_type="5", offset_site="p", output_dir=out,
        plot=False, exclude_start_codons=6, exclude_stop_codons=3,
        compute_phase_score=True, compute_fourier_spectrum=True,
        fourier_window_nt=99, fourier_render_plots=False, metagene_nt=90,
    )
    meta = json.loads((out / "by_length" / "periodicity.metadata.json").read_text())
    assert meta["exclude_start_codons"] == 6
    assert meta["exclude_stop_codons"] == 3
    assert meta["phase_score_enabled"] is True
    assert meta["fourier_spectrum_enabled"] is True
    assert meta["fourier_window_nt"] == 99
    assert meta["metagene_nt"] == 90
