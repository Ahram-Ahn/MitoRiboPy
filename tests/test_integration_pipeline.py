from __future__ import annotations

import itertools
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pandas as pd

from mitoribopy.analysis.offset_enrichment import create_csv_for_offset_enrichment
from mitoribopy.analysis.offset_selection import determine_p_site_offsets
from mitoribopy.analysis.translation_profile_analysis import run_translation_profile_analysis
from mitoribopy.data.transcript_annotations import human_annotation_df


def _all_codons() -> set[str]:
    return {"".join(parts) for parts in itertools.product("ATGC", repeat=3)}


def _human_annotation_with_codon_positions() -> pd.DataFrame:
    annotation = human_annotation_df.copy()
    annotation["start_codon"] = annotation["l_utr5"]
    annotation["stop_codon"] = annotation["l_tr"] - annotation["l_utr3"] - 3
    return annotation


def _build_tiny_offsets_bed_df() -> pd.DataFrame:
    # Majority reads are constructed to yield a deterministic RL30 P-site 5' offset of 13 at stop.
    rows: list[dict[str, int | str]] = []
    rows.extend({"chrom": "ND1", "start": 940, "end": 970, "read_length": 30} for _ in range(24))
    rows.extend({"chrom": "ND1", "start": 939, "end": 969, "read_length": 30} for _ in range(4))
    return pd.DataFrame(rows)


def _run_cli_smoke(
    *,
    fixture_dir: Path,
    output_dir: Path,
) -> subprocess.CompletedProcess[str]:
    input_dir = fixture_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)

    bed_lines = ["ND1\t940\t970" for _ in range(24)] + ["ND1\t939\t969" for _ in range(4)]
    (input_dir / "S1.bed").write_text("\n".join(bed_lines) + "\n", encoding="utf-8")

    # ND1 transcript length in built-in human annotation is 958 nt.
    fasta_path = fixture_dir / "tiny_human.fa"
    fasta_path.write_text(">ND1\n" + ("A" * 958) + "\n", encoding="utf-8")

    read_counts_path = fixture_dir / "read_counts.csv"
    read_counts_path.write_text("sample,read_count\nS1,28\n", encoding="utf-8")

    env = os.environ.copy()
    env["MPLCONFIGDIR"] = str(fixture_dir / "mplconfig")

    cmd = [
        sys.executable,
        "-m",
        "mitoribopy",
        "-s",
        "h",
        "-f",
        str(fasta_path),
        "--directory",
        str(input_dir),
        "-rpf",
        "30",
        "30",
        "-a",
        "stop",
        "--offset_type",
        "5",
        "--offset_site",
        "p",
        "--offset_pick_reference",
        "p_site",
        "--min_offset",
        "10",
        "--max_offset",
        "22",
        "--read_counts_file",
        str(read_counts_path),
        "--read_counts_sample_col",
        "sample",
        "--read_counts_reads_col",
        "read_count",
        "--rpm_norm_mode",
        "total",
        "--plot_format",
        "svg",
        "--output",
        str(output_dir),
        "-m",
    ]

    return subprocess.run(cmd, text=True, capture_output=True, check=False)


def test_offset_selection_psite_vs_asite_shifts_by_three_nt(tmp_path: Path) -> None:
    annotation = _human_annotation_with_codon_positions()
    bed_df = _build_tiny_offsets_bed_df()
    rpf_range = range(30, 31)

    p_summary_csv = tmp_path / "offset_stop_p.csv"
    _, p_offsets = create_csv_for_offset_enrichment(
        bed_df=bed_df,
        annotation_df=annotation,
        align_to="stop",
        rpf_range=rpf_range,
        output_csv=str(p_summary_csv),
        offset_limit=25,
        offset_site="p",
        codon_overlap_mode="full",
        strain="h",
    )
    assert p_offsets is not None

    p_selected = determine_p_site_offsets(
        offsets_df=p_offsets,
        align_to="stop",
        out_file=str(tmp_path / "p_site_offsets_stop_p.csv"),
        offset_min=10,
        offset_max=22,
        offset_site="p",
        selection_reference="p_site",
    )
    assert p_selected is not None

    a_summary_csv = tmp_path / "offset_stop_a.csv"
    _, a_offsets = create_csv_for_offset_enrichment(
        bed_df=bed_df,
        annotation_df=annotation,
        align_to="stop",
        rpf_range=rpf_range,
        output_csv=str(a_summary_csv),
        offset_limit=25,
        offset_site="a",
        codon_overlap_mode="full",
        strain="h",
    )
    assert a_offsets is not None

    a_selected = determine_p_site_offsets(
        offsets_df=a_offsets,
        align_to="stop",
        out_file=str(tmp_path / "p_site_offsets_stop_a.csv"),
        offset_min=10,
        offset_max=22,
        offset_site="a",
        selection_reference="p_site",
    )
    assert a_selected is not None

    merged = p_selected.merge(a_selected, on="Read Length", suffixes=("_p", "_a"))
    assert len(merged) == 1

    delta_5 = (
        merged["Most Enriched 5' Offset_a"] - merged["Most Enriched 5' Offset_p"]
    ).tolist()
    delta_3 = (
        merged["Most Enriched 3' Offset_a"] - merged["Most Enriched 3' Offset_p"]
    ).tolist()

    assert delta_5 == [3]
    assert delta_3 == [-3]


def test_offset_mask_excludes_near_anchor_bins_from_summary_and_selection(tmp_path: Path) -> None:
    annotation = _human_annotation_with_codon_positions()
    bed_df = pd.DataFrame(
        [
            {"chrom": "ND1", "start": 925, "end": 958, "read_length": 33},
            {"chrom": "ND1", "start": 925, "end": 960, "read_length": 35},
        ]
    )

    generated_summary, _ = create_csv_for_offset_enrichment(
        bed_df=bed_df,
        annotation_df=annotation,
        align_to="stop",
        rpf_range=range(33, 36),
        output_csv=str(tmp_path / "masked_offsets.csv"),
        offset_limit=10,
        offset_mask_nt=5,
        offset_site="p",
        codon_overlap_mode="full",
        strain="h",
    )

    assert generated_summary is not None
    row_33 = generated_summary[generated_summary["Read Length"] == 33].iloc[0]
    row_35 = generated_summary[generated_summary["Read Length"] == 35].iloc[0]
    assert row_33[4] == 0
    assert row_35[6] == 1

    offsets_df = pd.DataFrame(
        {
            "Read Length": [30] * 6,
            "5' Offset": [13, 13, 13, 13, 13, 13],
            "3' Offset": [4, 4, 4, 6, 6, 6],
        }
    )
    selected = determine_p_site_offsets(
        offsets_df=offsets_df,
        align_to="stop",
        out_file=str(tmp_path / "masked_selected.csv"),
        offset_min=1,
        offset_max=10,
        offset_mask_nt=5,
        offset_site="p",
        selection_reference="reported_site",
    )

    assert selected is not None
    assert selected["Most Enriched 3' Offset"].tolist() == [6]


def test_determine_p_site_offsets_supports_separate_five_and_three_ranges(tmp_path: Path) -> None:
    offsets_df = pd.DataFrame(
        {
            "Read Length": [30] * 6,
            "5' Offset": [8, 8, 8, 13, 13, 13],
            "3' Offset": [16, 16, 16, 16, 16, 16],
        }
    )

    selected = determine_p_site_offsets(
        offsets_df=offsets_df,
        align_to="stop",
        out_file=str(tmp_path / "separate_ranges.csv"),
        five_offset_min=12,
        five_offset_max=20,
        three_offset_min=15,
        three_offset_max=20,
        offset_site="p",
        selection_reference="reported_site",
    )

    assert selected is not None
    assert selected["Most Enriched 5' Offset"].tolist() == [13]
    assert selected["Most Enriched 3' Offset"].tolist() == [16]


def test_determine_p_site_offsets_respects_requested_output_range_for_asite_reporting(
    tmp_path: Path,
) -> None:
    offsets_df = pd.DataFrame(
        {
            "Read Length": [30] * 6,
            "5' Offset": [21, 21, 21, 18, 18, 18],
            "3' Offset": [15, 15, 15, 12, 12, 12],
        }
    )

    selected = determine_p_site_offsets(
        offsets_df=offsets_df,
        align_to="stop",
        out_file=str(tmp_path / "asite_ranges.csv"),
        five_offset_min=14,
        five_offset_max=20,
        three_offset_min=10,
        three_offset_max=20,
        offset_site="a",
        selection_reference="p_site",
    )

    assert selected is not None
    assert selected["Most Enriched 5' Offset"].tolist() == [18]
    assert selected["Most Enriched 3' Offset"].tolist() == [12]


def test_translation_profile_analysis_writes_transcript_level_asite_stop_codon_usage(
    tmp_path: Path,
) -> None:
    fasta_path = tmp_path / "tiny_tx.fa"
    tx_seq = ("A" * 9) + "ATG" + "AAA" + "AAA" + "AGA" + "TTT"
    fasta_path.write_text(f">TX1\n{tx_seq}\n", encoding="utf-8")

    annotation_df = pd.DataFrame(
        [
            {
                "transcript": "TX1",
                "start_codon": 9,
                "stop_codon": 18,
            }
        ]
    )
    filtered_bed_df = pd.DataFrame(
        [
            {
                "chrom": "TX1",
                "start": 3,
                "end": 33,
                "read_length": 30,
                "sample_name": "S1",
            }
        ]
    )
    offsets_df = pd.DataFrame(
        {
            "Read Length": [30],
            "Most Enriched 5' Offset": [13],
            "Most Enriched 3' Offset": [17],
        }
    )
    args = SimpleNamespace(
        codon_density_window=False,
        strain="h",
        cap_percentile=0.999,
        offset_site="p",
    )
    output_dir = tmp_path / "analysis"

    run_translation_profile_analysis(
        sample_dirs=["S1"],
        selected_offsets_df=offsets_df,
        offset_type="5",
        fasta_file=str(fasta_path),
        output_dir=str(output_dir),
        args=args,
        annotation_df=annotation_df,
        filtered_bed_df=filtered_bed_df,
        requested_sites=["p", "a"],
    )

    # v0.4.x flat layout: filenames consistently prefixed with the
    # site, regardless of which sites were requested.
    a_site_tx_csv = output_dir / "S1" / "codon_usage" / "a_site_codon_usage_TX1.csv"
    a_site_total_csv = output_dir / "S1" / "codon_usage" / "a_site_codon_usage_total.csv"

    assert a_site_tx_csv.exists()
    assert a_site_total_csv.exists()

    transcript_df = pd.read_csv(a_site_tx_csv)
    total_df = pd.read_csv(a_site_total_csv)

    tx_stop = transcript_df[(transcript_df["Codon"] == "AGA") & (transcript_df["AA"] == "*")]
    total_stop = total_df[(total_df["Codon"] == "AGA") & (total_df["AA"] == "*")]

    assert not tx_stop.empty
    assert not total_stop.empty
    assert tx_stop["Coverage"].iloc[0] > 0
    assert total_stop["Coverage"].iloc[0] >= tx_stop["Coverage"].iloc[0]


def test_end_to_end_cli_smoke_generates_codon_usage_outputs(tmp_path: Path) -> None:
    fixture_dir = tmp_path / "fixture"
    output_dir = tmp_path / "out"

    result = _run_cli_smoke(fixture_dir=fixture_dir, output_dir=output_dir)
    assert result.returncode == 0, (
        "CLI smoke run failed\n"
        f"STDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )

    for step_number in range(1, 8):
        assert f"[PIPELINE] Step {step_number}/7" in result.stdout

    # v0.4.x flat layout: offset CSVs land under
    # offset_diagnostics/csv/, plots under offset_diagnostics/plots/,
    # and translation_profile/coverage_profile_plots no longer use
    # p/ or a/ subfolders.
    offsets_csv = output_dir / "offset_diagnostics" / "csv" / "p_site_offsets_stop.csv"
    footprint_csv = (
        output_dir / "translation_profile" / "S1" / "footprint_density" / "ND1_footprint_density.csv"
    )
    codon_usage_csv = (
        output_dir / "translation_profile" / "S1" / "codon_usage" / "p_site_codon_usage_total.csv"
    )
    coverage_root = output_dir / "coverage_profile_plots"
    filtered_bed_dir = output_dir / "offset_diagnostics" / "filtered_bed"
    log_file = output_dir / "mitoribopy.log"

    assert offsets_csv.exists()
    assert footprint_csv.exists()
    assert codon_usage_csv.exists()
    # Read coverage outputs (default ON for both raw and rpm).
    assert (coverage_root / "read_coverage_rpm").exists()
    # Site-prefixed density plots, no p/ or a/ subdirs.
    assert (coverage_root / "p_site_density_rpm").exists()
    assert (coverage_root / "a_site_density_rpm").exists()
    # A-site translation_profile output sits in the same flat sample dir.
    assert (
        output_dir / "translation_profile" / "S1" / "codon_usage" / "a_site_codon_usage_total.csv"
    ).exists()
    assert not filtered_bed_dir.exists()
    assert list(coverage_root.rglob("*.svg")), "Expected at least one coverage-profile SVG output"
    # E_site column dropped from the CSV (v0.4.x cleanup); A_site comes
    # before P_site, no *_selected_depth tail column.
    fp_cols = pd.read_csv(footprint_csv, nrows=1).columns.tolist()
    assert fp_cols == ["Position", "Nucleotide", "A_site", "P_site"]
    assert log_file.exists()
    assert "[PIPELINE] Step 7/7" in log_file.read_text(encoding="utf-8")

    codon_df = pd.read_csv(codon_usage_csv)
    codon_col = next((col for col in codon_df.columns if col.lower() == "codon"), None)
    assert codon_col is not None
    assert len(codon_df) == 64
    assert set(codon_df[codon_col].astype(str).str.upper()) == _all_codons()
