from __future__ import annotations

import itertools
import os
from pathlib import Path
import subprocess
import sys

import pandas as pd

from mitoribopy.analysis.offset_enrichment import create_csv_for_offset_enrichment
from mitoribopy.analysis.offset_selection import determine_p_site_offsets
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


def test_end_to_end_cli_smoke_generates_codon_usage_outputs(tmp_path: Path) -> None:
    fixture_dir = tmp_path / "fixture"
    output_dir = tmp_path / "out"

    result = _run_cli_smoke(fixture_dir=fixture_dir, output_dir=output_dir)
    assert result.returncode == 0, (
        "CLI smoke run failed\n"
        f"STDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )

    offsets_csv = output_dir / "plots_and_csv" / "p_site_offsets_stop.csv"
    footprint_csv = output_dir / "S1" / "footprint_density" / "ND1_footprint_density.csv"
    codon_usage_csv = output_dir / "S1" / "codon_usage" / "codon_usage_total.csv"
    igv_plot_dir = output_dir / "igv_style_plots"

    assert offsets_csv.exists()
    assert footprint_csv.exists()
    assert codon_usage_csv.exists()
    assert igv_plot_dir.exists()
    assert list(igv_plot_dir.rglob("*.svg")), "Expected at least one IGV-style SVG output"

    codon_df = pd.read_csv(codon_usage_csv)
    codon_col = next((col for col in codon_df.columns if col.lower() == "codon"), None)
    assert codon_col is not None
    assert len(codon_df) == 64
    assert set(codon_df[codon_col].astype(str).str.upper()) == _all_codons()
