"""Unit tests for ``mitoribopy.io.bam_reader`` and the BAM-input path
through ``mitoribopy.io.bed_reader.process_bed_files``.

BAM fixtures are built in-memory via pysam so the tests run without any
bioinformatics tools installed.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pysam
import pytest

from mitoribopy.data import load_annotation_table
from mitoribopy.io.bam_reader import convert_bam_to_bed, prepare_bam_inputs
from mitoribopy.io.bed_reader import (
    compute_unfiltered_read_length_summary,
    process_bed_files,
)


# ---------- helpers ---------------------------------------------------------


def _make_bam(
    path: Path,
    records: list[dict],
    references: list[tuple[str, int]] | None = None,
) -> None:
    """Write a tiny BAM file from plain-dict records."""
    if references is None:
        references = [("ND1", 1000), ("COX1", 1000)]
    ref_names = [name for name, _length in references]
    ref_lengths = [length for _name, length in references]

    header = pysam.AlignmentHeader.from_references(ref_names, ref_lengths)

    with pysam.AlignmentFile(str(path), "wb", header=header) as sink:
        for rec in records:
            read = pysam.AlignedSegment(header)
            read.query_name = rec.get("name", "r")
            read.query_sequence = "A" * max(
                1, int(rec.get("end", 1) - rec.get("start", 0))
            )
            read.query_qualities = pysam.qualitystring_to_array(
                "I" * len(read.query_sequence)
            )
            ref = rec.get("ref")
            read.reference_id = ref_names.index(ref) if ref is not None else -1
            read.reference_start = int(rec.get("start", 0))
            read.cigarstring = f"{int(rec['end']) - int(rec['start'])}M"
            read.mapping_quality = int(rec.get("mapq", 30))

            flag = 0
            if rec.get("unmapped"):
                flag |= 0x4
                read.reference_id = -1
            if rec.get("reverse"):
                flag |= 0x10
            if rec.get("secondary"):
                flag |= 0x100
            if rec.get("supplementary"):
                flag |= 0x800
            read.flag = flag

            sink.write(read)


# ---------- convert_bam_to_bed ----------------------------------------------


def test_convert_bam_to_bed_writes_bed6_with_strand(tmp_path) -> None:
    bam = tmp_path / "sample.bam"
    _make_bam(
        bam,
        [
            {"name": "fwd", "ref": "ND1", "start": 10, "end": 40, "mapq": 30},
            {"name": "rev", "ref": "COX1", "start": 100, "end": 131, "mapq": 25, "reverse": True},
        ],
    )
    bed = tmp_path / "sample.bed"

    rows = convert_bam_to_bed(bam, bed, mapq_threshold=0)

    assert rows == 2
    lines = bed.read_text().splitlines()
    assert len(lines) == 2
    fwd_line, rev_line = (line.split("\t") for line in lines)
    assert fwd_line[5] == "+"
    assert rev_line[5] == "-"


def test_convert_bam_to_bed_applies_mapq_filter(tmp_path) -> None:
    bam = tmp_path / "sample.bam"
    _make_bam(
        bam,
        [
            {"name": "low", "ref": "ND1", "start": 0, "end": 30, "mapq": 5},
            {"name": "border", "ref": "ND1", "start": 50, "end": 80, "mapq": 10},
            {"name": "high", "ref": "COX1", "start": 0, "end": 30, "mapq": 42},
        ],
    )
    bed = tmp_path / "sample.bed"

    rows = convert_bam_to_bed(bam, bed, mapq_threshold=10)

    assert rows == 2
    names = sorted(line.split("\t")[3] for line in bed.read_text().splitlines())
    assert names == ["border", "high"]


# ---------- prepare_bam_inputs ----------------------------------------------


def test_prepare_bam_inputs_converts_every_bam(tmp_path) -> None:
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _make_bam(in_dir / "sampleA.bam", [{"name": "r1", "ref": "ND1", "start": 0, "end": 30}])
    _make_bam(in_dir / "sampleB.bam", [{"name": "r2", "ref": "COX1", "start": 0, "end": 30}])

    out = prepare_bam_inputs(in_dir, tmp_path / "converted")

    assert sorted(p.name for p in out) == ["sampleA.bed", "sampleB.bed"]
    for bed in out:
        assert bed.read_text().strip()


def test_prepare_bam_inputs_returns_empty_when_no_bams(tmp_path) -> None:
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    (in_dir / "sampleA.bed").write_text("ND1\t0\t30\n")

    out = prepare_bam_inputs(in_dir, tmp_path / "converted")

    assert out == []
    assert not (tmp_path / "converted").exists()  # not created when nothing to convert


def test_prepare_bam_inputs_skips_bam_when_same_name_bed_exists(tmp_path, caplog) -> None:
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _make_bam(in_dir / "sampleA.bam", [{"name": "r1", "ref": "ND1", "start": 0, "end": 30}])
    (in_dir / "sampleA.bed").write_text("ND1\t0\t30\n")

    out = prepare_bam_inputs(in_dir, tmp_path / "converted")

    # BAM was skipped because the explicit BED exists.
    assert out == []


# ---------- process_bed_files with BAM inputs -------------------------------


def _load_human_annotation_for_rpf() -> pd.DataFrame:
    return load_annotation_table(preset="h")


def test_process_bed_files_ingests_bam_only_directory(tmp_path) -> None:
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    # Reads at length 29 pass the default human RPF range (28-34).
    _make_bam(
        in_dir / "sampleA.bam",
        [{"name": f"r{i}", "ref": "MT-ND1", "start": i * 5, "end": i * 5 + 29}
         for i in range(12)],
        references=[("MT-ND1", 1000)],
    )

    out_dir = tmp_path / "out"
    out_dir.mkdir()

    df, samples = process_bed_files(
        input_dir=str(in_dir),
        output_dir=str(out_dir),
        organism="h",
        annotation_df=_load_human_annotation_for_rpf(),
        rpf_range=range(28, 35),
    )

    assert samples == ["sampleA"]
    assert not df.empty
    # The BAM converted BED lives under <output>/bam_converted/
    assert (out_dir / "bam_converted" / "sampleA.bed").exists()
    # Strand column is preserved end-to-end.
    assert "strand" in df.columns


def test_process_bed_files_ingests_mixed_bed_and_bam_directory(tmp_path) -> None:
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _make_bam(
        in_dir / "sampleA.bam",
        [{"name": f"r{i}", "ref": "MT-ND1", "start": i * 5, "end": i * 5 + 29}
         for i in range(8)],
        references=[("MT-ND1", 1000)],
    )
    (in_dir / "sampleB.bed").write_text(
        "\n".join(
            f"MT-COX1\t{i * 5}\t{i * 5 + 29}" for i in range(6)
        ) + "\n"
    )

    out_dir = tmp_path / "out"
    out_dir.mkdir()

    df, samples = process_bed_files(
        input_dir=str(in_dir),
        output_dir=str(out_dir),
        organism="h",
        annotation_df=_load_human_annotation_for_rpf(),
        rpf_range=range(28, 35),
    )

    assert sorted(samples) == ["sampleA", "sampleB"]
    # Both samples end up in the filtered frame.
    assert set(df["sample_name"].unique()) == {"sampleA", "sampleB"}


def test_process_bed_files_prefers_bed_over_bam_on_name_conflict(tmp_path) -> None:
    in_dir = tmp_path / "in"
    in_dir.mkdir()

    # Native BED with 5 reads of length 29.
    (in_dir / "sampleA.bed").write_text(
        "\n".join(
            f"MT-ND1\t{i * 5}\t{i * 5 + 29}" for i in range(5)
        ) + "\n"
    )
    # BAM with the SAME sample name but 10 reads. Should be skipped.
    _make_bam(
        in_dir / "sampleA.bam",
        [{"name": f"r{i}", "ref": "MT-ND1", "start": i * 5, "end": i * 5 + 29}
         for i in range(10)],
        references=[("MT-ND1", 1000)],
    )

    out_dir = tmp_path / "out"
    out_dir.mkdir()

    df, samples = process_bed_files(
        input_dir=str(in_dir),
        output_dir=str(out_dir),
        organism="h",
        annotation_df=_load_human_annotation_for_rpf(),
        rpf_range=range(28, 35),
    )

    # Only the 5 rows from the native BED should be ingested.
    assert samples == ["sampleA"]
    assert len(df) == 5
    # No BAM-converted file should have been written for sampleA.
    assert not (out_dir / "bam_converted" / "sampleA.bed").exists()


def test_process_bed_files_bam_mapq_filter_drops_low_quality_reads(tmp_path) -> None:
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _make_bam(
        in_dir / "sampleA.bam",
        [
            {"name": "keep", "ref": "MT-ND1", "start": 0, "end": 29, "mapq": 30},
            {"name": "drop_low_mapq", "ref": "MT-ND1", "start": 40, "end": 69, "mapq": 3},
        ],
        references=[("MT-ND1", 1000)],
    )

    out_dir = tmp_path / "out"
    out_dir.mkdir()

    df, _ = process_bed_files(
        input_dir=str(in_dir),
        output_dir=str(out_dir),
        organism="h",
        annotation_df=_load_human_annotation_for_rpf(),
        rpf_range=range(28, 35),
        bam_mapq=10,
    )

    assert len(df) == 1
    assert df.iloc[0]["name"] == "keep"


def test_unfiltered_read_length_qc_can_use_bam_converted_beds(tmp_path) -> None:
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    _make_bam(
        in_dir / "sampleA.bam",
        [{"name": f"r{i}", "ref": "MT-ND1", "start": 0, "end": 29 + i}
         for i in range(3)],
        references=[("MT-ND1", 1000)],
    )
    converted = prepare_bam_inputs(in_dir, tmp_path / "converted")
    out_csv = tmp_path / "unfiltered.csv"

    compute_unfiltered_read_length_summary(
        input_dir=str(in_dir),
        output_csv=str(out_csv),
        total_counts_map={},
        read_length_range=(15, 50),
        converted_bed_paths=converted,
    )

    df = pd.read_csv(out_csv)
    assert sorted(df["sample_name"].unique()) == ["sampleA"]
    assert set(df["read_length"]) == {29, 30, 31}


def test_process_bed_files_empty_dir_returns_empty_frame(tmp_path) -> None:
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    df, samples = process_bed_files(
        input_dir=str(in_dir),
        output_dir=str(out_dir),
        organism="h",
        annotation_df=_load_human_annotation_for_rpf(),
        rpf_range=range(28, 35),
    )

    assert df.empty
    assert samples == []
