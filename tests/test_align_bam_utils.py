"""Unit tests for ``mitoribopy.align.bam_utils``.

These build real (tiny) BAM files in-memory via pysam so we verify
behavior on actual htslib records rather than mocking pysam. Each test
creates a fresh temporary BAM in ``tmp_path`` so the tests are isolated
and deterministic.
"""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from mitoribopy.align.bam_utils import (
    bam_to_bed6,
    count_mapped_reads,
    filter_bam_mapq,
)


# ---------- helpers ---------------------------------------------------------


def _make_bam(
    path: Path,
    records: list[dict],
    references: list[tuple[str, int]] | None = None,
) -> None:
    """Write a tiny BAM file from a list of plain-dict records.

    Each record dict may contain: ``name`` (str), ``ref`` (str),
    ``start`` (int), ``end`` (int), ``mapq`` (int), ``reverse`` (bool),
    ``unmapped`` (bool), ``secondary`` (bool), ``supplementary`` (bool).
    """
    if references is None:
        references = [("ND1", 1000), ("COX1", 1000)]
    ref_names = [name for name, _length in references]
    ref_lengths = [length for _name, length in references]

    header = pysam.AlignmentHeader.from_references(ref_names, ref_lengths)

    with pysam.AlignmentFile(str(path), "wb", header=header) as sink:
        for rec in records:
            read = pysam.AlignedSegment(header)
            read.query_name = rec.get("name", "r")
            read.query_sequence = "A" * max(1, int(rec.get("end", 1) - rec.get("start", 0)))
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

    # pysam.sort is overkill for these tiny BAMs; they are already sorted
    # by the caller's order and we do not need coordinate order for the
    # unit tests.


# ---------- count_mapped_reads ----------------------------------------------


def test_count_mapped_reads_counts_only_primary_mapped(tmp_path) -> None:
    bam = tmp_path / "sample.bam"
    _make_bam(
        bam,
        [
            {"name": "r1", "ref": "ND1", "start": 100, "end": 130},
            {"name": "r2", "ref": "COX1", "start": 400, "end": 432},
            {"name": "r3_unmapped", "unmapped": True, "start": 0, "end": 20},
            {"name": "r4_secondary", "ref": "ND1", "start": 200, "end": 220, "secondary": True},
            {"name": "r5_supplementary", "ref": "COX1", "start": 500, "end": 520, "supplementary": True},
        ],
    )

    assert count_mapped_reads(bam) == 2


# ---------- filter_bam_mapq -------------------------------------------------


def test_filter_bam_mapq_keeps_at_or_above_threshold(tmp_path) -> None:
    src = tmp_path / "src.bam"
    _make_bam(
        src,
        [
            {"name": "low", "ref": "ND1", "start": 0, "end": 30, "mapq": 5},
            {"name": "border", "ref": "ND1", "start": 50, "end": 80, "mapq": 10},
            {"name": "high", "ref": "COX1", "start": 0, "end": 30, "mapq": 42},
        ],
    )
    dst = tmp_path / "dst.bam"

    kept = filter_bam_mapq(
        bam_in=src,
        bam_out=dst,
        mapq_threshold=10,
        indexer=lambda *a, **k: None,
    )

    assert kept == 2
    assert count_mapped_reads(dst) == 2
    with pysam.AlignmentFile(str(dst), "rb") as handle:
        names = sorted(read.query_name for read in handle)
    assert names == ["border", "high"]


def test_filter_bam_mapq_drops_unmapped_and_auxiliary_records(tmp_path) -> None:
    src = tmp_path / "src.bam"
    _make_bam(
        src,
        [
            {"name": "primary_high", "ref": "ND1", "start": 0, "end": 30, "mapq": 60},
            {"name": "unmapped_high", "unmapped": True, "start": 0, "end": 20, "mapq": 60},
            {"name": "secondary_high", "ref": "ND1", "start": 0, "end": 30, "mapq": 60, "secondary": True},
            {"name": "supplementary_high", "ref": "ND1", "start": 0, "end": 30, "mapq": 60, "supplementary": True},
        ],
    )
    dst = tmp_path / "dst.bam"

    kept = filter_bam_mapq(
        bam_in=src, bam_out=dst, mapq_threshold=0, indexer=lambda *a, **k: None
    )

    assert kept == 1


# ---------- bam_to_bed6 -----------------------------------------------------


def test_bam_to_bed6_emits_strand_plus_for_forward_reads(tmp_path) -> None:
    bam = tmp_path / "src.bam"
    _make_bam(
        bam,
        [
            {"name": "fwd", "ref": "ND1", "start": 10, "end": 40, "mapq": 30},
        ],
    )
    bed = tmp_path / "out.bed"

    rows = bam_to_bed6(bam_in=bam, bed_out=bed)

    assert rows == 1
    line = bed.read_text().strip()
    fields = line.split("\t")
    assert fields[0] == "ND1"
    assert fields[1] == "10"
    assert fields[2] == "40"
    assert fields[3] == "fwd"
    assert fields[4] == "30"
    assert fields[5] == "+"


def test_bam_to_bed6_emits_strand_minus_for_reverse_reads(tmp_path) -> None:
    bam = tmp_path / "src.bam"
    _make_bam(
        bam,
        [
            {"name": "rev", "ref": "ND1", "start": 500, "end": 530, "mapq": 40, "reverse": True},
        ],
    )
    bed = tmp_path / "out.bed"

    rows = bam_to_bed6(bam_in=bam, bed_out=bed)

    assert rows == 1
    line = bed.read_text().strip()
    fields = line.split("\t")
    assert fields[5] == "-"


def test_bam_to_bed6_clamps_mapq_score_into_valid_bed_range(tmp_path) -> None:
    # BED spec limits score to [0, 1000]. MAPQ is [0, 255] in practice;
    # the clamp is defensive and should never fire in real data, but we
    # verify the guard is in place.
    bam = tmp_path / "src.bam"
    _make_bam(
        bam,
        [
            {"name": "r", "ref": "ND1", "start": 0, "end": 30, "mapq": 200},
        ],
    )
    bed = tmp_path / "out.bed"
    bam_to_bed6(bam_in=bam, bed_out=bed)
    score = int(bed.read_text().strip().split("\t")[4])
    assert 0 <= score <= 1000


def test_bam_to_bed6_skips_unmapped_and_auxiliary(tmp_path) -> None:
    bam = tmp_path / "src.bam"
    _make_bam(
        bam,
        [
            {"name": "primary", "ref": "ND1", "start": 0, "end": 30, "mapq": 30},
            {"name": "unmapped", "unmapped": True, "start": 0, "end": 20},
            {"name": "secondary", "ref": "ND1", "start": 100, "end": 130, "mapq": 30, "secondary": True},
        ],
    )
    bed = tmp_path / "out.bed"

    rows = bam_to_bed6(bam_in=bam, bed_out=bed)

    assert rows == 1
    assert bed.read_text().strip().split("\t")[3] == "primary"


def test_bam_to_bed6_writes_one_line_per_primary_mapped_record(tmp_path) -> None:
    bam = tmp_path / "src.bam"
    _make_bam(
        bam,
        [
            {"name": "a", "ref": "ND1", "start": 0, "end": 30},
            {"name": "b", "ref": "ND1", "start": 50, "end": 80, "reverse": True},
            {"name": "c", "ref": "COX1", "start": 10, "end": 42},
        ],
    )
    bed = tmp_path / "out.bed"

    rows = bam_to_bed6(bam_in=bam, bed_out=bed)

    assert rows == 3
    bed_lines = bed.read_text().strip().splitlines()
    assert len(bed_lines) == 3
    for line in bed_lines:
        fields = line.split("\t")
        assert len(fields) == 6
