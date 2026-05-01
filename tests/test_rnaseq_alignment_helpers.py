"""Per-transcript counting + counts-matrix writers."""

from __future__ import annotations

from pathlib import Path

import pysam

from mitoribopy.align._types import ResolvedKit
from mitoribopy.rnaseq.alignment import (
    SampleAlignmentResult,
    count_per_transcript,
    write_counts_matrix,
    write_long_counts,
)
from mitoribopy.rnaseq.counts import load_ribo_counts


_HEADER = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [
        {"LN": 1000, "SN": "MT-ND1"},
        {"LN": 1500, "SN": "MT-CO1"},
    ],
}


def _make_record(
    name: str,
    *,
    ref_id: int = 0,
    pos: int = 100,
    flag: int = 0,
) -> pysam.AlignedSegment:
    rec = pysam.AlignedSegment()
    rec.query_name = name
    rec.flag = flag
    rec.reference_id = ref_id
    rec.reference_start = pos
    rec.mapping_quality = 60
    rec.cigarstring = "30M"
    rec.query_sequence = "A" * 30
    rec.query_qualities = pysam.qualitystring_to_array("I" * 30)
    return rec


def _build_bam(path: Path, records: list[pysam.AlignedSegment]) -> Path:
    sam_path = path.with_suffix(".sam")
    with pysam.AlignmentFile(str(sam_path), "wh", header=_HEADER) as handle:
        for rec in records:
            handle.write(rec)
    pysam.sort("-o", str(path), str(sam_path))
    pysam.index(str(path))
    sam_path.unlink()
    return path


def test_count_per_transcript_excludes_non_primary(tmp_path: Path) -> None:
    primary = _make_record("good", flag=0)  # primary, mapped
    unmapped = _make_record("unmapped", flag=4)  # unmapped
    secondary = _make_record("secondary", flag=256)  # secondary
    supplementary = _make_record("supp", flag=2048)  # supplementary
    bam = _build_bam(tmp_path / "se.bam", [primary, unmapped, secondary, supplementary])

    counts = count_per_transcript(bam)
    assert counts["MT-ND1"] == 1
    assert counts["MT-CO1"] == 0


def test_count_per_transcript_paired_counts_read1_only(tmp_path: Path) -> None:
    # Paired primary read 1 (flag 0x1 | 0x2 | 0x40 = 67), read 2 (0x1|0x2|0x80 = 131).
    r1 = _make_record("frag", flag=67, pos=100)
    r2 = _make_record("frag", flag=131, pos=400)
    bam = _build_bam(tmp_path / "pe.bam", [r1, r2])

    counts = count_per_transcript(bam)
    assert counts["MT-ND1"] == 1, "paired read1+read2 must count as ONE fragment"


def _result(name: str, counts: dict[str, int]) -> SampleAlignmentResult:
    return SampleAlignmentResult(
        sample=name,
        bam_path=Path("/tmp/dummy.bam"),
        counts=counts,
        paired=False,
        total_reads=sum(counts.values()),
        aligned_reads=sum(counts.values()),
        resolved_kit=ResolvedKit(
            kit="pretrimmed", adapter=None, umi_length=0, umi_position="5p"
        ),
    )


def test_write_counts_matrix_zero_fills(tmp_path: Path) -> None:
    a = _result("A", {"g1": 5, "g2": 3})
    b = _result("B", {"g2": 7, "g3": 2})
    out = tmp_path / "wide.tsv"
    write_counts_matrix([a, b], out)

    raw = out.read_text().splitlines()
    # P1.12: schema-version comment line precedes the column header.
    assert raw[0].startswith("# schema_version:")
    rows = [r for r in raw if not r.startswith("#")]
    assert rows[0] == "gene\tA\tB"
    body = {r.split("\t")[0]: r.split("\t")[1:] for r in rows[1:]}
    assert body["g1"] == ["5", "0"]
    assert body["g2"] == ["3", "7"]
    assert body["g3"] == ["0", "2"]


def test_write_long_counts_roundtrips_through_load_ribo_counts(tmp_path: Path) -> None:
    a = _result("A", {"g1": 5, "g2": 3})
    b = _result("B", {"g2": 7, "g3": 2})
    out = tmp_path / "long.tsv"
    write_long_counts([a, b], out)

    loaded = load_ribo_counts(out)
    assert loaded == {
        "g1": {"A": 5},
        "g2": {"A": 3, "B": 7},
        "g3": {"B": 2},
    }
