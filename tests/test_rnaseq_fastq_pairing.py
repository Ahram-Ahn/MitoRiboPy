"""Filename-token detection for SE / PE FASTQ inputs in rnaseq from-FASTQ mode."""

from __future__ import annotations

from pathlib import Path

import pytest

from mitoribopy.rnaseq.fastq_pairing import (
    FastqSample,
    detect_samples,
    enumerate_fastqs,
)


def _touch(p: Path) -> Path:
    p.write_bytes(b"")
    return p


def test_detect_pairs_for_every_supported_token(tmp_path: Path) -> None:
    cases = [
        ("bcl2fastq", "A_S1_L001_R1_001.fastq.gz", "A_S1_L001_R2_001.fastq.gz"),
        ("read12", "B_read1.fq.gz", "B_read2.fq.gz"),
        ("R12", "C_R1.fastq", "C_R2.fastq"),
        ("dot12", "D.1.fq.gz", "D.2.fq.gz"),
        ("us12", "E_1.fastq.gz", "E_2.fastq.gz"),
    ]
    paths = []
    for _label, r1, r2 in cases:
        paths.append(_touch(tmp_path / r1))
        paths.append(_touch(tmp_path / r2))

    samples = detect_samples(paths)
    by_name = {s.sample: s for s in samples}
    # Each pair collapses to one paired sample.
    assert set(by_name) == {"A_S1_L001", "B", "C", "D", "E"}
    for s in samples:
        assert s.paired, f"{s.sample!r} should be paired"
        assert s.r2 is not None


def test_lone_r1_is_single_end_with_stem_name(tmp_path: Path) -> None:
    p = _touch(tmp_path / "Sample42_R1.fastq.gz")
    [s] = detect_samples([p])
    assert s.sample == "Sample42"
    assert not s.paired
    assert s.r1 == p


def test_lone_r2_keeps_r2_in_sample_name(tmp_path: Path) -> None:
    p = _touch(tmp_path / "Sample42_R2.fastq.gz")
    [s] = detect_samples([p])
    assert s.sample == "Sample42_R2"
    assert not s.paired
    assert s.r1 == p


def test_no_mate_token_is_single_end(tmp_path: Path) -> None:
    p = _touch(tmp_path / "naked_sample.fastq.gz")
    [s] = detect_samples([p])
    assert s.sample == "naked_sample"
    assert not s.paired


def test_duplicate_pair_raises(tmp_path: Path) -> None:
    a = _touch(tmp_path / "X_R1.fastq.gz")
    b_dir = tmp_path / "alt"
    b_dir.mkdir()
    b = _touch(b_dir / "X_R1.fastq.gz")
    with pytest.raises(ValueError, match="Duplicate"):
        detect_samples([a, b])


def test_mixed_se_and_pe_in_one_call(tmp_path: Path) -> None:
    paths = [
        _touch(tmp_path / "pe_R1.fastq.gz"),
        _touch(tmp_path / "pe_R2.fastq.gz"),
        _touch(tmp_path / "se_only.fastq.gz"),
    ]
    samples = detect_samples(paths)
    by_name = {s.sample: s for s in samples}
    assert by_name["pe"].paired is True
    assert by_name["se_only"].paired is False


def test_enumerate_fastqs_directory_picks_up_all_suffixes(tmp_path: Path) -> None:
    d = tmp_path / "fqs"
    d.mkdir()
    for name in ("a.fq", "b.fq.gz", "c.fastq", "d.fastq.gz"):
        _touch(d / name)
    # Not a FASTQ — should be ignored when inside a directory.
    _touch(d / "README.md")
    out = enumerate_fastqs([d])
    assert sorted(p.name for p in out) == ["a.fq", "b.fq.gz", "c.fastq", "d.fastq.gz"]


def test_enumerate_fastqs_unknown_extension_raises(tmp_path: Path) -> None:
    bad = _touch(tmp_path / "weird.bam")
    with pytest.raises(ValueError, match="Unsupported FASTQ extension"):
        enumerate_fastqs([bad])


def test_enumerate_fastqs_dedupes_file_and_directory(tmp_path: Path) -> None:
    d = tmp_path / "fqs"
    d.mkdir()
    f = _touch(d / "x.fq.gz")
    out = enumerate_fastqs([d, f])
    # The same file shows up via both paths but must appear once.
    resolved = [p.resolve() for p in out]
    assert resolved.count(f.resolve()) == 1
