"""Tests for the parity-based FASTQ splitter used by auto-pseudo-replicate."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from mitoribopy.rnaseq.fastq_pairing import FastqSample
from mitoribopy.rnaseq.split_replicates import (
    split_sample_into_pseudo_replicates,
    stream_split_by_parity,
)


def _write_fastq(path: Path, n_records: int, *, prefix: str = "r") -> None:
    lines: list[str] = []
    for i in range(n_records):
        lines.append(f"@{prefix}{i}")
        lines.append("ACGTACGT")
        lines.append("+")
        lines.append("IIIIIIII")
    blob = "\n".join(lines) + "\n"
    if str(path).endswith(".gz"):
        with gzip.open(str(path), "wt") as h:
            h.write(blob)
    else:
        path.write_text(blob)


def _read_fastq_ids(path: Path) -> list[str]:
    open_fn = gzip.open if str(path).endswith(".gz") else open
    ids: list[str] = []
    with open_fn(str(path), "rt") as h:
        for i, line in enumerate(h):
            if i % 4 == 0:
                ids.append(line.strip().lstrip("@"))
    return ids


def test_stream_split_by_parity_alternates_records(tmp_path: Path) -> None:
    src = tmp_path / "src.fastq"
    _write_fastq(src, 10)
    a = tmp_path / "even.fastq"
    b = tmp_path / "odd.fastq"
    n = stream_split_by_parity(src, a, b)
    assert n == 10
    assert _read_fastq_ids(a) == ["r0", "r2", "r4", "r6", "r8"]
    assert _read_fastq_ids(b) == ["r1", "r3", "r5", "r7", "r9"]


def test_stream_split_by_parity_handles_odd_count(tmp_path: Path) -> None:
    src = tmp_path / "src.fastq"
    _write_fastq(src, 7)
    a = tmp_path / "even.fastq"
    b = tmp_path / "odd.fastq"
    stream_split_by_parity(src, a, b)
    # 7 records: even = {0,2,4,6}, odd = {1,3,5}
    assert _read_fastq_ids(a) == ["r0", "r2", "r4", "r6"]
    assert _read_fastq_ids(b) == ["r1", "r3", "r5"]


def test_stream_split_gzip_roundtrip(tmp_path: Path) -> None:
    src = tmp_path / "src.fq.gz"
    _write_fastq(src, 6)
    a = tmp_path / "even.fq.gz"
    b = tmp_path / "odd.fq.gz"
    stream_split_by_parity(src, a, b)
    assert _read_fastq_ids(a) == ["r0", "r2", "r4"]
    assert _read_fastq_ids(b) == ["r1", "r3", "r5"]


def test_stream_split_truncated_record_raises(tmp_path: Path) -> None:
    src = tmp_path / "bad.fastq"
    src.write_text("@r0\nACGT\n+\n")  # missing quality line
    with pytest.raises(ValueError, match="Truncated FASTQ record"):
        stream_split_by_parity(src, tmp_path / "a", tmp_path / "b")


def test_split_sample_se(tmp_path: Path) -> None:
    src = tmp_path / "ABC_R1.fq.gz"
    _write_fastq(src, 8)
    sample = FastqSample(sample="ABC", r1=src, r2=None)
    rep1, rep2 = split_sample_into_pseudo_replicates(sample, tmp_path / "split")

    assert rep1.sample == "ABC_rep1"
    assert rep2.sample == "ABC_rep2"
    assert rep1.r2 is None and rep2.r2 is None
    assert _read_fastq_ids(rep1.r1) == ["r0", "r2", "r4", "r6"]
    assert _read_fastq_ids(rep2.r1) == ["r1", "r3", "r5", "r7"]


def test_split_sample_pe_keeps_r1_r2_in_lockstep(tmp_path: Path) -> None:
    r1 = tmp_path / "ABC_R1.fq.gz"
    r2 = tmp_path / "ABC_R2.fq.gz"
    # R1 reads named r0..r5; R2 reads named s0..s5 — same record index
    # must end up in the same rep so pairs stay matched downstream.
    _write_fastq(r1, 6, prefix="r")
    _write_fastq(r2, 6, prefix="s")
    sample = FastqSample(sample="ABC", r1=r1, r2=r2)
    rep1, rep2 = split_sample_into_pseudo_replicates(sample, tmp_path / "split")

    rep1_r1_ids = _read_fastq_ids(rep1.r1)
    rep1_r2_ids = _read_fastq_ids(rep1.r2)
    rep2_r1_ids = _read_fastq_ids(rep2.r1)
    rep2_r2_ids = _read_fastq_ids(rep2.r2)

    # rep1: even indices 0, 2, 4 from each mate
    assert rep1_r1_ids == ["r0", "r2", "r4"]
    assert rep1_r2_ids == ["s0", "s2", "s4"]
    # rep2: odd indices 1, 3, 5
    assert rep2_r1_ids == ["r1", "r3", "r5"]
    assert rep2_r2_ids == ["s1", "s3", "s5"]
    # Pair invariant: R1[i] and R2[i] share the same numeric index
    for r1_id, r2_id in zip(rep1_r1_ids, rep1_r2_ids):
        assert r1_id[1:] == r2_id[1:]
