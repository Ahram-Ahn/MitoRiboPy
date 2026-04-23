"""Unit tests for ``mitoribopy.align.align``.

bowtie2, pysam.sort, and pysam.index are all mocked via injected
callables so these tests run with neither bowtie2 nor samtools installed
and without any real BAM machinery.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from mitoribopy.align._types import AlignResult
from mitoribopy.align.align import align_mt


STDERR_10k_500_UNALIGNED = """\
10000 reads; of these:
  10000 (100.00%) were unpaired; of these:
    500 (5.00%) aligned 0 times
    9000 (90.00%) aligned exactly 1 time
    500 (5.00%) aligned >1 times
95.00% overall alignment rate
"""


def _make_index(tmp_path: Path, prefix_name: str = "mt_idx") -> Path:
    prefix = tmp_path / prefix_name
    for suffix in (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"):
        (Path(str(prefix) + suffix)).write_text("dummy")
    return prefix


def _fake_runner(stderr: str = STDERR_10k_500_UNALIGNED, returncode: int = 0):
    recorded: list[list[str]] = []

    def runner(cmd, **kwargs):
        recorded.append(list(cmd))
        return SimpleNamespace(returncode=returncode, stdout="", stderr=stderr)

    return runner, recorded


def _touching_sorter(bam_out_path: list[Path]):
    def sorter(*args, **kwargs):
        # pysam.sort variadic signature: positional SAMtools-style flags.
        # We expect the command to include '-o <path>' - touch that path.
        args_list = list(args)
        if "-o" in args_list:
            target = Path(args_list[args_list.index("-o") + 1])
            target.write_bytes(b"")
            bam_out_path.append(target)

    return sorter


def _no_op_indexer(calls: list[str]):
    def indexer(path, *args, **kwargs):
        calls.append(str(path))

    return indexer


# ---------- happy path -------------------------------------------------------


def test_align_mt_builds_command_and_returns_counts(tmp_path) -> None:
    mt_index = _make_index(tmp_path)
    runner, recorded = _fake_runner()
    sort_targets: list[Path] = []
    index_calls: list[str] = []

    result = align_mt(
        fastq_in=tmp_path / "in.fq.gz",
        mt_index=mt_index,
        bam_out=tmp_path / "sample.bam",
        strandedness="forward",
        threads=4,
        seed=42,
        runner=runner,
        sorter=_touching_sorter(sort_targets),
        indexer=_no_op_indexer(index_calls),
    )

    assert isinstance(result, AlignResult)
    assert result.total_reads == 10000
    assert result.aligned == 9500
    assert result.bam_path == tmp_path / "sample.bam"

    assert len(recorded) == 1
    cmd = recorded[0]
    assert cmd[0] == "bowtie2"
    assert "-x" in cmd and cmd[cmd.index("-x") + 1] == str(mt_index)
    assert "-U" in cmd
    assert "--end-to-end" in cmd
    assert "--very-sensitive" in cmd
    assert "-L" in cmd and cmd[cmd.index("-L") + 1] == "18"
    assert "--no-unal" in cmd
    assert "--norc" in cmd  # strand=forward
    assert "-p" in cmd and cmd[cmd.index("-p") + 1] == "4"
    assert "--seed" in cmd and cmd[cmd.index("--seed") + 1] == "42"
    assert "-S" in cmd  # SAM output

    assert sort_targets == [tmp_path / "sample.bam"]
    assert index_calls == [str(tmp_path / "sample.bam")]


def test_align_mt_maps_reverse_strand(tmp_path) -> None:
    mt_index = _make_index(tmp_path)
    runner, recorded = _fake_runner()

    align_mt(
        fastq_in=tmp_path / "in.fq.gz",
        mt_index=mt_index,
        bam_out=tmp_path / "sample.bam",
        strandedness="reverse",
        runner=runner,
        sorter=_touching_sorter([]),
        indexer=_no_op_indexer([]),
    )

    cmd = recorded[0]
    assert "--nofw" in cmd
    assert "--norc" not in cmd


def test_align_mt_unstranded_omits_strand_flag(tmp_path) -> None:
    mt_index = _make_index(tmp_path)
    runner, recorded = _fake_runner()

    align_mt(
        fastq_in=tmp_path / "in.fq.gz",
        mt_index=mt_index,
        bam_out=tmp_path / "sample.bam",
        strandedness="unstranded",
        runner=runner,
        sorter=_touching_sorter([]),
        indexer=_no_op_indexer([]),
    )
    cmd = recorded[0]
    assert "--norc" not in cmd
    assert "--nofw" not in cmd


def test_align_mt_propagates_extra_flags(tmp_path) -> None:
    mt_index = _make_index(tmp_path)
    runner, recorded = _fake_runner()

    align_mt(
        fastq_in=tmp_path / "in.fq.gz",
        mt_index=mt_index,
        bam_out=tmp_path / "sample.bam",
        extra_flags=["-k", "5"],
        runner=runner,
        sorter=_touching_sorter([]),
        indexer=_no_op_indexer([]),
    )
    cmd = recorded[0]
    assert "-k" in cmd


# ---------- error paths -----------------------------------------------------


def test_align_mt_raises_when_mt_index_missing(tmp_path) -> None:
    with pytest.raises(FileNotFoundError) as exc:
        align_mt(
            fastq_in=tmp_path / "in.fq.gz",
            mt_index=tmp_path / "no_such_prefix",
            bam_out=tmp_path / "sample.bam",
            runner=lambda *a, **k: SimpleNamespace(returncode=0, stdout="", stderr=""),
            sorter=_touching_sorter([]),
            indexer=_no_op_indexer([]),
        )
    assert "bowtie2-build" in str(exc.value)


def test_align_mt_raises_on_bowtie2_failure(tmp_path) -> None:
    mt_index = _make_index(tmp_path)
    runner, _ = _fake_runner(stderr="(ERR): index corrupt", returncode=1)

    with pytest.raises(RuntimeError) as exc:
        align_mt(
            fastq_in=tmp_path / "in.fq.gz",
            mt_index=mt_index,
            bam_out=tmp_path / "sample.bam",
            runner=runner,
            sorter=_touching_sorter([]),
            indexer=_no_op_indexer([]),
        )
    assert "index corrupt" in str(exc.value)


def test_align_mt_raises_on_unparseable_stderr(tmp_path) -> None:
    mt_index = _make_index(tmp_path)
    runner, _ = _fake_runner(stderr="bowtie2 printed nothing useful\n")

    with pytest.raises(ValueError):
        align_mt(
            fastq_in=tmp_path / "in.fq.gz",
            mt_index=mt_index,
            bam_out=tmp_path / "sample.bam",
            runner=runner,
            sorter=_touching_sorter([]),
            indexer=_no_op_indexer([]),
        )
