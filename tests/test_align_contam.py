"""Unit tests for ``mitoribopy.align.contam``.

bowtie2 is mocked via an injected ``runner`` callable. The bowtie2-style
stderr summary is synthesized so parsing stays independent of whether
bowtie2 is installed.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from mitoribopy.align import contam
from mitoribopy.align._types import ContamResult
from mitoribopy.align.contam import (
    parse_bowtie2_stderr,
    subtract_contaminants,
)


STDERR_1000_950_UNALIGNED = """\
1000 reads; of these:
  1000 (100.00%) were unpaired; of these:
    950 (95.00%) aligned 0 times
    40 (4.00%) aligned exactly 1 time
    10 (1.00%) aligned >1 times
5.00% overall alignment rate
"""


# ---------- index validation -------------------------------------------------


def test_subtract_contaminants_raises_when_index_missing(tmp_path) -> None:
    with pytest.raises(FileNotFoundError) as exc:
        subtract_contaminants(
            fastq_in=tmp_path / "in.fq.gz",
            contam_index=tmp_path / "nonexistent_prefix",
            fastq_out_unaligned=tmp_path / "out.fq.gz",
            runner=lambda *a, **k: SimpleNamespace(returncode=0, stdout="", stderr=""),
        )
    assert "bowtie2-build" in str(exc.value)


def test_subtract_contaminants_accepts_small_index_suffix(tmp_path, monkeypatch) -> None:
    prefix = tmp_path / "contam_idx"
    (Path(str(prefix) + ".1.bt2")).write_text("dummy")
    for suffix in (".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"):
        (Path(str(prefix) + suffix)).write_text("dummy")

    def fake_runner(cmd, **kwargs):
        return SimpleNamespace(returncode=0, stdout="", stderr=STDERR_1000_950_UNALIGNED)

    result = subtract_contaminants(
        fastq_in=tmp_path / "in.fq.gz",
        contam_index=prefix,
        fastq_out_unaligned=tmp_path / "out.fq.gz",
        runner=fake_runner,
    )

    assert isinstance(result, ContamResult)
    assert result.total_reads == 1000
    assert result.aligned_to_contam == 50
    assert result.unaligned_reads == 950


def test_subtract_contaminants_accepts_large_index_suffix(tmp_path) -> None:
    prefix = tmp_path / "contam_idx_large"
    (Path(str(prefix) + ".1.bt2l")).write_text("dummy")

    def fake_runner(cmd, **kwargs):
        return SimpleNamespace(returncode=0, stdout="", stderr=STDERR_1000_950_UNALIGNED)

    result = subtract_contaminants(
        fastq_in=tmp_path / "in.fq.gz",
        contam_index=prefix,
        fastq_out_unaligned=tmp_path / "out.fq.gz",
        runner=fake_runner,
    )
    assert result.total_reads == 1000


# ---------- command construction --------------------------------------------


def test_subtract_contaminants_command_includes_strand_flag_for_forward(tmp_path) -> None:
    prefix = tmp_path / "contam_idx"
    (Path(str(prefix) + ".1.bt2")).write_text("dummy")

    recorded: list[list[str]] = []

    def fake_runner(cmd, **kwargs):
        recorded.append(list(cmd))
        return SimpleNamespace(returncode=0, stdout="", stderr=STDERR_1000_950_UNALIGNED)

    subtract_contaminants(
        fastq_in=tmp_path / "in.fq.gz",
        contam_index=prefix,
        fastq_out_unaligned=tmp_path / "out.fq.gz",
        strandedness="forward",
        threads=2,
        seed=42,
        runner=fake_runner,
    )

    assert len(recorded) == 1
    cmd = recorded[0]

    # Verify the invariants that matter for biology-correctness:
    assert cmd[0] == "bowtie2"
    assert "--end-to-end" in cmd
    assert "--very-sensitive" in cmd
    assert "--no-unal" in cmd
    assert "--un-gz" in cmd and str(tmp_path / "out.fq.gz") in cmd
    assert "-L" in cmd and cmd[cmd.index("-L") + 1] == "18"
    assert "-p" in cmd and cmd[cmd.index("-p") + 1] == "2"
    assert "--seed" in cmd and cmd[cmd.index("--seed") + 1] == "42"
    assert "--norc" in cmd  # strand=forward
    assert "--nofw" not in cmd


def test_subtract_contaminants_command_maps_reverse_strand(tmp_path) -> None:
    prefix = tmp_path / "contam_idx"
    (Path(str(prefix) + ".1.bt2")).write_text("dummy")

    recorded: list[list[str]] = []

    def fake_runner(cmd, **kwargs):
        recorded.append(list(cmd))
        return SimpleNamespace(returncode=0, stdout="", stderr=STDERR_1000_950_UNALIGNED)

    subtract_contaminants(
        fastq_in=tmp_path / "in.fq.gz",
        contam_index=prefix,
        fastq_out_unaligned=tmp_path / "out.fq.gz",
        strandedness="reverse",
        runner=fake_runner,
    )

    cmd = recorded[0]
    assert "--nofw" in cmd
    assert "--norc" not in cmd


def test_subtract_contaminants_command_omits_strand_flag_when_unstranded(tmp_path) -> None:
    prefix = tmp_path / "contam_idx"
    (Path(str(prefix) + ".1.bt2")).write_text("dummy")

    recorded: list[list[str]] = []

    def fake_runner(cmd, **kwargs):
        recorded.append(list(cmd))
        return SimpleNamespace(returncode=0, stdout="", stderr=STDERR_1000_950_UNALIGNED)

    subtract_contaminants(
        fastq_in=tmp_path / "in.fq.gz",
        contam_index=prefix,
        fastq_out_unaligned=tmp_path / "out.fq.gz",
        strandedness="unstranded",
        runner=fake_runner,
    )

    cmd = recorded[0]
    assert "--norc" not in cmd
    assert "--nofw" not in cmd


def test_subtract_contaminants_propagates_extra_flags(tmp_path) -> None:
    prefix = tmp_path / "contam_idx"
    (Path(str(prefix) + ".1.bt2")).write_text("dummy")

    recorded: list[list[str]] = []

    def fake_runner(cmd, **kwargs):
        recorded.append(list(cmd))
        return SimpleNamespace(returncode=0, stdout="", stderr=STDERR_1000_950_UNALIGNED)

    subtract_contaminants(
        fastq_in=tmp_path / "in.fq.gz",
        contam_index=prefix,
        fastq_out_unaligned=tmp_path / "out.fq.gz",
        extra_flags=["--mm", "-k", "3"],
        runner=fake_runner,
    )

    cmd = recorded[0]
    assert "--mm" in cmd
    assert "-k" in cmd


# ---------- stderr parsing --------------------------------------------------


def test_parse_bowtie2_stderr_extracts_total_and_unaligned() -> None:
    total, unaligned = parse_bowtie2_stderr(STDERR_1000_950_UNALIGNED)
    assert total == 1000
    assert unaligned == 950


def test_parse_bowtie2_stderr_raises_on_unexpected_format() -> None:
    with pytest.raises(ValueError):
        parse_bowtie2_stderr("nothing useful here")


# ---------- error propagation ------------------------------------------------


def test_subtract_contaminants_raises_runtime_error_on_nonzero_exit(tmp_path) -> None:
    prefix = tmp_path / "contam_idx"
    (Path(str(prefix) + ".1.bt2")).write_text("dummy")

    def fake_runner(cmd, **kwargs):
        return SimpleNamespace(returncode=1, stdout="", stderr="(ERR): unmatched index")

    with pytest.raises(RuntimeError) as exc:
        subtract_contaminants(
            fastq_in=tmp_path / "in.fq.gz",
            contam_index=prefix,
            fastq_out_unaligned=tmp_path / "out.fq.gz",
            runner=fake_runner,
        )
    assert "unmatched index" in str(exc.value)
