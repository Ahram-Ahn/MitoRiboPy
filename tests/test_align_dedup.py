"""Unit tests for ``mitoribopy.align.dedup``.

External tools (umi_tools) are mocked via injected runners and
shutil.which is stubbed so the module behaves as if they were present.
pysam.index is replaced with a recording no-op.

The picard ``mark-duplicates`` strategy was removed in v0.4.5 because
coordinate-only deduplication destroys codon-occupancy signal on
mt-Ribo-seq libraries (see docs/validation/taco1_ko_regression.md).
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from mitoribopy.align import dedup
from mitoribopy.align._types import DedupResult
from mitoribopy.align.dedup import (
    resolve_dedup_strategy,
    run_dedup,
    run_umi_tools_dedup,
    skip_dedup,
)
from mitoribopy.align.tool_check import ToolNotFoundError


# ---------- resolve_dedup_strategy ------------------------------------------


def test_resolve_auto_with_umis_picks_umi_tools() -> None:
    assert resolve_dedup_strategy("auto", umi_length=8) == "umi-tools"


def test_resolve_auto_without_umis_picks_skip() -> None:
    assert resolve_dedup_strategy("auto", umi_length=0) == "skip"


def test_resolve_umi_tools_without_umis_raises() -> None:
    with pytest.raises(ValueError) as exc:
        resolve_dedup_strategy("umi-tools", umi_length=0)
    assert "--umi-length > 0" in str(exc.value)


def test_resolve_umi_tools_with_umis_ok() -> None:
    assert resolve_dedup_strategy("umi-tools", umi_length=12) == "umi-tools"


def test_resolve_skip_is_always_allowed() -> None:
    assert resolve_dedup_strategy("skip", umi_length=0) == "skip"
    assert resolve_dedup_strategy("skip", umi_length=8) == "skip"


def test_resolve_unknown_strategy_raises() -> None:
    with pytest.raises(ValueError):
        resolve_dedup_strategy("does-not-exist", umi_length=0)  # type: ignore[arg-type]


# ---------- skip_dedup ------------------------------------------------------


def test_skip_dedup_creates_output_with_same_record_count(tmp_path, monkeypatch) -> None:
    bam_in = tmp_path / "in.bam"
    bam_in.write_bytes(b"dummy-bam")
    bam_out = tmp_path / "out.bam"

    monkeypatch.setattr(dedup.pysam, "index", lambda path, *a, **k: None)

    result = skip_dedup(
        bam_in=bam_in,
        bam_out=bam_out,
        counter=lambda p: 1234,
        indexer=lambda p, *a, **k: None,
    )

    assert isinstance(result, DedupResult)
    assert result.strategy == "skip"
    assert result.input_reads == result.output_reads == 1234
    assert bam_out.exists()
    assert bam_out.read_bytes() == b"dummy-bam"


# ---------- run_umi_tools_dedup ---------------------------------------------


def test_run_umi_tools_dedup_builds_expected_command(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(dedup.shutil, "which", lambda name: "/usr/local/bin/umi_tools")

    recorded: list[list[str]] = []

    def runner(cmd, **kwargs):
        recorded.append(list(cmd))
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    bam_in = tmp_path / "sorted.bam"
    bam_in.write_bytes(b"")

    bam_out = tmp_path / "deduped.bam"
    log_path = tmp_path / "umi.log"

    counts = iter([1000, 720])  # input then output

    def counter(_path):
        return next(counts)

    result = run_umi_tools_dedup(
        bam_in=bam_in,
        bam_out=bam_out,
        log_path=log_path,
        runner=runner,
        counter=counter,
        indexer=lambda *a, **k: None,
    )

    assert isinstance(result, DedupResult)
    assert result.strategy == "umi-tools"
    assert result.input_reads == 1000
    assert result.output_reads == 720

    cmd = recorded[0]
    assert cmd[:2] == ["umi_tools", "dedup"]
    assert f"--stdin={bam_in}" in cmd
    assert f"--stdout={bam_out}" in cmd
    assert "--method=unique" in cmd
    assert "--umi-separator=_" in cmd
    assert any(arg.startswith("--log=") for arg in cmd)


def test_run_umi_tools_dedup_raises_when_tool_missing(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(dedup.shutil, "which", lambda name: None)
    with pytest.raises(ToolNotFoundError) as exc:
        run_umi_tools_dedup(
            bam_in=tmp_path / "in.bam",
            bam_out=tmp_path / "out.bam",
            log_path=tmp_path / "log",
            runner=lambda *a, **k: SimpleNamespace(returncode=0, stdout="", stderr=""),
            counter=lambda _: 0,
            indexer=lambda *a, **k: None,
        )
    assert "umi_tools" in str(exc.value)
    assert "bioconda" in str(exc.value)


def test_run_umi_tools_dedup_raises_on_nonzero_exit(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(dedup.shutil, "which", lambda name: "/usr/local/bin/umi_tools")

    def runner(cmd, **kwargs):
        return SimpleNamespace(returncode=1, stdout="", stderr="bad umi pattern")

    with pytest.raises(RuntimeError) as exc:
        run_umi_tools_dedup(
            bam_in=tmp_path / "in.bam",
            bam_out=tmp_path / "out.bam",
            log_path=tmp_path / "log",
            runner=runner,
            counter=lambda _: 0,
            indexer=lambda *a, **k: None,
        )
    assert "bad umi pattern" in str(exc.value)


# ---------- run_dedup dispatcher --------------------------------------------


def test_run_dedup_auto_with_umis_dispatches_umi_tools(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(dedup.shutil, "which", lambda name: "/usr/local/bin/umi_tools")

    recorded: list[list[str]] = []

    def runner(cmd, **kwargs):
        recorded.append(list(cmd))
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    bam_in = tmp_path / "sorted.bam"
    bam_in.write_bytes(b"")
    counts = iter([1000, 720])

    result = run_dedup(
        strategy="auto",
        umi_length=8,
        bam_in=bam_in,
        bam_out=tmp_path / "out.bam",
        runner=runner,
        counter=lambda _: next(counts),
        indexer=lambda *a, **k: None,
    )

    assert result.strategy == "umi-tools"
    assert recorded[0][0] == "umi_tools"


def test_run_dedup_auto_without_umis_dispatches_skip(tmp_path, monkeypatch) -> None:
    bam_in = tmp_path / "sorted.bam"
    bam_in.write_bytes(b"data")
    recorded: list[list[str]] = []

    def runner(cmd, **kwargs):
        recorded.append(list(cmd))
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    result = run_dedup(
        strategy="auto",
        umi_length=0,
        bam_in=bam_in,
        bam_out=tmp_path / "out.bam",
        runner=runner,
        counter=lambda _: 1234,
        indexer=lambda *a, **k: None,
    )

    assert result.strategy == "skip"
    assert recorded == []  # skip never shells out
