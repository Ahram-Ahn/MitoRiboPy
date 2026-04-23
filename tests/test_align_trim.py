"""Unit tests for ``mitoribopy.align.trim``.

All external cutadapt invocations are mocked via an injected ``runner``
callable. JSON-log parsing is tested against a synthetic cutadapt log so
we do not depend on a specific cutadapt version being installed.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from mitoribopy.align import trim
from mitoribopy.align._types import CutadaptResult, ResolvedKit
from mitoribopy.align.trim import (
    parse_cutadapt_json,
    resolve_kit_settings,
    run_cutadapt,
)


# ---------- kit resolution ---------------------------------------------------


def test_resolve_kit_custom_without_adapter_raises() -> None:
    with pytest.raises(ValueError) as exc:
        resolve_kit_settings("custom")
    assert "--adapter" in str(exc.value)


def test_resolve_kit_custom_with_explicit_adapter_has_no_umi_by_default() -> None:
    resolved = resolve_kit_settings("custom", adapter="AAAA")
    assert resolved == ResolvedKit(
        kit="custom", adapter="AAAA", umi_length=0, umi_position="5p"
    )


def test_resolve_kit_truseq_smallrna_defaults() -> None:
    resolved = resolve_kit_settings("truseq_smallrna")
    assert resolved.adapter == "TGGAATTCTCGGGTGCCAAGG"
    assert resolved.umi_length == 0


def test_resolve_kit_nebnext_smallrna_defaults() -> None:
    resolved = resolve_kit_settings("nebnext_smallrna")
    assert resolved.adapter == "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    assert resolved.umi_length == 0


def test_resolve_kit_nebnext_ultra_umi_has_eight_nt_five_prime_umi() -> None:
    resolved = resolve_kit_settings("nebnext_ultra_umi")
    assert resolved.umi_length == 8
    assert resolved.umi_position == "5p"


def test_resolve_kit_qiaseq_has_twelve_nt_three_prime_umi() -> None:
    resolved = resolve_kit_settings("qiaseq_mirna")
    assert resolved.umi_length == 12
    assert resolved.umi_position == "3p"
    assert resolved.adapter == "AACTGTAGGCACCATCAAT"


def test_resolve_kit_explicit_adapter_overrides_preset() -> None:
    resolved = resolve_kit_settings("truseq_smallrna", adapter="CUSTOM_SEQ")
    assert resolved.adapter == "CUSTOM_SEQ"


def test_resolve_kit_explicit_umi_length_overrides_preset() -> None:
    resolved = resolve_kit_settings("truseq_smallrna", umi_length=6)
    assert resolved.umi_length == 6


def test_resolve_kit_rejects_unknown_preset_name() -> None:
    with pytest.raises(KeyError):
        resolve_kit_settings("not_a_kit")


def test_resolve_kit_rejects_negative_umi_length() -> None:
    with pytest.raises(ValueError):
        resolve_kit_settings("truseq_smallrna", umi_length=-1)


def test_resolve_kit_rejects_invalid_umi_position() -> None:
    with pytest.raises(ValueError):
        resolve_kit_settings(
            "custom", adapter="AAAA", umi_position="middle"  # type: ignore[arg-type]
        )


# ---------- JSON parsing -----------------------------------------------------


def test_parse_cutadapt_json_extracts_stage_counts(tmp_path) -> None:
    log = tmp_path / "sample.cutadapt.json"
    log.write_text(
        json.dumps(
            {
                "read_counts": {"input": 1000, "output": 875},
                "adapters_read1": [
                    {"total_matches": 800},
                    {"total_matches": 20},
                ],
            }
        )
    )

    counts = parse_cutadapt_json(log)

    assert counts == {
        "input_reads": 1000,
        "reads_with_adapter": 820,
        "reads_passing_filters": 875,
    }


def test_parse_cutadapt_json_handles_missing_keys(tmp_path) -> None:
    log = tmp_path / "empty.json"
    log.write_text(json.dumps({}))

    counts = parse_cutadapt_json(log)

    assert counts == {
        "input_reads": 0,
        "reads_with_adapter": 0,
        "reads_passing_filters": 0,
    }


# ---------- runner behavior: single pass (no UMI or 5' UMI) -----------------


def _fake_runner_factory():
    """Return ``(runner, recorded)`` — runner records every call."""
    recorded: list[list[str]] = []

    def runner(cmd, **kwargs):
        recorded.append(list(cmd))
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    return runner, recorded


def _write_log(path: Path, input_reads: int, adapter: int, output: int) -> None:
    path.write_text(
        json.dumps(
            {
                "read_counts": {"input": input_reads, "output": output},
                "adapters_read1": [{"total_matches": adapter}],
            }
        )
    )


def test_run_cutadapt_single_pass_truseq_builds_expected_command(tmp_path) -> None:
    runner, recorded = _fake_runner_factory()
    log_json = tmp_path / "sample.cutadapt.json"

    # Prepare the log before the (mocked) run so parse sees realistic data.
    _write_log(log_json, input_reads=1000, adapter=900, output=870)

    resolved = resolve_kit_settings("truseq_smallrna")
    result = run_cutadapt(
        fastq_in=tmp_path / "in.fq.gz",
        fastq_out=tmp_path / "out.fq.gz",
        resolved=resolved,
        threads=4,
        log_json=log_json,
        runner=runner,
    )

    assert len(recorded) == 1
    cmd = recorded[0]
    assert cmd[0] == "cutadapt"
    assert "-a" in cmd and cmd[cmd.index("-a") + 1] == "TGGAATTCTCGGGTGCCAAGG"
    assert "--minimum-length" in cmd and cmd[cmd.index("--minimum-length") + 1] == "15"
    assert "--maximum-length" in cmd and cmd[cmd.index("--maximum-length") + 1] == "45"
    assert "--cores" in cmd and cmd[cmd.index("--cores") + 1] == "4"
    assert "--json" in cmd
    assert cmd[-1] == str(tmp_path / "in.fq.gz")

    assert isinstance(result, CutadaptResult)
    assert result.input_reads == 1000
    assert result.reads_with_adapter == 900
    assert result.reads_passing_filters == 870


def test_run_cutadapt_single_pass_5p_umi_uses_cut_prefix_rename(tmp_path) -> None:
    runner, recorded = _fake_runner_factory()
    log_json = tmp_path / "sample.cutadapt.json"
    _write_log(log_json, 10, 9, 8)

    resolved = resolve_kit_settings("nebnext_ultra_umi")  # 8 nt 5' UMI
    run_cutadapt(
        fastq_in=tmp_path / "in.fq.gz",
        fastq_out=tmp_path / "out.fq.gz",
        resolved=resolved,
        log_json=log_json,
        runner=runner,
    )

    assert len(recorded) == 1
    cmd = recorded[0]
    assert "-u" in cmd
    ui = cmd.index("-u")
    assert cmd[ui + 1] == "8"  # 8 nt UMI cut from the 5' end
    rename_token = next(t for t in cmd if t.startswith("--rename="))
    assert "{cut_prefix}" in rename_token


def test_run_cutadapt_two_pass_3p_umi_issues_two_commands(tmp_path) -> None:
    runner, recorded = _fake_runner_factory()
    log_json = tmp_path / "sample.cutadapt.json"

    resolved = resolve_kit_settings("qiaseq_mirna")  # 12 nt 3' UMI

    # Prepare logs for both passes. Pass 1 log lives at <log>.pass1.json;
    # pass 2 at <log>.pass2.json (per the module's naming scheme).
    pass1_log = log_json.with_name(log_json.stem + ".pass1.json")
    pass2_log = log_json.with_name(log_json.stem + ".pass2.json")
    _write_log(pass1_log, input_reads=1000, adapter=900, output=900)
    _write_log(pass2_log, input_reads=900, adapter=0, output=888)

    result = run_cutadapt(
        fastq_in=tmp_path / "in.fq.gz",
        fastq_out=tmp_path / "out.fq.gz",
        resolved=resolved,
        log_json=log_json,
        runner=runner,
    )

    assert len(recorded) == 2
    pass1_cmd, pass2_cmd = recorded

    # Pass 1: trims adapter, no UMI handling (UMI handled in pass 2).
    assert "-a" in pass1_cmd
    assert pass1_cmd[pass1_cmd.index("-a") + 1] == "AACTGTAGGCACCATCAAT"
    # Pass 1 should NOT have cut_prefix / cut_suffix rename yet.
    assert not any(t.startswith("--rename=") for t in pass1_cmd)

    # Pass 2: extracts 3' UMI using -u -12 and cut_suffix rename.
    assert "-u" in pass2_cmd
    assert pass2_cmd[pass2_cmd.index("-u") + 1] == "-12"
    rename_token = next(t for t in pass2_cmd if t.startswith("--rename="))
    assert "{cut_suffix}" in rename_token

    # Input count comes from pass 1; output count from pass 2.
    assert result.input_reads == 1000
    assert result.reads_passing_filters == 888


def test_run_cutadapt_raises_when_subprocess_exits_nonzero(tmp_path) -> None:
    def runner(cmd, **kwargs):
        return SimpleNamespace(returncode=1, stdout="", stderr="adapter file bad")

    log_json = tmp_path / "sample.cutadapt.json"
    resolved = resolve_kit_settings("truseq_smallrna")
    with pytest.raises(RuntimeError) as exc:
        run_cutadapt(
            fastq_in=tmp_path / "in.fq.gz",
            fastq_out=tmp_path / "out.fq.gz",
            resolved=resolved,
            log_json=log_json,
            runner=runner,
        )
    assert "cutadapt failed" in str(exc.value)
    assert "adapter file bad" in str(exc.value)
