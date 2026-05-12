"""Tests for the MitoRiboPy subcommand CLI introduced in v0.3.0 Phase 2."""

from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest

from mitoribopy import cli
from mitoribopy.cli import common
from mitoribopy.config.runtime import load_user_config


# ---------- top-level dispatch -----------------------------------------------


def test_top_level_help_lists_all_subcommands(capsys) -> None:
    exit_code = cli.main(["--help"])
    captured = capsys.readouterr()
    assert exit_code == 0
    for subcommand in ("align", "rpf", "rnaseq", "all"):
        assert subcommand in captured.out


def test_top_level_help_does_not_import_matplotlib() -> None:
    probe = (
        "import sys; "
        "from mitoribopy import cli; "
        "rc = cli.main(['--help']); "
        "raise SystemExit(0 if rc == 0 and 'matplotlib' not in sys.modules "
        "and 'matplotlib.pyplot' not in sys.modules else 1)"
    )
    proc = subprocess.run(
        [sys.executable, "-c", probe],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr


def test_top_level_no_args_prints_help(capsys) -> None:
    exit_code = cli.main([])
    captured = capsys.readouterr()
    assert exit_code == 0
    assert "usage: mitoribopy" in captured.out


def test_top_level_version(capsys) -> None:
    from mitoribopy import __version__

    exit_code = cli.main(["--version"])
    captured = capsys.readouterr()
    assert exit_code == 0
    assert __version__ in captured.out


def test_unknown_subcommand_errors(capsys) -> None:
    exit_code = cli.main(["frobnicate"])
    captured = capsys.readouterr()
    assert exit_code == 2
    assert "unknown subcommand" in captured.err


# ---------- legacy backward-compat fallback (removed in v0.6.0) -------------


def test_legacy_no_subcommand_invocation_fails_with_clear_message(capsys) -> None:
    """Pre-v0.6.0 'mitoribopy -s h -f ref.fa ...' was silently routed to
    'mitoribopy rpf ...' with a deprecation warning. v0.6.0 removed that
    fallback for the publication freeze; the call must now exit with a
    non-zero status and a message that points at the right subcommand."""
    exit_code = cli.main(["-s", "h", "-f", "tiny.fa"])
    captured = capsys.readouterr()

    assert exit_code == 2
    assert "missing subcommand" in captured.err
    assert "mitoribopy rpf" in captured.err


# ---------- rpf subcommand ---------------------------------------------------


def test_rpf_subcommand_routes_to_pipeline_cli(monkeypatch) -> None:
    captured_argv: dict[str, list[str]] = {}

    def fake_run_pipeline_cli(argv):
        captured_argv["argv"] = list(argv)
        return 0

    monkeypatch.setattr(cli, "run_pipeline_cli", fake_run_pipeline_cli)

    exit_code = cli.main(["rpf", "-s", "h", "-f", "tiny.fa"])
    assert exit_code == 0
    assert captured_argv["argv"] == ["-s", "h", "-f", "tiny.fa"]


def test_rpf_subcommand_strips_common_args_before_pipeline(monkeypatch) -> None:
    captured_argv: dict[str, list[str]] = {}

    def fake_run_pipeline_cli(argv):
        captured_argv["argv"] = list(argv)
        return 0

    monkeypatch.setattr(cli, "run_pipeline_cli", fake_run_pipeline_cli)

    cli.main(
        [
            "rpf",
            "--threads",
            "2",
            "--log-level",
            "WARNING",
            "-s",
            "h",
            "-f",
            "tiny.fa",
        ]
    )

    assert captured_argv["argv"] == ["-s", "h", "-f", "tiny.fa"]
    assert os.environ.get("MITORIBOPY_THREADS") == "2"


def test_rpf_dry_run_prints_plan_and_exits_zero(capsys, monkeypatch) -> None:
    def boom(argv):
        raise AssertionError(
            "pipeline must not be invoked when --dry-run is set"
        )

    monkeypatch.setattr(cli, "run_pipeline_cli", boom)

    exit_code = cli.main(["rpf", "--dry-run"])
    captured = capsys.readouterr()

    assert exit_code == 0
    assert "dry-run" in captured.out
    assert "offset enrichment" in captured.out


# Subcommand stubs retired: 'align', 'rnaseq', and 'all' are now fully
# implemented and have their own dedicated test modules
# (tests/test_align_cli.py, tests/test_rnaseq_cli.py,
# tests/test_all_cli.py).


# ---------- common args: --threads, --log-level ------------------------------


def test_common_apply_sets_thread_env_vars(monkeypatch) -> None:
    for env_key in (
        "MITORIBOPY_THREADS",
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
    ):
        monkeypatch.delenv(env_key, raising=False)

    parser = __import__("argparse").ArgumentParser()
    common.add_common_arguments(parser)
    args = parser.parse_args(["--threads", "4"])
    common.apply_common_arguments(args)

    for env_key in (
        "MITORIBOPY_THREADS",
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
    ):
        assert os.environ[env_key] == "4"


def test_common_apply_rejects_non_positive_threads() -> None:
    parser = __import__("argparse").ArgumentParser()
    common.add_common_arguments(parser)
    args = parser.parse_args(["--threads", "0"])
    with pytest.raises(SystemExit):
        common.apply_common_arguments(args)


def test_common_apply_sets_logger_level(monkeypatch) -> None:
    import logging

    from mitoribopy.console import LOGGER

    monkeypatch.setattr(LOGGER, "level", logging.INFO, raising=False)

    parser = __import__("argparse").ArgumentParser()
    common.add_common_arguments(parser)
    args = parser.parse_args(["--log-level", "WARNING"])
    common.apply_common_arguments(args)

    assert LOGGER.level == logging.WARNING


# ---------- config file loading: JSON / YAML / TOML --------------------------


def test_load_user_config_reads_yaml(tmp_path) -> None:
    cfg_path = tmp_path / "cfg.yaml"
    cfg_path.write_text(
        "strain: h\nalign: stop\noffset_mask_nt: 7\n",
        encoding="utf-8",
    )

    cfg = load_user_config(str(cfg_path))

    assert cfg == {"strain": "h", "align": "stop", "offset_mask_nt": 7}


def test_load_user_config_reads_toml(tmp_path) -> None:
    cfg_path = tmp_path / "cfg.toml"
    cfg_path.write_text(
        'strain = "h"\nalign = "stop"\noffset_mask_nt = 7\n',
        encoding="utf-8",
    )

    cfg = load_user_config(str(cfg_path))

    assert cfg == {"strain": "h", "align": "stop", "offset_mask_nt": 7}


def test_load_user_config_still_reads_json(tmp_path) -> None:
    cfg_path = tmp_path / "cfg.json"
    cfg_path.write_text(
        json.dumps({"strain": "y", "rpm_norm_mode": "mt_mrna"}),
        encoding="utf-8",
    )

    cfg = load_user_config(str(cfg_path))

    assert cfg == {"strain": "y", "rpm_norm_mode": "mt_mrna"}


def test_common_load_config_file_yaml_round_trip(tmp_path) -> None:
    cfg_path = tmp_path / "cfg.yml"
    cfg_path.write_text(
        "offset_site: a\nmin_5_offset: 10\n",
        encoding="utf-8",
    )

    data = common.load_config_file(str(cfg_path))

    assert data == {"offset_site": "a", "min_5_offset": 10}
