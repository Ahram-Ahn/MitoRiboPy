"""Tests for the hash-validated resume guard (P0.5).

The orchestrator must NOT skip a stage on --resume when the user has
edited the config / sample sheet / reference between runs. Hash drift
detection is mandatory; --force-resume is the explicit override.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from mitoribopy import cli
from mitoribopy.cli import _resume_guard as guard


# ---------- sha256_of -------------------------------------------------------


def test_sha256_of_returns_none_for_missing_path(tmp_path: Path) -> None:
    assert guard.sha256_of(tmp_path / "nope.txt") is None
    assert guard.sha256_of(None) is None


def test_sha256_of_is_stable(tmp_path: Path) -> None:
    p = tmp_path / "x.bin"
    p.write_bytes(b"hello world\n")
    h1 = guard.sha256_of(p)
    h2 = guard.sha256_of(p)
    assert h1 == h2
    assert isinstance(h1, str) and len(h1) == 64


# ---------- force_resume_requested ------------------------------------------


def test_force_resume_cli_flag_wins(monkeypatch) -> None:
    monkeypatch.delenv(guard.FORCE_RESUME_ENV, raising=False)
    assert guard.force_resume_requested(cli_flag=True) is True
    assert guard.force_resume_requested(cli_flag=False) is False


def test_force_resume_env_var_fallback(monkeypatch) -> None:
    monkeypatch.setenv(guard.FORCE_RESUME_ENV, "1")
    assert guard.force_resume_requested(cli_flag=False) is True
    monkeypatch.setenv(guard.FORCE_RESUME_ENV, "true")
    assert guard.force_resume_requested(cli_flag=False) is True
    monkeypatch.setenv(guard.FORCE_RESUME_ENV, "0")
    assert guard.force_resume_requested(cli_flag=False) is False


# ---------- load_prior_manifest ---------------------------------------------


def test_load_prior_manifest_returns_none_when_absent(tmp_path: Path) -> None:
    assert guard.load_prior_manifest(tmp_path, "run_manifest.json") is None


def test_load_prior_manifest_returns_none_for_corrupt_json(tmp_path: Path) -> None:
    (tmp_path / "run_manifest.json").write_text("{not json}")
    assert guard.load_prior_manifest(tmp_path, "run_manifest.json") is None


def test_load_prior_manifest_loads_valid_manifest(tmp_path: Path) -> None:
    payload = {"schema_version": "1.0.0", "mitoribopy_version": "0.5.1"}
    (tmp_path / "run_manifest.json").write_text(json.dumps(payload))
    loaded = guard.load_prior_manifest(tmp_path, "run_manifest.json")
    assert loaded == payload


# ---------- validate_resume -------------------------------------------------


def _write(path: Path, content: bytes = b"hello") -> Path:
    path.write_bytes(content)
    return path


def test_validate_resume_no_prior_manifest_is_first_run() -> None:
    report = guard.validate_resume(
        prior_manifest=None,
        config_path=None,
        sample_sheet_path=None,
        reference_fasta=None,
        mitoribopy_version="0.5.1",
        manifest_schema_version="1.0.0",
    )
    assert report.no_prior_manifest is True
    assert report.ok is True


def test_validate_resume_all_match_returns_ok(tmp_path: Path) -> None:
    cfg = _write(tmp_path / "c.yaml", b"align: {}\n")
    sheet = _write(tmp_path / "samples.tsv", b"sample_id\n")
    ref = _write(tmp_path / "ref.fa", b">x\nA\n")
    prior = {
        "config_source_sha256": guard.sha256_of(cfg),
        "sample_sheet_sha256": guard.sha256_of(sheet),
        "reference_checksum": guard.sha256_of(ref),
        "mitoribopy_version": "0.5.1",
        "schema_version": "1.0.0",
    }
    report = guard.validate_resume(
        prior_manifest=prior,
        config_path=cfg,
        sample_sheet_path=sheet,
        reference_fasta=ref,
        mitoribopy_version="0.5.1",
        manifest_schema_version="1.0.0",
    )
    assert report.ok is True
    assert report.mismatches == ()


def test_validate_resume_config_drift_detected(tmp_path: Path) -> None:
    cfg = _write(tmp_path / "c.yaml", b"v1\n")
    prior = {"config_source_sha256": "deadbeef" * 8, "mitoribopy_version": "0.5.1"}
    report = guard.validate_resume(
        prior_manifest=prior,
        config_path=cfg,
        sample_sheet_path=None,
        reference_fasta=None,
        mitoribopy_version="0.5.1",
        manifest_schema_version="1.0.0",
    )
    assert not report.ok
    fields = [m.field for m in report.mismatches]
    assert "config_source_sha256" in fields


def test_validate_resume_sample_sheet_drift_detected(tmp_path: Path) -> None:
    sheet = _write(tmp_path / "s.tsv", b"new contents\n")
    prior = {"sample_sheet_sha256": "feedbeef" * 8, "mitoribopy_version": "0.5.1"}
    report = guard.validate_resume(
        prior_manifest=prior,
        config_path=None,
        sample_sheet_path=sheet,
        reference_fasta=None,
        mitoribopy_version="0.5.1",
        manifest_schema_version="1.0.0",
    )
    assert not report.ok
    assert any(m.field == "sample_sheet_sha256" for m in report.mismatches)


def test_validate_resume_reference_drift_detected(tmp_path: Path) -> None:
    ref = _write(tmp_path / "r.fa", b">x\nC\n")
    prior = {"reference_checksum": "cafef00d" * 8, "mitoribopy_version": "0.5.1"}
    report = guard.validate_resume(
        prior_manifest=prior,
        config_path=None,
        sample_sheet_path=None,
        reference_fasta=ref,
        mitoribopy_version="0.5.1",
        manifest_schema_version="1.0.0",
    )
    assert not report.ok
    assert any(m.field == "reference_checksum" for m in report.mismatches)


def test_validate_resume_version_drift_detected() -> None:
    prior = {"mitoribopy_version": "0.5.1", "schema_version": "1.0.0"}
    report = guard.validate_resume(
        prior_manifest=prior,
        config_path=None,
        sample_sheet_path=None,
        reference_fasta=None,
        mitoribopy_version="0.5.2",
        manifest_schema_version="1.0.0",
    )
    assert not report.ok
    assert any(m.field == "mitoribopy_version" for m in report.mismatches)


def test_validate_resume_skips_absent_prior_fields(tmp_path: Path) -> None:
    """An older manifest missing some hash fields must not produce
    spurious mismatches — we only flag fields the prior manifest set."""
    cfg = _write(tmp_path / "c.yaml")
    prior = {"mitoribopy_version": "0.5.1"}  # no SHA fields
    report = guard.validate_resume(
        prior_manifest=prior,
        config_path=cfg,
        sample_sheet_path=tmp_path / "missing.tsv",
        reference_fasta=tmp_path / "missing.fa",
        mitoribopy_version="0.5.1",
        manifest_schema_version="1.0.0",
    )
    assert report.ok is True


# ---------- end-to-end via mitoribopy all -----------------------------------


def _make_prior_run(tmp_path: Path, *, config_text: str) -> tuple[Path, Path]:
    """Build a fake prior `mitoribopy all` output dir with sentinels."""
    cfg = tmp_path / "pipeline.yaml"
    cfg.write_text(config_text)
    run_root = tmp_path / "results"
    (run_root / "align").mkdir(parents=True, exist_ok=True)
    (run_root / "rpf").mkdir(parents=True, exist_ok=True)
    (run_root / "align" / "read_counts.tsv").write_text("sample\n")
    (run_root / "align" / "run_settings.json").write_text("{}")
    (run_root / "rpf" / "rpf_counts.tsv").write_text("sample\tgene\tcount\n")
    (run_root / "rpf" / "run_settings.json").write_text("{}")
    return cfg, run_root


def test_all_resume_succeeds_when_config_unchanged(
    tmp_path: Path, monkeypatch
) -> None:
    """Sentinels exist + hashes match prior manifest -> resume skips both stages."""
    cfg_text = "align:\n  # adapter auto-detection (default)\nrpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
    cfg, run_root = _make_prior_run(tmp_path, config_text=cfg_text)

    # Fake prior manifest with the SHA of the current config.
    from mitoribopy import __version__
    from mitoribopy.cli.all_ import MANIFEST_SCHEMA_VERSION
    manifest = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "mitoribopy_version": __version__,
        "config_source_sha256": guard.sha256_of(cfg),
    }
    (run_root / "run_manifest.json").write_text(json.dumps(manifest))

    align_called = []
    rpf_called = []

    def fake_align(argv):
        align_called.append(argv)
        return 0

    def fake_rpf(argv):
        rpf_called.append(argv)
        return 0

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rpf as rpf_cli
    monkeypatch.setattr(align_cli, "run", fake_align)
    monkeypatch.setattr(rpf_cli, "run", fake_rpf)

    exit_code = cli.main(
        ["all", "--config", str(cfg), "--output", str(run_root), "--resume"]
    )
    assert exit_code == 0
    # Both stages were skipped (sentinels exist).
    assert align_called == []
    assert rpf_called == []


def test_all_resume_blocks_when_config_changed(
    tmp_path: Path, capsys
) -> None:
    """Edit the config after the prior run -> --resume must refuse to skip."""
    cfg_text = "align:\n  # adapter auto-detection (default)\nrpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
    cfg, run_root = _make_prior_run(tmp_path, config_text=cfg_text)

    from mitoribopy import __version__
    from mitoribopy.cli.all_ import MANIFEST_SCHEMA_VERSION
    manifest = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "mitoribopy_version": __version__,
        "config_source_sha256": "0" * 64,  # a deliberate mismatch
    }
    (run_root / "run_manifest.json").write_text(json.dumps(manifest))

    exit_code = cli.main(
        ["all", "--config", str(cfg), "--output", str(run_root), "--resume"]
    )
    assert exit_code == 2
    err = capsys.readouterr().err
    assert "config_source_sha256" in err
    assert "--force-resume" in err


def test_all_force_resume_bypasses_hash_guard(
    tmp_path: Path, capsys, monkeypatch
) -> None:
    cfg_text = "align:\n  # adapter auto-detection (default)\nrpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
    cfg, run_root = _make_prior_run(tmp_path, config_text=cfg_text)

    from mitoribopy import __version__
    from mitoribopy.cli.all_ import MANIFEST_SCHEMA_VERSION
    manifest = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "mitoribopy_version": __version__,
        "config_source_sha256": "0" * 64,
    }
    (run_root / "run_manifest.json").write_text(json.dumps(manifest))

    align_called = []
    rpf_called = []

    def fake_align(argv):
        align_called.append(argv)
        return 0

    def fake_rpf(argv):
        rpf_called.append(argv)
        return 0

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rpf as rpf_cli
    monkeypatch.setattr(align_cli, "run", fake_align)
    monkeypatch.setattr(rpf_cli, "run", fake_rpf)

    exit_code = cli.main(
        [
            "all",
            "--config", str(cfg),
            "--output", str(run_root),
            "--force-resume",
        ]
    )
    assert exit_code == 0
    err = capsys.readouterr().err
    # Warning printed but the stages were still skipped.
    assert "WARNING" in err
    assert "bypassing hash guard" in err
    assert align_called == []
    assert rpf_called == []


def test_all_force_resume_env_var_bypasses_hash_guard(
    tmp_path: Path, capsys, monkeypatch
) -> None:
    """MITORIBOPY_FORCE_RESUME=1 has the same effect as --force-resume."""
    cfg_text = "align:\n  # adapter auto-detection (default)\nrpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
    cfg, run_root = _make_prior_run(tmp_path, config_text=cfg_text)

    from mitoribopy import __version__
    from mitoribopy.cli.all_ import MANIFEST_SCHEMA_VERSION
    manifest = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "mitoribopy_version": __version__,
        "config_source_sha256": "0" * 64,
    }
    (run_root / "run_manifest.json").write_text(json.dumps(manifest))

    monkeypatch.setenv(guard.FORCE_RESUME_ENV, "1")

    align_called = []
    rpf_called = []

    def fake_align(argv):
        align_called.append(argv)
        return 0

    def fake_rpf(argv):
        rpf_called.append(argv)
        return 0

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rpf as rpf_cli
    monkeypatch.setattr(align_cli, "run", fake_align)
    monkeypatch.setattr(rpf_cli, "run", fake_rpf)

    exit_code = cli.main(
        ["all", "--config", str(cfg), "--output", str(run_root), "--resume"]
    )
    assert exit_code == 0
    # Sentinels still applied, env-var bypass let us through.
    assert align_called == []
    assert rpf_called == []


def test_all_resume_first_run_no_prior_manifest_is_ok(
    tmp_path: Path, monkeypatch
) -> None:
    """Sentinels exist but no prior manifest (e.g. interrupted older run);
    fall back to the historical sentinel-only behaviour."""
    cfg_text = "align:\n  # adapter auto-detection (default)\nrpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
    cfg, run_root = _make_prior_run(tmp_path, config_text=cfg_text)
    # Note: no run_manifest.json written.

    align_called = []
    rpf_called = []

    def fake_align(argv):
        align_called.append(argv)
        return 0

    def fake_rpf(argv):
        rpf_called.append(argv)
        return 0

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rpf as rpf_cli
    monkeypatch.setattr(align_cli, "run", fake_align)
    monkeypatch.setattr(rpf_cli, "run", fake_rpf)

    exit_code = cli.main(
        ["all", "--config", str(cfg), "--output", str(run_root), "--resume"]
    )
    assert exit_code == 0
    assert align_called == []
    assert rpf_called == []
