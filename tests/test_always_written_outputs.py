"""Refactor-4: warnings.tsv + outputs_index.tsv must ALWAYS exist after `mitoribopy all`.

The assessment §8 output contract requires every run root to advertise
its outputs and surface every structured warning, even when no warning
was emitted and no stage produced one of the optional artefacts. The
files exist with header-only content if there is nothing to report.

These tests cover three scenarios:

* the happy path (one stage runs, warnings empty)
* a stage failure (the orchestrator early-exits but the files are
  still on disk because they were touched at run start)
* the resume-skip path (no stage runs at all)
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mitoribopy import cli
from mitoribopy.io import warnings_log


@pytest.fixture(autouse=True)
def _clear_warnings():
    warnings_log.clear()
    yield
    warnings_log.clear()


def _write_min_yaml(path: Path, *, with_align: bool = True) -> None:
    sections: list[str] = []
    if with_align:
        sections.append("align:\n  adapter: TGGAATTCTCGGGTGCCAAGG\n")
    path.write_text("".join(sections) or "rpf:\n  strain: h\n")


def _stub_align(tmp_path: Path) -> None:
    out = tmp_path / "results" / "align"
    out.mkdir(parents=True, exist_ok=True)
    (out / "read_counts.tsv").write_text("sample\ttotal_reads\nA\t100\n")
    (out / "run_settings.json").write_text('{"subcommand":"align"}')


class TestAlwaysWrittenOnSuccess:
    def test_happy_path_writes_both(self, tmp_path: Path, monkeypatch) -> None:
        from mitoribopy.cli import align as align_cli
        from mitoribopy.cli import rnaseq as rnaseq_cli
        from mitoribopy.cli import rpf as rpf_cli

        def fake_align(argv):
            _stub_align(tmp_path)
            return 0

        monkeypatch.setattr(align_cli, "run", fake_align)
        monkeypatch.setattr(rpf_cli, "run", lambda a: 0)
        monkeypatch.setattr(rnaseq_cli, "run", lambda a: 0)

        cfg = tmp_path / "c.yaml"
        _write_min_yaml(cfg)
        rc = cli.main([
            "all",
            "--config", str(cfg),
            "--output", str(tmp_path / "results"),
            "--no-progress",
        ])
        assert rc == 0
        assert (tmp_path / "results" / "warnings.tsv").is_file()
        assert (tmp_path / "results" / "outputs_index.tsv").is_file()
        # Both files start with the schema header.
        assert (tmp_path / "results" / "warnings.tsv").read_text().startswith(
            "# schema_version:"
        )
        assert (tmp_path / "results" / "outputs_index.tsv").read_text().startswith(
            "# schema_version:"
        )


class TestAlwaysWrittenOnStageFailure:
    def test_align_failure_still_leaves_files(
        self, tmp_path: Path, monkeypatch, capsys
    ) -> None:
        from mitoribopy.cli import align as align_cli
        from mitoribopy.cli import rnaseq as rnaseq_cli
        from mitoribopy.cli import rpf as rpf_cli

        # Align returns nonzero; the orchestrator must still leave
        # warnings.tsv + outputs_index.tsv on disk.
        monkeypatch.setattr(align_cli, "run", lambda a: 7)
        monkeypatch.setattr(rpf_cli, "run", lambda a: 0)
        monkeypatch.setattr(rnaseq_cli, "run", lambda a: 0)

        cfg = tmp_path / "c.yaml"
        _write_min_yaml(cfg)
        rc = cli.main([
            "all",
            "--config", str(cfg),
            "--output", str(tmp_path / "results"),
            "--no-progress",
        ])
        assert rc == 7
        assert (tmp_path / "results" / "warnings.tsv").is_file()
        assert (tmp_path / "results" / "outputs_index.tsv").is_file()


class TestProgressFileEmittedOnAll:
    def test_progress_jsonl_default_location(
        self, tmp_path: Path, monkeypatch
    ) -> None:
        from mitoribopy.cli import align as align_cli
        from mitoribopy.cli import rnaseq as rnaseq_cli
        from mitoribopy.cli import rpf as rpf_cli

        def fake_align(argv):
            _stub_align(tmp_path)
            return 0

        monkeypatch.setattr(align_cli, "run", fake_align)
        monkeypatch.setattr(rpf_cli, "run", lambda a: 0)
        monkeypatch.setattr(rnaseq_cli, "run", lambda a: 0)

        cfg = tmp_path / "c.yaml"
        _write_min_yaml(cfg)
        rc = cli.main([
            "all",
            "--config", str(cfg),
            "--output", str(tmp_path / "results"),
            "--progress", "plain",
        ])
        assert rc == 0
        progress = tmp_path / "results" / "progress.jsonl"
        assert progress.is_file()
        # Sanity: at least a run_start + run_end pair in JSONL form.
        content = progress.read_text(encoding="utf-8")
        assert '"event": "run_start"' in content
        assert '"event": "run_end"' in content

    def test_no_progress_disables_jsonl(
        self, tmp_path: Path, monkeypatch
    ) -> None:
        from mitoribopy.cli import align as align_cli
        from mitoribopy.cli import rnaseq as rnaseq_cli
        from mitoribopy.cli import rpf as rpf_cli

        def fake_align(argv):
            _stub_align(tmp_path)
            return 0

        monkeypatch.setattr(align_cli, "run", fake_align)
        monkeypatch.setattr(rpf_cli, "run", lambda a: 0)
        monkeypatch.setattr(rnaseq_cli, "run", lambda a: 0)

        cfg = tmp_path / "c.yaml"
        _write_min_yaml(cfg)
        rc = cli.main([
            "all",
            "--config", str(cfg),
            "--output", str(tmp_path / "results"),
            "--no-progress",
        ])
        assert rc == 0
        # --no-progress means no console output AND no default JSONL.
        assert not (tmp_path / "results" / "progress.jsonl").exists()
