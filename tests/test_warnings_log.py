"""Tests for the structured warnings log (P1.11 / P1.12).

Every warning emitted via ``warnings_log.record`` must:

* be visible to the user via the existing console.log_warning channel
  (no behaviour regression for users who only watch stderr);
* be appended to the in-process record list so subsequent ``mitoribopy
  all`` orchestration can flush it;
* end up in ``warnings.tsv`` at the run root and in
  ``run_manifest.json`` under the ``warnings`` key.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from mitoribopy.io import warnings_log


@pytest.fixture(autouse=True)
def _isolate_warnings_log():
    """Clear the global record list before AND after every test so tests
    that exercise the same module don't pollute each other."""
    warnings_log.clear()
    yield
    warnings_log.clear()


# ---------- record + collected ---------------------------------------------


def test_record_appends_to_collected_list(capsys) -> None:
    rec = warnings_log.record(
        "RNASEQ",
        "boom",
        sample="S1",  # legacy alias still accepted
        code="UMI_INFERRED_NO_DECLARATION",
    )
    # Canonical attribute names (P5.8).
    assert rec.stage == "RNASEQ"
    assert rec.sample_id == "S1"
    # Backward-compat aliases still resolve.
    assert rec.component == "RNASEQ"
    assert rec.sample == "S1"
    assert rec.code == "UMI_INFERRED_NO_DECLARATION"
    assert rec.severity == "warn"
    assert rec.message == "boom"
    assert rec.suggested_action is None

    all_records = warnings_log.collected()
    assert len(all_records) == 1
    assert all_records[0] == rec


def test_record_mirrors_to_console_warning() -> None:
    """User-visible stderr behaviour preserved: log_warning is invoked.

    The mitoribopy LOGGER is built with a StreamHandler bound to the
    original sys.stdout at module-import time, so pytest's capture
    fixtures (which reassign sys.stdout) cannot intercept the line
    written by the handler. We instead assert that the LOGGER's
    `warning` method was called by attaching a one-shot handler.
    """
    import logging

    from mitoribopy.console import LOGGER

    captured: list[str] = []

    class _MemoryHandler(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:  # noqa: D401
            captured.append(self.format(record))

    handler = _MemoryHandler(level=logging.WARNING)
    handler.setFormatter(logging.Formatter("%(message)s"))
    LOGGER.addHandler(handler)
    try:
        warnings_log.record("RNASEQ", "boom", code="X")
    finally:
        LOGGER.removeHandler(handler)
    blob = " ".join(captured)
    assert "WARNING" in blob
    assert "RNASEQ" in blob
    assert "boom" in blob


def test_clear_resets_records() -> None:
    warnings_log.record("X", "msg")
    assert warnings_log.collected()
    warnings_log.clear()
    assert warnings_log.collected() == []


# ---------- flush_tsv -------------------------------------------------------


def test_flush_tsv_writes_header_and_one_row_per_record(tmp_path: Path) -> None:
    warnings_log.record(
        "RNASEQ",
        "boom",
        sample_id="S1",
        code="UMI",
        suggested_action="declare umi_length in samples.tsv",
    )
    warnings_log.record("ALIGN", "halt", code="DEDUP")
    out = warnings_log.flush_tsv(tmp_path / "warnings.tsv")
    assert out.exists()
    raw = out.read_text().splitlines()
    assert raw[0].startswith("# schema_version: ")
    rows = [r for r in raw if not r.startswith("#")]
    header = rows[0].split("\t")
    # P5.8: assessment §8 spec — no timestamp; suggested_action added.
    assert header == [
        "stage",
        "sample_id",
        "severity",
        "code",
        "message",
        "suggested_action",
    ]
    body = [r.split("\t") for r in rows[1:]]
    assert len(body) == 2
    assert body[0][0] == "RNASEQ"
    assert body[0][1] == "S1"
    assert body[0][3] == "UMI"
    assert body[0][5] == "declare umi_length in samples.tsv"
    assert body[1][0] == "ALIGN"
    assert body[1][5] == ""


def test_flush_tsv_creates_header_only_file_when_no_warnings(tmp_path: Path) -> None:
    out = warnings_log.flush_tsv(tmp_path / "warnings.tsv")
    raw = out.read_text().splitlines()
    rows = [r for r in raw if not r.startswith("#")]
    assert len(rows) == 1  # header only
    assert "stage" in rows[0]
    assert "sample_id" in rows[0]
    assert "suggested_action" in rows[0]


def test_flush_tsv_strips_tabs_and_newlines_from_messages(tmp_path: Path) -> None:
    """Messages with embedded tabs / newlines must not corrupt the TSV."""
    warnings_log.record("X", "first\tsecond\nthird")
    out = warnings_log.flush_tsv(tmp_path / "warnings.tsv")
    raw = out.read_text().splitlines()
    rows = [r for r in raw if not r.startswith("#")]
    # header + one data row
    assert len(rows) == 2
    cells = rows[1].split("\t")
    # P5.8: columns are stage, sample_id, severity, code, message,
    # suggested_action — message is index 4.
    assert "\n" not in cells[4]
    # tab-replacement: the original "first\tsecond\nthird" becomes
    # "first second third".
    assert cells[4] == "first second third"
    # Empty suggested_action still scrubbed cleanly.
    assert cells[5] == ""


# ---------- end-to-end via mitoribopy all -----------------------------------


def test_mitoribopy_all_writes_warnings_tsv_and_manifest(
    tmp_path: Path, monkeypatch
) -> None:
    """A real run must emit warnings.tsv at the run root and embed the
    same records in run_manifest.json's `warnings` field."""
    from mitoribopy import cli

    # Inject a warning through the real API mid-run (simulating a
    # downstream stage detecting an inferred UMI without declaration).
    def fake_align(argv):
        warnings_log.record(
            "ALIGN",
            "fake umi inference for test",
            sample="S1",
            code="UMI_INFERRED_NO_DECLARATION",
        )
        out = Path(tmp_path / "results" / "align")
        out.mkdir(parents=True, exist_ok=True)
        (out / "read_counts.tsv").write_text(
            "# schema_version: 1.0\nsample\n"
        )
        (out / "run_settings.json").write_text("{}")
        return 0

    def fake_rpf(argv):
        out = Path(tmp_path / "results" / "rpf")
        out.mkdir(parents=True, exist_ok=True)
        (out / "rpf_counts.tsv").write_text(
            "# schema_version: 1.0\nsample\tgene\tcount\n"
        )
        (out / "run_settings.json").write_text("{}")
        return 0

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rpf as rpf_cli
    monkeypatch.setattr(align_cli, "run", fake_align)
    monkeypatch.setattr(rpf_cli, "run", fake_rpf)

    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "align:\n  dedup_strategy: auto\n"
        "rpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
    )

    rc = cli.main([
        "all", "--config", str(cfg), "--output", str(tmp_path / "results"),
    ])
    assert rc == 0

    # File side: warnings.tsv exists with one record.
    warnings_path = tmp_path / "results" / "warnings.tsv"
    assert warnings_path.exists()
    rows = [
        r for r in warnings_path.read_text().splitlines() if not r.startswith("#")
    ]
    assert len(rows) >= 2  # header + at least one data row
    assert any("UMI_INFERRED_NO_DECLARATION" in r for r in rows)

    # Manifest side: same record(s) in run_manifest.json.
    manifest = json.loads(
        (tmp_path / "results" / "run_manifest.json").read_text()
    )
    warns = manifest.get("warnings")
    assert isinstance(warns, list)
    assert any(
        w.get("code") == "UMI_INFERRED_NO_DECLARATION" for w in warns
    )
