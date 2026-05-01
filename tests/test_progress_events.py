"""Tests for the event-driven progress system (refactor-4 / §3)."""

from __future__ import annotations

import io
import json
from pathlib import Path

import pytest

from mitoribopy.progress import (
    JsonlRenderer,
    PlainRenderer,
    ProgressManager,
    ProgressMode,
    RunEnd,
    RunStart,
    SampleStepEnd,
    StageEnd,
    StageStart,
    WarningEvent,
    resolve_mode,
)


# ---------------------------------------------------------------------------
# Mode resolution
# ---------------------------------------------------------------------------


class TestResolveMode:
    def test_auto_picks_bar_on_tty(self) -> None:
        class _Tty:
            def isatty(self) -> bool:
                return True

        assert resolve_mode("auto", stream=_Tty()) == ProgressMode.BAR

    def test_auto_picks_plain_on_non_tty(self) -> None:
        buf = io.StringIO()
        assert resolve_mode("auto", stream=buf) == ProgressMode.PLAIN

    def test_none_is_auto(self) -> None:
        buf = io.StringIO()
        assert resolve_mode(None, stream=buf) == ProgressMode.PLAIN

    def test_explicit_modes(self) -> None:
        buf = io.StringIO()
        for m in ("plain", "bar", "rich", "jsonl", "off"):
            assert resolve_mode(m, stream=buf) == ProgressMode(m)

    def test_unknown_falls_back_to_plain(self) -> None:
        buf = io.StringIO()
        assert resolve_mode("nonsense", stream=buf) == ProgressMode.PLAIN


# ---------------------------------------------------------------------------
# Event payload contracts
# ---------------------------------------------------------------------------


class TestEventToDict:
    def test_run_start_emits_event_and_timestamp(self) -> None:
        evt = RunStart(subcommand="all", n_samples=4)
        d = evt.to_dict()
        assert d["event"] == "run_start"
        assert d["subcommand"] == "all"
        assert d["n_samples"] == 4
        assert "timestamp" in d

    def test_optional_fields_dropped_when_none(self) -> None:
        evt = StageStart(stage="rpf")
        d = evt.to_dict()
        assert "n_samples" not in d
        assert d["event"] == "stage_start"

    def test_warning_event_round_trips_payload(self) -> None:
        evt = WarningEvent(
            stage="RNASEQ", sample_id="S1", code="X", message="m",
        )
        d = evt.to_dict()
        assert d["stage"] == "RNASEQ"
        assert d["sample_id"] == "S1"
        assert d["code"] == "X"


# ---------------------------------------------------------------------------
# ProgressManager — fan-out + history
# ---------------------------------------------------------------------------


class TestProgressManagerEmit:
    def test_null_manager_records_history(self) -> None:
        mgr = ProgressManager.null()
        mgr.run_start(subcommand="all")
        mgr.stage_start("align")
        mgr.stage_end("align", status="done", elapsed_seconds=1.0)
        mgr.run_end(status="done", elapsed_seconds=1.0)
        events = mgr.events()
        assert [type(e).__name__ for e in events] == [
            "RunStart", "StageStart", "StageEnd", "RunEnd",
        ]

    def test_emit_forwards_to_every_renderer(self) -> None:
        seen_a: list = []
        seen_b: list = []

        class _Capture:
            def __init__(self, target: list) -> None:
                self.target = target

            def handle(self, event) -> None:
                self.target.append(event)

            def close(self) -> None:
                pass

        mgr = ProgressManager([_Capture(seen_a), _Capture(seen_b)])
        mgr.stage_start("align")
        mgr.close()
        assert len(seen_a) == 1
        assert len(seen_b) == 1


# ---------------------------------------------------------------------------
# Renderers
# ---------------------------------------------------------------------------


class TestPlainRenderer:
    def test_one_line_per_event(self) -> None:
        buf = io.StringIO()
        renderer = PlainRenderer(buf)
        renderer.handle(RunStart(subcommand="all", n_samples=2))
        renderer.handle(StageEnd(stage="align", status="done", elapsed_seconds=12.4))
        renderer.close()
        lines = [
            line for line in buf.getvalue().splitlines() if line.strip()
        ]
        assert len(lines) == 2
        assert lines[0].startswith("[PROGRESS] event=run_start")
        assert "stage=align" in lines[1]
        assert "status=done" in lines[1]
        # Duration is human-formatted, not raw seconds.
        assert "elapsed_seconds=12.4s" in lines[1]


class TestJsonlRenderer:
    def test_writes_one_object_per_event(self, tmp_path: Path) -> None:
        path = tmp_path / "progress.jsonl"
        renderer = JsonlRenderer(path)
        renderer.handle(RunStart(subcommand="all"))
        renderer.handle(SampleStepEnd(
            stage="align", sample_id="S1", step="trim",
            elapsed_seconds=1.2, reads_in=100, reads_out=95,
        ))
        renderer.handle(RunEnd(status="done", elapsed_seconds=2.0))
        renderer.close()

        lines = path.read_text(encoding="utf-8").splitlines()
        assert len(lines) == 3
        first = json.loads(lines[0])
        assert first["event"] == "run_start"
        second = json.loads(lines[1])
        assert second["sample_id"] == "S1"
        assert second["reads_out"] == 95


class TestProgressManagerFromCli:
    def test_progress_file_attaches_jsonl(self, tmp_path: Path) -> None:
        target = tmp_path / "p.jsonl"
        mgr = ProgressManager.from_cli(
            mode="plain",
            progress_file=target,
            stream=io.StringIO(),
        )
        mgr.run_start(subcommand="all")
        mgr.run_end(status="done", elapsed_seconds=0.0)
        mgr.close()
        assert target.is_file()
        lines = target.read_text(encoding="utf-8").splitlines()
        assert len(lines) == 2
        assert json.loads(lines[0])["event"] == "run_start"

    def test_off_silences_console(self, tmp_path: Path) -> None:
        buf = io.StringIO()
        mgr = ProgressManager.from_cli(mode="off", stream=buf)
        mgr.run_start(subcommand="all")
        mgr.close()
        assert buf.getvalue() == ""
