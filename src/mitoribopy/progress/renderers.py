"""Pluggable renderers for ``ProgressEvent`` streams.

Three concrete renderers ship with the package:

* :class:`PlainRenderer` — one ``[PROGRESS] key=val ...`` line per event,
  HPC-safe (no ANSI, no rewrites). The default on non-TTY streams.
* :class:`TqdmRenderer` — compact tqdm sample counters per stage.
  Falls back to plain output when tqdm is missing or the stream is not
  a TTY. Per-event payload is also forwarded to plain so users keep the
  log trail in addition to the bar.
* :class:`JsonlRenderer` — one JSON object per line into a file (or
  passed file handle). Captures EVERY event verbatim — the most useful
  channel for downstream tooling.

A renderer must implement two methods:

* ``handle(event)`` — called by the manager for each event.
* ``close()`` — called by the manager on ``ProgressManager.close()``.

Renderers should never raise — failures are non-fatal and should be
swallowed (logged to stderr at most). The progress system is observability
only; it must not break a real run.
"""

from __future__ import annotations

import io
import json
import sys
from contextlib import suppress
from pathlib import Path
from typing import Protocol, TextIO

from .events import (
    ErrorEvent,
    OutputEvent,
    ProgressEvent,
    ResumeSkipEvent,
    RunEnd,
    RunStart,
    SampleEnd,
    SampleStart,
    SampleStepEnd,
    SampleStepStart,
    StageEnd,
    StageStart,
    WarningEvent,
)
from .timing import HAS_TQDM, format_duration


__all__ = [
    "JsonlRenderer",
    "PlainRenderer",
    "ProgressRenderer",
    "TqdmRenderer",
]


class ProgressRenderer(Protocol):
    """Minimal renderer protocol — implement ``handle`` and ``close``."""

    def handle(self, event: ProgressEvent) -> None: ...

    def close(self) -> None: ...


# ---------------------------------------------------------------------------
# PlainRenderer — HPC-safe one-line logs
# ---------------------------------------------------------------------------


def _fmt_kv(payload: dict) -> str:
    """Render ``{k: v}`` as ``k=v k=v`` with stable ordering.

    String values that contain spaces / equals signs / newlines are
    quoted; durations on the ``elapsed_seconds`` key are formatted
    through :func:`format_duration` for readability.
    """
    preferred_order = (
        "stage", "sample_id", "step", "status",
        "elapsed_seconds", "reads_in", "reads_out", "aligned",
        "code", "severity", "output_type", "schema_version",
        "path", "reason", "message", "suggested_action", "note",
    )
    keys = [k for k in preferred_order if k in payload]
    keys += [k for k in payload if k not in preferred_order
             and k not in {"event", "timestamp"}]

    parts: list[str] = []
    for key in keys:
        value = payload[key]
        if value is None:
            continue
        if key == "elapsed_seconds":
            value = format_duration(float(value))
        text = str(value)
        if any(c in text for c in (" ", "\t", "=", "\n")):
            text = '"' + text.replace('"', '\\"').replace("\n", " ") + '"'
        parts.append(f"{key}={text}")
    return " ".join(parts)


class PlainRenderer:
    """Emit one ``[PROGRESS] event=... key=val ...`` line per event.

    Designed for HPC log capture: stable, append-only, no in-place
    rewrites. Every event becomes one line; downstream `grep` works.
    """

    def __init__(self, stream: TextIO | None = None) -> None:
        self._stream = stream if stream is not None else sys.stderr

    def handle(self, event: ProgressEvent) -> None:
        payload = event.to_dict()
        line = f"[PROGRESS] event={payload.get('event', '?')} "
        line += _fmt_kv(payload)
        with suppress(Exception):
            print(line.rstrip(), file=self._stream, flush=True)

    def close(self) -> None:
        with suppress(Exception):
            self._stream.flush()


# ---------------------------------------------------------------------------
# TqdmRenderer — interactive progress bars
# ---------------------------------------------------------------------------


class TqdmRenderer:
    """Per-stage tqdm sample counter, with PlainRenderer fallback.

    Maintains one bar per active stage. Bars are advanced on
    ``SampleEnd`` events and closed on ``StageEnd``. Non-bar events
    (warnings, outputs, resume-skips) are forwarded to the wrapped
    ``PlainRenderer`` so users keep the structured log trail in
    addition to the bars.
    """

    def __init__(self, stream: TextIO | None = None) -> None:
        self._stream = stream if stream is not None else sys.stderr
        self._plain = PlainRenderer(self._stream)
        self._bars: dict[str, object] = {}
        self._enabled = HAS_TQDM and self._is_tty(self._stream)

    @staticmethod
    def _is_tty(stream: TextIO) -> bool:
        return bool(getattr(stream, "isatty", lambda: False)())

    def handle(self, event: ProgressEvent) -> None:
        if not self._enabled:
            self._plain.handle(event)
            return

        if isinstance(event, StageStart) and event.n_samples:
            self._open_bar(event.stage, total=event.n_samples)
            self._plain.handle(event)
            return
        if isinstance(event, SampleEnd):
            bar = self._bars.get(event.stage)
            if bar is not None:
                with suppress(Exception):
                    bar.set_postfix_str(event.sample_id, refresh=False)
                    bar.update(1)
            else:
                self._plain.handle(event)
            return
        if isinstance(event, StageEnd):
            bar = self._bars.pop(event.stage, None)
            if bar is not None:
                with suppress(Exception):
                    bar.close()
            self._plain.handle(event)
            return

        # Everything else (warnings, outputs, sample-step events, ...)
        # falls through to the structured plain log so reviewers still
        # see the trail when a bar is on screen.
        self._plain.handle(event)

    def close(self) -> None:
        for bar in list(self._bars.values()):
            with suppress(Exception):
                bar.close()
        self._bars.clear()
        self._plain.close()

    def _open_bar(self, stage: str, *, total: int) -> None:
        try:
            from tqdm.auto import tqdm  # type: ignore[import-not-found]
        except ImportError:
            self._enabled = False
            return
        bar = tqdm(
            total=total,
            desc=f"{stage} samples",
            unit="sample",
            dynamic_ncols=True,
            leave=True,
        )
        self._bars[stage] = bar


# ---------------------------------------------------------------------------
# JsonlRenderer — machine-readable progress sidecar
# ---------------------------------------------------------------------------


class JsonlRenderer:
    """Append every event to a JSONL stream (one JSON object per line).

    Accepts either a file path (the file is created / appended to) or
    an open text-mode handle. Failure to write is swallowed so progress
    rendering never breaks a real run.
    """

    def __init__(self, target: Path | str | TextIO) -> None:
        self._owns_handle = False
        if isinstance(target, (str, Path)):
            path = Path(target)
            path.parent.mkdir(parents=True, exist_ok=True)
            self._handle: TextIO = path.open("a", encoding="utf-8")
            self._owns_handle = True
        else:
            self._handle = target  # type: ignore[assignment]

    def handle(self, event: ProgressEvent) -> None:
        line = json.dumps(event.to_dict(), sort_keys=True)
        with suppress(Exception):
            self._handle.write(line + "\n")
            self._handle.flush()

    def close(self) -> None:
        if self._owns_handle:
            with suppress(Exception):
                self._handle.close()


# Internal: keep a sentinel renderer that drops every event, used when
# the user passes ``--no-progress`` / ``mode=off``. The manager still
# emits the events into the in-memory stream so downstream code that
# reads them via :py:meth:`ProgressManager.events` keeps working.


class _NullRenderer:
    def handle(self, event: ProgressEvent) -> None:  # noqa: D401
        return None

    def close(self) -> None:
        return None
