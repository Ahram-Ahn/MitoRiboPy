"""Central dispatcher for typed progress events.

A ``ProgressManager`` is created once per ``mitoribopy all`` invocation
(or once per standalone subcommand). Subcommands push events through a
small set of helper methods (``run_start``, ``stage_start``,
``sample_step_end``, ``warning``, ``output``, ...); the manager then
fans them out to one or more :class:`ProgressRenderer` instances.

The manager is intentionally cheap: when no renderer is attached every
helper is essentially a structured no-op, so existing call sites can
adopt it incrementally without paying for output they don't render.

Mode resolution
---------------

The ``--progress`` CLI flag accepts six values; :func:`resolve_mode`
maps any of them to a stable :class:`ProgressMode` enum after applying
the auto/no-progress logic:

* ``auto``  → ``bar`` on TTY, ``plain`` on non-TTY
* ``plain`` → :class:`PlainRenderer`
* ``bar``   → :class:`TqdmRenderer` (falls back to plain if tqdm missing)
* ``rich``  → :class:`TqdmRenderer` (rich support is a future addition)
* ``jsonl`` → :class:`JsonlRenderer` ONLY (machine-readable, no console)
* ``off``   → no console renderer (``progress.jsonl`` still attaches if
  ``--progress-file`` was set)
"""

from __future__ import annotations

import enum
import sys
from contextlib import suppress
from pathlib import Path
from typing import Iterable, TextIO

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
from .renderers import (
    JsonlRenderer,
    PlainRenderer,
    ProgressRenderer,
    TqdmRenderer,
    _NullRenderer,
)


__all__ = [
    "ProgressManager",
    "ProgressMode",
    "resolve_mode",
]


class ProgressMode(str, enum.Enum):
    AUTO = "auto"
    PLAIN = "plain"
    BAR = "bar"
    RICH = "rich"
    JSONL = "jsonl"
    OFF = "off"


_VALID_MODES: tuple[str, ...] = tuple(m.value for m in ProgressMode)


def _is_tty(stream: TextIO) -> bool:
    return bool(getattr(stream, "isatty", lambda: False)())


def resolve_mode(
    requested: str | None,
    *,
    stream: TextIO | None = None,
) -> ProgressMode:
    """Resolve a user-facing mode string to a concrete :class:`ProgressMode`.

    ``None`` is treated as ``auto``. Unknown values fall back to ``plain``
    so a typo never silently disables progress.
    """
    if requested is None or requested == ProgressMode.AUTO.value:
        target = stream if stream is not None else sys.stderr
        return ProgressMode.BAR if _is_tty(target) else ProgressMode.PLAIN
    if requested in _VALID_MODES:
        return ProgressMode(requested)
    return ProgressMode.PLAIN


class ProgressManager:
    """Central dispatcher for ``ProgressEvent`` instances.

    Construct via :meth:`from_cli` (preferred) or directly with an
    explicit list of renderers. Every emit helper builds the matching
    event, appends it to the in-memory stream, and forwards it to each
    attached renderer.

    Always close the manager via :meth:`close` (or use it as a context
    manager) so renderer file handles flush.
    """

    def __init__(
        self,
        renderers: Iterable[ProgressRenderer] | None = None,
        *,
        keep_history: bool = True,
    ) -> None:
        self._renderers: list[ProgressRenderer] = list(renderers or [])
        self._events: list[ProgressEvent] = [] if keep_history else []
        self._keep_history = keep_history
        self._closed = False

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    @classmethod
    def from_cli(
        cls,
        *,
        mode: str | None = None,
        progress_file: Path | str | None = None,
        stream: TextIO | None = None,
        keep_history: bool = True,
    ) -> "ProgressManager":
        """Build a manager from CLI flags.

        ``mode`` is one of the strings accepted by :func:`resolve_mode`.
        ``progress_file`` (the ``--progress-file`` flag) attaches a
        :class:`JsonlRenderer` regardless of ``mode`` so a user can
        run with ``--no-progress`` (mode=off) and still capture a
        machine-readable sidecar.
        """
        target_stream = stream if stream is not None else sys.stderr
        resolved = resolve_mode(mode, stream=target_stream)

        renderers: list[ProgressRenderer] = []
        if resolved == ProgressMode.OFF:
            renderers.append(_NullRenderer())
        elif resolved == ProgressMode.JSONL:
            # JSONL mode means: only the JSONL sidecar, no console.
            # If the user passed ``mode=jsonl`` without ``--progress-file``,
            # write to stderr as JSONL so a tee captures everything.
            if progress_file is not None:
                renderers.append(JsonlRenderer(progress_file))
            else:
                renderers.append(JsonlRenderer(target_stream))
        elif resolved == ProgressMode.PLAIN:
            renderers.append(PlainRenderer(target_stream))
        else:  # BAR, RICH (rich currently degrades to bar)
            renderers.append(TqdmRenderer(target_stream))

        # An explicit ``--progress-file`` always attaches a JSONL
        # renderer alongside the console one (unless it was already
        # added above for ``mode=jsonl``).
        if progress_file is not None and resolved != ProgressMode.JSONL:
            renderers.append(JsonlRenderer(progress_file))

        return cls(renderers, keep_history=keep_history)

    @classmethod
    def null(cls) -> "ProgressManager":
        """A manager with no renderers — useful for tests + library use."""
        return cls([], keep_history=True)

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def __enter__(self) -> "ProgressManager":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def close(self) -> None:
        if self._closed:
            return
        for renderer in self._renderers:
            with suppress(Exception):
                renderer.close()
        self._closed = True

    # ------------------------------------------------------------------
    # Event emission
    # ------------------------------------------------------------------

    def emit(self, event: ProgressEvent) -> None:
        """Forward *event* to every attached renderer."""
        if self._keep_history:
            self._events.append(event)
        for renderer in self._renderers:
            with suppress(Exception):
                renderer.handle(event)

    def events(self) -> list[ProgressEvent]:
        """Return a shallow copy of every event seen so far."""
        return list(self._events)

    # ---- ergonomic helpers (one per concrete event) ------------------

    def run_start(
        self,
        *,
        subcommand: str = "all",
        n_samples: int | None = None,
        config_source: str | None = None,
    ) -> None:
        self.emit(RunStart(
            subcommand=subcommand,
            n_samples=n_samples,
            config_source=config_source,
        ))

    def run_end(
        self,
        *,
        status: str = "done",
        elapsed_seconds: float = 0.0,
        n_warnings: int = 0,
        n_errors: int = 0,
    ) -> None:
        self.emit(RunEnd(
            status=status,
            elapsed_seconds=elapsed_seconds,
            n_warnings=n_warnings,
            n_errors=n_errors,
        ))

    def stage_start(self, stage: str, *, n_samples: int | None = None) -> None:
        self.emit(StageStart(stage=stage, n_samples=n_samples))

    def stage_end(
        self,
        stage: str,
        *,
        status: str = "done",
        elapsed_seconds: float = 0.0,
        reason: str | None = None,
    ) -> None:
        self.emit(StageEnd(
            stage=stage,
            status=status,
            elapsed_seconds=elapsed_seconds,
            reason=reason,
        ))

    def sample_start(self, stage: str, sample_id: str) -> None:
        self.emit(SampleStart(stage=stage, sample_id=sample_id))

    def sample_end(
        self,
        stage: str,
        sample_id: str,
        *,
        status: str = "done",
        elapsed_seconds: float = 0.0,
    ) -> None:
        self.emit(SampleEnd(
            stage=stage,
            sample_id=sample_id,
            status=status,
            elapsed_seconds=elapsed_seconds,
        ))

    def sample_step_start(self, stage: str, sample_id: str, step: str) -> None:
        self.emit(SampleStepStart(stage=stage, sample_id=sample_id, step=step))

    def sample_step_end(
        self,
        stage: str,
        sample_id: str,
        step: str,
        *,
        status: str = "done",
        elapsed_seconds: float = 0.0,
        reads_in: int | None = None,
        reads_out: int | None = None,
        aligned: int | None = None,
        note: str | None = None,
    ) -> None:
        self.emit(SampleStepEnd(
            stage=stage,
            sample_id=sample_id,
            step=step,
            status=status,
            elapsed_seconds=elapsed_seconds,
            reads_in=reads_in,
            reads_out=reads_out,
            aligned=aligned,
            note=note,
        ))

    def warning(
        self,
        *,
        stage: str,
        message: str,
        sample_id: str | None = None,
        code: str | None = None,
        severity: str = "warn",
        suggested_action: str | None = None,
    ) -> None:
        self.emit(WarningEvent(
            stage=stage,
            sample_id=sample_id,
            severity=severity,
            code=code,
            message=message,
            suggested_action=suggested_action,
        ))

    def error(
        self,
        *,
        stage: str,
        message: str,
        sample_id: str | None = None,
        code: str | None = None,
    ) -> None:
        self.emit(ErrorEvent(
            stage=stage,
            sample_id=sample_id,
            code=code,
            message=message,
        ))

    def output(
        self,
        *,
        stage: str,
        output_type: str,
        path: str,
        schema_version: str | None = None,
    ) -> None:
        self.emit(OutputEvent(
            stage=stage,
            output_type=output_type,
            path=path,
            schema_version=schema_version,
        ))

    def resume_skip(self, *, stage: str, reason: str) -> None:
        self.emit(ResumeSkipEvent(stage=stage, reason=reason))
