"""Typed event vocabulary for the MitoRiboPy progress manager.

All events are frozen dataclasses with a stable ``event`` discriminator
string and a ``to_dict()`` method that returns a plain JSON-serialisable
dict (no datetime objects, no enums) suitable for streaming to
``progress.jsonl``. The canonical key set is:

* ``event``           — discriminator (``run_start``, ``stage_end``, ...)
* ``timestamp``       — ISO-8601 UTC, second resolution
* ``stage`` / ``sample_id`` / ``step`` — present where applicable
* event-specific keys (``status``, ``elapsed_seconds``, ``code``, ...)

Adding a new event:

1. Subclass :class:`ProgressEvent` (frozen dataclass) and set ``EVENT``.
2. Override ``to_dict()`` only if you need extra keys beyond the
   declared dataclass fields.
3. Register the new event in the renderer of interest if it should
   surface in plain / tqdm output (the JSONL renderer captures every
   event automatically).
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields
from datetime import datetime, timezone


__all__ = [
    "ErrorEvent",
    "OutputEvent",
    "ProgressEvent",
    "ResumeSkipEvent",
    "RunEnd",
    "RunStart",
    "SampleEnd",
    "SampleStart",
    "SampleStepEnd",
    "SampleStepStart",
    "StageEnd",
    "StageStart",
    "WarningEvent",
]


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@dataclass(frozen=True)
class ProgressEvent:
    """Base class for every event surfaced through ``ProgressManager``.

    Concrete subclasses set the ``EVENT`` class attribute to a stable
    discriminator string and add typed fields. Use :func:`to_dict` to
    serialise; do not subclass ``__post_init__`` for I/O — the manager
    handles dispatch.
    """

    EVENT: str = field(default="", init=False, repr=False)
    timestamp: str = field(default_factory=_utc_now_iso)

    def to_dict(self) -> dict:
        out: dict = {"event": self.EVENT, "timestamp": self.timestamp}
        for f in fields(self):
            if f.name == "timestamp":
                continue
            value = getattr(self, f.name)
            if value is None:
                continue
            out[f.name] = value
        return out


@dataclass(frozen=True)
class RunStart(ProgressEvent):
    """Emitted once at the very start of a ``mitoribopy all`` invocation."""

    EVENT: str = field(default="run_start", init=False, repr=False)
    subcommand: str = "all"
    n_samples: int | None = None
    config_source: str | None = None


@dataclass(frozen=True)
class RunEnd(ProgressEvent):
    """Emitted once at the end of a run, regardless of success / failure."""

    EVENT: str = field(default="run_end", init=False, repr=False)
    status: str = "done"  # done | error
    elapsed_seconds: float = 0.0
    n_warnings: int = 0
    n_errors: int = 0


@dataclass(frozen=True)
class StageStart(ProgressEvent):
    EVENT: str = field(default="stage_start", init=False, repr=False)
    stage: str = ""  # align | rpf | rnaseq | summarize | validate-figures
    n_samples: int | None = None


@dataclass(frozen=True)
class StageEnd(ProgressEvent):
    EVENT: str = field(default="stage_end", init=False, repr=False)
    stage: str = ""
    status: str = "done"  # done | skipped | error
    elapsed_seconds: float = 0.0
    reason: str | None = None  # populated when status=skipped|error


@dataclass(frozen=True)
class SampleStart(ProgressEvent):
    EVENT: str = field(default="sample_start", init=False, repr=False)
    stage: str = ""
    sample_id: str = ""


@dataclass(frozen=True)
class SampleEnd(ProgressEvent):
    EVENT: str = field(default="sample_end", init=False, repr=False)
    stage: str = ""
    sample_id: str = ""
    status: str = "done"  # done | error
    elapsed_seconds: float = 0.0


@dataclass(frozen=True)
class SampleStepStart(ProgressEvent):
    EVENT: str = field(default="sample_step_start", init=False, repr=False)
    stage: str = ""
    sample_id: str = ""
    step: str = ""  # trim | mt_align | dedup | offset_selection | ...


@dataclass(frozen=True)
class SampleStepEnd(ProgressEvent):
    EVENT: str = field(default="sample_step_end", init=False, repr=False)
    stage: str = ""
    sample_id: str = ""
    step: str = ""
    status: str = "done"  # done | error
    elapsed_seconds: float = 0.0
    # Optional payload for steps that report read counts / yields.
    reads_in: int | None = None
    reads_out: int | None = None
    aligned: int | None = None
    note: str | None = None


@dataclass(frozen=True)
class WarningEvent(ProgressEvent):
    """Mirror of one row in ``warnings.tsv`` so it streams in real time."""

    EVENT: str = field(default="warning", init=False, repr=False)
    stage: str = ""
    sample_id: str | None = None
    severity: str = "warn"  # warn | error | info
    code: str | None = None
    message: str = ""
    suggested_action: str | None = None


@dataclass(frozen=True)
class ErrorEvent(ProgressEvent):
    """A hard error that aborted a stage or sample step."""

    EVENT: str = field(default="error", init=False, repr=False)
    stage: str = ""
    sample_id: str | None = None
    code: str | None = None
    message: str = ""


@dataclass(frozen=True)
class OutputEvent(ProgressEvent):
    """A meaningful artefact (TSV, plot, sidecar) was just written."""

    EVENT: str = field(default="output", init=False, repr=False)
    stage: str = ""
    output_type: str = ""  # rpf_counts | te_table | offset_diagnostics | ...
    path: str = ""
    schema_version: str | None = None


@dataclass(frozen=True)
class ResumeSkipEvent(ProgressEvent):
    """A stage was skipped because ``--resume`` matched a prior good run."""

    EVENT: str = field(default="resume_skip", init=False, repr=False)
    stage: str = ""
    reason: str = ""  # human-readable why-it-was-skipped
