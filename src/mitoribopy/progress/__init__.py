"""Per-stage timing + event-driven progress for MitoRiboPy.

Two coexisting layers:

1. **Timing primitives** (:mod:`.timing`) — wall-clock ``Stopwatch`` /
   ``StageTimings`` accumulators and the optional tqdm-backed
   ``SampleCounter``. Re-exported at the package root for back-compat
   with call sites that imported them directly.
2. **Event-driven progress** (:mod:`.events`, :mod:`.renderers`,
   :mod:`.manager`) — a small typed event vocabulary
   (``RunStart``, ``StageStart``, ``SampleStepEnd``, ``WarningEvent``,
   ``OutputEvent``, ``ResumeSkipEvent``, ...) emitted through a single
   ``ProgressManager`` and rendered by one or more pluggable renderers
   (plain log lines, tqdm bars, JSONL stream).

The split keeps existing call sites working while letting the
orchestrator publish a richer set of progress events that map cleanly
to HPC-safe one-line logs and a machine-readable ``progress.jsonl``
sidecar.
"""

from __future__ import annotations

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
from .manager import ProgressManager, ProgressMode, resolve_mode
from .renderers import (
    JsonlRenderer,
    PlainRenderer,
    ProgressRenderer,
    TqdmRenderer,
)
from .timing import (
    HAS_TQDM,
    SampleCounter,
    StageTimings,
    Stopwatch,
    format_duration,
    progress,
    render_step_timeline,
    render_summary_lines,
    stage_timer,
)

__all__ = [
    # Timing primitives (back-compat surface).
    "HAS_TQDM",
    "SampleCounter",
    "StageTimings",
    "Stopwatch",
    "format_duration",
    "progress",
    "render_step_timeline",
    "render_summary_lines",
    "stage_timer",
    # Event-driven progress.
    "ErrorEvent",
    "JsonlRenderer",
    "OutputEvent",
    "PlainRenderer",
    "ProgressEvent",
    "ProgressManager",
    "ProgressMode",
    "ProgressRenderer",
    "ResumeSkipEvent",
    "RunEnd",
    "RunStart",
    "SampleEnd",
    "SampleStart",
    "SampleStepEnd",
    "SampleStepStart",
    "StageEnd",
    "StageStart",
    "TqdmRenderer",
    "WarningEvent",
    "resolve_mode",
]
