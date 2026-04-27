"""Per-stage timing and compact CLI progress helpers.

Two roles:

1. Time the per-sample stages (cutadapt, bowtie2, dedup, ...) so the
   user can see how long each step took at a glance, both on the
   console and in the per-run ``mitoribopy.log``.
2. Provide an optional thin tqdm wrapper for compact "N of M done"
   counters, with a graceful no-op fallback when tqdm is not
   installed (the package's hard deps do not include tqdm).

Usage shape mirrors the lightweight ``progress.py`` style: a
``Stopwatch`` context manager, a ``StageTimings`` accumulator, and a
``progress(...)`` iterable wrapper.
"""

from __future__ import annotations

import threading
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Iterable, Iterator, TypeVar

T = TypeVar("T")

try:  # pragma: no cover - tqdm is optional
    from tqdm.auto import tqdm as _tqdm

    HAS_TQDM = True
except ImportError:  # pragma: no cover
    HAS_TQDM = False
    _tqdm = None


# ---------------------------------------------------------------------------
# Duration formatting
# ---------------------------------------------------------------------------


def format_duration(seconds: float) -> str:
    """Render a wall-clock duration in a compact, human-friendly form.

    Examples::

        format_duration(0.0034)  -> '<1ms'
        format_duration(0.42)    -> '420ms'
        format_duration(12.345)  -> '12.3s'
        format_duration(83.0)    -> '1m 23s'
        format_duration(3725.0)  -> '1h 02m 05s'

    The format aims for a fixed-width-friendly representation so that
    columns of timings line up when printed in a table.
    """
    if seconds < 0:
        seconds = 0.0
    if seconds < 0.001:
        return "<1ms"
    if seconds < 1.0:
        return f"{int(round(seconds * 1000))}ms"
    if seconds < 60.0:
        return f"{seconds:.1f}s"
    if seconds < 3600.0:
        minutes = int(seconds // 60)
        rem = int(round(seconds - minutes * 60))
        if rem == 60:  # rounding overflow
            minutes += 1
            rem = 0
        return f"{minutes}m {rem:02d}s"
    hours = int(seconds // 3600)
    rem_seconds = seconds - hours * 3600
    minutes = int(rem_seconds // 60)
    secs = int(round(rem_seconds - minutes * 60))
    if secs == 60:
        minutes += 1
        secs = 0
    if minutes == 60:
        hours += 1
        minutes = 0
    return f"{hours}h {minutes:02d}m {secs:02d}s"


# ---------------------------------------------------------------------------
# Stopwatch
# ---------------------------------------------------------------------------


@dataclass
class Stopwatch:
    """Context-manager-friendly wall-clock timer.

    The elapsed duration in seconds is exposed via ``.seconds`` after
    the ``with`` block exits, and via ``.elapsed`` for the live value
    while still inside the block. Uses ``time.perf_counter`` so the
    measurements are unaffected by wall-clock adjustments.

    Example::

        with Stopwatch() as sw:
            do_work()
        print(f"took {sw.seconds:.2f}s")
    """

    _start: float | None = None
    _end: float | None = None

    def __enter__(self) -> "Stopwatch":
        self._start = time.perf_counter()
        self._end = None
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self._end = time.perf_counter()

    @property
    def elapsed(self) -> float:
        """Live elapsed seconds (works inside or after the with-block)."""
        if self._start is None:
            return 0.0
        end = self._end if self._end is not None else time.perf_counter()
        return end - self._start

    @property
    def seconds(self) -> float:
        """Final elapsed seconds; valid only after ``__exit__``."""
        if self._start is None or self._end is None:
            return 0.0
        return self._end - self._start


# ---------------------------------------------------------------------------
# Per-sample, per-stage accumulator
# ---------------------------------------------------------------------------


@dataclass
class StageTimings:
    """Thread-safe per-sample, per-stage timing accumulator.

    Designed for the parallel-sample align flow where multiple worker
    threads record stage durations concurrently. Records are keyed by
    ``(sample, stage)``; a stable ``stage_order`` tracks first-seen
    insertion order so summary tables match the actual pipeline order
    instead of dictionary iteration order.

    Use ``record(sample, stage, seconds)`` to log one stage; use
    ``aggregate()`` to build a per-stage summary (total / mean / max)
    over all samples for the end-of-run table.
    """

    _per_sample: dict[str, dict[str, float]] = field(default_factory=dict)
    _stage_order: list[str] = field(default_factory=list)
    _lock: threading.Lock = field(default_factory=threading.Lock)

    def record(self, sample: str, stage: str, seconds: float) -> None:
        with self._lock:
            self._per_sample.setdefault(sample, {})[stage] = float(seconds)
            if stage not in self._stage_order:
                self._stage_order.append(stage)

    def total_for(self, sample: str) -> float:
        with self._lock:
            return float(sum(self._per_sample.get(sample, {}).values()))

    def stage_order(self) -> list[str]:
        with self._lock:
            return list(self._stage_order)

    def per_sample_view(self) -> dict[str, dict[str, float]]:
        with self._lock:
            return {
                sample: dict(stages) for sample, stages in self._per_sample.items()
            }

    def aggregate(self) -> list[tuple[str, float, float, float, int]]:
        """Per-stage summary: list of ``(stage, total, mean, max, n)`` tuples."""
        with self._lock:
            stages = list(self._stage_order)
            samples = list(self._per_sample.keys())
            rows: list[tuple[str, float, float, float, int]] = []
            for stage in stages:
                values = [
                    self._per_sample[s][stage]
                    for s in samples
                    if stage in self._per_sample[s]
                ]
                n = len(values)
                if n == 0:
                    continue
                total = sum(values)
                mean = total / n
                rows.append((stage, total, mean, max(values), n))
            return rows


# ---------------------------------------------------------------------------
# Summary rendering
# ---------------------------------------------------------------------------


def render_summary_lines(
    timings: StageTimings,
    *,
    wall_seconds: float | None = None,
    samples_total: int | None = None,
) -> list[str]:
    """Render a compact per-stage summary as a list of log-ready lines.

    Each returned string is one log line (no trailing newlines, no
    component prefix). The caller decides which logger / component to
    push them through; see :func:`mitoribopy.console.log_info`.

    The table is padded so columns line up in a fixed-width font.
    """
    rows = timings.aggregate()
    samples_present = len(timings.per_sample_view())
    sample_count = samples_total if samples_total is not None else samples_present

    if not rows:
        return [f"Timing summary: no stages recorded ({sample_count} sample(s))."]

    stage_w = max(len("stage"), max(len(stage) for stage, *_ in rows))
    headers = ("stage", "total", "mean", "max")
    col_widths = [stage_w, 10, 10, 10]

    def fmt_row(values: tuple[str, str, str, str]) -> str:
        return "  ".join(
            value.ljust(col_widths[i]) if i == 0 else value.rjust(col_widths[i])
            for i, value in enumerate(values)
        )

    out: list[str] = []
    label = (
        f"Timing summary ({sample_count} sample(s), {len(rows)} stage(s)):"
    )
    out.append(label)
    out.append(fmt_row(headers))
    for stage, total, mean, mx, _n in rows:
        out.append(
            fmt_row(
                (
                    stage,
                    format_duration(total),
                    format_duration(mean),
                    format_duration(mx),
                )
            )
        )
    if wall_seconds is not None:
        out.append(f"wall: {format_duration(wall_seconds)}")
    return out


# ---------------------------------------------------------------------------
# Optional tqdm progress bar
# ---------------------------------------------------------------------------


def progress(
    it: Iterable[T],
    *,
    desc: str | None = None,
    total: int | None = None,
    leave: bool = True,
    unit: str = "it",
    disable: bool = False,
) -> Iterator[T]:
    """Wrap an iterable with a tqdm progress bar.

    Falls back to plain iteration when tqdm is unavailable or when
    ``disable=True``. Designed for serial loops; for parallel
    ThreadPoolExecutor-driven flows, prefer :class:`SampleCounter`.
    """
    if disable or not HAS_TQDM:
        return iter(it)
    return _tqdm(  # type: ignore[no-any-return]
        it,
        desc=desc,
        total=total,
        leave=leave,
        unit=unit,
        dynamic_ncols=True,
    )


class SampleCounter:
    """Compact "samples done" counter for parallel-sample flows.

    A thin, thread-safe wrapper around tqdm that updates a single
    line as workers finish samples. When tqdm is unavailable or
    ``disable=True``, the calls become no-ops; the per-sample log
    lines emitted by the caller already convey the same information.

    Use as a context manager so the bar closes cleanly on success or
    exception::

        with SampleCounter(total=N, desc="ALIGN samples") as bar:
            for sample in completed_samples():
                bar.advance(sample)
    """

    def __init__(
        self,
        *,
        total: int,
        desc: str = "samples",
        disable: bool = False,
    ) -> None:
        self._total = max(0, int(total))
        self._desc = desc
        self._disabled = disable or not HAS_TQDM or self._total == 0
        self._bar = None
        self._lock = threading.Lock()

    def __enter__(self) -> "SampleCounter":
        if not self._disabled:
            self._bar = _tqdm(  # type: ignore[misc]
                total=self._total,
                desc=self._desc,
                unit="sample",
                dynamic_ncols=True,
                leave=True,
            )
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if self._bar is not None:
            self._bar.close()
            self._bar = None

    def advance(self, sample: str | None = None) -> None:
        with self._lock:
            if self._bar is None:
                return
            if sample:
                self._bar.set_postfix_str(sample, refresh=False)
            self._bar.update(1)


@contextmanager
def stage_timer(
    timings: StageTimings | None,
    sample: str,
    stage: str,
) -> Iterator[Stopwatch]:
    """Time a stage and (optionally) record into a ``StageTimings``.

    Example::

        with stage_timer(timings, sample, "trim") as sw:
            run_cutadapt(...)
        log(f"trim {format_duration(sw.seconds)}")
    """
    sw = Stopwatch()
    sw.__enter__()
    try:
        yield sw
    finally:
        sw.__exit__(None, None, None)
        if timings is not None:
            timings.record(sample, stage, sw.seconds)
