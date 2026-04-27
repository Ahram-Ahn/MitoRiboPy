"""Unit tests for ``mitoribopy.progress``."""

from __future__ import annotations

import threading
import time

import pytest

from mitoribopy.progress import (
    StageTimings,
    Stopwatch,
    format_duration,
    render_summary_lines,
    stage_timer,
)


# ---------------------------------------------------------------------------
# format_duration
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "seconds, expected",
    [
        (-1.0, "<1ms"),
        (0.0, "<1ms"),
        (0.0005, "<1ms"),
        (0.05, "50ms"),
        (0.5, "500ms"),
        (1.0, "1.0s"),
        (12.345, "12.3s"),
        (59.4, "59.4s"),  # under 60 -> seconds with one decimal
        (60.0, "1m 00s"),
        (83.0, "1m 23s"),
        (3599.0, "59m 59s"),  # below 1h boundary, still mm format
        (3600.0, "1h 00m 00s"),
        (3725.0, "1h 02m 05s"),
        (7325.4, "2h 02m 05s"),
    ],
)
def test_format_duration_buckets(seconds, expected) -> None:
    assert format_duration(seconds) == expected


# ---------------------------------------------------------------------------
# Stopwatch
# ---------------------------------------------------------------------------


def test_stopwatch_records_positive_elapsed() -> None:
    with Stopwatch() as sw:
        time.sleep(0.01)
    assert sw.seconds > 0
    assert sw.seconds < 1.0  # sanity: not measuring whole seconds


def test_stopwatch_elapsed_is_zero_before_enter() -> None:
    sw = Stopwatch()
    assert sw.seconds == 0.0
    assert sw.elapsed == 0.0


def test_stopwatch_elapsed_is_live_inside_block() -> None:
    sw = Stopwatch()
    sw.__enter__()
    time.sleep(0.005)
    mid = sw.elapsed
    assert mid > 0
    sw.__exit__(None, None, None)
    assert sw.seconds >= mid


# ---------------------------------------------------------------------------
# StageTimings
# ---------------------------------------------------------------------------


def test_stage_timings_record_and_view() -> None:
    t = StageTimings()
    t.record("s1", "trim", 1.0)
    t.record("s1", "mt-align", 5.0)
    t.record("s2", "trim", 2.0)
    view = t.per_sample_view()
    assert view == {"s1": {"trim": 1.0, "mt-align": 5.0}, "s2": {"trim": 2.0}}


def test_stage_timings_preserves_first_seen_order() -> None:
    t = StageTimings()
    t.record("s1", "trim", 1.0)
    t.record("s1", "mt-align", 2.0)
    t.record("s2", "dedup", 3.0)
    # 'dedup' was added AFTER trim/mt-align via s2; should appear last.
    assert t.stage_order() == ["trim", "mt-align", "dedup"]


def test_stage_timings_aggregate_computes_total_mean_max() -> None:
    t = StageTimings()
    t.record("s1", "trim", 10.0)
    t.record("s2", "trim", 30.0)
    t.record("s1", "mt-align", 100.0)
    rows = t.aggregate()
    by_stage = {r[0]: r for r in rows}
    assert by_stage["trim"] == ("trim", 40.0, 20.0, 30.0, 2)
    # mt-align only seen on s1
    assert by_stage["mt-align"] == ("mt-align", 100.0, 100.0, 100.0, 1)


def test_stage_timings_total_for_sample() -> None:
    t = StageTimings()
    t.record("s1", "a", 1.5)
    t.record("s1", "b", 2.5)
    t.record("s2", "a", 7.0)
    assert t.total_for("s1") == pytest.approx(4.0)
    assert t.total_for("s2") == pytest.approx(7.0)
    assert t.total_for("missing") == 0.0


def test_stage_timings_is_thread_safe() -> None:
    """Multiple writers should not corrupt the per-sample dict."""
    t = StageTimings()
    n_threads = 8
    n_per_thread = 200

    def worker(idx: int) -> None:
        for j in range(n_per_thread):
            t.record(f"s{idx}", f"stage{j % 3}", float(j))

    threads = [threading.Thread(target=worker, args=(i,)) for i in range(n_threads)]
    for th in threads:
        th.start()
    for th in threads:
        th.join()

    view = t.per_sample_view()
    assert set(view.keys()) == {f"s{i}" for i in range(n_threads)}
    # Each sample saw stage0..stage2 with the LAST value written for that
    # stage (overwrite semantics; not an accumulation).
    for sample, stages in view.items():
        assert set(stages.keys()) == {"stage0", "stage1", "stage2"}


# ---------------------------------------------------------------------------
# stage_timer context manager
# ---------------------------------------------------------------------------


def test_stage_timer_records_into_timings() -> None:
    t = StageTimings()
    with stage_timer(t, "sX", "trim") as sw:
        time.sleep(0.005)
    view = t.per_sample_view()
    assert "sX" in view and "trim" in view["sX"]
    assert view["sX"]["trim"] > 0
    assert sw.seconds == view["sX"]["trim"]


def test_stage_timer_with_none_timings_still_times() -> None:
    """Passing timings=None must not crash; the Stopwatch still works."""
    with stage_timer(None, "sX", "trim") as sw:
        time.sleep(0.005)
    assert sw.seconds > 0


def test_stage_timer_records_even_on_exception() -> None:
    t = StageTimings()
    with pytest.raises(RuntimeError):
        with stage_timer(t, "sX", "trim") as sw:
            time.sleep(0.005)
            raise RuntimeError("boom")
    # The exception must not skip the recording.
    assert sw.seconds > 0
    assert t.per_sample_view()["sX"]["trim"] == sw.seconds


# ---------------------------------------------------------------------------
# render_summary_lines
# ---------------------------------------------------------------------------


def test_render_summary_lines_includes_header_and_rows() -> None:
    t = StageTimings()
    t.record("s1", "trim", 12.3)
    t.record("s2", "trim", 14.7)
    t.record("s1", "mt-align", 45.0)
    lines = render_summary_lines(t, wall_seconds=60.0)
    # First line is the title.
    assert "Timing summary" in lines[0]
    assert "2 sample(s)" in lines[0]
    # The header row has the column names.
    assert "stage" in lines[1] and "total" in lines[1] and "max" in lines[1]
    # One row per stage.
    body = "\n".join(lines)
    assert "trim" in body
    assert "mt-align" in body
    # Wall-time line at the end.
    assert lines[-1].startswith("wall: ")


def test_render_summary_lines_handles_empty() -> None:
    t = StageTimings()
    lines = render_summary_lines(t, samples_total=0)
    assert any("no stages recorded" in line for line in lines)
