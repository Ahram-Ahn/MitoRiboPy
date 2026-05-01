"""Tests for :mod:`mitoribopy.pipeline.resource_plan`.

The plan is a pure function over (n_samples, threads, parallel,
memory). It is the new default scheduler for ``mitoribopy align``;
these tests pin its behaviour so a future refactor cannot silently
change the worker count a publication run gets.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from mitoribopy.pipeline.resource_plan import (
    ResourcePlan,
    plan_parallelism,
    write_resource_plan,
)


def test_auto_default_uses_all_threads_when_enough_samples() -> None:
    plan = plan_parallelism(
        n_samples=8,
        requested_threads=16,
        requested_parallel="auto",
        memory_gb=None,
    )
    # min(8 samples, 16/2) = 8 parallel; 16/8 = 2 threads/worker.
    assert plan.parallel_samples == 8
    assert plan.per_sample_threads == 2
    assert plan.single_sample_mode is False
    assert plan.reason.startswith("auto")


def test_auto_caps_to_n_samples_when_threads_abundant() -> None:
    plan = plan_parallelism(
        n_samples=2,
        requested_threads=64,
        requested_parallel="auto",
    )
    assert plan.parallel_samples == 2
    assert plan.per_sample_threads == 32


def test_single_sample_mode_forces_serial() -> None:
    plan = plan_parallelism(
        n_samples=10,
        requested_threads=32,
        requested_parallel=1,
    )
    assert plan.parallel_samples == 1
    assert plan.per_sample_threads == 32
    assert plan.single_sample_mode is True


def test_parallel_string_alias_single() -> None:
    plan = plan_parallelism(
        n_samples=4,
        requested_threads=8,
        requested_parallel="single",
    )
    assert plan.parallel_samples == 1
    assert plan.single_sample_mode is True


def test_explicit_int_overrides_auto() -> None:
    plan = plan_parallelism(
        n_samples=10,
        requested_threads=32,
        requested_parallel=4,
    )
    assert plan.parallel_samples == 4
    assert plan.per_sample_threads == 8


def test_memory_constraint_caps_parallelism() -> None:
    plan = plan_parallelism(
        n_samples=20,
        requested_threads=32,
        requested_parallel="auto",
        memory_gb=8,  # 8 GiB total / 4 GiB per-sample = 2 workers
    )
    assert plan.parallel_samples == 2
    assert plan.memory_constraint_active is True


def test_threads_floor_at_one_when_t_lt_n() -> None:
    plan = plan_parallelism(
        n_samples=4,
        requested_threads=2,
        requested_parallel=4,
    )
    # T=2 / N=4 = 0 -> floor to 1
    assert plan.per_sample_threads == 1


def test_write_resource_plan_round_trips(tmp_path: Path) -> None:
    plan = plan_parallelism(
        n_samples=3, requested_threads=6, requested_parallel="auto",
    )
    out = write_resource_plan(plan, tmp_path)
    assert out == tmp_path / "resource_plan.json"
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["parallel_samples"] == plan.parallel_samples
    assert payload["per_sample_threads"] == plan.per_sample_threads
    assert payload["reason"] == plan.reason


def test_zero_or_negative_n_samples_clamps_to_one() -> None:
    plan = plan_parallelism(
        n_samples=0,
        requested_threads=4,
        requested_parallel="auto",
    )
    assert plan.parallel_samples >= 1
    assert plan.per_sample_threads >= 1
