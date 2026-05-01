"""Resource planning for parallel sample execution.

Replaces ad-hoc ``--max-parallel-samples 1`` defaults with an explicit,
serialisable :class:`ResourcePlan`. The plan is written to disk
(``resource_plan.json``) so a reviewer can audit how many samples ran
concurrently and how many threads each tool worker received without
re-reading the run banner.

Design contract
---------------

* Default (``parallel_samples="auto"``) maximises throughput safely:
  ``min(n_samples, total_threads // min_threads_per_sample)`` and, when
  ``memory_gb`` is set, ``total_memory_gb // estimated_memory_per_sample_gb``.
* ``parallel_samples=1`` (or the alias ``"single"``) forces serial.
* ``threads="auto"`` resolves to ``os.cpu_count()`` at plan time.
* Per-sample thread budget is ``max(1, total_threads // parallel_samples)``,
  matching the existing :func:`_per_worker_threads` semantics in the
  align CLI.

Recorded fields are stable across releases — see ``to_dict()``. Every
new field must be additive (downstream summarisers must tolerate
missing keys).
"""

from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


__all__ = [
    "ResourcePlan",
    "plan_parallelism",
    "write_resource_plan",
]


# Reasonable defaults for an mt-Ribo-seq align worker. Per-sample memory
# is dominated by bowtie2 + a sorted BAM held in memory during dedup;
# 4 GiB is a safe overestimate for typical metazoan mt-transcriptomes.
_MIN_THREADS_PER_SAMPLE_DEFAULT = 2
_EST_MEMORY_PER_SAMPLE_GB_DEFAULT = 4.0


@dataclass(frozen=True)
class ResourcePlan:
    """Concrete plan for running ``n_samples`` samples in parallel."""

    n_samples: int
    total_threads: int
    total_memory_gb: float | None
    parallel_samples: int
    per_sample_threads: int
    min_threads_per_sample: int
    estimated_memory_per_sample_gb: float
    requested_threads: int | str
    requested_parallel: int | str
    memory_constraint_active: bool
    single_sample_mode: bool
    reason: str

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        # Stringify "auto" without losing type info for ints.
        d["requested_threads"] = (
            self.requested_threads
            if isinstance(self.requested_threads, str)
            else int(self.requested_threads)
        )
        d["requested_parallel"] = (
            self.requested_parallel
            if isinstance(self.requested_parallel, str)
            else int(self.requested_parallel)
        )
        return d


def _coerce_threads(requested: int | str | None) -> int:
    if requested is None or (
        isinstance(requested, str) and requested.lower() in {"auto", ""}
    ):
        return max(1, os.cpu_count() or 1)
    return max(1, int(requested))


def _coerce_memory(memory_gb: float | str | None) -> float | None:
    if memory_gb is None:
        return None
    if isinstance(memory_gb, str):
        s = memory_gb.strip().lower()
        if s in {"", "auto", "none"}:
            return None
        return float(s)
    return float(memory_gb)


def plan_parallelism(
    n_samples: int,
    *,
    requested_threads: int | str | None = "auto",
    requested_parallel: int | str | None = "auto",
    memory_gb: float | str | None = "auto",
    min_threads_per_sample: int = _MIN_THREADS_PER_SAMPLE_DEFAULT,
    estimated_memory_per_sample_gb: float = _EST_MEMORY_PER_SAMPLE_GB_DEFAULT,
) -> ResourcePlan:
    """Plan how many samples to run concurrently and at what thread budget.

    Implements the policy described in the module docstring. Returns a
    frozen :class:`ResourcePlan`; never raises on edge cases (empty
    sample list, T < N, memory unset) — the plan always converges to
    ``parallel_samples >= 1`` and ``per_sample_threads >= 1``.
    """
    n_samples = max(1, int(n_samples))
    total_threads = _coerce_threads(requested_threads)
    total_memory = _coerce_memory(memory_gb)

    rp_raw = requested_parallel
    rp_norm: str | int
    if rp_raw is None:
        rp_norm = "auto"
    elif isinstance(rp_raw, str):
        rp_norm = rp_raw.lower().strip() or "auto"
    else:
        rp_norm = int(rp_raw)

    single_sample = (
        rp_norm == 1 or (isinstance(rp_norm, str) and rp_norm in {"1", "single"})
    )

    memory_constraint_active = False
    reason: str
    if single_sample:
        parallel = 1
        reason = "user override: single-sample mode"
    elif rp_norm == "auto":
        by_cpu = max(1, total_threads // max(1, min_threads_per_sample))
        if total_memory is not None and estimated_memory_per_sample_gb > 0:
            by_mem = max(1, int(total_memory // estimated_memory_per_sample_gb))
            parallel = min(n_samples, by_cpu, by_mem)
            memory_constraint_active = parallel == by_mem and by_mem < min(
                n_samples, by_cpu
            )
            reason = (
                "auto: min(n_samples, threads/min_per_sample, "
                "memory_gb/est_per_sample)"
            )
        else:
            parallel = min(n_samples, by_cpu)
            reason = "auto: min(n_samples, threads/min_per_sample)"
    else:
        parallel = max(1, int(rp_norm))
        reason = "user override: explicit parallel_samples"

    parallel = max(1, min(parallel, n_samples))
    per_sample = max(1, total_threads // parallel)

    return ResourcePlan(
        n_samples=n_samples,
        total_threads=total_threads,
        total_memory_gb=total_memory,
        parallel_samples=parallel,
        per_sample_threads=per_sample,
        min_threads_per_sample=min_threads_per_sample,
        estimated_memory_per_sample_gb=estimated_memory_per_sample_gb,
        requested_threads=(
            requested_threads if requested_threads is not None else "auto"
        ),
        requested_parallel=(
            requested_parallel if requested_parallel is not None else "auto"
        ),
        memory_constraint_active=memory_constraint_active,
        single_sample_mode=parallel == 1,
        reason=reason,
    )


def write_resource_plan(plan: ResourcePlan, output_dir: str | Path) -> Path:
    """Persist a plan to ``<output_dir>/resource_plan.json``.

    Returns the written path so callers can log it.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "resource_plan.json"
    path.write_text(
        json.dumps(plan.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path
