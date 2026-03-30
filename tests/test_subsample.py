from __future__ import annotations

from pathlib import Path

import pytest

from mitoribopy.tools.subsample import reservoir_sample_lines


def _write_lines(path: Path, n: int) -> None:
    path.write_text("".join(f"line_{idx}\n" for idx in range(n)), encoding="utf-8")


def test_reservoir_sample_is_deterministic_for_seed(tmp_path) -> None:
    input_path = tmp_path / "small.bed"
    _write_lines(input_path, 50)

    first = reservoir_sample_lines(input_path, n=10, seed=7)
    second = reservoir_sample_lines(input_path, n=10, seed=7)

    assert len(first) == 10
    assert first == second


def test_reservoir_sample_n_larger_than_file_returns_all_lines(tmp_path) -> None:
    input_path = tmp_path / "tiny.bed"
    _write_lines(input_path, 3)

    sampled = reservoir_sample_lines(input_path, n=10, seed=1)

    assert len(sampled) == 3


def test_reservoir_sample_rejects_nonpositive_n(tmp_path) -> None:
    input_path = tmp_path / "tiny.bed"
    _write_lines(input_path, 3)

    with pytest.raises(ValueError):
        reservoir_sample_lines(input_path, n=0, seed=1)
