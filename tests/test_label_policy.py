"""Tests for ``mitoribopy.plotting.label_policy`` (refactor-4 / §5)."""

from __future__ import annotations

from pathlib import Path

import pytest

from mitoribopy.plotting.label_policy import (
    decide_labels,
    write_label_candidates,
)


class TestDecideLabels:
    def test_small_n_labels_all(self) -> None:
        decision = decide_labels(list(range(13)))
        assert decision.policy == "all"
        assert decision.n_labelled == 13
        assert decision.indices == tuple(range(13))

    def test_medium_n_labels_top_15_by_score(self) -> None:
        n = 50
        scores = [1.0] * n
        scores[7] = 99.0
        scores[42] = 80.0
        decision = decide_labels(list(range(n)), scores=scores)
        assert decision.policy == "top_n_small"
        assert decision.n_labelled == 15
        assert 7 in decision.indices
        assert 42 in decision.indices

    def test_large_n_labels_top_10(self) -> None:
        n = 1000
        scores = [float(i) for i in range(n)]
        decision = decide_labels(list(range(n)), scores=scores)
        assert decision.policy == "top_n_large"
        assert decision.n_labelled == 10
        # Top 10 by score = highest indices.
        assert decision.indices == tuple(range(990, 1000))

    def test_explicit_all_overrides_size(self) -> None:
        decision = decide_labels(list(range(500)), annotate="all")
        assert decision.policy == "all"
        assert decision.n_labelled == 500

    def test_explicit_none(self) -> None:
        decision = decide_labels(list(range(13)), annotate="none")
        assert decision.policy == "none"
        assert decision.n_labelled == 0

    def test_explicit_int_labels_top_n(self) -> None:
        scores = [float(i) for i in range(100)]
        decision = decide_labels(list(range(100)), scores=scores, annotate=5)
        assert decision.policy == "top_5"
        assert decision.n_labelled == 5
        assert decision.indices == (95, 96, 97, 98, 99)

    def test_indices_returned_in_input_order(self) -> None:
        n = 30
        scores = [0.0] * n
        scores[0] = 100.0
        scores[29] = 99.0
        decision = decide_labels(list(range(n)), scores=scores)
        assert decision.indices == tuple(sorted(decision.indices))

    def test_empty_input(self) -> None:
        decision = decide_labels([])
        assert decision.policy == "all"
        assert decision.n_labelled == 0


class TestWriteLabelCandidates:
    def test_writes_descending_score_order(self, tmp_path: Path) -> None:
        items = [("a", 1.0), ("b", 5.5), ("c", 3.2), ("d", None)]
        out = write_label_candidates(tmp_path / "cand.tsv", items=items)
        lines = out.read_text(encoding="utf-8").splitlines()
        assert lines[0] == "label\tscore"
        # b > c > a > d (None last).
        assert lines[1].startswith("b\t")
        assert lines[2].startswith("c\t")
        assert lines[3].startswith("a\t")
        assert lines[4] == "d\t"
