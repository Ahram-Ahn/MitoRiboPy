"""Tests for the v0.9.0 W_OFFSET_FALLBACK warning (was silent before)."""

from __future__ import annotations

import pandas as pd
import pytest

from mitoribopy.analysis.periodicity import (
    _reset_offset_fallback_reporter,
    _resolve_offsets_for_sample,
)
from mitoribopy.errors import E_OFFSET_FALLBACK_USED
from mitoribopy.io import warnings_log


@pytest.fixture(autouse=True)
def _isolated_warnings_log():
    warnings_log.clear()
    _reset_offset_fallback_reporter()
    yield
    warnings_log.clear()
    _reset_offset_fallback_reporter()


def _offset_table(values: dict[int, int]) -> pd.DataFrame:
    return pd.DataFrame({
        "Read Length": list(values.keys()),
        "Most Enriched 5' Offset": list(values.values()),
        "Most Enriched 3' Offset": [-v for v in values.values()],
    })


def test_fallback_warning_emitted_when_sample_missing_from_per_sample_table() -> None:
    per_sample = {"WT_R1": _offset_table({29: 12, 30: 12})}
    combined = _offset_table({29: 13, 30: 13})

    out = _resolve_offsets_for_sample(per_sample, combined, "KO_R1")
    assert out is combined

    records = warnings_log.collected()
    codes = [r.code for r in records]
    assert E_OFFSET_FALLBACK_USED in codes
    rec = next(r for r in records if r.code == E_OFFSET_FALLBACK_USED)
    assert rec.sample_id == "KO_R1"
    assert rec.severity == "warn"
    assert "falling back" in rec.message.lower()


def test_no_warning_when_per_sample_table_is_none() -> None:
    """Intentional combined-mode run (no per-sample dict at all) is silent."""
    combined = _offset_table({29: 13})
    out = _resolve_offsets_for_sample(None, combined, "WT_R1")
    assert out is combined
    assert not [r for r in warnings_log.collected() if r.code == E_OFFSET_FALLBACK_USED]


def test_no_warning_when_sample_present_in_per_sample_table() -> None:
    per_sample = {"WT_R1": _offset_table({29: 12})}
    combined = _offset_table({29: 13})
    out = _resolve_offsets_for_sample(per_sample, combined, "WT_R1")
    assert out is per_sample["WT_R1"]
    assert not [r for r in warnings_log.collected() if r.code == E_OFFSET_FALLBACK_USED]


def test_warning_dedups_per_stage_sample_pair() -> None:
    per_sample = {"WT_R1": _offset_table({29: 12})}
    combined = _offset_table({29: 13})

    for _ in range(5):
        _resolve_offsets_for_sample(per_sample, combined, "KO_R1", stage="PERIODICITY")

    fallback_records = [
        r for r in warnings_log.collected() if r.code == E_OFFSET_FALLBACK_USED
    ]
    assert len(fallback_records) == 1


def test_warning_can_be_suppressed_with_record_warning_false() -> None:
    per_sample = {"WT_R1": _offset_table({29: 12})}
    combined = _offset_table({29: 13})
    _resolve_offsets_for_sample(
        per_sample, combined, "KO_R1", record_warning=False,
    )
    assert not [r for r in warnings_log.collected() if r.code == E_OFFSET_FALLBACK_USED]
