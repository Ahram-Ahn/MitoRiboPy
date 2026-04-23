"""Unit tests for ``mitoribopy.align.adapter_detect``.

The ``opener`` parameter is the dependency-injection point: every test
feeds an in-memory ``StringIO`` so we never touch the filesystem and
never depend on a real gzip codec being available.
"""

from __future__ import annotations

import io
from collections.abc import Callable
from pathlib import Path
from typing import IO

import pytest

from mitoribopy.align.adapter_detect import (
    DetectionResult,
    detect_adapter,
    format_per_kit_rates,
)


# Reference adapter sequences (from KIT_PRESETS — kept inline so the test
# is self-documenting and would fail loudly if a preset's adapter ever
# changes silently).
TRUSEQ_ADAPTER = "TGGAATTCTCGGGTGCCAAGG"
NEB_ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
QIASEQ_ADAPTER = "AACTGTAGGCACCATCAAT"


def _make_fastq(seqs: list[str]) -> str:
    """Build a multi-record FASTQ string from a list of sequence lines."""
    parts = []
    for index, seq in enumerate(seqs):
        parts.append(f"@read{index}\n{seq}\n+\n{'I' * len(seq)}\n")
    return "".join(parts)


def _opener_for(text: str) -> Callable[[Path], IO[str]]:
    """Return an ``opener(path)`` callable that always yields ``text``."""

    def open_text(_path: Path) -> IO[str]:
        return io.StringIO(text)

    return open_text


# ---------- happy-path detection per kit ------------------------------------


def test_detect_truseq_when_all_reads_have_truseq_adapter() -> None:
    seqs = ["ACGTACGT" + TRUSEQ_ADAPTER + "AACTCCAGTCAC"] * 100
    result = detect_adapter(
        Path("dummy.fq"), opener=_opener_for(_make_fastq(seqs))
    )
    assert result.preset_name == "truseq_smallrna"
    assert result.match_rate == pytest.approx(1.0)
    assert result.n_reads_scanned == 100
    assert not result.ambiguous


def test_detect_neb_smallrna_with_ambiguous_flag_for_neb_family_share() -> None:
    seqs = ["ACGTACGT" + NEB_ADAPTER + "AAAAAA"] * 100
    result = detect_adapter(
        Path("dummy.fq"), opener=_opener_for(_make_fastq(seqs))
    )
    # Both nebnext_smallrna and nebnext_ultra_umi share this adapter;
    # alphabetical tie-break picks nebnext_smallrna and ambiguous=True
    # so the user is warned to re-confirm whether they actually have UMIs.
    assert result.preset_name == "nebnext_smallrna"
    assert result.match_rate == pytest.approx(1.0)
    assert result.ambiguous is True
    assert result.per_kit_rates["nebnext_ultra_umi"] == pytest.approx(1.0)


def test_detect_qiaseq_mirna() -> None:
    seqs = ["ACGTACGT" + QIASEQ_ADAPTER + "GGG"] * 100
    result = detect_adapter(
        Path("dummy.fq"), opener=_opener_for(_make_fastq(seqs))
    )
    assert result.preset_name == "qiaseq_mirna"
    assert result.match_rate == pytest.approx(1.0)
    assert not result.ambiguous


# ---------- negative + threshold cases --------------------------------------


def test_detect_returns_none_when_no_adapter_matches() -> None:
    # Random-looking ACGT 24-mers that contain no adapter prefix.
    seqs = ["ACGTACGTACGTACGTACGTACGT"] * 100
    result = detect_adapter(
        Path("dummy.fq"), opener=_opener_for(_make_fastq(seqs))
    )
    assert result.preset_name is None
    assert all(rate < 0.05 for rate in result.per_kit_rates.values())


def test_detect_returns_none_when_match_rate_below_threshold() -> None:
    # 20% TruSeq, 80% adapter-free => below default min_match_rate=0.30
    seqs = ["ACGT" + TRUSEQ_ADAPTER] * 20 + ["ACGTACGTACGTACGTACGT"] * 80
    result = detect_adapter(
        Path("dummy.fq"),
        opener=_opener_for(_make_fastq(seqs)),
    )
    assert result.preset_name is None
    assert result.per_kit_rates["truseq_smallrna"] == pytest.approx(0.20)


def test_detect_passes_at_50_percent_with_default_threshold() -> None:
    # 50% TruSeq, 50% adapter-free => above default min_match_rate=0.30
    seqs = ["ACGT" + TRUSEQ_ADAPTER] * 50 + ["ACGTACGTACGTACGTACGT"] * 50
    result = detect_adapter(
        Path("dummy.fq"), opener=_opener_for(_make_fastq(seqs))
    )
    assert result.preset_name == "truseq_smallrna"
    assert result.match_rate == pytest.approx(0.5)


def test_detect_respects_n_reads_cap() -> None:
    seqs = ["ACGT" + TRUSEQ_ADAPTER] * 1000
    result = detect_adapter(
        Path("dummy.fq"),
        opener=_opener_for(_make_fastq(seqs)),
        n_reads=50,
    )
    assert result.n_reads_scanned == 50


def test_detect_handles_empty_fastq() -> None:
    result = detect_adapter(Path("dummy.fq"), opener=_opener_for(""))
    assert result.preset_name is None
    assert result.n_reads_scanned == 0
    assert result.match_rate == 0.0


def test_detect_handles_truncated_final_record() -> None:
    # Record 0 is intact; record 1 is missing the +/quality lines.
    raw = (
        f"@read0\n{TRUSEQ_ADAPTER}AAAA\n+\n{'I' * (len(TRUSEQ_ADAPTER) + 4)}\n"
        f"@read1\n{TRUSEQ_ADAPTER}AAAA\n"
    )
    result = detect_adapter(Path("dummy.fq"), opener=_opener_for(raw))
    # Both records' sequence lines should have been read; iterator stops
    # cleanly when the underlying stream runs out mid-record.
    assert result.n_reads_scanned >= 1
    assert result.preset_name == "truseq_smallrna"


# ---------- end-of-window / determinism behaviour ---------------------------


def test_detect_is_deterministic_across_runs() -> None:
    seqs = ["ACGT" + NEB_ADAPTER] * 10
    text = _make_fastq(seqs)
    a = detect_adapter(Path("dummy.fq"), opener=_opener_for(text))
    b = detect_adapter(Path("dummy.fq"), opener=_opener_for(text))
    assert a == b


def test_detect_short_adapter_skipped_from_candidates() -> None:
    # A min_match_len longer than every adapter would skip every kit;
    # detection should return None without raising.
    seqs = ["ACGT" + TRUSEQ_ADAPTER + "AAA"] * 20
    result = detect_adapter(
        Path("dummy.fq"),
        opener=_opener_for(_make_fastq(seqs)),
        min_match_len=999,
    )
    assert result.preset_name is None
    # No candidates means per_kit_rates is empty.
    assert result.per_kit_rates == {}


# ---------- format_per_kit_rates ---------------------------------------------


def test_format_per_kit_rates_sorts_alphabetically_and_uses_percents() -> None:
    rates = {"truseq_smallrna": 0.941, "nebnext_smallrna": 0.012}
    formatted = format_per_kit_rates(rates)
    # Sorted alphabetically => nebnext first, truseq second.
    assert formatted == "nebnext_smallrna=1.2%, truseq_smallrna=94.1%"


def test_format_per_kit_rates_handles_empty_dict() -> None:
    assert format_per_kit_rates({}) == ""


# ---------- DetectionResult dataclass ----------------------------------------


def test_detection_result_is_frozen_dataclass() -> None:
    result = DetectionResult(
        preset_name="truseq_smallrna",
        match_rate=0.94,
        per_kit_rates={"truseq_smallrna": 0.94},
        n_reads_scanned=5000,
        ambiguous=False,
    )
    with pytest.raises(Exception):
        # FrozenInstanceError on Python 3.10+; just assert any error.
        result.preset_name = "qiaseq_mirna"  # type: ignore[misc]
