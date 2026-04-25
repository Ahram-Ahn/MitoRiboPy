"""Tests for the per-sample .sample_done resume contract in mitoribopy align.

After a successful sample, the orchestrator writes
``<output>/.sample_done/<sample>.json`` with the full SampleCounts
payload. On a subsequent invocation with ``--resume`` those samples
are reloaded instead of re-run; samples without a marker are
processed normally and their counts merged into the same
``read_counts.tsv``.
"""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import pytest

from mitoribopy.align._types import SampleCounts
from mitoribopy.cli import align as align_cli


def _counts(sample: str, total: int = 1000) -> SampleCounts:
    return SampleCounts(
        sample=sample,
        total_reads=total,
        post_trim=int(total * 0.85),
        rrna_aligned=int(total * 0.35),
        post_rrna_filter=int(total * 0.50),
        mt_aligned=int(total * 0.48),
        unaligned_to_mt=int(total * 0.02),
        mt_aligned_after_mapq=int(total * 0.47),
        mt_aligned_after_dedup=int(total * 0.47),
    )


def test_write_then_load_round_trip(tmp_path) -> None:
    counts = _counts("sampleA")
    align_cli._write_sample_done(tmp_path, counts)
    loaded = align_cli._load_sample_done(tmp_path, "sampleA")
    assert loaded == counts


def test_load_returns_none_when_marker_missing(tmp_path) -> None:
    assert align_cli._load_sample_done(tmp_path, "missing") is None


def test_load_returns_none_when_marker_is_corrupt(tmp_path) -> None:
    """A truncated JSON marker must NOT silently feed a wrong row into
    the aggregated read_counts.tsv -- the resume helper returns None
    so the orchestrator re-runs the sample."""
    bad = align_cli._sample_done_path(tmp_path, "sampleA")
    bad.parent.mkdir(parents=True, exist_ok=True)
    bad.write_text("{not valid json,,,")
    assert align_cli._load_sample_done(tmp_path, "sampleA") is None


def test_load_returns_none_on_schema_drift(tmp_path) -> None:
    """A marker from an older release whose schema no longer matches
    SampleCounts must be ignored, not crash the pipeline."""
    bad = align_cli._sample_done_path(tmp_path, "sampleA")
    bad.parent.mkdir(parents=True, exist_ok=True)
    bad.write_text(json.dumps({"sample": "sampleA", "obsolete_field": 1}))
    assert align_cli._load_sample_done(tmp_path, "sampleA") is None


def test_write_is_atomic_via_rename(tmp_path) -> None:
    """The writer must rename a .tmp into place so a kill mid-write
    leaves either the previous good marker or nothing -- never a
    half-flushed file that the loader would parse as JSON garbage."""
    counts = _counts("sampleA", total=2000)
    align_cli._write_sample_done(tmp_path, counts)
    target = align_cli._sample_done_path(tmp_path, "sampleA")
    assert target.is_file()
    assert not target.with_suffix(target.suffix + ".tmp").exists()
    payload = json.loads(target.read_text(encoding="utf-8"))
    assert payload == asdict(counts)


def test_resume_skips_completed_sample_in_per_sample_loop(
    tmp_path, monkeypatch
) -> None:
    """End-to-end-ish: pre-seed a .sample_done marker for sampleA,
    then run the per-sample loop with two samples and --resume.
    sampleA must be skipped (its cached counts reused); sampleB must
    invoke _process_one_sample exactly once."""
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    sample_done = out_dir / ".sample_done"

    # Pre-seed sampleA as already complete.
    cached_a = _counts("sampleA", total=500)
    align_cli._write_sample_done(sample_done, cached_a)

    # Stub out _process_one_sample so we do not need real BAMs.
    invocations: list[str] = []

    def fake_process(*, output_dir, args, resolution):
        invocations.append(resolution.sample)
        return _counts(resolution.sample, total=999)

    monkeypatch.setattr(align_cli, "_process_one_sample", fake_process)

    # Mimic the per-sample loop body.
    from mitoribopy.align.sample_resolve import SampleResolution
    from mitoribopy.align._types import ResolvedKit

    def make_res(name: str) -> SampleResolution:
        kit = ResolvedKit(
            kit="illumina_smallrna",
            adapter="TGGAATTCTCGGGTGCCAAGG",
            umi_length=0,
            umi_position="5p",
        )
        return SampleResolution(
            sample=name,
            fastq=tmp_path / f"{name}.fq.gz",
            kit=kit,
            dedup_strategy="skip",
            detected_kit=None,
            detection_match_rate=0.0,
            detection_ambiguous=False,
            source="explicit_off",
        )

    resolutions = [make_res("sampleA"), make_res("sampleB")]

    rows: list[SampleCounts] = []
    for resolution in resolutions:
        cached = align_cli._load_sample_done(sample_done, resolution.sample)
        if cached is not None:
            rows.append(cached)
            continue
        counts = align_cli._process_one_sample(
            output_dir=out_dir, args=None, resolution=resolution
        )
        align_cli._write_sample_done(sample_done, counts)
        rows.append(counts)

    # sampleA was reused from cache; only sampleB invoked the real
    # processing path.
    assert invocations == ["sampleB"]
    assert rows[0] == cached_a            # cached value preserved
    assert rows[1].sample == "sampleB"
    assert rows[1].total_reads == 999

    # sampleB now has its own .done marker for the next resume.
    assert align_cli._load_sample_done(sample_done, "sampleB") is not None
