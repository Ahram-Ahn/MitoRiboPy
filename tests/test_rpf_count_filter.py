"""Unit tests for the read-length auto-filter (``_apply_rpf_count_filter``).

The filter is data-driven: it drops read-length bins whose total count
across all samples is below ``--rpf_min_count_frac`` x the most-enriched
length's count. This shrinks ``context.rpf_range`` and re-filters
``context.filtered_bed_df`` so the downstream offset-enrichment step
only sees signal-bearing bins.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from mitoribopy.pipeline.context import PipelineContext
from mitoribopy.pipeline.steps import _apply_rpf_count_filter


def _make_context(read_lengths: list[int], counts: list[int], frac: float) -> PipelineContext:
    """Build a minimal PipelineContext with a synthetic filtered_bed_df."""
    rows = []
    for rl, n in zip(read_lengths, counts):
        for i in range(n):
            rows.append(
                {
                    "chrom": "tx",
                    "start": i,
                    "end": i + rl,
                    "read_length": rl,
                    "sample_name": "S1",
                }
            )
    df = pd.DataFrame(rows)
    args = SimpleNamespace(rpf_min_count_frac=frac)
    return PipelineContext(
        args=args,
        base_output_dir=Path("/tmp/unused"),
        plot_output_dir=Path("/tmp/unused/plots"),
        csv_dir=Path("/tmp/unused/plots/csv"),
        plot_subdir=Path("/tmp/unused/plots/plots"),
        annotation_df=pd.DataFrame(),
        resolved_codon_table={},
        resolved_start_codons=[],
        rpf_range=list(read_lengths),
        unfiltered_read_length_range=(min(read_lengths), max(read_lengths)),
        filtered_bed_df=df,
        sample_dirs=["S1"],
    )


def _capture_status() -> tuple[list[str], callable]:
    captured: list[str] = []
    return captured, captured.append


def test_filter_drops_bins_below_threshold() -> None:
    # Most enriched length (30) has 1000 reads; threshold at 20% = 200.
    # Read length 27 has 100 (drop), 28 has 250 (keep), 29 has 800 (keep),
    # 30 has 1000 (keep), 31 has 50 (drop).
    ctx = _make_context(
        read_lengths=[27, 28, 29, 30, 31],
        counts=[100, 250, 800, 1000, 50],
        frac=0.20,
    )
    log, emit = _capture_status()

    _apply_rpf_count_filter(ctx, emit)

    assert ctx.rpf_range == [28, 29, 30]
    assert sorted(ctx.filtered_bed_df["read_length"].unique().tolist()) == [28, 29, 30]
    assert any("kept [28, 29, 30]" in line for line in log)
    assert any("dropped [27, 31]" in line for line in log)
    assert any("20%" in line for line in log)


def test_filter_disabled_when_frac_zero() -> None:
    ctx = _make_context(
        read_lengths=[27, 30],
        counts=[10, 1000],
        frac=0.0,
    )
    log, emit = _capture_status()

    _apply_rpf_count_filter(ctx, emit)

    assert ctx.rpf_range == [27, 30]
    assert log == []
    # The DataFrame is untouched.
    assert sorted(ctx.filtered_bed_df["read_length"].unique().tolist()) == [27, 30]


def test_filter_quiet_when_nothing_pruned() -> None:
    # All bins are equally populated => nothing to drop.
    ctx = _make_context(
        read_lengths=[28, 29, 30],
        counts=[500, 500, 500],
        frac=0.20,
    )
    log, emit = _capture_status()

    _apply_rpf_count_filter(ctx, emit)

    assert ctx.rpf_range == [28, 29, 30]
    assert log == []  # quiet when there is nothing to report


def test_filter_handles_empty_dataframe() -> None:
    df = pd.DataFrame(columns=["chrom", "start", "end", "read_length", "sample_name"])
    args = SimpleNamespace(rpf_min_count_frac=0.20)
    ctx = PipelineContext(
        args=args,
        base_output_dir=Path("/tmp/unused"),
        plot_output_dir=Path("/tmp/unused/plots"),
        csv_dir=Path("/tmp/unused/plots/csv"),
        plot_subdir=Path("/tmp/unused/plots/plots"),
        annotation_df=pd.DataFrame(),
        resolved_codon_table={},
        resolved_start_codons=[],
        rpf_range=[28, 29, 30],
        unfiltered_read_length_range=(28, 30),
        filtered_bed_df=df,
        sample_dirs=[],
    )
    log, emit = _capture_status()

    _apply_rpf_count_filter(ctx, emit)

    assert ctx.rpf_range == [28, 29, 30]
    assert log == []


def test_filter_threshold_is_inclusive() -> None:
    # A bin exactly at the threshold (20% of max) is kept.
    ctx = _make_context(
        read_lengths=[28, 30],
        counts=[200, 1000],  # 200 == 0.20 * 1000
        frac=0.20,
    )
    log, emit = _capture_status()

    _apply_rpf_count_filter(ctx, emit)

    assert ctx.rpf_range == [28, 30]  # both retained
    assert log == []  # nothing pruned, quiet


# ---------- Phase 2.3: rpf_counts.metadata.json sidecar --------------------


def _ctx_with_args(read_lengths, counts, *, frac=0.20, fasta=None) -> PipelineContext:
    """Like _make_context but with a richer args namespace covering the
    fields the metadata writer reads (offset_*, footprint_class, strain,
    fasta)."""
    rows = []
    for rl, n in zip(read_lengths, counts):
        for i in range(n):
            rows.append({
                "chrom": "tx",
                "start": i,
                "end": i + rl,
                "read_length": rl,
                "sample_name": "S1",
            })
    df = pd.DataFrame(rows)
    args = SimpleNamespace(
        rpf_min_count_frac=frac,
        offset_mode="per_sample",
        offset_type="5",
        offset_site="p",
        footprint_class="monosome",
        strain="h.sapiens",
        fasta=str(fasta) if fasta else None,
    )
    return PipelineContext(
        args=args,
        base_output_dir=Path("/tmp/unused"),
        plot_output_dir=Path("/tmp/unused/plots"),
        csv_dir=Path("/tmp/unused/plots/csv"),
        plot_subdir=Path("/tmp/unused/plots/plots"),
        annotation_df=pd.DataFrame(),
        resolved_codon_table={},
        resolved_start_codons=[],
        rpf_range=list(read_lengths),
        unfiltered_read_length_range=(min(read_lengths), max(read_lengths)),
        filtered_bed_df=df,
        sample_dirs=["S1"],
    )


def test_metadata_sidecar_records_filter_decision_when_pruned(tmp_path: Path) -> None:
    """When the read-length auto-filter drops bins, the sidecar must
    record both the requested range and the survivors with the
    threshold rule + value spelled out."""
    import json

    from mitoribopy.pipeline.steps import write_rpf_counts_metadata

    ctx = _ctx_with_args(
        read_lengths=[27, 28, 29, 30, 31],
        counts=[100, 250, 800, 1000, 50],
        frac=0.20,
    )
    _apply_rpf_count_filter(ctx, lambda _: None)
    out = tmp_path / "rpf_counts.metadata.json"
    write_rpf_counts_metadata(context=ctx, output_path=out)

    payload = json.loads(out.read_text())
    assert payload["schema_version"] == "1.0.0"
    assert payload["subcommand"] == "rpf"
    assert payload["counts_path"] == "rpf_counts.tsv"

    rlf = payload["read_length_filter"]
    assert rlf["requested_rpf_range"] == [27, 28, 29, 30, 31]
    assert rlf["retained_lengths"] == [28, 29, 30]
    assert rlf["dropped_lengths"] == [27, 31]
    assert rlf["threshold_rule"] == "dominant_fraction"
    assert rlf["threshold_value"] == pytest.approx(0.20)
    # Observed counts cover EVERY observed length, including dropped ones.
    obs = rlf["observed_counts_by_length"]
    assert obs == {"27": 100, "28": 250, "29": 800, "30": 1000, "31": 50}

    # Offset block is faithfully threaded from args.
    assert payload["offset_mode"] == "per_sample"
    assert payload["offset_type"] == "5"
    assert payload["offset_site"] == "p"
    # Footprint class + strain.
    assert payload["footprint_class"] == "monosome"
    assert payload["strain"] == "h.sapiens"
    # Reference fields default to None when no FASTA was passed.
    assert payload["reference_fasta"] is None
    assert payload["reference_checksum"] is None
    # Aggregate row counts.
    assert payload["n_samples"] == 1
    assert payload["n_genes"] == 1
    assert payload["total_reads"] == 250 + 800 + 1000  # only retained reads


def test_metadata_sidecar_records_disabled_filter() -> None:
    """When --rpf_min_count_frac is 0, threshold_rule must be 'disabled'
    and dropped_lengths must be empty even though observed_counts is
    populated."""
    import json
    import tempfile

    from mitoribopy.pipeline.steps import write_rpf_counts_metadata

    ctx = _ctx_with_args(
        read_lengths=[28, 30],
        counts=[10, 1000],
        frac=0.0,
    )
    _apply_rpf_count_filter(ctx, lambda _: None)
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "rpf_counts.metadata.json"
        write_rpf_counts_metadata(context=ctx, output_path=out)
        payload = json.loads(out.read_text())
    rlf = payload["read_length_filter"]
    assert rlf["threshold_rule"] == "disabled"
    assert rlf["threshold_value"] == 0.0
    assert rlf["dropped_lengths"] == []
    assert rlf["retained_lengths"] == [28, 30]
    assert rlf["observed_counts_by_length"] == {"28": 10, "30": 1000}


def test_metadata_sidecar_hashes_real_fasta(tmp_path: Path) -> None:
    """When args.fasta points at a real file, the SHA-256 must land in
    reference_checksum so a downstream rnaseq run can verify it
    matches without re-reading the file."""
    import hashlib
    import json

    from mitoribopy.pipeline.steps import write_rpf_counts_metadata

    fa = tmp_path / "tx.fa"
    fa.write_text(">ND1\nATGGCC\n", encoding="utf-8")
    expected = hashlib.sha256(fa.read_bytes()).hexdigest()

    ctx = _ctx_with_args(
        read_lengths=[30],
        counts=[100],
        frac=0.0,
        fasta=fa,
    )
    _apply_rpf_count_filter(ctx, lambda _: None)
    out = tmp_path / "meta.json"
    write_rpf_counts_metadata(context=ctx, output_path=out)
    payload = json.loads(out.read_text())
    assert payload["reference_fasta"] == str(fa)
    assert payload["reference_checksum"] == expected
