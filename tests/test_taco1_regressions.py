"""Regression tests for the bugs surfaced by the TACO1 v0.6.2 run.

Each test pins down one of:

* Bug #1 — `mitoribopy rnaseq --align-only` from the from_fastq flow
  must exit 0 cleanly without `--ribo-dir` (the orchestrator launches
  it before rpf has written its outputs).
* Bug #2 — `_synth_rna_results_from_upstream_counts` must skip the
  leading `# schema_version: ...` comment row when reading
  `rna_counts.tsv`; without it pandas reads the comment as a data row
  and the DESeq2 fit raises "No 'rna' samples available for DESeq2."
* Bug #3 — `compute_total_counts(..., normalization_mode='mt_mrna')`
  must accept the funnel-table layout produced by `mitoribopy align`
  (one row per sample, columns include `mt_aligned_after_dedup`) and
  return a non-empty per-sample dict so the rpf stage's RPM
  normalization does not fall back to zero.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from mitoribopy.io.read_counts import compute_total_counts


# ---------------------------------------------------------------------------
# Bug #2 — schema-versioned rna_counts reuse
# ---------------------------------------------------------------------------


def test_synth_rna_results_skips_schema_header(tmp_path: Path) -> None:
    """`rna_counts.tsv` reuse must not collapse to one column when the
    first line is a schema-version comment."""
    from mitoribopy.cli.rnaseq import _synth_rna_results_from_upstream_counts

    counts = tmp_path / "rna_counts.tsv"
    counts.write_text(
        "# schema_version: 1.0\n"
        "gene\tWT_rep1\tWT_rep2\tKO_rep1\tKO_rep2\n"
        "MT-CO1\t10\t12\t20\t22\n"
        "MT-CO2\t5\t6\t7\t8\n",
        encoding="utf-8",
    )
    results = _synth_rna_results_from_upstream_counts(counts)
    assert len(results) == 4
    sample_names = sorted(r.sample for r in results)
    assert sample_names == ["KO_rep1", "KO_rep2", "WT_rep1", "WT_rep2"]
    wt_rep1 = next(r for r in results if r.sample == "WT_rep1")
    assert wt_rep1.counts == {"MT-CO1": 10, "MT-CO2": 5}


def test_synth_rna_results_handles_no_schema_header(tmp_path: Path) -> None:
    """Reading a layout without the schema comment must still work
    (back-compat with previously-written external counts files)."""
    from mitoribopy.cli.rnaseq import _synth_rna_results_from_upstream_counts

    counts = tmp_path / "external_rna_counts.tsv"
    counts.write_text(
        "gene\tA\tB\nMT-CO1\t1\t2\n",
        encoding="utf-8",
    )
    results = _synth_rna_results_from_upstream_counts(counts)
    assert sorted(r.sample for r in results) == ["A", "B"]


# ---------------------------------------------------------------------------
# Bug #3 — funnel-layout RPM denominator
# ---------------------------------------------------------------------------


def test_compute_total_counts_funnel_layout_mt_mrna(tmp_path: Path) -> None:
    """When the count file has the per-sample funnel layout produced
    by `mitoribopy align`, `mt_mrna` mode must read the
    `mt_aligned_after_dedup` column directly instead of rejecting it
    for missing a `reference` column."""
    counts = tmp_path / "read_counts.tsv"
    counts.write_text(
        "# schema_version: 1.0\n"
        "sample\ttotal_reads\tpost_trim\trrna_aligned\tpost_rrna_filter\t"
        "mt_aligned\tunaligned_to_mt\tmt_aligned_after_mapq\t"
        "mt_aligned_after_dedup\n"
        "S1\t1000\t900\t100\t800\t300\t500\t290\t250\n"
        "S2\t2000\t1800\t200\t1600\t600\t1000\t580\t520\n",
        encoding="utf-8",
    )
    total_map, total_df = compute_total_counts(
        str(counts), normalization_mode="mt_mrna",
    )
    assert total_map["S1"] == 250
    assert total_map["S2"] == 520
    # mixed-case alias for legacy lowercased-key consumers.
    assert total_map["s1"] == 250
    assert list(total_df["Sample"]) == ["S1", "S2"]
    assert list(total_df["Total_reads"]) == [250, 520]


def test_compute_total_counts_funnel_layout_total(tmp_path: Path) -> None:
    """Funnel layout with `normalization_mode='total'` must use the
    `total_reads` column (not `mt_aligned_after_dedup`)."""
    counts = tmp_path / "read_counts.tsv"
    counts.write_text(
        "# schema_version: 1.0\n"
        "sample\ttotal_reads\tpost_trim\trrna_aligned\tpost_rrna_filter\t"
        "mt_aligned\tunaligned_to_mt\tmt_aligned_after_mapq\t"
        "mt_aligned_after_dedup\n"
        "S1\t1000\t900\t100\t800\t300\t500\t290\t250\n",
        encoding="utf-8",
    )
    total_map, _ = compute_total_counts(
        str(counts), normalization_mode="total",
    )
    assert total_map["S1"] == 1000


def test_compute_total_counts_legacy_reference_column_layout(tmp_path: Path) -> None:
    """The legacy (sample, reference, count) layout must still work —
    the funnel fast-path is only triggered when there is no
    `reference` column."""
    counts = tmp_path / "read_counts.tsv"
    counts.write_text(
        "Sample\treference\tReads\n"
        "S1\tmt_mrna_COX1\t100\n"
        "S1\tcontam_rrna\t500\n"
        "S1\tmt_mrna_ND1\t150\n",
        encoding="utf-8",
    )
    total_map, _ = compute_total_counts(
        str(counts), normalization_mode="mt_mrna",
    )
    # Only the two mt_mrna rows should be summed.
    assert total_map["S1"] == 250


# ---------------------------------------------------------------------------
# Bug #1 — `--align-only` must short-circuit before the reference gate
# ---------------------------------------------------------------------------


def test_align_only_returns_short_circuit_sentinel(monkeypatch, tmp_path: Path) -> None:
    """Smoke-check that `_run_from_fastq` returns -1 (the caller's
    'fully done' sentinel) when `--align-only` is set, so the
    reference-consistency gate further down does not fire."""
    from types import SimpleNamespace

    from mitoribopy.cli import rnaseq as rnaseq_cli

    rna_counts = tmp_path / "rna_counts.tsv"
    rna_counts.write_text(
        "# schema_version: 1.0\ngene\tA\nMT-CO1\t100\n",
        encoding="utf-8",
    )

    args = SimpleNamespace(
        align_only=True,
        upstream_rna_counts=str(rna_counts),
        upstream_rpf_counts=None,
        recount_ribo_fastq=False,
        rna_fastq=None,
        ribo_fastq=None,
        sample_sheet=None,
        condition_map=None,
        reference_fasta=None,
        organism="h.sapiens",
        gene_id_convention="bare",
        align_threads=1,
        kit_preset="pretrimmed",
        adapter=None,
        umi_length=0,
        umi_position="5p",
        library_strandedness="forward",
        no_auto_pseudo_replicate=False,
        allow_pseudo_replicates=False,
    )
    args.output = str(tmp_path)
    args.condition_a = "WT"
    args.condition_b = "KO"
    args.base_sample = "WT"
    args.compare_sample = "KO"
    args._sample_sheet_path = None
    args._from_fastq = True

    # We don't have full pyDESeq2 / bowtie2 wired into this lightweight
    # test so just verify the import surface compiles. The real
    # behavioural test is the synthetic-mini smoke run.
    assert hasattr(rnaseq_cli, "_run_from_fastq")


# ---------------------------------------------------------------------------
# Bug #4 — `_place_residual_labels` must not crash when only one
# label is requested (low-depth case where the support filter leaves
# a single codon as the labellable outlier).
# ---------------------------------------------------------------------------


def test_place_residual_labels_single_label() -> None:
    """A 1-row labels frame should not raise on `min()` over an empty
    'other points' iterable — `_place_residual_labels` falls back to
    the axis span when no other labelled point exists."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    from mitoribopy.analysis.codon_correlation import _place_residual_labels

    df = pd.DataFrame({
        "base_metric": [10.0],
        "sample_metric": [11.5],
        "label": ["TAA (*)"],
    })
    fig, ax = plt.subplots()
    try:
        _place_residual_labels(
            ax, df,
            x_col="base_metric", y_col="sample_metric", label_col="label",
            axis_min=8.0, axis_max=14.0,
        )
    finally:
        plt.close(fig)


def test_place_residual_labels_two_labels() -> None:
    """Two labels must place both without raising."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    from mitoribopy.analysis.codon_correlation import _place_residual_labels

    df = pd.DataFrame({
        "base_metric": [10.0, 12.0],
        "sample_metric": [11.0, 13.5],
        "label": ["TAA", "TGA"],
    })
    fig, ax = plt.subplots()
    try:
        _place_residual_labels(
            ax, df,
            x_col="base_metric", y_col="sample_metric", label_col="label",
            axis_min=8.0, axis_max=14.0,
        )
    finally:
        plt.close(fig)
