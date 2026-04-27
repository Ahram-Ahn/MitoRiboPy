"""Tests for the from-FASTQ DE-analysis helpers (no pyDESeq2 required)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mitoribopy.align._types import ResolvedKit
from mitoribopy.rnaseq.alignment import SampleAlignmentResult
from mitoribopy.rnaseq.de_analysis import (
    build_sample_sheet,
    deseq2_to_de_table,
    write_de_table_tsv,
)
from mitoribopy.rnaseq.de_loader import load_de_table


def _result(name: str, counts: dict[str, int]) -> SampleAlignmentResult:
    return SampleAlignmentResult(
        sample=name,
        bam_path=Path("/tmp/dummy.bam"),
        counts=counts,
        paired=False,
        total_reads=sum(counts.values()),
        aligned_reads=sum(counts.values()),
        resolved_kit=ResolvedKit(
            kit="pretrimmed", adapter=None, umi_length=0, umi_position="5p"
        ),
    )


def test_build_sample_sheet_aligns_samples_in_order() -> None:
    rna = [_result("R1", {"g1": 10, "g2": 20}), _result("R2", {"g1": 12, "g2": 22})]
    ribo = [_result("B1", {"g1": 5}), _result("B2", {"g1": 7})]
    cmap = {"R1": "ctrl", "R2": "treated", "B1": "ctrl", "B2": "treated"}
    samples, counts_df, metadata_df = build_sample_sheet(rna, ribo, cmap)

    assert samples == ["R1", "R2", "B1", "B2"]
    assert list(counts_df.index) == samples
    assert list(metadata_df.index) == samples
    assert metadata_df.loc["R1", "assay"] == "rna"
    assert metadata_df.loc["B1", "assay"] == "ribo"
    assert metadata_df.loc["R2", "condition"] == "treated"
    # Zero-fill for genes that one assay never saw:
    assert counts_df.loc["B1", "g2"] == 0
    assert counts_df.loc["R1", "g1"] == 10


def test_build_sample_sheet_missing_condition_raises() -> None:
    rna = [_result("R1", {"g1": 1}), _result("R2", {"g1": 2})]
    ribo: list[SampleAlignmentResult] = []
    cmap = {"R1": "ctrl"}  # R2 missing
    with pytest.raises(ValueError, match="missing entries"):
        build_sample_sheet(rna, ribo, cmap)


def test_build_sample_sheet_single_level_condition_raises() -> None:
    rna = [_result("R1", {"g1": 1}), _result("R2", {"g1": 2})]
    ribo: list[SampleAlignmentResult] = []
    cmap = {"R1": "ctrl", "R2": "ctrl"}
    with pytest.raises(ValueError, match="fewer than two"):
        build_sample_sheet(rna, ribo, cmap)


def test_deseq2_to_de_table_converts_nan_and_inf_to_none() -> None:
    df = pd.DataFrame(
        {
            "baseMean": [100.0, 50.0, np.nan],
            "log2FoldChange": [0.5, np.inf, -np.inf],
            "padj": [0.01, np.nan, 1.0],
        },
        index=["g1", "g2", "g3"],
    )
    table = deseq2_to_de_table(df)
    rows = {r["gene_id"]: r for r in table.rows}
    assert table.format == "deseq2"

    assert rows["g1"]["log2fc"] == 0.5
    assert rows["g1"]["padj"] == 0.01
    assert rows["g1"]["basemean"] == 100.0

    assert rows["g2"]["log2fc"] is None  # +Inf -> None
    assert rows["g2"]["padj"] is None    # NaN -> None
    assert rows["g3"]["basemean"] is None
    assert rows["g3"]["log2fc"] is None  # -Inf -> None


def test_write_de_table_tsv_uses_canonical_deseq2_header(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "baseMean": [100.0, 50.0],
            "log2FoldChange": [0.5, -0.3],
            "padj": [0.01, 0.05],
        },
        index=["MT-ND1", "MT-CO1"],
    )
    table = deseq2_to_de_table(df)
    out = tmp_path / "de.tsv"
    write_de_table_tsv(table, out)

    lines = out.read_text().splitlines()
    assert lines[0] == "gene_id\tbaseMean\tlog2FoldChange\tpadj"

    # And it must roundtrip through the existing loader with auto-detection.
    reloaded = load_de_table(out)
    assert reloaded.format == "deseq2"
    by_id = {r["gene_id"]: r for r in reloaded.rows}
    assert math.isclose(by_id["MT-ND1"]["log2fc"], 0.5, rel_tol=1e-3)
    assert math.isclose(by_id["MT-CO1"]["padj"], 0.05, rel_tol=1e-3)
