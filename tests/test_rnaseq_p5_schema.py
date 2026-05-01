"""Tests for the §5 TE / dTE schema upgrade.

Covers:

* The expanded :class:`TeRow` / :class:`DTeRow` field set and the
  back-compat property aliases.
* The canonical note codes (``publication_grade``,
  ``exploratory_mt_only``, ``insufficient_replicates``,
  ``missing_rna_gene``, ``pseudo_replicate_no_statistics``).
* The new ``rna_only`` mode in :data:`RNASEQ_MODES`.
* The TSV writer column orders.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mitoribopy.rnaseq._types import (
    CANONICAL_NOTE_CODES,
    DTeRow,
    TeRow,
    NOTE_EXPLORATORY_MT_ONLY,
    NOTE_INSUFFICIENT_REPLICATES,
    NOTE_MISSING_RNA_GENE,
    NOTE_PSEUDO_REPLICATE_NO_STATISTICS,
    NOTE_PUBLICATION_GRADE,
)
from mitoribopy.rnaseq.te import compute_delta_te, compute_te


def _de_table(rows):
    from mitoribopy.rnaseq._types import DeColumnMap, DeTable

    cm = DeColumnMap(
        gene_id="gene_id", log2fc="log2FoldChange", padj="padj", basemean="baseMean"
    )
    canonical = [
        {
            "gene_id": r["gene_id"],
            "log2fc": r.get("log2fc"),
            "padj": r.get("padj"),
            "basemean": r.get("basemean"),
        }
        for r in rows
    ]
    return DeTable(format="deseq2", column_map=cm, rows=canonical)


class TestRnaseqModes:
    def test_rna_only_mode_registered(self) -> None:
        from mitoribopy.cli.rnaseq import RNASEQ_MODES

        assert "rna_only" in RNASEQ_MODES
        # The publication-grade and exploratory entries must still be there.
        assert "de_table" in RNASEQ_MODES
        assert "from_fastq" in RNASEQ_MODES
        assert "none" in RNASEQ_MODES


class TestNoteCodes:
    def test_canonical_codes_match_spec(self) -> None:
        # Every code in the §5 spec MUST be exported.
        spec = {
            "publication_grade",
            "exploratory_mt_only",
            "pseudo_replicate_no_statistics",
            "insufficient_replicates",
            "missing_rna_gene",
            "missing_rpf_gene",
            "zero_count_in_condition",
            "gene_id_unmatched",
        }
        assert set(CANONICAL_NOTE_CODES) == spec

    def test_de_table_mode_tags_publication_grade(self) -> None:
        ribo = {"MT-ND1": {"A1": 100, "B1": 150}}
        de = _de_table([
            {"gene_id": "MT-ND1", "log2fc": 0.5, "padj": 0.01, "basemean": 1000.0}
        ])
        rows = compute_delta_te(
            ribo,
            de,
            condition_map={"A1": "A", "B1": "B"},
            base_condition="A",
            compare_condition="B",
            mode="de_table",
        )
        assert rows[0].note == NOTE_PUBLICATION_GRADE
        # ...and a singleton-replicate run with mode=de_table cannot
        # claim publication grade — it should still tag insufficient_replicates.
        ribo_singleton = {"MT-ND1": {"A": 100}}
        rows = compute_delta_te(ribo_singleton, de, mode="de_table")
        assert rows[0].note == NOTE_INSUFFICIENT_REPLICATES

    def test_from_fastq_mode_tags_exploratory(self) -> None:
        rows = compute_te(
            {"MT-ND1": {"A": 100}},
            {"MT-ND1": 1000.0},
            mode="from_fastq",
        )
        assert rows[0].note == NOTE_EXPLORATORY_MT_ONLY

    def test_pseudo_replicate_overrides_mode(self) -> None:
        rows = compute_te(
            {"MT-ND1": {"A": 100}},
            {"MT-ND1": 1000.0},
            mode="from_fastq",
            pseudo_replicate=True,
        )
        assert rows[0].note == NOTE_PSEUDO_REPLICATE_NO_STATISTICS

    def test_missing_rna_gene_replaces_legacy_note(self) -> None:
        ribo = {"MT-NOVEL": {"A": 10, "B": 20}}
        de = _de_table([])
        rows = compute_delta_te(ribo, de)
        assert rows[0].note == NOTE_MISSING_RNA_GENE


class TestTeRowFieldSet:
    def test_terow_has_all_spec_columns(self) -> None:
        r = TeRow(
            sample_id="A",
            gene="MT-ND1",
            rpf_count=10,
            rna_abundance=5.0,
            te=2.0,
            log2_te=1.0,
            condition="WT",
            assay="ribo",
            note=NOTE_PUBLICATION_GRADE,
        )
        for field in (
            "sample_id",
            "condition",
            "assay",
            "gene",
            "rpf_count",
            "rna_abundance",
            "te",
            "log2_te",
            "note",
        ):
            assert hasattr(r, field), f"TeRow missing spec column: {field}"

    def test_terow_back_compat_aliases(self) -> None:
        r = TeRow(
            sample_id="A",
            gene="MT-ND1",
            rpf_count=10,
            rna_abundance=5.0,
            te=2.0,
        )
        assert r.sample == "A"
        assert r.mrna_abundance == 5.0


class TestDTeRowFieldSet:
    def test_dterow_has_all_spec_columns(self) -> None:
        r = DTeRow(
            gene="MT-ND1",
            base_condition="WT",
            compare_condition="KO",
            mrna_log2fc=0.5,
            rpf_log2fc=1.0,
            delta_te_log2=0.5,
            padj_mrna=0.01,
            padj_rpf=0.02,
            padj_delta_te=0.03,
            method="external_de_table",
            note=NOTE_PUBLICATION_GRADE,
        )
        for field in (
            "gene",
            "base_condition",
            "compare_condition",
            "mrna_log2fc",
            "rpf_log2fc",
            "delta_te_log2",
            "padj_mrna",
            "padj_rpf",
            "padj_delta_te",
            "method",
            "note",
        ):
            assert hasattr(r, field), f"DTeRow missing spec column: {field}"

    def test_dterow_padj_alias_resolves_to_padj_mrna(self) -> None:
        r = DTeRow(gene="x", padj_mrna=0.04)
        assert r.padj == 0.04


class TestTsvWriters:
    def test_te_tsv_columns_match_spec(self, tmp_path: Path) -> None:
        from mitoribopy.cli.rnaseq import _write_te_table

        rows = compute_te(
            {"MT-ND1": {"A": 100}},
            {"MT-ND1": 1000.0},
            condition_map={"A": "WT"},
            mode="de_table",
        )
        path = tmp_path / "te.tsv"
        _write_te_table(rows, path)
        lines = [
            l for l in path.read_text().splitlines() if not l.startswith("#")
        ]
        header = lines[0].split("\t")
        assert header == [
            "sample_id",
            "condition",
            "assay",
            "gene",
            "rpf_count",
            "rna_abundance",
            "te",
            "log2_te",
            "note",
        ]
        cells = lines[1].split("\t")
        # spot check that condition + note made it into the row
        assert cells[1] == "WT"
        assert cells[8] == NOTE_PUBLICATION_GRADE

    def test_delta_te_tsv_columns_match_spec(self, tmp_path: Path) -> None:
        from mitoribopy.cli.rnaseq import _write_delta_te_table

        rows = compute_delta_te(
            {"MT-ND1": {"A1": 100, "B1": 200}},
            _de_table(
                [
                    {
                        "gene_id": "MT-ND1",
                        "log2fc": 0.5,
                        "padj": 0.01,
                        "basemean": 100.0,
                    }
                ]
            ),
            condition_map={"A1": "A", "B1": "B"},
            base_condition="A",
            compare_condition="B",
            mode="de_table",
        )
        path = tmp_path / "delta_te.tsv"
        _write_delta_te_table(rows, path)
        lines = [
            l for l in path.read_text().splitlines() if not l.startswith("#")
        ]
        header = lines[0].split("\t")
        assert header == [
            "gene",
            "base_condition",
            "compare_condition",
            "mrna_log2fc",
            "rpf_log2fc",
            "delta_te_log2",
            "padj_mrna",
            "padj_rpf",
            "padj_delta_te",
            "method",
            "note",
        ]
