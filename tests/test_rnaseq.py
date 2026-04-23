"""Unit tests for ``mitoribopy.rnaseq`` (Phase 5)."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from mitoribopy.rnaseq import (
    ReferenceMismatchError,
    compute_delta_te,
    compute_reference_checksum,
    compute_te,
    detect_de_format,
    load_de_table,
    load_ribo_counts,
    match_mt_mrnas,
    normalize_gene_id,
    verify_reference_consistency,
)
from mitoribopy.rnaseq._types import DeColumnMap


# ---------- gene_ids --------------------------------------------------------


def test_match_mt_mrnas_human_hgnc_full_13() -> None:
    de_ids = [
        "MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8", "MT-ATP6",
        "MT-CO3", "MT-ND3", "MT-ND4L", "MT-ND4", "MT-ND5", "MT-ND6", "MT-CYB",
    ]
    out = match_mt_mrnas(de_ids, "hgnc", "h")
    assert len(out["matched"]) == 13
    assert out["missing"] == []


def test_match_mt_mrnas_bare_case_insensitive() -> None:
    de_ids = ["nd1", "Nd2", "CO1"]
    out = match_mt_mrnas(de_ids, "bare", "h")
    assert "ND1" in out["matched"]
    assert "ND2" in out["matched"]
    assert "CO1" in out["matched"]


def test_match_mt_mrnas_reports_missing() -> None:
    de_ids = ["MT-ND1"]
    out = match_mt_mrnas(de_ids, "hgnc", "h")
    assert out["matched"] == ["MT-ND1"]
    assert "MT-CYB" in out["missing"]


def test_match_mt_mrnas_yeast_bare() -> None:
    de_ids = ["COX1", "COX2", "COB", "ATP6"]
    out = match_mt_mrnas(de_ids, "bare", "y")
    assert "COX1" in out["matched"]
    assert "ATP6" in out["matched"]


def test_normalize_gene_id_strips_and_uppercases() -> None:
    assert normalize_gene_id("  mt-nd1  ") == "MT-ND1"


# ---------- de_loader -------------------------------------------------------


def test_detect_de_format_deseq2() -> None:
    header = ["gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    fmt, mapping = detect_de_format(header)
    assert fmt == "deseq2"
    assert mapping.log2fc == "log2FoldChange"
    assert mapping.padj == "padj"
    assert mapping.basemean == "baseMean"


def test_detect_de_format_xtail() -> None:
    header = ["gene_id", "mRNA_log2FC_A", "log2FC_TE_final", "pvalue.adjust"]
    fmt, mapping = detect_de_format(header)
    assert fmt == "xtail"
    assert mapping.log2fc == "log2FC_TE_final"
    assert mapping.padj == "pvalue.adjust"


def test_detect_de_format_anota2seq() -> None:
    header = ["gene_id", "slopeTranslation", "pAdjustedRvmRvm", "TotalExpression"]
    fmt, mapping = detect_de_format(header)
    assert fmt == "anota2seq"
    assert mapping.log2fc == "slopeTranslation"


def test_detect_de_format_custom_when_unknown() -> None:
    header = ["gene_id", "weird_col_a", "weird_col_b"]
    fmt, mapping = detect_de_format(header)
    assert fmt == "custom"
    assert mapping.log2fc is None


def test_load_de_table_deseq2(tmp_path) -> None:
    path = tmp_path / "de.tsv"
    path.write_text(
        "gene_id\tbaseMean\tlog2FoldChange\tpadj\n"
        "MT-ND1\t1000.0\t1.2\t0.01\n"
        "MT-ND2\t500.0\t-0.5\t0.05\n"
    )
    table = load_de_table(path)
    assert table.format == "deseq2"
    assert len(table.rows) == 2
    assert table.rows[0]["gene_id"] == "MT-ND1"
    assert table.rows[0]["log2fc"] == pytest.approx(1.2)
    assert table.rows[0]["padj"] == pytest.approx(0.01)
    assert table.rows[0]["basemean"] == pytest.approx(1000.0)


def test_load_de_table_handles_NA_values(tmp_path) -> None:
    path = tmp_path / "de.tsv"
    path.write_text(
        "gene_id\tbaseMean\tlog2FoldChange\tpadj\n"
        "MT-ND1\tNA\t1.2\tNA\n"
    )
    table = load_de_table(path)
    assert table.rows[0]["basemean"] is None
    assert table.rows[0]["padj"] is None
    assert table.rows[0]["log2fc"] == pytest.approx(1.2)


def test_load_de_table_custom_column_map(tmp_path) -> None:
    path = tmp_path / "de.csv"
    path.write_text(
        "my_gene,my_l2fc,my_p\n"
        "MT-ND1,0.5,0.01\n"
    )
    mapping = DeColumnMap(gene_id="my_gene", log2fc="my_l2fc", padj="my_p", basemean=None)
    table = load_de_table(path, column_map=mapping)
    assert len(table.rows) == 1
    assert table.rows[0]["log2fc"] == pytest.approx(0.5)
    assert table.rows[0]["basemean"] is None


# ---------- reference_gate --------------------------------------------------


def test_compute_reference_checksum_matches_hashlib(tmp_path) -> None:
    path = tmp_path / "ref.fa"
    payload = b">chrM\nACGTACGT\n"
    path.write_bytes(payload)
    expected = hashlib.sha256(payload).hexdigest()
    assert compute_reference_checksum(path) == expected


def test_verify_reference_consistency_matches(tmp_path) -> None:
    ribo_dir = tmp_path / "rpf_out"
    ribo_dir.mkdir()
    ref = tmp_path / "ref.fa"
    ref.write_bytes(b">m\nAAAA\n")
    digest = compute_reference_checksum(ref)

    (ribo_dir / "run_settings.json").write_text(
        json.dumps({"reference_checksum": digest})
    )

    result = verify_reference_consistency(ribo_dir=ribo_dir, reference_path=ref)
    assert result.lower() == digest.lower()


def test_verify_reference_consistency_raises_on_mismatch(tmp_path) -> None:
    ribo_dir = tmp_path / "rpf_out"
    ribo_dir.mkdir()
    ref = tmp_path / "ref.fa"
    ref.write_bytes(b">m\nAAAA\n")

    (ribo_dir / "run_settings.json").write_text(
        json.dumps({"reference_checksum": "a" * 64})
    )

    with pytest.raises(ReferenceMismatchError):
        verify_reference_consistency(ribo_dir=ribo_dir, reference_path=ref)


def test_verify_reference_consistency_raises_when_no_recorded_digest(tmp_path) -> None:
    ribo_dir = tmp_path / "rpf_out"
    ribo_dir.mkdir()
    ref = tmp_path / "ref.fa"
    ref.write_bytes(b">m\nAAAA\n")

    with pytest.raises(ReferenceMismatchError):
        verify_reference_consistency(ribo_dir=ribo_dir, reference_path=ref)


def test_verify_reference_consistency_requires_exactly_one_of_two(tmp_path) -> None:
    ribo_dir = tmp_path / "rpf_out"
    ribo_dir.mkdir()
    with pytest.raises(ValueError):
        verify_reference_consistency(ribo_dir=ribo_dir)
    with pytest.raises(ValueError):
        verify_reference_consistency(
            ribo_dir=ribo_dir,
            reference_path=tmp_path / "a",
            reference_checksum="abc",
        )


# ---------- counts ---------------------------------------------------------


def test_load_ribo_counts(tmp_path) -> None:
    path = tmp_path / "rpf_counts.tsv"
    path.write_text(
        "sample\tgene\tcount\n"
        "A\tMT-ND1\t100\n"
        "B\tMT-ND1\t80\n"
        "A\tMT-CO1\t200\n"
    )
    counts = load_ribo_counts(path)
    assert counts["MT-ND1"] == {"A": 100, "B": 80}
    assert counts["MT-CO1"] == {"A": 200}


def test_load_ribo_counts_raises_on_missing_columns(tmp_path) -> None:
    path = tmp_path / "bad.tsv"
    path.write_text("sample\tvalue\nA\t100\n")
    with pytest.raises(ValueError):
        load_ribo_counts(path)


# ---------- TE / delta-TE --------------------------------------------------


def test_compute_te_produces_rows_for_each_sample_gene_pair() -> None:
    ribo = {"MT-ND1": {"A": 100, "B": 80}, "MT-CO1": {"A": 200}}
    mrna = {"MT-ND1": 1000.0, "MT-CO1": 2000.0}
    rows = compute_te(ribo, mrna)
    assert len(rows) == 3
    by_key = {(r.sample, r.gene): r for r in rows}
    assert by_key[("A", "MT-ND1")].te == pytest.approx((100 + 0.5) / (1000 + 0.5))


def test_compute_te_skips_genes_missing_from_de() -> None:
    ribo = {"MT-ND1": {"A": 100}, "MT-ND99": {"A": 50}}
    mrna = {"MT-ND1": 1000.0}
    rows = compute_te(ribo, mrna)
    assert [r.gene for r in rows] == ["MT-ND1"]


def _synthetic_de_table(rows: list[dict]):
    from mitoribopy.rnaseq._types import DeTable

    column_map = DeColumnMap(
        gene_id="gene_id", log2fc="log2FoldChange", padj="padj", basemean="baseMean"
    )
    canonical = [
        {"gene_id": r["gene_id"], "log2fc": r.get("log2fc"),
         "padj": r.get("padj"), "basemean": r.get("basemean")}
        for r in rows
    ]
    return DeTable(format="deseq2", column_map=column_map, rows=canonical)


def test_compute_delta_te_single_replicate_records_note() -> None:
    ribo = {"MT-ND1": {"A": 100}}
    de = _synthetic_de_table([
        {"gene_id": "MT-ND1", "log2fc": 0.5, "padj": 0.01, "basemean": 1000.0},
    ])
    dte = compute_delta_te(ribo, de)
    assert len(dte) == 1
    assert dte[0].note == "single_replicate_no_statistics"
    assert dte[0].rpf_log2fc is None
    assert dte[0].mrna_log2fc == pytest.approx(0.5)


def test_compute_delta_te_with_replicates() -> None:
    ribo = {"MT-ND1": {"A1": 100, "A2": 110, "B1": 200, "B2": 220}}
    de = _synthetic_de_table([
        {"gene_id": "MT-ND1", "log2fc": 0.5, "padj": 0.01, "basemean": 1000.0},
    ])
    dte = compute_delta_te(
        ribo,
        de,
        condition_map={"A1": "A", "A2": "A", "B1": "B", "B2": "B"},
        condition_a="A",
        condition_b="B",
    )
    assert dte[0].rpf_log2fc is not None
    # B/A ~= 2 -> log2 ~= 1
    assert dte[0].rpf_log2fc == pytest.approx(1.0, abs=0.05)
    # delta_te = rpf_log2fc - mrna_log2fc
    assert dte[0].delta_te_log2 == pytest.approx(
        dte[0].rpf_log2fc - dte[0].mrna_log2fc, abs=1e-6
    )


def test_compute_delta_te_flags_gene_missing_from_de() -> None:
    ribo = {"MT-NOVEL": {"A": 10, "B": 20}}
    de = _synthetic_de_table([])
    dte = compute_delta_te(ribo, de)
    assert dte[0].note == "missing_from_de_table"
    assert dte[0].mrna_log2fc is None
