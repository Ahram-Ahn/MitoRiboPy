"""End-to-end test of the ``mitoribopy rnaseq`` CLI orchestrator."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from mitoribopy import cli
from mitoribopy.rnaseq import compute_reference_checksum


def _make_ribo_dir(tmp_path: Path, reference_hash: str) -> Path:
    """Build a minimal rpf-like output directory."""
    ribo_dir = tmp_path / "rpf_out"
    ribo_dir.mkdir()
    (ribo_dir / "run_settings.json").write_text(
        json.dumps({"subcommand": "rpf", "reference_checksum": reference_hash})
    )
    (ribo_dir / "rpf_counts.tsv").write_text(
        "sample\tgene\tcount\n"
        "A1\tMT-ND1\t100\n"
        "A2\tMT-ND1\t110\n"
        "B1\tMT-ND1\t200\n"
        "B2\tMT-ND1\t220\n"
        "A1\tMT-CO1\t500\n"
        "B1\tMT-CO1\t550\n"
    )
    return ribo_dir


def _make_de_table(tmp_path: Path) -> Path:
    path = tmp_path / "de.tsv"
    path.write_text(
        "gene_id\tbaseMean\tlog2FoldChange\tpadj\n"
        "MT-ND1\t1000.0\t0.5\t0.01\n"
        "MT-CO1\t2000.0\t-0.3\t0.05\n"
    )
    return path


def test_rnaseq_dry_run_prints_plan_and_exits_zero(capsys, tmp_path) -> None:
    exit_code = cli.main([
        "rnaseq",
        "--dry-run",
        "--gene-id-convention",
        "hgnc",
    ])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "dry-run" in out
    assert "hgnc" in out


def test_rnaseq_errors_listing_every_missing_required_arg(capsys) -> None:
    """No flags → default flow (from raw FASTQ) is selected; the error
    message must list every from-FASTQ required arg plus a HINT pointing
    at the alternative --de-table flow.
    """
    exit_code = cli.main(["rnaseq"])
    err = capsys.readouterr().err
    assert exit_code == 2
    for expected_flag in (
        "--rna-fastq",
        "--reference-fasta",
        "--gene-id-convention",
        "--output",
        "--condition-map",
        "--condition-a",
        "--condition-b",
    ):
        assert expected_flag in err
    # The HINT line tells users about the alternative flow.
    assert "--de-table" in err


def test_rnaseq_errors_listing_de_table_flow_required_args(capsys, tmp_path) -> None:
    """When --de-table IS set, the missing-args list switches to the
    DE-table-flow required args (no --rna-fastq, --reference-fasta etc.).
    """
    de = tmp_path / "de.tsv"
    de.write_text("gene_id\tlog2FoldChange\tpadj\tbaseMean\n")
    exit_code = cli.main(["rnaseq", "--de-table", str(de)])
    err = capsys.readouterr().err
    assert exit_code == 2
    for expected_flag in (
        "--gene-id-convention",
        "--ribo-dir or --ribo-counts",
        "--output",
        "--reference-gtf or --reference-checksum",
    ):
        assert expected_flag in err
    # The from-FASTQ-only flags should NOT appear.
    assert "--rna-fastq" not in err
    assert "--reference-fasta" not in err


def test_rnaseq_end_to_end_produces_te_and_delta_te(tmp_path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_bytes(b">MT-ND1\nACGTACGTACGT\n")
    ref_hash = compute_reference_checksum(ref)

    ribo_dir = _make_ribo_dir(tmp_path, ref_hash)
    de_table = _make_de_table(tmp_path)

    condition_map = tmp_path / "conditions.tsv"
    condition_map.write_text(
        "sample\tcondition\n"
        "A1\tA\nA2\tA\nB1\tB\nB2\tB\n"
    )

    out_dir = tmp_path / "rnaseq_out"

    exit_code = cli.main([
        "rnaseq",
        "--de-table", str(de_table),
        "--gene-id-convention", "hgnc",
        "--ribo-dir", str(ribo_dir),
        "--reference-gtf", str(ref),
        "--condition-map", str(condition_map),
        "--condition-a", "A",
        "--condition-b", "B",
        "--output", str(out_dir),
    ])
    assert exit_code == 0

    # Output files present
    assert (out_dir / "te.tsv").is_file()
    assert (out_dir / "delta_te.tsv").is_file()
    assert (out_dir / "plots" / "mrna_vs_rpf.png").is_file()
    assert (out_dir / "plots" / "delta_te_volcano.png").is_file()
    assert (out_dir / "run_settings.json").is_file()

    # TE table has a row per (sample, gene) present in both counts and DE.
    # P1.12: a `# schema_version: X.Y` comment line precedes the column header.
    raw_te = (out_dir / "te.tsv").read_text().splitlines()
    assert raw_te[0].startswith("# schema_version:")
    te_lines = [line for line in raw_te if not line.startswith("#")]
    # P5.5: te.tsv schema 2.0 — assessment §5 column set.
    assert te_lines[0].split("\t") == [
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
    assert len(te_lines) - 1 >= 6  # 4 ND1 rows + 2 CO1 rows

    # delta_te table has per-gene rows with replicate-based RPF log2FC
    raw_dte = (out_dir / "delta_te.tsv").read_text().splitlines()
    assert raw_dte[0].startswith("# schema_version:")
    dte_lines = [line for line in raw_dte if not line.startswith("#")]
    header = dte_lines[0].split("\t")
    nd1_row = [l.split("\t") for l in dte_lines[1:] if l.startswith("MT-ND1")][0]
    nd1 = dict(zip(header, nd1_row))
    # ND1: A mean ~105, B mean ~210 -> ratio ~2 -> rpf_log2fc ~1
    assert nd1["rpf_log2fc"]
    assert abs(float(nd1["rpf_log2fc"]) - 1.0) < 0.1
    # delta_te = rpf_log2fc - mrna_log2fc
    assert nd1["delta_te_log2"]
    assert abs(float(nd1["delta_te_log2"]) - (float(nd1["rpf_log2fc"]) - 0.5)) < 1e-6


def test_rnaseq_de_table_aliases_are_remapped_to_rpf_gene_names(tmp_path, capsys) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_bytes(b">ND1\nACGT\n>COX1\nACGT\n>CYTB\nACGT\n")
    ref_hash = compute_reference_checksum(ref)

    ribo_dir = tmp_path / "rpf_out"
    ribo_dir.mkdir()
    (ribo_dir / "run_settings.json").write_text(
        json.dumps({"subcommand": "rpf", "reference_checksum": ref_hash})
    )
    (ribo_dir / "rpf_counts.tsv").write_text(
        "sample\tgene\tcount\n"
        "A\tND1\t100\n"
        "A\tCOX1\t200\n"
        "A\tCYTB\t300\n"
    )
    de_table = tmp_path / "de.tsv"
    de_table.write_text(
        "gene_id\tbaseMean\tlog2FoldChange\tpadj\n"
        "MT-ND1\t1000\t0.1\t0.5\n"
        "MT-CO1\t2000\t0.2\t0.4\n"
        "MT-CYB\t3000\t0.3\t0.3\n"
    )
    out_dir = tmp_path / "rnaseq_out"

    exit_code = cli.main([
        "rnaseq",
        "--de-table", str(de_table),
        "--gene-id-convention", "hgnc",
        "--ribo-dir", str(ribo_dir),
        "--reference-gtf", str(ref),
        "--output", str(out_dir),
    ])

    assert exit_code == 0
    err = capsys.readouterr().err
    assert "MT-CO1->COX1" in err
    te_text = (out_dir / "te.tsv").read_text()
    assert "\tCOX1\t" in te_text
    assert "\tCYTB\t" in te_text
    assert "\tMT-CO1\t" not in te_text
    settings = json.loads((out_dir / "run_settings.json").read_text())
    assert settings["de_to_ribo_gene_map"]["de_to_ribo"]["MT-CYB"] == "CYTB"


def test_rnaseq_reference_mismatch_hard_fails(tmp_path, capsys) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_bytes(b">MT-ND1\nACGTACGT\n")
    ribo_dir = _make_ribo_dir(tmp_path, reference_hash="a" * 64)  # wrong hash
    de_table = _make_de_table(tmp_path)
    out_dir = tmp_path / "rnaseq_out"

    exit_code = cli.main([
        "rnaseq",
        "--de-table", str(de_table),
        "--gene-id-convention", "hgnc",
        "--ribo-dir", str(ribo_dir),
        "--reference-gtf", str(ref),
        "--output", str(out_dir),
    ])

    assert exit_code == 2
    err = capsys.readouterr().err
    assert "MISMATCH" in err
    assert not (out_dir / "te.tsv").exists()


def test_rnaseq_warns_on_missing_mt_mrnas_but_still_runs(tmp_path, capsys) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_bytes(b">MT-ND1\nACGTACGT\n")
    ref_hash = compute_reference_checksum(ref)

    ribo_dir = _make_ribo_dir(tmp_path, ref_hash)
    # DE table only has 2 out of 13 expected mt-mRNAs
    de_table = _make_de_table(tmp_path)
    out_dir = tmp_path / "rnaseq_out"

    exit_code = cli.main([
        "rnaseq",
        "--de-table", str(de_table),
        "--gene-id-convention", "hgnc",
        "--ribo-dir", str(ribo_dir),
        "--reference-gtf", str(ref),
        "--output", str(out_dir),
    ])
    assert exit_code == 0
    err = capsys.readouterr().err
    assert "WARNING" in err
    assert "not found" in err
