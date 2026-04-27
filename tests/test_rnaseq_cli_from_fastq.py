"""CLI-level guards for the rnaseq from-FASTQ mode (mutual exclusion + missing args)."""

from __future__ import annotations

from pathlib import Path

from mitoribopy import cli


def test_de_table_and_rna_fastq_are_mutually_exclusive(capsys, tmp_path: Path) -> None:
    fq = tmp_path / "x.fq.gz"
    fq.write_bytes(b"")
    de = tmp_path / "de.tsv"
    de.write_text("gene_id\tlog2FoldChange\tpadj\tbaseMean\n")
    exit_code = cli.main([
        "rnaseq",
        "--de-table", str(de),
        "--rna-fastq", str(fq),
    ])
    assert exit_code == 2
    err = capsys.readouterr().err
    assert "--de-table" in err
    assert "--rna-fastq" in err


def test_missing_required_from_fastq_args_lists_all_missing(capsys, tmp_path: Path) -> None:
    fq = tmp_path / "x.fq.gz"
    fq.write_bytes(b"")
    exit_code = cli.main(["rnaseq", "--rna-fastq", str(fq)])
    assert exit_code == 2
    err = capsys.readouterr().err
    for flag in (
        "--reference-fasta",
        "--gene-id-convention",
        "--output",
        "--condition-map",
        "--condition-a",
        "--condition-b",
    ):
        assert flag in err, f"expected {flag!r} in error message; got:\n{err}"
