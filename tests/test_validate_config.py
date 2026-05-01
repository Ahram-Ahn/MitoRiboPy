"""Tests for the `mitoribopy validate-config` subcommand (P1.9)."""

from __future__ import annotations

from pathlib import Path

import pytest

from mitoribopy import cli


def _write(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    return path


# ---------- happy paths -----------------------------------------------------


def test_validate_config_ok_for_minimal_canonical_config(
    tmp_path: Path, capsys
) -> None:
    cfg = _write(
        tmp_path / "c.yaml",
        "rpf:\n  strain: h.sapiens\n  fasta: /tmp/tx.fa\n",
    )
    rc = cli.main(["validate-config", str(cfg), "--no-path-checks"])
    assert rc == 0
    err = capsys.readouterr().err
    assert "OK" in err


def test_validate_config_warns_on_legacy_keys_but_passes_by_default(
    tmp_path: Path, capsys
) -> None:
    cfg = _write(
        tmp_path / "c.yaml",
        "rpf:\n  strain: h\n  merge_density: true\n  fasta: /tmp/tx.fa\n",
    )
    rc = cli.main(["validate-config", str(cfg), "--no-path-checks"])
    assert rc == 0
    err = capsys.readouterr().err
    assert "WARNING" in err
    assert "merge_density" in err
    assert "OK" in err


def test_validate_config_strict_mode_fails_on_legacy_keys(
    tmp_path: Path, capsys
) -> None:
    cfg = _write(
        tmp_path / "c.yaml",
        "rpf:\n  strain: h\n  fasta: /tmp/tx.fa\n",
    )
    rc = cli.main([
        "validate-config", str(cfg), "--no-path-checks", "--strict"
    ])
    assert rc == 2
    err = capsys.readouterr().err
    assert "WARNING" in err


# ---------- structural errors -----------------------------------------------


def test_validate_config_unknown_top_level_key(
    tmp_path: Path, capsys
) -> None:
    cfg = _write(
        tmp_path / "c.yaml",
        "mystery: 42\nrpf:\n  strain: h.sapiens\n  fasta: /tmp/tx.fa\n",
    )
    rc = cli.main(["validate-config", str(cfg), "--no-path-checks"])
    assert rc == 2
    err = capsys.readouterr().err
    assert "unknown top-level key" in err
    assert "mystery" in err


def test_validate_config_mutually_exclusive_rnaseq_inputs(
    tmp_path: Path, capsys
) -> None:
    cfg = _write(
        tmp_path / "c.yaml",
        "rnaseq:\n"
        "  de_table: de.tsv\n"
        "  rna_fastq:\n    - rna1.fq.gz\n",
    )
    rc = cli.main(["validate-config", str(cfg), "--no-path-checks"])
    assert rc == 2
    err = capsys.readouterr().err
    assert "mutually exclusive" in err


def test_validate_config_sample_sheet_conflict(
    tmp_path: Path, capsys
) -> None:
    sheet = _write(
        tmp_path / "samples.tsv",
        "sample_id\tassay\tcondition\tfastq_1\n"
        "WT\tribo\tWT\twt.fq.gz\n",
    )
    cfg = _write(
        tmp_path / "c.yaml",
        f"samples:\n  table: {sheet}\n"
        "align:\n  fastq_dir: input/\n",
    )
    rc = cli.main(["validate-config", str(cfg), "--no-path-checks"])
    assert rc == 2
    err = capsys.readouterr().err
    assert "samples" in err
    assert "fastq_dir" in err


def test_validate_config_explicit_mode_input_mismatch(
    tmp_path: Path, capsys
) -> None:
    cfg = _write(
        tmp_path / "c.yaml",
        "rnaseq:\n"
        "  rnaseq_mode: de_table\n"
        "  rna_fastq:\n    - rna1.fq.gz\n",
    )
    rc = cli.main(["validate-config", str(cfg), "--no-path-checks"])
    assert rc == 2
    err = capsys.readouterr().err
    assert "rnaseq" in err
    assert "de_table" in err


# ---------- path checks -----------------------------------------------------


def test_validate_config_missing_fasta_path_errors(
    tmp_path: Path, capsys
) -> None:
    cfg = _write(
        tmp_path / "c.yaml",
        f"rpf:\n  strain: h.sapiens\n  fasta: {tmp_path}/missing.fa\n",
    )
    rc = cli.main(["validate-config", str(cfg)])
    assert rc == 2
    err = capsys.readouterr().err
    assert "missing.fa" in err
    assert "rpf.fasta" in err


def test_validate_config_present_fasta_path_passes(
    tmp_path: Path, capsys
) -> None:
    fasta = _write(tmp_path / "tx.fa", ">x\nACGT\n")
    cfg = _write(
        tmp_path / "c.yaml",
        f"rpf:\n  strain: h.sapiens\n  fasta: {fasta}\n",
    )
    rc = cli.main(["validate-config", str(cfg)])
    assert rc == 0


def test_validate_config_no_path_checks_flag_skips_existence(
    tmp_path: Path, capsys
) -> None:
    cfg = _write(
        tmp_path / "c.yaml",
        "rpf:\n  strain: h.sapiens\n  fasta: /nope/missing.fa\n",
    )
    rc = cli.main(["validate-config", str(cfg), "--no-path-checks"])
    assert rc == 0


def test_validate_config_bowtie2_index_prefix_missing_sidecar(
    tmp_path: Path, capsys
) -> None:
    cfg = _write(
        tmp_path / "c.yaml",
        f"align:\n  contam_index: {tmp_path}/contam\n  mt_index: {tmp_path}/mt\n",
    )
    rc = cli.main(["validate-config", str(cfg)])
    assert rc == 2
    err = capsys.readouterr().err
    assert ".1.bt2" in err


def test_validate_config_bowtie2_index_prefix_with_sidecar_passes(
    tmp_path: Path, capsys
) -> None:
    (tmp_path / "contam.1.bt2").write_bytes(b"")
    (tmp_path / "mt.1.bt2").write_bytes(b"")
    cfg = _write(
        tmp_path / "c.yaml",
        f"align:\n  contam_index: {tmp_path}/contam\n  mt_index: {tmp_path}/mt\n",
    )
    rc = cli.main(["validate-config", str(cfg)])
    assert rc == 0


# ---------- sample sheet validation -----------------------------------------


def test_validate_config_invalid_sample_sheet_surfaces_loader_error(
    tmp_path: Path, capsys
) -> None:
    sheet = _write(
        tmp_path / "samples.tsv",
        "sample_id\tassay\n"  # missing required columns
        "WT\tribo\n",
    )
    cfg = _write(
        tmp_path / "c.yaml",
        f"samples:\n  table: {sheet}\n",
    )
    rc = cli.main(["validate-config", str(cfg)])
    assert rc == 2
    err = capsys.readouterr().err
    assert "samples" in err
    assert "missing required column" in err


# ---------- top-level help --------------------------------------------------


def test_validate_config_appears_in_top_level_help(capsys) -> None:
    rc = cli.main(["--help"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "validate-config" in out
