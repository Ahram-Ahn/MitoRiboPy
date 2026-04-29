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


def test_yaml_config_values_reach_args(tmp_path: Path, capsys) -> None:
    """`mitoribopy rnaseq --config x.yaml` must fold YAML keys into
    argparse defaults so the user does not have to repeat every flag
    on the CLI. Regression for v0.5.0 where --config was registered
    but never loaded for the rnaseq subcommand.

    We assert that a YAML setting `gene_id_convention: bare` plus
    `condition_a` / `condition_b` is enough to clear those required-
    arg checks (we still trigger an error from the missing
    --rna-fastq / --reference-fasta, but the gene-id / condition flags
    must NOT appear in the missing list).
    """
    rna = tmp_path / "rna.fq.gz"
    rna.write_bytes(b"")
    cmap = tmp_path / "conditions.tsv"
    cmap.write_text("sample\tcondition\nA\tWT\nB\tKO\n")
    yaml_path = tmp_path / "rnaseq.yaml"
    yaml_path.write_text(
        "gene_id_convention: bare\n"
        f"condition_map: {cmap}\n"
        "condition_a: WT\n"
        "condition_b: KO\n"
    )
    # Run with only --config (and one --rna-fastq so we are in default
    # flow). The other YAML-supplied keys should satisfy the missing-
    # args check.
    exit_code = cli.main([
        "rnaseq",
        "--config", str(yaml_path),
        "--rna-fastq", str(rna),
    ])
    assert exit_code == 2
    err = capsys.readouterr().err
    # Still missing because we did not set them anywhere:
    assert "--reference-fasta" in err
    assert "--output" in err
    # MUST NOT be in the missing list — the YAML supplied them:
    assert "--gene-id-convention" not in err
    assert "--condition-map" not in err
    assert "--condition-a" not in err
    assert "--condition-b" not in err


def test_yaml_config_unknown_keys_warn_but_do_not_crash(
    tmp_path: Path, capsys
) -> None:
    yaml_path = tmp_path / "rnaseq.yaml"
    yaml_path.write_text(
        "gene_id_convention: bare\n"
        "this_key_does_not_exist: hello\n"
    )
    exit_code = cli.main(["rnaseq", "--config", str(yaml_path)])
    # The point of the test: unknown keys are tolerated (warned via
    # the package logger), not a hard error. We exit 2 only because
    # other required args are missing.
    assert exit_code == 2


def test_base_sample_alias_satisfies_required_flags(
    tmp_path: Path, capsys
) -> None:
    """``--base-sample`` / ``--compare-sample`` populate ``condition_a``
    / ``condition_b`` so the missing-flag check no longer trips for the
    legacy spellings. Adding ``--rna-fastq`` keeps us in the default
    flow; we expect the run to fail later (no real reference / output)
    but the missing-args block must NOT mention the condition flags
    when their aliases supplied the values."""
    fq = tmp_path / "x.fq.gz"
    fq.write_bytes(b"")
    cmap = tmp_path / "conditions.tsv"
    cmap.write_text("sample\tcondition\nA\tWT\nB\tKO\n")
    exit_code = cli.main([
        "rnaseq",
        "--rna-fastq", str(fq),
        "--gene-id-convention", "bare",
        "--condition-map", str(cmap),
        "--base-sample", "WT",
        "--compare-sample", "KO",
    ])
    assert exit_code == 2
    err = capsys.readouterr().err
    # We still trip on --reference-fasta / --output (those were never
    # provided), but the alias-satisfied condition flags must NOT be
    # listed as missing.
    assert "--reference-fasta" in err
    assert "--output" in err
    assert "--condition-a" not in err
    assert "--condition-b" not in err


def test_base_sample_conflicts_with_condition_a(tmp_path: Path, capsys) -> None:
    """When the alias and the legacy form disagree, exit with code 2 and
    a clear error pointing at both flags."""
    fq = tmp_path / "x.fq.gz"
    fq.write_bytes(b"")
    cmap = tmp_path / "conditions.tsv"
    cmap.write_text("sample\tcondition\nA\tWT\nB\tKO\n")
    exit_code = cli.main([
        "rnaseq",
        "--rna-fastq", str(fq),
        "--gene-id-convention", "bare",
        "--condition-map", str(cmap),
        "--condition-a", "WT",
        "--base-sample", "RESC",
    ])
    assert exit_code == 2
    err = capsys.readouterr().err
    assert "--base-sample" in err
    assert "--condition-a" in err


def test_yaml_base_sample_key_satisfies_condition_a(
    tmp_path: Path, capsys
) -> None:
    """``base_sample:`` in YAML must populate ``condition_a`` (via the
    alias-reconciliation step) so the missing-arg check doesn't fire."""
    rna = tmp_path / "rna.fq.gz"
    rna.write_bytes(b"")
    cmap = tmp_path / "conditions.tsv"
    cmap.write_text("sample\tcondition\nA\tWT\nB\tKO\n")
    yaml_path = tmp_path / "rnaseq.yaml"
    yaml_path.write_text(
        "gene_id_convention: bare\n"
        f"condition_map: {cmap}\n"
        "base_sample: WT\n"
        "compare_sample: KO\n"
    )
    exit_code = cli.main([
        "rnaseq",
        "--config", str(yaml_path),
        "--rna-fastq", str(rna),
    ])
    assert exit_code == 2
    err = capsys.readouterr().err
    assert "--reference-fasta" in err
    assert "--output" in err
    assert "--condition-a" not in err
    assert "--condition-b" not in err
