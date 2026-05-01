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


# ---------- pseudo-replicate gate (Phase 1.2 of the refactor) ----------------


def test_conditions_with_one_sample_helper() -> None:
    from mitoribopy.cli.rnaseq import _conditions_with_one_sample
    from mitoribopy.rnaseq.fastq_pairing import FastqSample
    from pathlib import Path as _P

    samples = [
        FastqSample(sample="WT_a", r1=_P("WT_a.fq.gz")),
        FastqSample(sample="WT_b", r1=_P("WT_b.fq.gz")),
        FastqSample(sample="KO_a", r1=_P("KO_a.fq.gz")),
        FastqSample(sample="UNUSED", r1=_P("UNUSED.fq.gz")),
    ]
    cmap = {"WT_a": "WT", "WT_b": "WT", "KO_a": "KO"}
    # KO has only one sample; WT has two; UNUSED is not in the map.
    assert _conditions_with_one_sample(samples, cmap) == ["KO"]


def test_exploratory_sidecar_calls_out_outputs_to_avoid(tmp_path) -> None:
    from mitoribopy.cli.rnaseq import _write_exploratory_sidecar

    out = tmp_path / "rnaseq"
    out.mkdir()
    _write_exploratory_sidecar(out, ["WT", "KO"])
    text = (out / "EXPLORATORY.md").read_text()
    assert "EXPLORATORY" in text
    assert "WT" in text and "KO" in text
    # Make sure the sidecar names the dangerous outputs explicitly.
    for token in ("padj", "p-value", "significant"):
        assert token in text.lower()


def _patch_fastq_layer(monkeypatch, rna_conds: dict, ribo_conds: dict | None = None):
    """Mock enumerate_fastqs / detect_samples so _run_from_fastq can reach
    the singleton gate without doing any real I/O.

    ``rna_conds`` is a dict mapping ``sample_name -> condition``; one
    FastqSample is fabricated per entry. ``ribo_conds`` does the same for
    Ribo-seq, or ``None`` to leave the Ribo side empty.
    """
    from mitoribopy.cli import rnaseq as rnaseq_cli
    from mitoribopy.rnaseq import fastq_pairing as fp

    samples_rna = [
        fp.FastqSample(sample=name, r1=Path(f"{name}.fq.gz"))
        for name in rna_conds
    ]
    samples_ribo = [
        fp.FastqSample(sample=name, r1=Path(f"{name}.fq.gz"))
        for name in (ribo_conds or {})
    ]

    def fake_enumerate(_inputs):
        return [Path("dummy.fq.gz")]

    def fake_detect(_paths):
        # Distinguish RNA vs Ribo by checking the call order via a flag.
        # We use a list-pop trick: first call returns RNA, second returns Ribo.
        return fake_detect._queue.pop(0)

    fake_detect._queue = [samples_rna, samples_ribo]
    monkeypatch.setattr(fp, "enumerate_fastqs", fake_enumerate)
    monkeypatch.setattr(fp, "detect_samples", fake_detect)
    return rnaseq_cli


def test_pseudo_replicate_gate_fails_fast_without_optin(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """An n=1 condition without the explicit opt-in flag must exit 2
    BEFORE bowtie2 / pyDESeq2 can run."""
    rnaseq_cli = _patch_fastq_layer(
        monkeypatch,
        rna_conds={"WT_a": "WT", "KO_a": "KO"},  # both n=1
    )

    cmap = tmp_path / "conditions.tsv"
    cmap.write_text("sample\tcondition\nWT_a\tWT\nKO_a\tKO\n")

    import argparse as _ap
    args = _ap.Namespace(
        rna_fastq=["dummy.fq.gz"],
        ribo_fastq=None,
        reference_fasta=str(tmp_path / "ref.fa"),
        bowtie2_index=None,
        workdir=str(tmp_path / "work"),
        align_threads=1,
        condition_map=str(cmap),
        condition_a="WT",
        condition_b="KO",
        allow_pseudo_replicates=False,
        no_auto_pseudo_replicate=False,
    )
    rc = rnaseq_cli._run_from_fastq(args, tmp_path / "out")
    assert rc == 2
    err = capsys.readouterr().err
    assert "only 1 sample" in err
    assert "--allow-pseudo-replicates-for-demo-not-publication" in err
    # Both singleton conditions should be reported.
    assert "'WT'" in err and "'KO'" in err


def test_no_auto_pseudo_replicate_is_deprecated_no_op(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """The pre-v0.5.2 opt-out flag is now a no-op that warns and falls
    through to the new safe default (i.e. n=1 still hard-fails)."""
    rnaseq_cli = _patch_fastq_layer(
        monkeypatch,
        rna_conds={"WT_a": "WT", "KO_a": "KO"},
    )

    cmap = tmp_path / "conditions.tsv"
    cmap.write_text("sample\tcondition\nWT_a\tWT\nKO_a\tKO\n")

    import argparse as _ap
    args = _ap.Namespace(
        rna_fastq=["dummy.fq.gz"],
        ribo_fastq=None,
        reference_fasta=str(tmp_path / "ref.fa"),
        bowtie2_index=None,
        workdir=str(tmp_path / "work"),
        align_threads=1,
        condition_map=str(cmap),
        condition_a="WT",
        condition_b="KO",
        allow_pseudo_replicates=False,
        no_auto_pseudo_replicate=True,  # legacy opt-out
    )
    rc = rnaseq_cli._run_from_fastq(args, tmp_path / "out")
    err = capsys.readouterr().err
    # Deprecation warning emitted for the legacy flag…
    assert "DEPRECATED" in err
    # …but the run still hard-fails because n=1 is now the safe default.
    assert rc == 2
    assert "--allow-pseudo-replicates-for-demo-not-publication" in err


def test_allow_pseudo_replicates_flag_is_accepted_by_parser() -> None:
    """The new opt-in flag must be parseable without error, and must set
    args.allow_pseudo_replicates to True."""
    from mitoribopy.cli.rnaseq import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "--rna-fastq", "x.fq.gz",
        "--reference-fasta", "ref.fa",
        "--gene-id-convention", "bare",
        "--output", "out",
        "--condition-map", "cmap.tsv",
        "--condition-a", "WT",
        "--condition-b", "KO",
        "--allow-pseudo-replicates-for-demo-not-publication",
    ])
    assert args.allow_pseudo_replicates is True
