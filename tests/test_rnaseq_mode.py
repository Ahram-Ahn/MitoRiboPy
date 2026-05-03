"""Tests for the explicit ``--rnaseq-mode`` flag and resolver.

P0.3 of the publication-readiness refactor: the rnaseq stage now has
an explicit ``mode`` of ``de_table | from_fastq | none``. The resolver
must:

* infer the mode from supplied inputs when the flag is omitted
  (back-compat with v0.5.x configs);
* reject ambiguous configs that supply inputs for both flows;
* reject explicit modes that conflict with the supplied inputs;
* surface the resolved mode in ``run_settings.json`` and in the
  ``--print-canonical-config`` output of ``mitoribopy all``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pytest

from mitoribopy import cli
from mitoribopy.cli import all_ as all_cli
from mitoribopy.cli.rnaseq import (
    RNASEQ_MODES,
    _normalize_mode,
    _resolve_rnaseq_mode,
)


def _ns(**kwargs) -> argparse.Namespace:
    """Build a Namespace covering every attribute the resolver reads."""
    defaults = dict(
        rnaseq_mode=None,
        de_table=None,
        rna_fastq=None,
        sample_sheet=None,
    )
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


# ---------- _normalize_mode -------------------------------------------------


def test_normalize_mode_handles_hyphen_and_case() -> None:
    assert _normalize_mode("from-fastq") == "from_fastq"
    assert _normalize_mode("DE_TABLE") == "de_table"
    assert _normalize_mode(None) is None


# ---------- _resolve_rnaseq_mode (inference) --------------------------------


def test_inferred_mode_is_de_table_when_only_de_table_present() -> None:
    mode, err = _resolve_rnaseq_mode(_ns(de_table="de.tsv"))
    assert mode == "de_table"
    assert err is None


def test_inferred_mode_is_from_fastq_for_rna_fastq() -> None:
    mode, err = _resolve_rnaseq_mode(_ns(rna_fastq=["a.fq.gz"]))
    assert mode == "from_fastq"
    assert err is None


def test_inferred_mode_is_from_fastq_for_sample_sheet() -> None:
    mode, err = _resolve_rnaseq_mode(_ns(sample_sheet="samples.tsv"))
    assert mode == "from_fastq"
    assert err is None


def test_inferred_mode_is_none_when_no_inputs() -> None:
    mode, err = _resolve_rnaseq_mode(_ns())
    assert mode == "none"
    assert err is None


def test_ambiguous_no_mode_with_both_inputs_errors() -> None:
    """Both flow inputs without an explicit mode is a hard error."""
    mode, err = _resolve_rnaseq_mode(_ns(de_table="de.tsv", rna_fastq=["a.fq.gz"]))
    assert mode is None
    assert err is not None
    assert "mutually exclusive" in err


# ---------- _resolve_rnaseq_mode (explicit) ---------------------------------


def test_explicit_de_table_with_de_table_input_ok() -> None:
    mode, err = _resolve_rnaseq_mode(
        _ns(rnaseq_mode="de_table", de_table="de.tsv")
    )
    assert mode == "de_table"
    assert err is None


def test_explicit_from_fastq_with_rna_fastq_ok() -> None:
    mode, err = _resolve_rnaseq_mode(
        _ns(rnaseq_mode="from_fastq", rna_fastq=["a.fq.gz"])
    )
    assert mode == "from_fastq"
    assert err is None


def test_explicit_from_fastq_accepts_hyphenated_form() -> None:
    """`--rnaseq-mode from-fastq` is normalized to `from_fastq`."""
    mode, err = _resolve_rnaseq_mode(
        _ns(rnaseq_mode="from-fastq", rna_fastq=["a.fq.gz"])
    )
    assert mode == "from_fastq"
    assert err is None


def test_explicit_de_table_with_rna_fastq_conflict_errors() -> None:
    mode, err = _resolve_rnaseq_mode(
        _ns(rnaseq_mode="de_table", rna_fastq=["a.fq.gz"])
    )
    assert mode is None
    assert err is not None
    assert "de_table" in err and "from_fastq" in err


def test_explicit_from_fastq_with_de_table_conflict_errors() -> None:
    mode, err = _resolve_rnaseq_mode(
        _ns(rnaseq_mode="from_fastq", de_table="de.tsv")
    )
    assert mode is None
    assert err is not None


def test_explicit_de_table_without_de_table_input_errors() -> None:
    mode, err = _resolve_rnaseq_mode(_ns(rnaseq_mode="de_table"))
    assert mode is None
    assert err is not None
    assert "--de-table" in err


def test_explicit_from_fastq_without_inputs_errors() -> None:
    mode, err = _resolve_rnaseq_mode(_ns(rnaseq_mode="from_fastq"))
    assert mode is None
    assert err is not None
    assert "--rna-fastq" in err or "--sample-sheet" in err


def test_explicit_unknown_mode_errors() -> None:
    """argparse.choices=... already rejects this at parse time, but the
    resolver is also defence-in-depth for the YAML path which bypasses
    argparse via ``rnaseq.mode``."""
    mode, err = _resolve_rnaseq_mode(
        _ns(rnaseq_mode="bogus", rna_fastq=["a.fq.gz"])
    )
    assert mode is None
    assert err is not None
    assert "bogus" in err


# ---------- _resolve_rnaseq_mode_from_config (YAML path) --------------------


def test_yaml_mode_alias_is_accepted() -> None:
    """A user can write `rnaseq.mode:` instead of `rnaseq.rnaseq_mode:`."""
    mode, err = all_cli._resolve_rnaseq_mode_from_config(
        {"mode": "de_table", "de_table": "de.tsv"}
    )
    assert mode == "de_table"
    assert err is None


def test_yaml_mode_inference_from_de_table() -> None:
    mode, err = all_cli._resolve_rnaseq_mode_from_config({"de_table": "de.tsv"})
    assert mode == "de_table"
    assert err is None


def test_yaml_mode_inference_from_rna_fastq() -> None:
    mode, err = all_cli._resolve_rnaseq_mode_from_config(
        {"rna_fastq": ["a.fq.gz"]}
    )
    assert mode == "from_fastq"
    assert err is None


def test_yaml_mode_ambiguous_inputs_errors() -> None:
    mode, err = all_cli._resolve_rnaseq_mode_from_config(
        {"de_table": "de.tsv", "rna_fastq": ["a.fq.gz"]}
    )
    assert err is not None
    assert "mutually exclusive" in err


def test_yaml_mode_explicit_de_table_with_fastq_input_errors() -> None:
    mode, err = all_cli._resolve_rnaseq_mode_from_config(
        {"rnaseq_mode": "de_table", "rna_fastq": ["a.fq.gz"]}
    )
    assert err is not None
    assert "de_table" in err


# ---------- end-to-end: --print-canonical-config records mode ---------------


def test_print_canonical_config_records_resolved_mode(tmp_path, capsys) -> None:
    """`--print-canonical-config` should pin the resolved mode into the
    rnaseq section so reviewers can see what mode actually drove the run."""
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\tassay\tcondition\tfastq_1\n"
        "WT_RNA_1\trna\tWT\trna/WT.fq.gz\n"
        "WT_Ribo_1\tribo\tWT\tribo/WT.fq.gz\n"
        "KO_RNA_1\trna\tKO\trna/KO.fq.gz\n"
        "KO_Ribo_1\tribo\tKO\tribo/KO.fq.gz\n"
    )
    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        f"samples:\n  table: {sheet}\n"
        "align:\n  # adapter auto-detection (default)\n"
        "rpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
        "rnaseq:\n  gene_id_convention: bare\n"
        "  condition_a: WT\n  condition_b: KO\n"
    )
    exit_code = cli.main([
        "all",
        "--print-canonical-config",
        "--config", str(cfg),
        "--output", str(tmp_path / "results"),
    ])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "rnaseq_mode" in out
    assert "from_fastq" in out


def test_yaml_explicit_mode_de_table_pinned_in_canonical_config(
    tmp_path, capsys
) -> None:
    """When the user writes `mode: de_table` we must echo the canonical
    spelling `rnaseq_mode: de_table` back so config diffs are stable."""
    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "rpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
        "rnaseq:\n"
        "  mode: de_table\n"
        "  de_table: de.tsv\n"
        "  gene_id_convention: hgnc\n"
        "  reference_gtf: ref.fa\n"
    )
    exit_code = cli.main([
        "all",
        "--print-canonical-config",
        "--config", str(cfg),
        "--output", str(tmp_path / "results"),
    ])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "rnaseq_mode: de_table" in out


def test_exploratory_banner_emits_for_from_fastq_runs(
    tmp_path, capsys, monkeypatch
) -> None:
    """The exploratory banner must fire on every from_fastq run so users
    are reminded that the mt-only pyDESeq2 path is not publication-grade."""
    from mitoribopy.cli import rnaseq as rnaseq_cli

    def fake_run_from_fastq(args, output_dir):
        return -1  # short-circuit before any heavy work

    monkeypatch.setattr(rnaseq_cli, "_run_from_fastq", fake_run_from_fastq)

    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\tassay\tcondition\tfastq_1\n"
        "WT_RNA_1\trna\tWT\trna/WT.fq.gz\n"
        "KO_RNA_1\trna\tKO\trna/KO.fq.gz\n"
    )
    exit_code = cli.main([
        "rnaseq",
        "--sample-sheet", str(sheet),
        "--reference-fasta", str(tmp_path / "ref.fa"),
        "--gene-id-convention", "bare",
        "--output", str(tmp_path / "out"),
        "--condition-a", "WT",
        "--condition-b", "KO",
    ])
    assert exit_code == 0
    err = capsys.readouterr().err
    assert "exploratory" in err.lower()
    assert "from_fastq" in err


def test_explicit_rnaseq_mode_argv_passes_through_to_subcommand(
    tmp_path, monkeypatch
) -> None:
    """`mitoribopy all` should serialize the resolved mode as
    `--rnaseq-mode <mode>` in the argv handed to the rnaseq subcommand."""
    captured: dict[str, list[str]] = {}

    def fake_align(argv):
        out = Path(tmp_path / "results" / "align")
        out.mkdir(parents=True, exist_ok=True)
        (out / "read_counts.tsv").write_text("sample\n")
        (out / "run_settings.json").write_text("{}")
        return 0

    def fake_rpf(argv):
        out = Path(tmp_path / "results" / "rpf")
        out.mkdir(parents=True, exist_ok=True)
        (out / "rpf_counts.tsv").write_text("sample\tgene\tcount\n")
        (out / "run_settings.json").write_text("{}")
        return 0

    def fake_rnaseq(argv):
        captured["rnaseq"] = list(argv)
        out = Path(tmp_path / "results" / "rnaseq")
        out.mkdir(parents=True, exist_ok=True)
        (out / "delta_te.tsv").write_text("gene\n")
        (out / "run_settings.json").write_text("{}")
        return 0

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rnaseq as rnaseq_cli
    from mitoribopy.cli import rpf as rpf_cli
    monkeypatch.setattr(align_cli, "run", fake_align)
    monkeypatch.setattr(rpf_cli, "run", fake_rpf)
    monkeypatch.setattr(rnaseq_cli, "run", fake_rnaseq)

    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "align:\n  adapter: TGGAATTCTCGGGTGCCAAGG\n"
        "rpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
        "rnaseq:\n"
        "  rna_fastq:\n    - rna1.fq.gz\n"
        "  reference_fasta: /tmp/ref.fa\n"
        "  gene_id_convention: hgnc\n"
        "  condition_a: WT\n  condition_b: KO\n"
    )
    exit_code = cli.main([
        "all",
        "--config", str(cfg),
        "--output", str(tmp_path / "results"),
    ])
    assert exit_code == 0
    rnaseq_argv = captured["rnaseq"]
    assert "--rnaseq-mode" in rnaseq_argv
    idx = rnaseq_argv.index("--rnaseq-mode")
    assert rnaseq_argv[idx + 1] == "from_fastq"
