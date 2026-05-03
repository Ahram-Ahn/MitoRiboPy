"""Unit tests for ``mitoribopy all`` (Phase 6 orchestrator)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from mitoribopy import cli
from mitoribopy.cli import all_ as all_cli


# ---------- dict -> argv serializer -----------------------------------------


def test_dict_to_argv_handles_scalars_bools_and_lists() -> None:
    argv = all_cli._dict_to_argv(
        {
            "adapter": "TGGAATTCTCGGGTGCCAAGG",
            "umi_length": 0,
            "resume": True,
            "skip_align": False,
            "mt_mrna_substring_patterns": ["mt_genome", "mt-mrna"],
            "optional": None,
        }
    )
    # Flag-name transform: underscores -> dashes.
    assert "--adapter" in argv and argv[argv.index("--adapter") + 1] == "TGGAATTCTCGGGTGCCAAGG"
    assert "--umi-length" in argv and argv[argv.index("--umi-length") + 1] == "0"
    # Bool True -> bare flag
    assert "--resume" in argv
    # Bool False -> omitted entirely
    assert "--skip-align" not in argv
    # None values -> omitted
    assert "--optional" not in argv
    # Lists -> flag followed by space-separated values
    idx = argv.index("--mt-mrna-substring-patterns")
    assert argv[idx + 1] == "mt_genome"
    assert argv[idx + 2] == "mt-mrna"


def test_dict_to_argv_underscore_style_preserves_underscores() -> None:
    """rpf's argparse declares --offset_type etc. with underscores, so the
    orchestrator must NOT hyphen-convert when serializing the rpf section."""
    argv = all_cli._dict_to_argv(
        {
            "offset_type": "5",
            "min_5_offset": 10,
            "mt_mrna_substring_patterns": ["mt_genome", "mt-mrna"],
            "codon_density_window": True,
        },
        flag_style="underscore",
    )
    assert "--offset_type" in argv and "--offset-type" not in argv
    assert "--min_5_offset" in argv and "--min-5-offset" not in argv
    assert "--mt_mrna_substring_patterns" in argv
    assert "--codon_density_window" in argv


def test_dict_to_argv_rejects_unknown_flag_style() -> None:
    with pytest.raises(ValueError):
        all_cli._dict_to_argv({"x": 1}, flag_style="weird")


def test_dict_to_argv_round_trip_preserves_list_and_bool() -> None:
    """Serializing and re-parsing preserves list-valued and boolean flags."""
    import argparse

    section = {
        "rpf": [29, 34],
        "mt_mrna_substring_patterns": ["mt_genome", "mt-mrna"],
        "codon_density_window": True,
        "structure_density": False,
    }
    argv = all_cli._dict_to_argv(section, flag_style="hyphen")
    parser = argparse.ArgumentParser()
    parser.add_argument("--rpf", nargs="+", type=int)
    parser.add_argument("--mt-mrna-substring-patterns", nargs="+")
    parser.add_argument("--codon-density-window", action="store_true")
    parser.add_argument("--structure-density", action="store_true")
    ns = parser.parse_args(argv)
    assert ns.rpf == [29, 34]
    assert ns.mt_mrna_substring_patterns == ["mt_genome", "mt-mrna"]
    assert ns.codon_density_window is True
    assert ns.structure_density is False  # bool False was dropped entirely


def test_dict_to_argv_repeats_append_style_flags() -> None:
    argv = all_cli._dict_to_argv(
        {"fastq": ["sample_a.fq.gz", "sample_b.fq.gz"]},
        repeat_flags={"fastq"},
    )
    assert argv == [
        "--fastq",
        "sample_a.fq.gz",
        "--fastq",
        "sample_b.fq.gz",
    ]


def test_normalize_align_inputs_promotes_string_to_fastq_dir() -> None:
    cfg = {"fastq": "input_data/", "dedup_strategy": "auto"}
    out = all_cli._normalize_align_inputs(cfg)
    assert out["fastq_dir"] == "input_data/"
    # Original `fastq` is rewritten to None so _dict_to_argv drops it.
    assert out["fastq"] is None
    # Caller's dict is not mutated.
    assert cfg["fastq"] == "input_data/"


def test_normalize_align_inputs_rejects_legacy_kit_preset() -> None:
    """v0.7.1: align.kit_preset is a hard error with a migration message."""
    cfg = {"fastq": "input_data/", "kit_preset": "auto"}
    with pytest.raises(ValueError, match="kit_preset"):
        all_cli._normalize_align_inputs(cfg)


def test_normalize_align_inputs_keeps_explicit_list() -> None:
    cfg = {"fastq": ["a.fq.gz", "b.fq.gz"]}
    out = all_cli._normalize_align_inputs(cfg)
    assert out["fastq"] == ["a.fq.gz", "b.fq.gz"]
    assert "fastq_dir" not in out


def test_normalize_align_inputs_does_not_clobber_existing_fastq_dir() -> None:
    cfg = {"fastq": "should_be_ignored/", "fastq_dir": "explicit_dir/"}
    out = all_cli._normalize_align_inputs(cfg)
    assert out["fastq_dir"] == "explicit_dir/"
    assert out["fastq"] is None


def test_all_dry_run_polymorphic_fastq_string_becomes_fastq_dir(tmp_path, capsys) -> None:
    cfg = tmp_path / "p.yaml"
    cfg.write_text(
        "align:\n"
        "  dedup_strategy: auto\n"
        "  fastq: input_data/\n"
    )
    exit_code = cli.main([
        "all",
        "--dry-run",
        "--config", str(cfg),
        "--output", str(tmp_path / "out"),
    ])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "--fastq-dir input_data/" in out
    # The polymorphic string MUST NOT appear as a --fastq value.
    assert "--fastq input_data/" not in out


def test_dict_to_argv_supports_legacy_flag_overrides() -> None:
    argv = all_cli._dict_to_argv(
        {"rpf": [29, 34], "plot_format": "svg"},
        flag_style="underscore",
        flag_overrides={"rpf": "-rpf"},
    )
    assert argv == ["-rpf", "29", "34", "--plot_format", "svg"]


# ---------- CLI help / dry-run ---------------------------------------------


def test_all_help_mentions_three_stages(capsys) -> None:
    with pytest.raises(SystemExit) as exc:
        cli.main(["all", "--help"])
    assert exc.value.code == 0
    out = capsys.readouterr().out
    assert "mitoribopy all" in out
    assert "--resume" in out
    assert "--skip-align" in out
    assert "--skip-rpf" in out
    assert "--skip-rnaseq" in out
    assert "--manifest" in out
    assert "--show-stage-help" in out


def test_all_dry_run_without_config_prints_plan_and_exits_zero(capsys) -> None:
    exit_code = cli.main(["all", "--dry-run"])
    assert exit_code == 0
    assert "dry-run" in capsys.readouterr().out


def test_all_dry_run_with_config_lists_per_stage_argv(tmp_path, capsys) -> None:
    cfg = tmp_path / "pipeline.yaml"
    cfg.write_text(
        "align:\n"
        "  adapter: TGGAATTCTCGGGTGCCAAGG\n"
        "  fastq_dir: fastqs/\n"
        "rpf:\n"
        "  strain: h\n"
        "rnaseq:\n"
        "  de_table: de.tsv\n"
        "  gene_id_convention: hgnc\n"
    )

    exit_code = cli.main([
        "all",
        "--dry-run",
        "--config", str(cfg),
        "--output", str(tmp_path / "results"),
    ])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "align: " in out
    assert "--adapter TGGAATTCTCGGGTGCCAAGG" in out
    assert "rpf: " in out
    assert "rnaseq: " in out
    assert "--gene-id-convention hgnc" in out
    # As of v0.6.0 the rpf parser accepts hyphenated flags as the
    # canonical public form, so the orchestrator emits --read-counts-file.
    # The underscore form (--read_counts_file) is still parsable as a
    # legacy alias but is no longer surfaced through the orchestrator.
    assert "--read-counts-file" in out


def test_all_show_stage_help_prints_full_stage_help(capsys) -> None:
    exit_code = cli.main(["all", "--show-stage-help", "rpf"])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "Run the MitoRiboPy Ribo-seq analysis stage" in out
    # Hyphenated form is canonical as of v0.6.0; verify the help surfaces it.
    assert "--unfiltered-read-length-range" in out


def test_all_print_config_template_exits_zero_with_all_three_sections(capsys) -> None:
    exit_code = cli.main(["all", "--print-config-template"])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "align:" in out
    assert "rpf:" in out
    assert "rnaseq:" in out.lower()  # rnaseq is commented out, but the header appears
    # Template is valid YAML (the active parts, at least).
    import yaml

    yaml_text = out
    parsed = yaml.safe_load(yaml_text)
    assert "align" in parsed and "rpf" in parsed


def test_all_dry_run_rpf_section_uses_hyphenated_flags(tmp_path, capsys) -> None:
    """As of v0.6.0 the rpf parser accepts hyphenated flags as the
    canonical public form, so the dry-run plan must emit them. The
    underscored aliases still parse but are no longer surfaced through
    the orchestrator."""
    cfg = tmp_path / "p.yaml"
    cfg.write_text(
        "rpf:\n"
        "  strain: h\n"
        "  offset_type: '5'\n"
        "  offset_site: p\n"
        "  min_5_offset: 10\n"
        "  codon_density_window: true\n"
    )
    exit_code = cli.main([
        "all",
        "--dry-run",
        "--config", str(cfg),
        "--output", str(tmp_path / "out"),
    ])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "--offset-type 5" in out
    assert "--offset_type 5" not in out
    assert "-rpf" not in out  # absent because this config did not set rpf
    assert "--min-5-offset 10" in out
    assert "--codon-density-window" in out


# ---------- required args ---------------------------------------------------


def test_all_requires_config(capsys) -> None:
    exit_code = cli.main(["all", "--output", "/tmp/x"])
    assert exit_code == 2
    assert "--config" in capsys.readouterr().err


def test_all_requires_output_when_not_dry_run(tmp_path, capsys) -> None:
    cfg = tmp_path / "c.yaml"
    cfg.write_text("align:\n  adapter: TGGAATTCTCGGGTGCCAAGG\n")
    exit_code = cli.main(["all", "--config", str(cfg)])
    assert exit_code == 2
    assert "--output" in capsys.readouterr().err


# ---------- orchestration end-to-end (mocked subcommands) -------------------


def test_all_runs_align_then_rpf_then_rnaseq(tmp_path, monkeypatch) -> None:
    calls: list[str] = []

    def fake_align_run(argv):
        calls.append("align")
        # Simulate align writing its output files.
        out = Path(tmp_path / "results" / "align")
        out.mkdir(parents=True, exist_ok=True)
        (out / "read_counts.tsv").write_text("sample\ttotal_reads\nA\t100\n")
        (out / "run_settings.json").write_text('{"subcommand":"align"}')
        return 0

    def fake_rpf_run(argv):
        calls.append("rpf")
        out = Path(tmp_path / "results" / "rpf")
        out.mkdir(parents=True, exist_ok=True)
        (out / "rpf_counts.tsv").write_text("sample\tgene\tcount\nA\tMT-ND1\t50\n")
        (out / "run_settings.json").write_text(
            '{"subcommand":"rpf","reference_checksum":"abc123"}'
        )
        return 0

    def fake_rnaseq_run(argv):
        calls.append("rnaseq")
        out = Path(tmp_path / "results" / "rnaseq")
        out.mkdir(parents=True, exist_ok=True)
        (out / "delta_te.tsv").write_text("gene\tmrna_log2fc\nMT-ND1\t0.5\n")
        (out / "run_settings.json").write_text('{"subcommand":"rnaseq"}')
        return 0

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rnaseq as rnaseq_cli
    from mitoribopy.cli import rpf as rpf_cli

    monkeypatch.setattr(align_cli, "run", fake_align_run)
    monkeypatch.setattr(rpf_cli, "run", fake_rpf_run)
    monkeypatch.setattr(rnaseq_cli, "run", fake_rnaseq_run)

    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "align:\n  adapter: TGGAATTCTCGGGTGCCAAGG\n"
        "rpf:\n  strain: h\n"
        "rnaseq:\n  de_table: de.tsv\n  gene_id_convention: hgnc\n"
    )

    exit_code = cli.main([
        "all",
        "--config", str(cfg),
        "--output", str(tmp_path / "results"),
    ])
    assert exit_code == 0
    assert calls == ["align", "rpf", "rnaseq"]

    manifest = json.loads((tmp_path / "results" / "run_manifest.json").read_text())
    assert manifest["subcommand"] == "all"
    # Schema-v1 stages block:
    from mitoribopy.cli.all_ import MANIFEST_SCHEMA_VERSION as _MV
    assert manifest["schema_version"] == _MV
    assert manifest["stages"]["align"]["status"] == "completed"
    assert manifest["stages"]["rpf"]["status"] == "completed"
    assert manifest["stages"]["rnaseq"]["status"] == "completed"
    # Per-stage runtime is recorded for completed stages.
    assert "runtime_seconds" in manifest["stages"]["align"]
    assert manifest["rpf"]["reference_checksum"] == "abc123"
    # Reference checksum promoted to top level for easy consumption.
    assert manifest["reference_checksum"] == "abc123"


def test_all_resume_skips_stages_with_existing_outputs(tmp_path, monkeypatch) -> None:
    calls: list[str] = []

    def align_should_not_run(argv):
        calls.append("align")
        return 0

    def rpf_ran(argv):
        calls.append("rpf")
        out = Path(tmp_path / "results" / "rpf")
        out.mkdir(parents=True, exist_ok=True)
        (out / "rpf_counts.tsv").write_text("sample\tgene\tcount\n")
        return 0

    # Pre-create align's sentinel output so --resume skips it.
    results = tmp_path / "results"
    (results / "align").mkdir(parents=True, exist_ok=True)
    (results / "align" / "read_counts.tsv").write_text("sample\n")

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rpf as rpf_cli
    monkeypatch.setattr(align_cli, "run", align_should_not_run)
    monkeypatch.setattr(rpf_cli, "run", rpf_ran)

    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "align:\n  adapter: TGGAATTCTCGGGTGCCAAGG\n"
        "rpf:\n  strain: h\n"
    )

    exit_code = cli.main([
        "all",
        "--config", str(cfg),
        "--output", str(results),
        "--resume",
    ])
    assert exit_code == 0
    assert "align" not in calls  # skipped by --resume
    assert "rpf" in calls

    manifest = json.loads((results / "run_manifest.json").read_text())
    assert manifest["stages"]["align"]["status"] == "skipped"
    assert "resume" in manifest["stages"]["align"]["reason"]
    assert manifest["stages"]["rpf"]["status"] == "completed"


def test_all_propagates_nonzero_exit_from_align(tmp_path, monkeypatch, capsys) -> None:
    from mitoribopy.cli import align as align_cli
    monkeypatch.setattr(align_cli, "run", lambda argv: 2)

    cfg = tmp_path / "c.yaml"
    cfg.write_text("align:\n  adapter: TGGAATTCTCGGGTGCCAAGG\n")
    exit_code = cli.main([
        "all",
        "--config", str(cfg),
        "--output", str(tmp_path / "results"),
    ])
    assert exit_code == 2
    assert "align stage failed" in capsys.readouterr().err


def test_all_skips_rnaseq_when_de_table_absent(tmp_path, monkeypatch) -> None:
    calls: list[str] = []

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rnaseq as rnaseq_cli
    from mitoribopy.cli import rpf as rpf_cli

    def stub_align(argv):
        calls.append("align")
        return 0

    def stub_rpf(argv):
        calls.append("rpf")
        return 0

    def stub_rnaseq(argv):  # should not be called
        calls.append("rnaseq")
        return 0

    monkeypatch.setattr(align_cli, "run", stub_align)
    monkeypatch.setattr(rpf_cli, "run", stub_rpf)
    monkeypatch.setattr(rnaseq_cli, "run", stub_rnaseq)

    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "align:\n  adapter: TGGAATTCTCGGGTGCCAAGG\n"
        "rpf:\n  strain: h\n"
        # rnaseq section is present but has no de_table -> skipped.
        "rnaseq:\n  gene_id_convention: hgnc\n"
    )
    results = tmp_path / "results"
    exit_code = cli.main([
        "all",
        "--config", str(cfg),
        "--output", str(results),
    ])
    assert exit_code == 0
    assert "rnaseq" not in calls
    manifest = json.loads((results / "run_manifest.json").read_text())
    assert manifest["stages"]["rnaseq"]["status"] == "skipped"


def test_all_respects_skip_rpf_flag(tmp_path, monkeypatch) -> None:
    calls: list[str] = []

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rpf as rpf_cli

    monkeypatch.setattr(align_cli, "run", lambda argv: calls.append("align") or 0)
    monkeypatch.setattr(rpf_cli, "run", lambda argv: calls.append("rpf") or 0)

    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "align:\n  adapter: TGGAATTCTCGGGTGCCAAGG\n"
        "rpf:\n  strain: h\n"
    )
    exit_code = cli.main([
        "all",
        "--config", str(cfg),
        "--output", str(tmp_path / "results"),
        "--skip-rpf",
    ])
    assert exit_code == 0
    assert "align" in calls
    assert "rpf" not in calls


def test_all_runs_rnaseq_in_from_fastq_mode(tmp_path, monkeypatch) -> None:
    """rnaseq must execute when the section configures the from-FASTQ
    flow (rna_fastq), even without a de_table."""
    calls: list[str] = []
    captured_argv: dict[str, list[str]] = {}

    def fake_align_run(argv):
        calls.append("align")
        out = Path(tmp_path / "results" / "align")
        out.mkdir(parents=True, exist_ok=True)
        (out / "read_counts.tsv").write_text("sample\n")
        (out / "run_settings.json").write_text('{"subcommand":"align"}')
        return 0

    def fake_rpf_run(argv):
        calls.append("rpf")
        out = Path(tmp_path / "results" / "rpf")
        out.mkdir(parents=True, exist_ok=True)
        (out / "rpf_counts.tsv").write_text("sample\tgene\tcount\n")
        (out / "run_settings.json").write_text(
            '{"subcommand":"rpf","reference_checksum":"abc"}'
        )
        return 0

    def fake_rnaseq_run(argv):
        argv_list = list(argv)
        # The orchestrator now runs the rnaseq stage TWICE in from-FASTQ
        # mode: once with --align-only (parallel to the Ribo align stage,
        # to produce rna_counts.tsv) and once at the end for DE/TE.
        # Differentiate by inspecting the argv.
        if "--align-only" in argv_list:
            calls.append("rnaseq_prealign")
            captured_argv["rnaseq_prealign"] = argv_list
            out = Path(tmp_path / "results" / "rnaseq")
            out.mkdir(parents=True, exist_ok=True)
            (out / "rna_counts.tsv").write_text("gene\tS1\n")
        else:
            calls.append("rnaseq")
            captured_argv["rnaseq"] = argv_list
            out = Path(tmp_path / "results" / "rnaseq")
            out.mkdir(parents=True, exist_ok=True)
            (out / "delta_te.tsv").write_text("gene\n")
            (out / "run_settings.json").write_text('{"subcommand":"rnaseq"}')
        return 0

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rnaseq as rnaseq_cli
    from mitoribopy.cli import rpf as rpf_cli

    monkeypatch.setattr(align_cli, "run", fake_align_run)
    monkeypatch.setattr(rpf_cli, "run", fake_rpf_run)
    monkeypatch.setattr(rnaseq_cli, "run", fake_rnaseq_run)

    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "align:\n  adapter: TGGAATTCTCGGGTGCCAAGG\n"
        "rpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
        "rnaseq:\n"
        "  rna_fastq:\n"
        "    - rna1.fq.gz\n"
        "    - rna2.fq.gz\n"
        "  ribo_fastq:\n"
        "    - ribo1.fq.gz\n"
        "  gene_id_convention: hgnc\n"
        "  condition_map: cond.tsv\n"
    )
    results = tmp_path / "results"
    exit_code = cli.main([
        "all",
        "--config", str(cfg),
        "--output", str(results),
    ])
    assert exit_code == 0
    # The new ordering is:
    #   1. rnaseq --align-only (kicked off in parallel with align)
    #   2. align (Ribo)
    #   3. rpf
    #   4. rnaseq DE/TE
    # Step 1 may finish before or after step 2 depending on scheduling;
    # step 1's stage_start always precedes step 2's though, and the
    # orchestrator joins on step 1 before step 3 starts.
    assert "rnaseq_prealign" in calls
    # rpf must come AFTER both align and rnaseq_prealign.
    assert calls.index("rpf") > calls.index("align")
    assert calls.index("rpf") > calls.index("rnaseq_prealign")
    # Final DE rnaseq comes last.
    assert calls[-1] == "rnaseq"

    rnaseq_argv = captured_argv["rnaseq"]
    # rnaseq.output auto-wired from --output.
    assert "--output" in rnaseq_argv
    assert rnaseq_argv[rnaseq_argv.index("--output") + 1] == str(results / "rnaseq")
    # reference_fasta auto-wired from rpf.fasta when not explicitly set.
    assert "--reference-fasta" in rnaseq_argv
    assert rnaseq_argv[rnaseq_argv.index("--reference-fasta") + 1] == "/tmp/tx.fa"
    # ribo-dir is NOT wired in from-FASTQ mode (uses --ribo-fastq instead).
    assert "--ribo-dir" not in rnaseq_argv
    # The DE invocation gets --upstream-rna-counts pointing at the
    # matrix the align-only worker just wrote.
    assert "--upstream-rna-counts" in rnaseq_argv

    # The align-only invocation must NOT carry the upstream-rna-counts
    # flag (otherwise it would short-circuit and produce no counts).
    align_only_argv = captured_argv["rnaseq_prealign"]
    assert "--align-only" in align_only_argv
    assert "--upstream-rna-counts" not in align_only_argv

    manifest = json.loads((results / "run_manifest.json").read_text())
    assert manifest["stages"]["rnaseq"]["status"] == "completed"


def test_all_rejects_rnaseq_with_both_de_table_and_rna_fastq(
    tmp_path, monkeypatch, capsys
) -> None:
    """The two rnaseq flows are mutually exclusive; fail fast at the
    orchestrator instead of letting the rnaseq CLI error mid-run."""
    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "rnaseq:\n"
        "  de_table: de.tsv\n"
        "  rna_fastq:\n    - rna1.fq.gz\n"
        "  gene_id_convention: hgnc\n"
    )
    exit_code = cli.main([
        "all",
        "--config", str(cfg),
        "--output", str(tmp_path / "results"),
    ])
    assert exit_code == 2
    err = capsys.readouterr().err
    assert "mutually exclusive" in err


def test_manifest_carries_schema_v1_metadata(tmp_path, monkeypatch) -> None:
    """run_manifest.json v1.0.0 must record schema_version, command,
    config_source + sha256, config_canonical, git_commit, tools, and
    a stages map with per-stage status + runtime."""
    def fake_align(argv):
        out = Path(tmp_path / "results" / "align")
        out.mkdir(parents=True, exist_ok=True)
        (out / "read_counts.tsv").write_text("sample\n")
        (out / "run_settings.json").write_text(
            '{"subcommand":"align","cutadapt_version":"4.9"}'
        )
        return 0

    def fake_rpf(argv):
        out = Path(tmp_path / "results" / "rpf")
        out.mkdir(parents=True, exist_ok=True)
        (out / "rpf_counts.tsv").write_text("sample\tgene\tcount\n")
        (out / "run_settings.json").write_text(
            '{"subcommand":"rpf","reference_checksum":"abc123"}'
        )
        return 0

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rpf as rpf_cli
    monkeypatch.setattr(align_cli, "run", fake_align)
    monkeypatch.setattr(rpf_cli, "run", fake_rpf)

    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "align:\n  adapter: TGGAATTCTCGGGTGCCAAGG\n"
        "rpf:\n  strain: h\n"
    )
    results = tmp_path / "results"
    exit_code = cli.main([
        "all", "--config", str(cfg), "--output", str(results),
    ])
    assert exit_code == 0

    manifest = json.loads((results / "run_manifest.json").read_text())
    from mitoribopy.cli.all_ import MANIFEST_SCHEMA_VERSION as _MV
    assert manifest["schema_version"] == _MV
    assert manifest["command"].startswith("mitoribopy all --config")
    assert manifest["config_source"] == str(cfg)
    # Hash of the YAML the user wrote (stable for unchanged inputs).
    assert manifest["config_source_sha256"] is not None
    assert len(manifest["config_source_sha256"]) == 64
    # No sample sheet was set; both fields stay null but present.
    assert manifest["sample_sheet"] is None
    assert manifest["sample_sheet_sha256"] is None
    # The merged + auto-wired config is embedded so a reviewer can see
    # what was actually executed.
    assert manifest["config_canonical"]["align"]["output"] == str(results / "align")
    assert manifest["config_canonical"]["rpf"]["output"] == str(results / "rpf")
    # Stages reshape (align + rpf completed; rnaseq not configured).
    assert manifest["stages"]["align"]["status"] == "completed"
    assert manifest["stages"]["rpf"]["status"] == "completed"
    assert manifest["stages"]["rnaseq"]["status"] == "not_configured"
    assert "runtime_seconds" in manifest["stages"]["align"]
    # tools map lifts cutadapt_version from align's run_settings.json
    # and always carries python + mitoribopy.
    tools = manifest["tools"]
    assert tools["cutadapt"] == "4.9"
    assert "python" in tools and "mitoribopy" in tools
    # Reference checksum is still promoted to the top level.
    assert manifest["reference_checksum"] == "abc123"
    # Warnings is a placeholder list.
    assert manifest["warnings"] == []


def test_manifest_distinguishes_skipped_vs_not_configured(
    tmp_path, monkeypatch
) -> None:
    """A stage actively skipped (--skip-<stage>) shows status='skipped'
    with a reason; a stage with no config section shows
    status='not_configured' and carries no reason."""
    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rpf as rpf_cli

    def fake_align(argv):
        out = Path(tmp_path / "results" / "align")
        out.mkdir(parents=True, exist_ok=True)
        (out / "read_counts.tsv").write_text("sample\n")
        (out / "run_settings.json").write_text("{}")
        return 0

    def fake_rpf(argv):
        # Should NOT run because we --skip-rpf below.
        raise AssertionError("rpf must be skipped")

    monkeypatch.setattr(align_cli, "run", fake_align)
    monkeypatch.setattr(rpf_cli, "run", fake_rpf)

    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "align:\n  adapter: TGGAATTCTCGGGTGCCAAGG\n"
        "rpf:\n  strain: h\n"
        # rnaseq section absent -> not_configured.
    )
    results = tmp_path / "results"
    exit_code = cli.main([
        "all", "--config", str(cfg), "--output", str(results),
        "--skip-rpf",
    ])
    assert exit_code == 0

    manifest = json.loads((results / "run_manifest.json").read_text())
    # rpf was actively skipped via --skip-rpf -> skipped + reason.
    assert manifest["stages"]["rpf"]["status"] == "skipped"
    assert manifest["stages"]["rpf"]["reason"] == "--skip-rpf flag set"
    # rnaseq has no config section -> not_configured + no reason field.
    assert manifest["stages"]["rnaseq"]["status"] == "not_configured"
    assert "reason" not in manifest["stages"]["rnaseq"]


def test_manifest_records_sample_sheet_sha256(tmp_path, monkeypatch) -> None:
    """When a top-level sample sheet drives the run, both its path and
    its SHA256 land in the manifest so a reviewer can verify the inputs
    have not changed since the run."""
    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rpf as rpf_cli
    monkeypatch.setattr(
        align_cli, "run",
        lambda argv: (
            Path(tmp_path / "results" / "align").mkdir(parents=True, exist_ok=True),
            (tmp_path / "results" / "align" / "read_counts.tsv").write_text("x"),
            (tmp_path / "results" / "align" / "run_settings.json").write_text("{}"),
        ) and 0,
    )
    monkeypatch.setattr(
        rpf_cli, "run",
        lambda argv: (
            Path(tmp_path / "results" / "rpf").mkdir(parents=True, exist_ok=True),
            (tmp_path / "results" / "rpf" / "rpf_counts.tsv").write_text("x"),
            (tmp_path / "results" / "rpf" / "run_settings.json").write_text("{}"),
        ) and 0,
    )

    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\tassay\tcondition\tfastq_1\n"
        "WT_Ribo_1\tribo\tWT\tribo/WT.fq.gz\n"
        "WT_RNA_1\trna\tWT\trna/WT.fq.gz\n"
    )
    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        f"samples:\n  table: {sheet}\n"
        "align:\n  dedup_strategy: auto\n"
        "rpf:\n  strain: h\n"
    )
    results = tmp_path / "results"
    exit_code = cli.main([
        "all", "--config", str(cfg), "--output", str(results),
    ])
    assert exit_code == 0

    manifest = json.loads((results / "run_manifest.json").read_text())
    assert manifest["sample_sheet"] == str(sheet)
    assert manifest["sample_sheet_sha256"] is not None
    assert len(manifest["sample_sheet_sha256"]) == 64


def test_print_canonical_config_outputs_merged_yaml(tmp_path, capsys) -> None:
    """--print-canonical-config must emit the same blob the manifest
    would embed under config_canonical, without running any stage."""
    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "align:\n  adapter: TGGAATTCTCGGGTGCCAAGG\n"
        "rpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
    )
    exit_code = cli.main([
        "all",
        "--config", str(cfg),
        "--output", str(tmp_path / "results"),
        "--print-canonical-config",
    ])
    assert exit_code == 0
    out = capsys.readouterr().out
    # The auto-wired stage outputs must appear in the printed canonical
    # form (proves --print-canonical-config ran the auto-wiring step,
    # not just dumped the input verbatim).
    assert "align" in out and "rpf" in out
    assert str(tmp_path / "results" / "align") in out
    assert str(tmp_path / "results" / "rpf") in out


def test_print_canonical_config_requires_config_path(tmp_path, capsys) -> None:
    exit_code = cli.main([
        "all", "--print-canonical-config",
    ])
    assert exit_code == 2
    err = capsys.readouterr().err
    assert "--config" in err


def test_all_top_level_samples_drives_align_and_rnaseq(
    tmp_path, monkeypatch
) -> None:
    """A top-level `samples:` block in the YAML should auto-wire
    align.fastq + align.sample_overrides + rnaseq.sample_sheet from the
    unified sheet, replacing the need for per-stage input flags."""
    captured: dict[str, list[str]] = {}

    def fake_align(argv):
        captured["align"] = list(argv)
        out = Path(tmp_path / "results" / "align")
        out.mkdir(parents=True, exist_ok=True)
        (out / "read_counts.tsv").write_text("sample\n")
        (out / "run_settings.json").write_text('{"subcommand":"align"}')
        return 0

    def fake_rpf(argv):
        captured["rpf"] = list(argv)
        out = Path(tmp_path / "results" / "rpf")
        out.mkdir(parents=True, exist_ok=True)
        (out / "rpf_counts.tsv").write_text("sample\tgene\tcount\n")
        (out / "run_settings.json").write_text('{"subcommand":"rpf"}')
        return 0

    def fake_rnaseq(argv):
        captured["rnaseq"] = list(argv)
        out = Path(tmp_path / "results" / "rnaseq")
        out.mkdir(parents=True, exist_ok=True)
        (out / "delta_te.tsv").write_text("gene\n")
        (out / "run_settings.json").write_text('{"subcommand":"rnaseq"}')
        return 0

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rnaseq as rnaseq_cli
    from mitoribopy.cli import rpf as rpf_cli
    monkeypatch.setattr(align_cli, "run", fake_align)
    monkeypatch.setattr(rpf_cli, "run", fake_rpf)
    monkeypatch.setattr(rnaseq_cli, "run", fake_rnaseq)

    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\tassay\tcondition\tfastq_1\tadapter\tumi_length\n"
        "WT_Ribo_1\tribo\tWT\tribo/WT.fq.gz\tAGATCGGAAGAGCACACGTCTGAACTCCAGTCA\t8\n"
        "KO_Ribo_1\tribo\tKO\tribo/KO.fq.gz\tAGATCGGAAGAGCACACGTCTGAACTCCAGTCA\t8\n"
        "WT_RNA_1\trna\tWT\trna/WT.fq.gz\t\t\n"
        "KO_RNA_1\trna\tKO\trna/KO.fq.gz\t\t\n"
    )
    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        f"samples:\n  table: {sheet}\n"
        "align:\n  dedup_strategy: auto\n"
        "rpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
        "rnaseq:\n  gene_id_convention: bare\n"
        "  condition_a: WT\n  condition_b: KO\n"
    )
    results = tmp_path / "results"
    exit_code = cli.main([
        "all", "--config", str(cfg), "--output", str(results),
    ])
    assert exit_code == 0

    align_argv = captured["align"]
    # The align serializer uses repeat_flags={"fastq"}, so --fastq is
    # emitted once per FASTQ. Pull the repeated values back out.
    fastqs = [
        align_argv[i + 1]
        for i, tok in enumerate(align_argv)
        if tok == "--fastq" and i + 1 < len(align_argv)
    ]
    assert fastqs == ["ribo/WT.fq.gz", "ribo/KO.fq.gz"]
    # Per-sample overrides materialised under <run_root>/align/.
    assert "--sample-overrides" in align_argv
    overrides_path = Path(
        align_argv[align_argv.index("--sample-overrides") + 1]
    )
    assert overrides_path.exists()
    body = overrides_path.read_text()
    assert "WT_Ribo_1\tAGATCGGAAGAGCACACGTCTGAACTCCAGTCA\t" in body
    assert "umi_length" in body.splitlines()[0]

    rnaseq_argv = captured["rnaseq"]
    # The sheet path threads through to the rnaseq stage.
    assert "--sample-sheet" in rnaseq_argv
    assert rnaseq_argv[rnaseq_argv.index("--sample-sheet") + 1] == str(sheet)


def test_all_top_level_samples_rejects_align_fastq_conflict(
    tmp_path, capsys
) -> None:
    """Declaring `samples:` AND `align.fastq:` at the same time must
    fail loudly so the user does not get a silent shadow."""
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\tassay\tcondition\tfastq_1\n"
        "WT_Ribo_1\tribo\tWT\tribo/WT.fq.gz\n"
        "WT_RNA_1\trna\tWT\trna/WT.fq.gz\n"
    )
    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        f"samples:\n  table: {sheet}\n"
        "align:\n  fastq_dir: input/\n"
    )
    exit_code = cli.main([
        "all", "--config", str(cfg), "--output", str(tmp_path / "out"),
    ])
    assert exit_code == 2
    err = capsys.readouterr().err
    assert "samples" in err and "align.fastq" in err


def test_all_top_level_samples_shorthand_string_form(tmp_path, monkeypatch) -> None:
    """`samples: path/to/sheet.tsv` as a bare string should also work."""
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\tassay\tcondition\tfastq_1\n"
        "WT_Ribo_1\tribo\tWT\tribo/WT.fq.gz\n"
        "WT_RNA_1\trna\tWT\trna/WT.fq.gz\n"
    )
    captured: dict[str, list[str]] = {}

    def fake_align(argv):
        captured["align"] = list(argv)
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

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rpf as rpf_cli
    monkeypatch.setattr(align_cli, "run", fake_align)
    monkeypatch.setattr(rpf_cli, "run", fake_rpf)

    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        f"samples: {sheet}\n"
        "align:\n  dedup_strategy: auto\n"
        "rpf:\n  strain: h\n"
    )
    exit_code = cli.main([
        "all", "--config", str(cfg), "--output", str(tmp_path / "results"),
    ])
    assert exit_code == 0
    align_argv = captured["align"]
    assert "--fastq" in align_argv
    # repeated --fastq emission (one per file) — only one Ribo row in this sheet.
    fastqs = [
        align_argv[i + 1]
        for i, tok in enumerate(align_argv)
        if tok == "--fastq" and i + 1 < len(align_argv)
    ]
    assert fastqs == ["ribo/WT.fq.gz"]


def test_all_keeps_explicit_rnaseq_reference_override(tmp_path, monkeypatch) -> None:
    """User-specified rnaseq.reference_fasta must NOT be clobbered by
    rpf.fasta auto-wiring (RNA-seq may use a different transcriptome)."""
    captured_argv: dict[str, list[str]] = {}

    def fake_align_run(argv):
        out = Path(tmp_path / "results" / "align")
        out.mkdir(parents=True, exist_ok=True)
        (out / "read_counts.tsv").write_text("sample\n")
        (out / "run_settings.json").write_text('{"subcommand":"align"}')
        return 0

    def fake_rpf_run(argv):
        out = Path(tmp_path / "results" / "rpf")
        out.mkdir(parents=True, exist_ok=True)
        (out / "rpf_counts.tsv").write_text("sample\tgene\tcount\n")
        (out / "run_settings.json").write_text('{"subcommand":"rpf"}')
        return 0

    def fake_rnaseq_run(argv):
        captured_argv["rnaseq"] = list(argv)
        out = Path(tmp_path / "results" / "rnaseq")
        out.mkdir(parents=True, exist_ok=True)
        (out / "delta_te.tsv").write_text("gene\n")
        (out / "run_settings.json").write_text('{"subcommand":"rnaseq"}')
        return 0

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rnaseq as rnaseq_cli
    from mitoribopy.cli import rpf as rpf_cli

    monkeypatch.setattr(align_cli, "run", fake_align_run)
    monkeypatch.setattr(rpf_cli, "run", fake_rpf_run)
    monkeypatch.setattr(rnaseq_cli, "run", fake_rnaseq_run)

    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "align:\n  adapter: TGGAATTCTCGGGTGCCAAGG\n"
        "rpf:\n  strain: h\n  fasta: /tmp/rpf.fa\n"
        "rnaseq:\n"
        "  rna_fastq:\n    - rna1.fq.gz\n"
        "  reference_fasta: /tmp/rnaseq_specific.fa\n"
        "  gene_id_convention: hgnc\n"
    )
    exit_code = cli.main([
        "all",
        "--config", str(cfg),
        "--output", str(tmp_path / "results"),
    ])
    assert exit_code == 0
    rnaseq_argv = captured_argv["rnaseq"]
    assert "--reference-fasta" in rnaseq_argv
    idx = rnaseq_argv.index("--reference-fasta")
    # rpf.fasta did NOT clobber the explicit rnaseq.reference_fasta.
    assert rnaseq_argv[idx + 1] == "/tmp/rnaseq_specific.fa"
