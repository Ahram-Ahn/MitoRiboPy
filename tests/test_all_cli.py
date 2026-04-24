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
            "kit_preset": "truseq_smallrna",
            "umi_length": 0,
            "resume": True,
            "skip_align": False,
            "mrna_ref_patterns": ["mt_genome", "mt-mrna"],
            "optional": None,
        }
    )
    # Flag-name transform: underscores -> dashes.
    assert "--kit-preset" in argv and argv[argv.index("--kit-preset") + 1] == "truseq_smallrna"
    assert "--umi-length" in argv and argv[argv.index("--umi-length") + 1] == "0"
    # Bool True -> bare flag
    assert "--resume" in argv
    # Bool False -> omitted entirely
    assert "--skip-align" not in argv
    # None values -> omitted
    assert "--optional" not in argv
    # Lists -> flag followed by space-separated values
    idx = argv.index("--mrna-ref-patterns")
    assert argv[idx + 1] == "mt_genome"
    assert argv[idx + 2] == "mt-mrna"


def test_dict_to_argv_underscore_style_preserves_underscores() -> None:
    """rpf's argparse declares --offset_type etc. with underscores, so the
    orchestrator must NOT hyphen-convert when serializing the rpf section."""
    argv = all_cli._dict_to_argv(
        {
            "offset_type": "5",
            "min_5_offset": 10,
            "mrna_ref_patterns": ["mt_genome", "mt-mrna"],
            "merge_density": True,
        },
        flag_style="underscore",
    )
    assert "--offset_type" in argv and "--offset-type" not in argv
    assert "--min_5_offset" in argv and "--min-5-offset" not in argv
    assert "--mrna_ref_patterns" in argv
    assert "--merge_density" in argv


def test_dict_to_argv_rejects_unknown_flag_style() -> None:
    with pytest.raises(ValueError):
        all_cli._dict_to_argv({"x": 1}, flag_style="weird")


def test_dict_to_argv_round_trip_preserves_list_and_bool() -> None:
    """Serializing and re-parsing preserves list-valued and boolean flags."""
    import argparse

    section = {
        "rpf": [29, 34],
        "mrna_ref_patterns": ["mt_genome", "mt-mrna"],
        "merge_density": True,
        "structure_density": False,
    }
    argv = all_cli._dict_to_argv(section, flag_style="hyphen")
    parser = argparse.ArgumentParser()
    parser.add_argument("--rpf", nargs="+", type=int)
    parser.add_argument("--mrna-ref-patterns", nargs="+")
    parser.add_argument("--merge-density", action="store_true")
    parser.add_argument("--structure-density", action="store_true")
    ns = parser.parse_args(argv)
    assert ns.rpf == [29, 34]
    assert ns.mrna_ref_patterns == ["mt_genome", "mt-mrna"]
    assert ns.merge_density is True
    assert ns.structure_density is False  # bool False was dropped entirely


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
        "  kit_preset: truseq_smallrna\n"
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
    assert "--kit-preset truseq_smallrna" in out
    assert "rpf: " in out
    assert "rnaseq: " in out
    assert "--gene-id-convention hgnc" in out
    # rpf uses underscored flag style because pipeline.runner declares
    # --read_counts_file that way. The orchestrator must emit it verbatim.
    assert "--read_counts_file" in out


def test_all_show_stage_help_prints_full_stage_help(capsys) -> None:
    exit_code = cli.main(["all", "--show-stage-help", "rpf"])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "Run the standalone MitoRiboPy Ribo-seq pipeline." in out
    assert "--unfiltered_read_length_range" in out


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


def test_all_dry_run_rpf_section_uses_underscored_flags(tmp_path, capsys) -> None:
    """Regression: the rpf parser declares --offset_type with underscores,
    so the dry-run plan must emit underscored flags for rpf."""
    cfg = tmp_path / "p.yaml"
    cfg.write_text(
        "rpf:\n"
        "  strain: h\n"
        "  offset_type: '5'\n"
        "  offset_site: p\n"
        "  min_5_offset: 10\n"
        "  merge_density: true\n"
    )
    exit_code = cli.main([
        "all",
        "--dry-run",
        "--config", str(cfg),
        "--output", str(tmp_path / "out"),
    ])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "--offset_type 5" in out
    assert "--offset-type" not in out
    assert "--min_5_offset 10" in out
    assert "--merge_density" in out


# ---------- required args ---------------------------------------------------


def test_all_requires_config(capsys) -> None:
    exit_code = cli.main(["all", "--output", "/tmp/x"])
    assert exit_code == 2
    assert "--config" in capsys.readouterr().err


def test_all_requires_output_when_not_dry_run(tmp_path, capsys) -> None:
    cfg = tmp_path / "c.yaml"
    cfg.write_text("align:\n  kit_preset: truseq_smallrna\n")
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
        "align:\n  kit_preset: truseq_smallrna\n"
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
    assert manifest["stages_run"] == ["align", "rpf", "rnaseq"]
    assert manifest["stages_skipped"] == []
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
        "align:\n  kit_preset: truseq_smallrna\n"
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
    assert "align" in manifest["stages_skipped"]
    assert "rpf" in manifest["stages_run"]


def test_all_propagates_nonzero_exit_from_align(tmp_path, monkeypatch, capsys) -> None:
    from mitoribopy.cli import align as align_cli
    monkeypatch.setattr(align_cli, "run", lambda argv: 2)

    cfg = tmp_path / "c.yaml"
    cfg.write_text("align:\n  kit_preset: truseq_smallrna\n")
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
        "align:\n  kit_preset: truseq_smallrna\n"
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
    assert "rnaseq" in manifest["stages_skipped"]


def test_all_respects_skip_rpf_flag(tmp_path, monkeypatch) -> None:
    calls: list[str] = []

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rpf as rpf_cli

    monkeypatch.setattr(align_cli, "run", lambda argv: calls.append("align") or 0)
    monkeypatch.setattr(rpf_cli, "run", lambda argv: calls.append("rpf") or 0)

    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "align:\n  kit_preset: truseq_smallrna\n"
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
