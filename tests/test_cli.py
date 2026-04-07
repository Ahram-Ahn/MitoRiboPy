import pytest

import mitoribopy.cli as cli
from mitoribopy.cli import _normalize_args
from mitoribopy.config import DEFAULT_CONFIG
from mitoribopy.pipeline import runner


def test_normalize_args_strips_run_and_separator() -> None:
    assert _normalize_args(["run", "--", "-s", "h"]) == ["-s", "h"]


def test_normalize_args_leaves_standard_argv() -> None:
    assert _normalize_args(["-s", "y", "--align", "start"]) == ["-s", "y", "--align", "start"]


def test_main_delegates_to_package_runner(monkeypatch) -> None:
    captured = {}

    def fake_run_pipeline_cli(argv):
        captured["argv"] = argv
        return 17

    monkeypatch.setattr(cli, "run_pipeline_cli", fake_run_pipeline_cli)

    assert cli.main(["run", "--", "-s", "h", "-f", "tiny.fa"]) == 17
    assert captured["argv"] == ["-s", "h", "-f", "tiny.fa"]


def test_build_parser_help_groups_and_key_guidance() -> None:
    help_text = runner.build_parser(dict(DEFAULT_CONFIG)).format_help()

    assert help_text.startswith("usage: mitoribopy")
    assert "Core Inputs:" in help_text
    assert "Offset Enrichment and Selection:" in help_text
    assert "Read-Count Normalization:" in help_text
    assert "Outputs and Plotting:" in help_text
    assert "Optional Modules:" in help_text
    assert "Which read end defines downstream site placement" in help_text
    assert "How the best offset is chosen from enrichment tables" in help_text
    assert "sum only rows whose reference matches --mrna_ref_patterns" in help_text
    assert "--structure_density" in help_text
    assert "--annotation_file" in help_text
    assert "--codon_table_name" in help_text
    assert "--start_codons" in help_text
    assert "--atp8_atp6_baseline" in help_text
    assert "--varna" not in help_text
    assert "[default: current working directory]" in help_text
    assert "[default: y: 37-41, h: 28-34]" in help_text
    assert "[default: same as --min_offset]" in help_text
    assert "[default: total]" in help_text


def test_run_pipeline_cli_help_exits_before_plot_setup(monkeypatch) -> None:
    def fail_configure_plot_defaults() -> None:
        raise AssertionError("configure_plot_defaults should not run for --help")

    monkeypatch.setattr(runner, "configure_plot_defaults", fail_configure_plot_defaults)

    with pytest.raises(SystemExit) as exc:
        runner.run_pipeline_cli(["--help"])

    assert exc.value.code == 0


def test_custom_strain_requires_annotation_and_rpf() -> None:
    with pytest.raises(SystemExit) as exc:
        runner.parse_pipeline_args(
            [
                "-s",
                "custom",
                "-f",
                "tiny.fa",
                "--directory",
                "beds",
                "--codon_table_name",
                "standard",
            ]
        )

    assert exc.value.code == 2
