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
    assert "--unfiltered_read_length_range" in help_text
    assert "--atp8_atp6_baseline" in help_text
    assert "--varna" not in help_text
    assert "[default: current working directory]" in help_text
    # -rpf's default_display is footprint-class-aware now.
    assert "[default: monosome h/vm: 28-34, y/ym: 37-41" in help_text
    assert "--footprint_class" in help_text
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


def test_vm_strain_requires_annotation_and_rpf(capsys) -> None:
    """vm (any vertebrate mito) uses the built-in codon table but the
    user still has to pass --annotation_file and an explicit -rpf range."""
    with pytest.raises(SystemExit) as exc:
        runner.parse_pipeline_args(
            ["-s", "vm", "-f", "tiny.fa", "--directory", "beds"]
        )
    assert exc.value.code == 2
    err = capsys.readouterr().err
    assert "--strain vm requires --annotation_file" in err


def test_ym_strain_requires_explicit_rpf(tmp_path, capsys) -> None:
    """ym can pick up the built-in yeast_mitochondrial codon table, but
    without --rpf we cannot guess its RPF range."""
    ann = tmp_path / "ann.csv"
    ann.write_text("transcript,l_tr,l_utr5,l_utr3\nGENE,300,0,0\n")
    with pytest.raises(SystemExit) as exc:
        runner.parse_pipeline_args(
            [
                "-s", "ym",
                "-f", "tiny.fa",
                "--directory", "beds",
                "--annotation_file", str(ann),
            ]
        )
    assert exc.value.code == 2
    err = capsys.readouterr().err
    assert "--strain ym requires an explicit -rpf" in err


def test_footprint_class_disome_widens_unfiltered_default(tmp_path) -> None:
    """Default footprint_class=monosome keeps [15, 50]; disome widens to [40, 110]
    unless the user overrides --unfiltered_read_length_range."""
    args = runner.parse_pipeline_args(
        ["-s", "h", "-f", "tiny.fa", "--directory", "beds",
         "--footprint_class", "disome"]
    )
    assert list(args.unfiltered_read_length_range) == [40, 110]


def test_footprint_class_user_override_wins() -> None:
    """An explicit --unfiltered_read_length_range must not be overwritten by the
    footprint_class default injection."""
    args = runner.parse_pipeline_args(
        ["-s", "h", "-f", "tiny.fa", "--directory", "beds",
         "--footprint_class", "disome",
         "--unfiltered_read_length_range", "50", "120"]
    )
    assert list(args.unfiltered_read_length_range) == [50, 120]


def test_resolve_rpf_range_disome_default_differs_from_monosome() -> None:
    """Check the RPF defaults actually differ by footprint class."""
    from mitoribopy.config.runtime import resolve_rpf_range

    mono = resolve_rpf_range("h", None, footprint_class="monosome")
    di = resolve_rpf_range("h", None, footprint_class="disome")
    assert mono == list(range(28, 35))
    assert di == list(range(60, 91))


def test_resolve_rpf_range_custom_strain_without_rpf_raises() -> None:
    """custom / vm / ym without -rpf should raise a clear ValueError
    (the CLI layer catches this and surfaces it as a parser error)."""
    from mitoribopy.config.runtime import resolve_rpf_range

    with pytest.raises(ValueError) as exc:
        resolve_rpf_range("custom", None, footprint_class="monosome")
    assert "Cannot default RPF range" in str(exc.value)
