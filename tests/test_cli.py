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
    assert "Which coordinate space the best offset is chosen in" in help_text
    assert "sum only rows whose reference matches" in help_text
    assert "--mt_mrna_substring_patterns" in help_text
    assert "--structure_density" in help_text
    assert "--annotation_file" in help_text
    assert "--codon_table_name" in help_text
    assert "--start_codons" in help_text
    assert "--unfiltered_read_length_range" in help_text
    assert "--atp8_atp6_baseline" in help_text
    assert "--varna" not in help_text
    assert "[default: current working directory]" in help_text
    # -rpf's default_display is footprint-class-aware.
    assert "monosome h.sapiens: 28-34, s.cerevisiae: 37-41" in help_text
    assert "disome   h.sapiens: 50-70, s.cerevisiae: 60-90" in help_text
    assert "short    h.sapiens / s.cerevisiae: 16-24" in help_text
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


def test_strain_h_alias_canonicalises_to_h_sapiens(capsys) -> None:
    """The deprecated short form -s h is rewritten to h.sapiens with one
    DEPRECATED line on stderr; the rest of the pipeline only sees the
    canonical name."""
    args = runner.parse_pipeline_args(
        ["-s", "h", "-f", "tiny.fa", "--directory", "beds"]
    )
    assert args.strain == "h.sapiens"
    err = capsys.readouterr().err
    assert "DEPRECATED" in err
    assert "h -> h.sapiens" in err


def test_strain_y_alias_canonicalises_to_s_cerevisiae(capsys) -> None:
    args = runner.parse_pipeline_args(
        ["-s", "y", "-f", "tiny.fa", "--directory", "beds"]
    )
    assert args.strain == "s.cerevisiae"
    err = capsys.readouterr().err
    assert "DEPRECATED" in err
    assert "y -> s.cerevisiae" in err


def test_canonical_strain_does_not_warn(capsys) -> None:
    args = runner.parse_pipeline_args(
        ["-s", "h.sapiens", "-f", "tiny.fa", "--directory", "beds"]
    )
    assert args.strain == "h.sapiens"
    err = capsys.readouterr().err
    assert "DEPRECATED" not in err


def test_footprint_class_disome_widens_unfiltered_default(tmp_path) -> None:
    """Default footprint_class=monosome keeps [15, 50]; disome widens to [40, 100]
    unless the user overrides --unfiltered_read_length_range."""
    args = runner.parse_pipeline_args(
        ["-s", "h.sapiens", "-f", "tiny.fa", "--directory", "beds",
         "--footprint_class", "disome"]
    )
    assert list(args.unfiltered_read_length_range) == [40, 100]


def test_footprint_class_short_uses_truncated_window(tmp_path) -> None:
    """short footprint class targets RNase truncation products (16-24 nt)
    and widens the unfiltered QC window down to [10, 30] so the short
    tail is visible in the read-length heatmap."""
    args = runner.parse_pipeline_args(
        ["-s", "h.sapiens", "-f", "tiny.fa", "--directory", "beds",
         "--footprint_class", "short"]
    )
    assert list(args.unfiltered_read_length_range) == [10, 30]
    assert args.rpf is None  # rpf is resolved lazily downstream from footprint_class


def test_footprint_class_user_override_wins() -> None:
    """An explicit --unfiltered_read_length_range must not be overwritten by the
    footprint_class default injection."""
    args = runner.parse_pipeline_args(
        ["-s", "h.sapiens", "-f", "tiny.fa", "--directory", "beds",
         "--footprint_class", "disome",
         "--unfiltered_read_length_range", "50", "120"]
    )
    assert list(args.unfiltered_read_length_range) == [50, 120]


def test_resolve_rpf_range_disome_default_differs_from_monosome() -> None:
    """Check the RPF defaults actually differ by footprint class."""
    from mitoribopy.config.runtime import resolve_rpf_range

    mono = resolve_rpf_range("h.sapiens", None, footprint_class="monosome")
    di = resolve_rpf_range("h.sapiens", None, footprint_class="disome")
    short = resolve_rpf_range("h.sapiens", None, footprint_class="short")
    assert mono == list(range(28, 35))
    assert di == list(range(50, 71))
    assert short == list(range(16, 25))


def test_resolve_rpf_range_yeast_disome_window() -> None:
    from mitoribopy.config.runtime import resolve_rpf_range

    assert resolve_rpf_range("s.cerevisiae", None, footprint_class="disome") == list(
        range(60, 91)
    )


def test_resolve_rpf_range_strain_alias_resolves_via_canonical() -> None:
    """The deprecated short strain aliases must still resolve through
    the same FOOTPRINT_CLASS_DEFAULTS table as the canonical names."""
    from mitoribopy.config.runtime import resolve_rpf_range

    assert resolve_rpf_range("h", None, footprint_class="monosome") == list(
        range(28, 35)
    )
    assert resolve_rpf_range("y", None, footprint_class="monosome") == list(
        range(37, 42)
    )


def test_resolve_rpf_range_custom_strain_without_rpf_raises() -> None:
    """`custom` strain without -rpf should raise a clear ValueError
    (the CLI layer catches this and surfaces it as a parser error)."""
    from mitoribopy.config.runtime import resolve_rpf_range

    with pytest.raises(ValueError) as exc:
        resolve_rpf_range("custom", None, footprint_class="monosome")
    assert "Cannot default the RPF range" in str(exc.value)
