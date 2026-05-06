"""End-to-end tests for ``mitoribopy validate-figures``."""

from __future__ import annotations

import struct
import zlib
from pathlib import Path

import pytest

from mitoribopy.cli import validate_figures as cli_vf
from mitoribopy.io import warnings_log
from mitoribopy.plotting.figure_validator import write_plot_metadata


def _make_png(path: Path, *, dpi: int) -> None:
    px_per_metre = int(round(dpi / 0.0254))
    sig = b"\x89PNG\r\n\x1a\n"

    def chunk(t: bytes, d: bytes) -> bytes:
        return (
            struct.pack(">I", len(d))
            + t + d
            + struct.pack(">I", zlib.crc32(t + d) & 0xFFFFFFFF)
        )

    path.write_bytes(
        sig
        + chunk(b"IHDR", struct.pack(">IIBBBBB", 1, 1, 8, 0, 0, 0, 0))
        + chunk(b"pHYs", struct.pack(">IIB", px_per_metre, px_per_metre, 1))
        + chunk(b"IDAT", zlib.compress(b"\x00\xff"))
        + chunk(b"IEND", b"")
    )


@pytest.fixture(autouse=True)
def _isolate_warnings():
    warnings_log.clear()
    yield
    warnings_log.clear()


class TestValidateFiguresCli:
    def test_clean_run_exits_zero(self, tmp_path: Path) -> None:
        plot = tmp_path / "rnaseq" / "v.png"
        plot.parent.mkdir(parents=True, exist_ok=True)
        _make_png(plot, dpi=300)
        write_plot_metadata(
            plot, plot_type="de_volcano", stage="rnaseq",
            source_data="rnaseq/de_table.tsv",
            n_points_expected=3, n_points_drawn=3,
            label_overlap_count=0,
        )
        rc = cli_vf.run([str(tmp_path)])
        assert rc == 0
        assert (tmp_path / "figure_qc.tsv").is_file()

    def test_warn_only_exits_one(self, tmp_path: Path) -> None:
        plot = tmp_path / "rnaseq" / "v.png"
        plot.parent.mkdir(parents=True, exist_ok=True)
        # Low DPI is a warn in default mode.
        _make_png(plot, dpi=100)
        write_plot_metadata(
            plot,
            plot_type="de_volcano",
            stage="rnaseq",
            n_points_expected=3,
            n_points_drawn=3,
        )
        rc = cli_vf.run([str(tmp_path)])
        assert rc == 1

    def test_fail_exits_two(self, tmp_path: Path) -> None:
        plot = tmp_path / "rnaseq" / "v.png"
        plot.parent.mkdir(parents=True, exist_ok=True)
        _make_png(plot, dpi=300)
        write_plot_metadata(
            plot, plot_type="de_volcano", stage="rnaseq",
            n_points_expected=10, n_points_drawn=9,
        )
        rc = cli_vf.run([str(tmp_path)])
        assert rc == 2

    def test_strict_upgrades_warn(self, tmp_path: Path) -> None:
        plot = tmp_path / "rnaseq" / "v.png"
        plot.parent.mkdir(parents=True, exist_ok=True)
        _make_png(plot, dpi=100)  # below 300
        write_plot_metadata(
            plot, plot_type="de_volcano", stage="rnaseq",
            n_points_expected=3, n_points_drawn=3,
        )
        rc_default = cli_vf.run([str(tmp_path)])
        assert rc_default == 1
        rc_strict = cli_vf.run([str(tmp_path), "--strict"])
        assert rc_strict == 2

    def test_warnings_are_mirrored_into_warnings_log(self, tmp_path: Path) -> None:
        plot = tmp_path / "rnaseq" / "v.png"
        plot.parent.mkdir(parents=True, exist_ok=True)
        _make_png(plot, dpi=300)
        write_plot_metadata(
            plot, plot_type="de_volcano", stage="rnaseq",
            n_points_expected=10, n_points_drawn=9,
        )
        cli_vf.run([str(tmp_path)])
        records = warnings_log.collected()
        assert any(
            r.stage == "VALIDATE_FIGURES" and r.code == "FIGURE_QC_FAIL"
            for r in records
        )

    def test_missing_run_root_exits_two(self, tmp_path: Path) -> None:
        rc = cli_vf.run([str(tmp_path / "no_such_dir")])
        assert rc == 2

    def test_custom_out_path(self, tmp_path: Path) -> None:
        plot = tmp_path / "rnaseq" / "v.png"
        plot.parent.mkdir(parents=True, exist_ok=True)
        _make_png(plot, dpi=300)
        write_plot_metadata(
            plot, plot_type="de_volcano", stage="rnaseq",
            n_points_expected=3, n_points_drawn=3,
        )
        custom = tmp_path / "elsewhere" / "fqc.tsv"
        cli_vf.run([str(tmp_path), "--out", str(custom)])
        assert custom.is_file()


class TestSubcommandRegistration:
    def test_validate_figures_is_registered(self) -> None:
        from mitoribopy.cli import _SUBCOMMANDS, _SUBCOMMAND_SUMMARIES

        assert "validate-figures" in _SUBCOMMANDS
        names = [name for name, _ in _SUBCOMMAND_SUMMARIES]
        assert "validate-figures" in names
        # The README claim is "seven utility subcommands"; verify the
        # CLI has at least seven non-pipeline subcommands so the README
        # never drifts back to a stale count.
        utilities = {
            "migrate-config", "validate-config", "validate-reference",
            "validate-figures", "summarize", "benchmark",
        }
        assert utilities.issubset(set(_SUBCOMMANDS))
