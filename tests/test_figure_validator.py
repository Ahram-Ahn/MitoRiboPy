"""Tests for ``mitoribopy.plotting.figure_validator`` (refactor-4 / §4 §6)."""

from __future__ import annotations

import json
import struct
import zlib
from pathlib import Path

import pytest

from mitoribopy.plotting.figure_validator import (
    build_figure_qc_rows,
    discover_plots,
    load_plot_metadata,
    metadata_sidecar_path,
    validate_figures,
    write_figure_qc,
    write_plot_metadata,
)


# ---------------------------------------------------------------------------
# Helpers — build minimal but valid PNG / SVG fixtures without matplotlib
# ---------------------------------------------------------------------------


def _png_chunk(chunk_type: bytes, data: bytes) -> bytes:
    return (
        struct.pack(">I", len(data))
        + chunk_type
        + data
        + struct.pack(">I", zlib.crc32(chunk_type + data) & 0xFFFFFFFF)
    )


def _make_png(path: Path, *, dpi: int) -> None:
    """Write a 1×1 white PNG with the requested dpi in its pHYs chunk."""
    px_per_metre = int(round(dpi / 0.0254))
    sig = b"\x89PNG\r\n\x1a\n"
    ihdr = _png_chunk(
        b"IHDR",
        struct.pack(">IIBBBBB", 1, 1, 8, 0, 0, 0, 0),
    )
    phys = _png_chunk(
        b"pHYs",
        struct.pack(">IIB", px_per_metre, px_per_metre, 1),
    )
    raw = b"\x00\xff"  # 1 row, filter byte + one white pixel
    idat = _png_chunk(b"IDAT", zlib.compress(raw))
    iend = _png_chunk(b"IEND", b"")
    path.write_bytes(sig + ihdr + phys + idat + iend)


def _make_svg(path: Path, *, editable_text: bool) -> None:
    body = (
        '<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10">'
        + ('<text x="0" y="10">label</text>' if editable_text else '<path d="M0 0"/>')
        + "</svg>"
    )
    path.write_text(body, encoding="utf-8")


# ---------------------------------------------------------------------------
# write_plot_metadata + load_plot_metadata
# ---------------------------------------------------------------------------


class TestPlotMetadata:
    def test_writes_sibling_json(self, tmp_path: Path) -> None:
        plot = tmp_path / "rnaseq" / "delta_te_volcano.png"
        plot.parent.mkdir(parents=True, exist_ok=True)
        _make_png(plot, dpi=300)
        sidecar = write_plot_metadata(
            plot,
            plot_type="delta_te_volcano",
            stage="rnaseq",
            source_data="rnaseq/delta_te.tsv",
            n_points_expected=13,
            n_points_drawn=13,
            label_overlap_count=0,
        )
        assert sidecar == plot.with_suffix(".metadata.json")
        assert sidecar.is_file()
        data = json.loads(sidecar.read_text())
        assert data["plot_type"] == "delta_te_volcano"
        assert data["stage"] == "rnaseq"
        assert data["n_points_expected"] == 13
        assert data["n_points_drawn"] == 13
        assert data["palette"] == "Okabe-Ito"
        assert data["dpi"] == 300

    def test_load_returns_none_when_missing(self, tmp_path: Path) -> None:
        assert load_plot_metadata(tmp_path / "no_such.png") is None

    def test_metadata_path_is_pure_function(self) -> None:
        p = Path("/run/rnaseq/foo.svg")
        assert metadata_sidecar_path(p) == Path("/run/rnaseq/foo.metadata.json")


# ---------------------------------------------------------------------------
# discover_plots
# ---------------------------------------------------------------------------


class TestDiscoverPlots:
    def test_finds_png_and_svg_recursively(self, tmp_path: Path) -> None:
        a = tmp_path / "rnaseq" / "a.png"
        b = tmp_path / "rpf" / "diag" / "b.svg"
        for p in (a, b):
            p.parent.mkdir(parents=True, exist_ok=True)
        _make_png(a, dpi=300)
        _make_svg(b, editable_text=True)
        found = discover_plots(tmp_path)
        assert a in found
        assert b in found

    def test_dedups_png_svg_pair_by_stem(self, tmp_path: Path) -> None:
        png = tmp_path / "rnaseq" / "v.png"
        svg = png.with_suffix(".svg")
        png.parent.mkdir(parents=True, exist_ok=True)
        _make_png(png, dpi=300)
        _make_svg(svg, editable_text=True)
        found = discover_plots(tmp_path)
        # PNG is preferred when both exist for the same stem.
        assert found == [png]


# ---------------------------------------------------------------------------
# validate_figures — score per plot
# ---------------------------------------------------------------------------


class TestValidateFigures:
    def test_clean_run_is_pass(self, tmp_path: Path) -> None:
        plot = tmp_path / "rnaseq" / "delta_te_volcano.png"
        plot.parent.mkdir(parents=True, exist_ok=True)
        _make_png(plot, dpi=300)
        _make_svg(plot.with_suffix(".svg"), editable_text=True)
        write_plot_metadata(
            plot, plot_type="delta_te_volcano", stage="rnaseq",
            source_data="rnaseq/delta_te.tsv",
            n_points_expected=13, n_points_drawn=13,
            label_overlap_count=0, label_outside_axes_count=0,
            n_labels=13,
        )
        records = validate_figures(tmp_path)
        assert len(records) == 1
        assert records[0].status == "pass"
        assert records[0].plot_type == "delta_te_volcano"
        assert records[0].stage == "rnaseq"

    def test_point_count_mismatch_is_fail(self, tmp_path: Path) -> None:
        plot = tmp_path / "rnaseq" / "v.png"
        plot.parent.mkdir(parents=True, exist_ok=True)
        _make_png(plot, dpi=300)
        write_plot_metadata(
            plot, plot_type="de_volcano", stage="rnaseq",
            n_points_expected=13, n_points_drawn=12,
        )
        records = validate_figures(tmp_path)
        assert records[0].status == "fail"
        assert any("point-count mismatch" in w for w in records[0].warnings)

    def test_label_overlap_is_fail(self, tmp_path: Path) -> None:
        plot = tmp_path / "rnaseq" / "v.png"
        plot.parent.mkdir(parents=True, exist_ok=True)
        _make_png(plot, dpi=300)
        write_plot_metadata(
            plot, plot_type="de_volcano", stage="rnaseq",
            label_overlap_count=2,
        )
        records = validate_figures(tmp_path)
        assert records[0].status == "fail"

    def test_low_dpi_is_warn_in_default_mode(self, tmp_path: Path) -> None:
        plot = tmp_path / "rnaseq" / "v.png"
        plot.parent.mkdir(parents=True, exist_ok=True)
        _make_png(plot, dpi=100)
        write_plot_metadata(
            plot, plot_type="de_volcano", stage="rnaseq",
            n_points_expected=5, n_points_drawn=5,
        )
        records = validate_figures(tmp_path)
        assert records[0].status == "warn"

    def test_strict_upgrades_low_dpi_to_fail(self, tmp_path: Path) -> None:
        plot = tmp_path / "rnaseq" / "v.png"
        plot.parent.mkdir(parents=True, exist_ok=True)
        _make_png(plot, dpi=100)
        write_plot_metadata(
            plot, plot_type="de_volcano", stage="rnaseq",
            n_points_expected=5, n_points_drawn=5,
        )
        records = validate_figures(tmp_path, strict=True)
        assert records[0].status == "fail"

    def test_svg_path_text_is_warn(self, tmp_path: Path) -> None:
        plot = tmp_path / "rnaseq" / "v.png"
        plot.parent.mkdir(parents=True, exist_ok=True)
        _make_png(plot, dpi=300)
        _make_svg(plot.with_suffix(".svg"), editable_text=False)
        write_plot_metadata(
            plot, plot_type="de_volcano", stage="rnaseq",
            n_points_expected=5, n_points_drawn=5,
        )
        records = validate_figures(tmp_path)
        # SVG text non-editable is a warn.
        assert records[0].status == "warn"
        assert records[0].svg_text_editable is False

    def test_missing_sidecar_is_warn(self, tmp_path: Path) -> None:
        plot = tmp_path / "rpf" / "no_sidecar.png"
        plot.parent.mkdir(parents=True, exist_ok=True)
        _make_png(plot, dpi=300)
        records = validate_figures(tmp_path)
        assert records[0].status == "warn"
        assert any("sidecar" in w for w in records[0].warnings)


# ---------------------------------------------------------------------------
# write_figure_qc — TSV contract
# ---------------------------------------------------------------------------


class TestWriteFigureQc:
    def test_header_includes_schema_version(self, tmp_path: Path) -> None:
        plot = tmp_path / "rnaseq" / "v.png"
        plot.parent.mkdir(parents=True, exist_ok=True)
        _make_png(plot, dpi=300)
        write_plot_metadata(
            plot, plot_type="de_volcano", stage="rnaseq",
            n_points_expected=3, n_points_drawn=3,
        )
        records = validate_figures(tmp_path)
        rows = build_figure_qc_rows(records, run_root=tmp_path)
        out = write_figure_qc(rows, tmp_path / "figure_qc.tsv")
        text = out.read_text(encoding="utf-8")
        assert text.startswith("# schema_version: 1.0\n")
        # Path is relative to run root.
        assert "rnaseq/v.png" in text
        # Column header order is stable.
        header = text.splitlines()[1]
        assert header.split("\t")[0] == "plot_path"
        assert "status" in header

    def test_empty_records_writes_header_only(self, tmp_path: Path) -> None:
        out = write_figure_qc([], tmp_path / "figure_qc.tsv")
        text = out.read_text(encoding="utf-8")
        assert text.startswith("# schema_version: 1.0\n")
        assert text.count("\n") == 2  # schema + header lines
