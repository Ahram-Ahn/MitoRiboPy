"""Tests for ``mitoribopy.io.outputs_index`` (§8 output contracts)."""

from __future__ import annotations

from pathlib import Path

import pytest

from mitoribopy.io.outputs_index import (
    build_outputs_index_rows,
    write_outputs_index,
)


def _seed(run_root: Path) -> None:
    (run_root / "align").mkdir(parents=True, exist_ok=True)
    (run_root / "align" / "read_counts.tsv").write_text("# stub\n")
    (run_root / "align" / "kit_resolution.tsv").write_text("# stub\n")
    (run_root / "rpf").mkdir(parents=True, exist_ok=True)
    (run_root / "rpf" / "rpf_counts.tsv").write_text("# stub\n")
    (run_root / "rnaseq").mkdir(parents=True, exist_ok=True)
    (run_root / "rnaseq" / "te.tsv").write_text("# stub\n")
    (run_root / "warnings.tsv").write_text("# stub\n")
    (run_root / "SUMMARY.md").write_text("stub\n")


class TestBuildRows:
    def test_skips_missing_outputs(self, tmp_path: Path) -> None:
        # Empty run root → no rows.
        rows = build_outputs_index_rows(tmp_path)
        assert rows == []

    def test_only_existing_outputs_included(self, tmp_path: Path) -> None:
        _seed(tmp_path)
        rows = build_outputs_index_rows(tmp_path)
        types = {r["output_type"] for r in rows}
        assert "read_counts" in types
        assert "kit_resolution" in types
        assert "rpf_counts" in types
        assert "te_table" in types
        assert "warnings" in types
        assert "summary_md" in types
        # delta_te.tsv was not seeded so the row must be absent.
        assert "delta_te_table" not in types

    def test_columns_present_per_row(self, tmp_path: Path) -> None:
        _seed(tmp_path)
        rows = build_outputs_index_rows(tmp_path)
        required = {
            "output_type",
            "stage",
            "path",
            "description",
            "recommended_for",
            "schema_version",
        }
        for row in rows:
            assert set(row.keys()) == required

    def test_schema_version_filled_for_known_tsv(self, tmp_path: Path) -> None:
        _seed(tmp_path)
        rows = build_outputs_index_rows(tmp_path)
        rc = next(r for r in rows if r["output_type"] == "read_counts")
        assert rc["schema_version"]  # not empty
        # SUMMARY.md is not a registered TSV → schema_version is empty.
        smd = next(r for r in rows if r["output_type"] == "summary_md")
        assert smd["schema_version"] == ""


class TestWriteOutputsIndex:
    def test_writes_header_and_schema_comment(self, tmp_path: Path) -> None:
        _seed(tmp_path)
        path = write_outputs_index(tmp_path)
        assert path.exists()
        lines = path.read_text().splitlines()
        assert lines[0].startswith("# schema_version: ")
        header = lines[1].split("\t")
        assert header == [
            "output_type",
            "stage",
            "path",
            "description",
            "recommended_for",
            "schema_version",
        ]

    def test_empty_run_root_still_writes_header(self, tmp_path: Path) -> None:
        path = write_outputs_index(tmp_path)
        assert path.exists()
        lines = path.read_text().splitlines()
        # comment + header, no data rows
        assert len(lines) == 2

    def test_paths_are_relative_posix(self, tmp_path: Path) -> None:
        _seed(tmp_path)
        write_outputs_index(tmp_path)
        text = (tmp_path / "outputs_index.tsv").read_text()
        # No absolute paths from our base dir leaked in.
        assert str(tmp_path) not in text
        # The kit_resolution row uses POSIX-style relative path.
        assert "align/kit_resolution.tsv" in text
