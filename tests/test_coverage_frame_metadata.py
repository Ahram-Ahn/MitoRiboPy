"""Tests for the coverage-plot frame metadata sidecar (publication-readiness)."""

from __future__ import annotations

import json
from pathlib import Path

from mitoribopy.plotting.coverage_profile_plots import (
    _FRAME_LABELS,
    _write_coverage_frame_metadata,
)


def test_p_site_metadata_uses_p_site_formula(tmp_path: Path) -> None:
    out = _write_coverage_frame_metadata(
        frame_dir=tmp_path / "p_site_density_rpm_frame",
        site="p",
        normalization="RPM",
        offset_type="5",
        offset_site="p",
        included_read_lengths=[28, 30, 32],
    )
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["site"] == "P-site"
    assert payload["frame_formula"] == "(P_site_nt - CDS_start_nt) % 3"
    assert payload["included_read_lengths"] == [28, 30, 32]
    assert payload["normalization"] == "RPM"


def test_a_site_metadata_uses_a_site_formula(tmp_path: Path) -> None:
    out = _write_coverage_frame_metadata(
        frame_dir=tmp_path / "a_site_density_raw_frame",
        site="a",
        normalization="raw_counts",
        offset_type="5",
        offset_site="p",
    )
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["site"] == "A-site"
    assert "A-site = P-site + 3 nt" in payload["frame_formula"]


def test_metadata_includes_publication_friendly_labels(tmp_path: Path) -> None:
    out = _write_coverage_frame_metadata(
        frame_dir=tmp_path / "p_site_density_rpm_frame",
        site="p",
        normalization="RPM",
        offset_type="5",
        offset_site="p",
    )
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["frame_labels"] == list(_FRAME_LABELS)
    # The labels must be unambiguous about the coordinate system.
    assert payload["frame_labels"][0].startswith("Frame 0")
    assert "annotated CDS frame" in payload["frame_labels"][0]
