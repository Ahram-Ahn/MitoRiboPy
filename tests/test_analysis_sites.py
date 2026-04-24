"""Tests for the --analysis_sites flag and per-site downstream layout."""

from __future__ import annotations

from pathlib import Path

from mitoribopy import cli as _cli  # noqa: F401  pre-import to break import cycle
from mitoribopy.config.runtime import DEFAULT_CONFIG
from mitoribopy.pipeline.runner import build_parser


def test_analysis_sites_default_is_both() -> None:
    parser = build_parser(dict(DEFAULT_CONFIG))
    ns = parser.parse_args(["-f", "ref.fa"])
    assert ns.analysis_sites == "both"


def test_analysis_sites_accepts_explicit_p_or_a() -> None:
    parser = build_parser(dict(DEFAULT_CONFIG))
    p_args = parser.parse_args(["-f", "ref.fa", "--analysis_sites", "p"])
    a_args = parser.parse_args(["-f", "ref.fa", "--analysis_sites", "a"])
    assert p_args.analysis_sites == "p"
    assert a_args.analysis_sites == "a"


def test_offset_mode_default_is_per_sample() -> None:
    parser = build_parser(dict(DEFAULT_CONFIG))
    ns = parser.parse_args(["-f", "ref.fa"])
    assert ns.offset_mode == "per_sample"


def test_offset_site_help_clarifies_relationship_with_analysis_sites() -> None:
    parser = build_parser(dict(DEFAULT_CONFIG))
    help_text = parser.format_help()
    # Both flags should be present and mutually referenced in the help so
    # users can find the right knob without trial and error.
    assert "--analysis_sites" in help_text
    assert "--offset_site" in help_text
    assert "analysis_sites" in help_text  # cross-reference in offset_site help


def test_offset_mode_help_lists_per_sample_and_combined() -> None:
    parser = build_parser(dict(DEFAULT_CONFIG))
    help_text = parser.format_help()
    assert "per_sample" in help_text
    assert "combined" in help_text
    assert "offset_drift" in help_text  # surface the diagnostic plot name
