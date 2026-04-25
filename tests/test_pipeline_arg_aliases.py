"""Tests for backward-compatible aliases on renamed pipeline CLI flags.

The renames in v0.4.x kept the old flag names as deprecated argparse
aliases so existing YAML configs and shell scripts keep working. These
tests pin that contract so a future refactor cannot silently drop an
alias and break a user's invocation.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mitoribopy.config.runtime import (
    DEPRECATED_CONFIG_KEY_ALIASES,
    load_user_config,
)
from mitoribopy.pipeline.runner import parse_pipeline_args


# ---------- CLI aliases -----------------------------------------------------


def _base_argv(tmp_path: Path) -> list[str]:
    """Smallest valid argv for parse_pipeline_args (no --config required)."""
    return [
        "-s", "h.sapiens",
        "-f", str(tmp_path / "fasta.fa"),
        "-d", str(tmp_path),
    ]


def test_legacy_merge_density_flag_still_parses(tmp_path) -> None:
    args = parse_pipeline_args(_base_argv(tmp_path) + ["--merge_density"])
    assert args.codon_density_window is True
    # Old attribute name doesn't exist on the namespace any more.
    assert not hasattr(args, "merge_density")


def test_legacy_mrna_ref_patterns_flag_still_parses(tmp_path) -> None:
    args = parse_pipeline_args(
        _base_argv(tmp_path) + ["--mrna_ref_patterns", "foo", "bar"]
    )
    assert args.mt_mrna_substring_patterns == ["foo", "bar"]
    assert not hasattr(args, "mrna_ref_patterns")


def test_legacy_selected_site_canonicalises_to_reported_site(
    tmp_path, capsys
) -> None:
    args = parse_pipeline_args(
        _base_argv(tmp_path)
        + ["--offset_pick_reference", "selected_site"]
    )
    assert args.offset_pick_reference == "reported_site"
    err = capsys.readouterr().err
    assert "DEPRECATED" in err
    assert "selected_site" in err


def test_canonical_reported_site_does_not_warn(tmp_path, capsys) -> None:
    args = parse_pipeline_args(
        _base_argv(tmp_path)
        + ["--offset_pick_reference", "reported_site"]
    )
    assert args.offset_pick_reference == "reported_site"
    err = capsys.readouterr().err
    assert "DEPRECATED" not in err


# ---------- YAML key aliases ------------------------------------------------


def test_legacy_yaml_keys_canonicalise(tmp_path) -> None:
    cfg = tmp_path / "legacy.yaml"
    cfg.write_text(
        "merge_density: true\n"
        "mrna_ref_patterns: [aa, bb]\n"
    )
    loaded = load_user_config(str(cfg))
    assert loaded["codon_density_window"] is True
    assert loaded["mt_mrna_substring_patterns"] == ["aa", "bb"]
    # Legacy keys are not present in the canonical dict.
    assert "merge_density" not in loaded
    assert "mrna_ref_patterns" not in loaded


def test_legacy_alias_map_pins_expected_keys() -> None:
    """Pin the contract of which legacy keys still resolve."""
    assert DEPRECATED_CONFIG_KEY_ALIASES == {
        "merge_density": "codon_density_window",
        "mrna_ref_patterns": "mt_mrna_substring_patterns",
    }


def test_canonical_yaml_keys_win_over_legacy_when_both_present(tmp_path) -> None:
    cfg = tmp_path / "both.yaml"
    cfg.write_text(
        "codon_density_window: true\n"
        "merge_density: false\n"
    )
    loaded = load_user_config(str(cfg))
    assert loaded["codon_density_window"] is True
    assert "merge_density" not in loaded
