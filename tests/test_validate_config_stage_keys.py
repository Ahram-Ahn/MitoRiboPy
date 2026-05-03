"""Tests for the stage-level unknown-key check added in v0.8.0.

The audit before publication called out that ``validate-config`` only
caught unknown keys at the *top* level. A typo nested inside a stage
section (``rpf.offset_pick_refernce``) silently fell through to default
behaviour, which is dangerous on long HPC runs that back a manuscript.

These tests cover:

* the stage allowlist is non-empty for every stage we validate,
* the unknown-key reporter returns ``did you mean`` suggestions,
* ``validate-config`` surfaces stage-level unknowns as warnings in
  normal mode and as exit-2 errors under ``--strict``.
"""

from __future__ import annotations

import io
import sys
from pathlib import Path
from textwrap import dedent

import pytest

from mitoribopy.cli import validate_config as vc
from mitoribopy.config import stage_keys


@pytest.mark.parametrize(
    "stage", ["align", "rpf", "rnaseq", "execution", "periodicity"]
)
def test_stage_allowlist_is_non_empty(stage: str) -> None:
    keys = stage_keys.allowed_keys_for(stage)
    assert keys, f"Stage {stage!r} has an empty allowlist; nothing would parse."


def test_typo_in_rpf_is_flagged_with_suggestion() -> None:
    section = {"offset_pick_refernce": "p_site"}
    unknowns = stage_keys.find_unknown_keys(section, "rpf")
    assert len(unknowns) == 1
    key, suggestions = unknowns[0]
    assert key == "offset_pick_refernce"
    assert "offset_pick_reference" in suggestions


def test_known_rpf_keys_pass() -> None:
    section = {
        "offset_pick_reference": "p_site",
        "footprint_class": "monosome",
        "rpf": [29, 34],
    }
    assert stage_keys.find_unknown_keys(section, "rpf") == []


def test_known_align_keys_pass() -> None:
    section = {
        "kit_preset": "auto",
        "adapter_detection": "auto",
        "dedup_strategy": "auto",
        "library_strandedness": "forward",
        "samples": [],  # YAML-only convention; allowed.
        "fastq": "data/",  # YAML-only convention; allowed.
    }
    assert stage_keys.find_unknown_keys(section, "align") == []


def test_typo_in_execution_is_flagged() -> None:
    section = {"thread": 8}  # missing 's'
    unknowns = stage_keys.find_unknown_keys(section, "execution")
    assert len(unknowns) == 1
    key, suggestions = unknowns[0]
    assert key == "thread"
    assert "threads" in suggestions


def _write_config(tmp_path: Path, body: str) -> Path:
    cfg = tmp_path / "pipeline_config.yaml"
    cfg.write_text(dedent(body))
    return cfg


def _run_validator(argv: list[str]) -> tuple[int, str]:
    captured = io.StringIO()
    real_stderr = sys.stderr
    sys.stderr = captured
    try:
        rc = vc.run(argv)
    finally:
        sys.stderr = real_stderr
    return rc, captured.getvalue()


def test_validate_config_warns_on_stage_typo_in_normal_mode(
    tmp_path: Path,
) -> None:
    cfg = _write_config(
        tmp_path,
        """
        rpf:
          offset_pick_refernce: p_site
          footprint_class: monosome
        """,
    )
    rc, stderr = _run_validator([str(cfg), "--no-path-checks"])
    assert "rpf.offset_pick_refernce" in stderr
    assert "Did you mean 'offset_pick_reference'" in stderr
    # Normal mode: warnings do not flip the exit code.
    assert rc == 0


def test_validate_config_fails_on_stage_typo_in_strict_mode(
    tmp_path: Path,
) -> None:
    cfg = _write_config(
        tmp_path,
        """
        rpf:
          offset_pick_refernce: p_site
          footprint_class: monosome
        """,
    )
    rc, stderr = _run_validator(
        [str(cfg), "--no-path-checks", "--strict"]
    )
    assert "ERROR" in stderr
    assert "rpf.offset_pick_refernce" in stderr
    assert rc == 2


def test_validate_config_passes_clean_publication_config(
    tmp_path: Path,
) -> None:
    cfg = _write_config(
        tmp_path,
        """
        execution:
          threads: 8
        rpf:
          footprint_class: monosome
          offset_pick_reference: p_site
        periodicity:
          enabled: true
        """,
    )
    rc, stderr = _run_validator(
        [str(cfg), "--no-path-checks", "--strict"]
    )
    assert rc == 0, stderr
