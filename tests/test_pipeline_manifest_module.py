"""Tests for the v0.7.0 pipeline.manifest module.

The manifest builder was extracted from cli/all_.py in v0.7.0. This
test locks in:

* Backwards-compat: ``MANIFEST_SCHEMA_VERSION`` is still importable from
  ``mitoribopy.cli.all_`` (the original location).
* The new ``mitoribopy.pipeline.manifest`` module exposes the helpers
  with public names (``sha256_of``, ``git_commit``, ``yaml_dump``,
  ``build_stages_block``, ``lift_tool_versions``, ``write_manifest``).
* ``build_stages_block`` produces the right status enum values
  (``completed``, ``skipped``, ``not_configured``).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest


def test_manifest_schema_version_re_exports_from_cli_all_() -> None:
    from mitoribopy.cli.all_ import MANIFEST_SCHEMA_VERSION as cli_version
    from mitoribopy.pipeline.manifest import (
        MANIFEST_SCHEMA_VERSION as pipeline_version,
    )
    assert cli_version == pipeline_version


def test_pipeline_manifest_public_surface() -> None:
    from mitoribopy.pipeline import manifest as m

    expected = {
        "MANIFEST_SCHEMA_VERSION",
        "W_MANIFEST_SCHEMA_DRIFT_CODE",
        "sha256_of",
        "read_stage_settings",
        "git_commit",
        "yaml_dump",
        "build_stages_block",
        "lift_tool_versions",
        "write_manifest",
    }
    assert expected.issubset(set(dir(m)))
    assert expected.issubset(set(m.__all__))


def test_sha256_of_handles_none_and_missing_paths(tmp_path: Path) -> None:
    from mitoribopy.pipeline.manifest import sha256_of

    assert sha256_of(None) is None
    assert sha256_of(tmp_path / "does_not_exist") is None


def test_sha256_of_is_deterministic_for_real_file(tmp_path: Path) -> None:
    from mitoribopy.pipeline.manifest import sha256_of

    p = tmp_path / "x.txt"
    p.write_bytes(b"hello mitoribopy v0.7.0")
    digest = sha256_of(p)
    assert isinstance(digest, str) and len(digest) == 64
    # Stable across two reads.
    assert digest == sha256_of(p)


def test_build_stages_block_uses_canonical_status_strings() -> None:
    from mitoribopy.pipeline.manifest import build_stages_block

    out = build_stages_block(
        stages_run=["align", "rpf"],
        stages_skipped=["rnaseq"],
        runtimes={"align": 12.345, "rpf": 4.5},
        skip_reasons={"rnaseq": "no rnaseq section"},
    )
    assert out["align"]["status"] == "completed"
    assert out["rpf"]["status"] == "completed"
    assert out["rnaseq"]["status"] == "skipped"
    assert out["rnaseq"]["reason"] == "no rnaseq section"
    # Runtime is rounded to 3 decimals.
    assert out["align"]["runtime_seconds"] == 12.345


def test_build_stages_block_marks_unconfigured_stages() -> None:
    from mitoribopy.pipeline.manifest import build_stages_block

    out = build_stages_block(
        stages_run=["align"], stages_skipped=[],
        runtimes={"align": 1.0}, skip_reasons={},
    )
    assert out["rpf"]["status"] == "not_configured"
    assert out["rnaseq"]["status"] == "not_configured"


def test_lift_tool_versions_includes_python_and_mitoribopy() -> None:
    from mitoribopy.pipeline.manifest import lift_tool_versions

    out = lift_tool_versions(None, {"tools": {"cutadapt": "5.0"}}, None)
    assert "python" in out
    assert "mitoribopy" in out
    assert out["cutadapt"] == "5.0"


def test_yaml_dump_falls_back_to_json_when_yaml_unavailable() -> None:
    """The dumped string must round-trip through json.loads when PyYAML
    is missing — the lean install path."""
    from mitoribopy.pipeline.manifest import yaml_dump

    payload = {"a": 1, "b": [2, 3]}
    text = yaml_dump(payload)
    # If PyYAML is present we get yaml; otherwise we get JSON. Either
    # way the YAML parser (which accepts JSON) should round-trip.
    try:
        import yaml  # type: ignore[import-not-found]

        assert yaml.safe_load(text) == payload
    except ImportError:
        assert json.loads(text) == payload


def test_w_manifest_schema_drift_code_constant_is_stable() -> None:
    """The warning code is part of the public CLI contract — it shows
    up in warnings.tsv and run_manifest.json. Renaming it would break
    downstream parsers, so the constant value is itself a fixture."""
    from mitoribopy.pipeline.manifest import W_MANIFEST_SCHEMA_DRIFT_CODE

    assert W_MANIFEST_SCHEMA_DRIFT_CODE == "W_MANIFEST_SCHEMA_DRIFT"
