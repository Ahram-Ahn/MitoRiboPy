"""Tests for the v0.9.0 run_manifest.json schema.

The schema is intentionally minimal: it pins the publication-grade
reproducibility surface (mitoribopy_version, command, config + sample-
sheet hashes, stages, output_schemas, outputs, tools, warnings,
runtime + platform fields) and leaves stage-specific blocks open. A
silent rename or removal of any of those required fields breaks
downstream reproducibility tooling, so the schema is the canonical
guard.

The fallback validator (used when ``jsonschema`` isn't installed)
covers the required-key + top-level type checks. The full draft
2020-12 validator covers patterns and nested structures.
"""

from __future__ import annotations

import json

from mitoribopy.cli.all_ import MANIFEST_SCHEMA_VERSION
from mitoribopy.data import load_run_manifest_schema, validate_run_manifest


def _minimal_valid_manifest() -> dict:
    """Smallest manifest that satisfies the v0.9.0 schema."""
    return {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "subcommand": "all",
        "mitoribopy_version": "0.9.0",
        "git_commit": None,
        "command": "mitoribopy all --config foo.yaml --output run/",
        "config_source": "foo.yaml",
        "config_source_sha256": "0" * 64,
        "config_canonical": {"align": {}, "rpf": {}},
        "sample_sheet": None,
        "sample_sheet_sha256": None,
        "reference_checksum": None,
        "stages": {
            "align": {"status": "completed", "runtime_seconds": 12.3},
            "rpf":   {"status": "completed", "runtime_seconds": 8.1},
            "rnaseq":{"status": "skipped", "reason": "no rnaseq section"},
        },
        "align": {},
        "rpf": {},
        "rnaseq": None,
        "tools": {"cutadapt": "5.0", "bowtie2": "2.5.4"},
        "output_schemas": {"read_counts.tsv": "1.0"},
        "warnings": [],
        "outputs": [],
        "runtime_seconds": 25.7,
        "platform": "Linux-6.5.0-x86_64",
        "python_version": "3.12.0",
        "resource_plan": None,
    }


# ---------------------------------------------------------------------------
# Schema invariants
# ---------------------------------------------------------------------------


def test_schema_loads_and_declares_required_keys() -> None:
    schema = load_run_manifest_schema()
    assert schema["title"] == "MitoRiboPy run_manifest.json"
    assert "schema_version" in schema["required"]
    assert "mitoribopy_version" in schema["required"]
    assert "stages" in schema["required"]
    assert "output_schemas" in schema["required"]


def test_schema_id_matches_canonical_repo_path() -> None:
    schema = load_run_manifest_schema()
    assert schema["$id"].endswith("run_manifest.schema.json")


# ---------------------------------------------------------------------------
# Validation: positive case
# ---------------------------------------------------------------------------


def test_validate_run_manifest_accepts_minimal_valid_manifest() -> None:
    manifest = _minimal_valid_manifest()
    errors = validate_run_manifest(manifest)
    assert errors == [], f"unexpected validation errors: {errors}"


# ---------------------------------------------------------------------------
# Validation: negative cases — every required key must be flagged
# ---------------------------------------------------------------------------


def test_validate_rejects_missing_mitoribopy_version() -> None:
    manifest = _minimal_valid_manifest()
    del manifest["mitoribopy_version"]
    errors = validate_run_manifest(manifest)
    assert any("mitoribopy_version" in e for e in errors), errors


def test_validate_rejects_missing_command() -> None:
    manifest = _minimal_valid_manifest()
    del manifest["command"]
    errors = validate_run_manifest(manifest)
    assert any("command" in e for e in errors), errors


def test_validate_rejects_missing_stages() -> None:
    manifest = _minimal_valid_manifest()
    del manifest["stages"]
    errors = validate_run_manifest(manifest)
    assert any("stages" in e for e in errors), errors


def test_validate_rejects_wrong_runtime_type() -> None:
    manifest = _minimal_valid_manifest()
    manifest["runtime_seconds"] = "twenty-five seconds"  # str, not number
    errors = validate_run_manifest(manifest)
    assert any("runtime_seconds" in e for e in errors), errors


def test_validate_rejects_missing_output_schemas() -> None:
    manifest = _minimal_valid_manifest()
    del manifest["output_schemas"]
    errors = validate_run_manifest(manifest)
    assert any("output_schemas" in e for e in errors), errors


# ---------------------------------------------------------------------------
# JSON-serialisability invariant
# ---------------------------------------------------------------------------


def test_minimal_manifest_round_trips_through_json() -> None:
    """The schema-test fixture must itself round-trip through json.dumps
    so a real manifest that adds a non-serialisable field would fail
    fast in CI."""
    manifest = _minimal_valid_manifest()
    blob = json.dumps(manifest)
    parsed = json.loads(blob)
    assert validate_run_manifest(parsed) == []
