"""Parity test between the runtime output registries and the user docs.

The runtime is the single source of truth for what files MitoRiboPy
emits. The user-facing reference at ``docs/reference/output_schema.md``
is hand-curated; without this test there is no mechanical guarantee
that adding a new output (or bumping a schema version) updates the
doc.

Asserted invariants:

* every ``relative_path`` in
  ``mitoribopy.io.outputs_index._KNOWN_OUTPUTS`` is mentioned in the
  doc (either verbatim or, for stage prefixes shared by sibling
  outputs like the codon-correlation `{p_site,a_site}` parallel
  files, the doc may reference the brace-expansion form);
* every TSV name in
  ``mitoribopy.io.schema_versions.OUTPUT_SCHEMA_VERSIONS`` appears in
  the doc with its current schema version (e.g. `read_counts.tsv`
  with `schema 1.0`);
* every TSV file mentioned in the doc maps back to either a registry
  entry or the OUTPUT_SCHEMA_VERSIONS map (catches stale doc rows
  for files the runtime no longer writes).

The test is cheap (one file read per registry walk) and runs in the
default pytest pass.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from mitoribopy.io.outputs_index import _KNOWN_OUTPUTS
from mitoribopy.io.schema_versions import OUTPUT_SCHEMA_VERSIONS


_DOC_PATH = (
    Path(__file__).resolve().parent.parent
    / "docs"
    / "reference"
    / "output_schema.md"
)


def _doc_text() -> str:
    return _DOC_PATH.read_text(encoding="utf-8")


def _doc_mentions_path(text: str, relative_path: str) -> bool:
    """Return True when *relative_path* is referenced by the doc.

    Accepts a few canonicalised forms so that paired sibling files
    (codon_correlation/{p_site,a_site}/...) can be documented with
    a single brace-expansion mention rather than two near-identical
    rows. The text is rendered Markdown; we only need substring
    presence — the doc convention is to backtick the path.
    """
    if relative_path in text:
        return True
    # Brace-expansion: rpf/codon_correlation/p_site/foo.tsv vs
    # rpf/codon_correlation/{p_site,a_site}/foo.tsv.
    parts = relative_path.split("/")
    for i in range(1, len(parts) - 1):
        peer = parts[i]
        if peer in {"p_site", "a_site"}:
            brace = "/".join(
                parts[:i] + ["{p_site,a_site}"] + parts[i + 1:]
            )
            if brace in text:
                return True
    return False


def test_every_registry_path_is_documented() -> None:
    """Every ``_KNOWN_OUTPUTS`` entry must appear in the doc."""
    text = _doc_text()
    missing = sorted(
        desc.relative_path
        for desc in _KNOWN_OUTPUTS
        if not _doc_mentions_path(text, desc.relative_path)
    )
    assert not missing, (
        "These outputs are emitted by `_KNOWN_OUTPUTS` but missing "
        f"from {_DOC_PATH.relative_to(_DOC_PATH.parents[2])}: {missing}. "
        "Add a row so users can interpret the file when they see "
        "it under their run root."
    )


def test_every_schema_versioned_tsv_appears_in_doc() -> None:
    """Every entry in ``OUTPUT_SCHEMA_VERSIONS`` must be mentioned in
    the doc paired with its current schema version (e.g.
    ``read_counts.tsv`` documented as ``schema 1.0``).
    """
    text = _doc_text()
    missing: list[str] = []
    wrong_version: list[str] = []
    for filename, version in OUTPUT_SCHEMA_VERSIONS.items():
        if filename not in text:
            missing.append(filename)
            continue
        # Anchor pattern: the doc convention is "(schema 1.0)" or
        # "(schema 1.1)" right after the filename heading.
        # Look for the version token anywhere in the same line as
        # the filename to keep the check tolerant of formatting.
        pattern = re.compile(
            re.escape(filename)
            + r"[^\n]{0,200}schema "
            + re.escape(version)
        )
        if not pattern.search(text):
            wrong_version.append(f"{filename} (expected schema {version})")
    assert not missing, (
        "These OUTPUT_SCHEMA_VERSIONS entries are not mentioned in "
        f"{_DOC_PATH.name}: {missing}. Add a row."
    )
    assert not wrong_version, (
        "These OUTPUT_SCHEMA_VERSIONS entries are mentioned in the "
        "doc but with a stale version annotation; bump the doc to "
        f"match the runtime: {wrong_version}."
    )


_DOC_FILE_PATTERN = re.compile(r"`([A-Za-z0-9_./-]+\.tsv)`")


def test_every_documented_tsv_resolves_to_a_known_output() -> None:
    """Every ``*.tsv`` filename mentioned in backticks in the doc must
    map to either a registry path or an OUTPUT_SCHEMA_VERSIONS entry."""
    text = _doc_text()
    documented = set(_DOC_FILE_PATTERN.findall(text))
    known: set[str] = set(OUTPUT_SCHEMA_VERSIONS.keys())
    for desc in _KNOWN_OUTPUTS:
        if desc.relative_path.endswith(".tsv"):
            known.add(desc.relative_path)
            known.add(desc.relative_path.split("/")[-1])
    extras = sorted(documented - known)
    assert not extras, (
        "These .tsv filenames appear in the doc but are not produced "
        "by any registry entry or OUTPUT_SCHEMA_VERSIONS row: "
        f"{extras}. Either add the runtime entry or remove the stale "
        "doc reference."
    )


@pytest.mark.parametrize("entry", _KNOWN_OUTPUTS, ids=lambda e: e.output_type)
def test_each_output_descriptor_has_description(entry) -> None:
    """An output without a description is half-useless to a user."""
    assert entry.description.strip(), (
        f"{entry.output_type}: empty description in _KNOWN_OUTPUTS"
    )
    assert entry.recommended_for in {
        "downstream-scripting",
        "reviewer-spot-check",
    }, (
        f"{entry.output_type}: recommended_for={entry.recommended_for!r} "
        "is not one of the documented audiences."
    )
