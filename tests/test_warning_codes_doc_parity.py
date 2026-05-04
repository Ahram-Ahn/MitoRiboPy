"""Parity test between the WARNING_CODES registry and the user docs.

The runtime registry (``src/mitoribopy/io/warning_codes.py``) is the
single source of truth for every warning / error emitted by
``warnings_log.record(code=...)``. The user-facing reference page at
``docs/reference/warning_codes.md`` is hand-curated; without this
test there is no mechanical guarantee that adding a new code (or
removing an old one) updates the doc.

Asserted invariants:

* every code in ``WARNING_CODES`` appears verbatim in
  ``docs/reference/warning_codes.md``;
* every code mentioned in the doc maps back to a registry entry
  (catches typos and codes the doc keeps after the registry drops
  them).

The test is cheap (one file read + one frozen-tuple iteration) and
runs in the default pytest pass.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from mitoribopy.io.warning_codes import WARNING_CODES


_DOC_PATH = (
    Path(__file__).resolve().parent.parent
    / "docs"
    / "reference"
    / "warning_codes.md"
)

# Match every code-shaped token that looks like the registry's
# canonical form: prefix `E_` or `W_` followed by SCREAMING_SNAKE.
# Restricting to that shape means casual prose mentioning unrelated
# uppercase tokens does not pollute the parity check.
_CODE_PATTERN = re.compile(r"\b([EW]_[A-Z][A-Z0-9_]+)\b")


def _doc_codes() -> set[str]:
    text = _DOC_PATH.read_text(encoding="utf-8")
    return set(_CODE_PATTERN.findall(text))


def _registry_codes() -> set[str]:
    return {entry.code for entry in WARNING_CODES}


def test_every_registry_code_is_documented() -> None:
    """Every WARNING_CODES entry must appear in the user-facing doc."""
    missing = sorted(_registry_codes() - _doc_codes())
    assert not missing, (
        "These registry codes are missing from "
        f"{_DOC_PATH.relative_to(_DOC_PATH.parents[2])}: {missing}. "
        "Add a row under the matching section so users can interpret "
        "the warning when it shows up in `warnings.tsv`."
    )


def test_every_documented_code_exists_in_registry() -> None:
    """Every code mentioned in the doc must map to a registry entry."""
    extras = sorted(_doc_codes() - _registry_codes())
    assert not extras, (
        "These codes appear in the doc but not in WARNING_CODES: "
        f"{extras}. Either add the code to "
        "src/mitoribopy/io/warning_codes.py, or remove the stale "
        "doc row."
    )


@pytest.mark.parametrize("entry", WARNING_CODES, ids=lambda e: e.code)
def test_each_registry_entry_has_remediation_text(entry) -> None:
    """A code without remediation guidance is half-useless to a user."""
    assert entry.summary.strip(), f"{entry.code}: empty summary"
    assert entry.remediation.strip(), f"{entry.code}: empty remediation"
