"""§9 P0 — every YAML template under ``examples/templates`` parses cleanly.

This is a minimal "smoke" test that complements
``tests/test_docs_execution.py`` (which also verifies key→flag
mapping). The contract checked here is narrower:

* the file loads via ``mitoribopy.cli.common.load_config_file``;
* the result canonicalises through
  :func:`mitoribopy.config.canonicalize_config` without producing
  any ``warn``-level :class:`ConfigChange` (i.e. no unknown
  top-level sections, no legacy keys that ought to have been
  rewritten earlier).

Templates land in this test simply by being added to the directory:
the loop below globs ``examples/templates/*.yaml`` and runs both
checks against each.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mitoribopy.cli.common import load_config_file
from mitoribopy.config import canonicalize_config


REPO_ROOT = Path(__file__).resolve().parents[1]
TEMPLATES_DIR = REPO_ROOT / "examples" / "templates"


def _template_paths() -> list[Path]:
    return sorted(TEMPLATES_DIR.glob("*.yaml"))


@pytest.mark.parametrize("template", _template_paths(), ids=lambda p: p.name)
def test_template_parses_via_load_config_file(template: Path) -> None:
    cfg = load_config_file(str(template))
    assert isinstance(cfg, dict)
    # Templates are intentionally exhaustive; the result must be
    # non-empty so the test catches "empty file slipped in".
    assert cfg, f"template {template.name} parsed to an empty dict"


@pytest.mark.parametrize("template", _template_paths(), ids=lambda p: p.name)
def test_template_canonicalises_without_warnings(template: Path) -> None:
    cfg = load_config_file(str(template))
    result = canonicalize_config(cfg)
    assert not result.has_warnings, (
        f"template {template.name} emitted canonicalisation warnings: "
        + "; ".join(result.warning_messages())
    )


def test_template_directory_is_non_empty() -> None:
    """A safety net: a refactor that accidentally moves the templates
    directory should fail this test loudly rather than silently
    parametrising over zero cases."""
    assert _template_paths(), (
        f"no *.yaml templates found under {TEMPLATES_DIR} — has the "
        "directory moved?"
    )
