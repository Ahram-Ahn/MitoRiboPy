"""Cross-source version consistency check.

Catches the failure mode where pyproject.toml, the package fallback,
the README, the generated CLI reference, and example template headers
disagree about which release of MitoRiboPy is current.

A bioinformatics package whose user-facing surfaces all advertise
different versions looks abandoned. This test fails CI before such a
mismatch can ship.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

if sys.version_info >= (3, 11):  # pragma: no cover - python-version branch
    import tomllib
else:  # pragma: no cover - python-version branch
    import tomli as tomllib


ROOT = Path(__file__).resolve().parent.parent


def _project_version() -> str:
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text())
    return pyproject["project"]["version"]


def test_pyproject_matches_package_fallback() -> None:
    version = _project_version()
    init_text = (ROOT / "src" / "mitoribopy" / "__init__.py").read_text()
    match = re.search(r'__version__\s*=\s*"([^"]+)"', init_text)
    assert match is not None, "src/mitoribopy/__init__.py has no __version__ fallback"
    assert match.group(1) == version, (
        f"pyproject.toml version {version!r} does not match "
        f"src/mitoribopy/__init__.py fallback {match.group(1)!r}"
    )


def test_cli_reference_matches_pyproject() -> None:
    version = _project_version()
    text = (ROOT / "docs" / "reference" / "cli.md").read_text()
    needle = f"Generated against MitoRiboPy v{version}."
    assert needle in text, (
        f"docs/reference/cli.md does not mention v{version}; "
        f"regenerate with `python docs/generate_cli_reference.py`."
    )


def test_readme_pins_to_current_version() -> None:
    version = _project_version()
    text = (ROOT / "README.md").read_text()
    assert f"v{version}" in text or f">={version}" in text, (
        f"README.md does not advertise v{version}; the install / pin "
        "instructions are out of date."
    )


def test_template_headers_advertise_current_version() -> None:
    version = _project_version()
    template_dir = ROOT / "examples" / "templates"
    needle = f"# Compatible with: MitoRiboPy {version}+"
    stale: list[str] = []
    for template in sorted(template_dir.glob("*")):
        if not template.is_file():
            continue
        if template.suffix not in {".yaml", ".yml", ".sh"}:
            continue
        text = template.read_text()
        if "Compatible with: MitoRiboPy" not in text:
            continue
        if needle not in text:
            stale.append(template.name)
    assert not stale, (
        f"Templates do not advertise compatibility with MitoRiboPy {version}+: "
        f"{', '.join(stale)}"
    )


def test_changelog_has_current_version_section() -> None:
    version = _project_version()
    text = (ROOT / "CHANGELOG.md").read_text()
    pattern = re.compile(rf"^##\s*\[{re.escape(version)}\]", re.MULTILINE)
    assert pattern.search(text), (
        f"CHANGELOG.md has no [{version}] section. Move the relevant "
        "Unreleased entries into a versioned release header before tagging."
    )


@pytest.mark.parametrize(
    "path",
    [
        Path("README.md"),
        Path("docs") / "reference" / "cli.md",
    ],
)
def test_no_stale_version_pins(path: Path) -> None:
    """Reject install / verify snippets that still pin to a prior release."""

    version = _project_version()
    text = (ROOT / path).read_text()
    # Look for `>=X.Y.Z` install pins. Anything not equal to the current
    # version is a stale pin (manuscript-reproducibility pins should
    # quote a git tag, e.g. `git checkout vX.Y.Z`, not a pip pin).
    pins = set(re.findall(r">=\s*(\d+\.\d+\.\d+)", text))
    pins.discard(version)
    assert not pins, (
        f"{path} still pins to prior MitoRiboPy versions: {sorted(pins)}. "
        f"Update install instructions to >={version}."
    )
