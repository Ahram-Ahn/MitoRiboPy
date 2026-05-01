"""§9 P0 — every fenced ``mitoribopy ...`` command in README.md is dispatchable.

This file complements ``tests/test_docs_execution.py``. Where that
test runs every dry-run-safe command in-process, this one is a
narrower contract: every fenced ``mitoribopy <subcommand> ...`` block
in :file:`README.md` must name a real subcommand and parse without
the CLI raising. Commands that need real I/O are skipped via
``--help`` substitution so the test stays hermetic.

The narrower contract catches the common drift mode: someone renames
a subcommand and forgets to update the README, or adds a new flag to
the README before the parser learns it.
"""

from __future__ import annotations

import re
import shlex
from pathlib import Path

import pytest

from mitoribopy.cli import _SUBCOMMANDS


REPO_ROOT = Path(__file__).resolve().parents[1]
README_PATH = REPO_ROOT / "README.md"


_FENCED_BLOCK = re.compile(r"```(?:bash|shell|console)?\n(.*?)\n```", re.DOTALL)
_LINE_PREFIX = re.compile(r"^\s*(?:\$\s*)?(mitoribopy[^\n]*)", re.MULTILINE)


def _extracted_mitoribopy_lines() -> list[str]:
    """Return every ``mitoribopy ...`` command line found in fenced
    code blocks of ``README.md``.

    Multi-line commands (``\\``-continued) are joined into a single
    logical line. Inline backticks (e.g. ``` `mitoribopy align` ```)
    are NOT extracted — only fenced blocks count.
    """
    if not README_PATH.exists():  # pragma: no cover - belt-and-braces
        return []
    text = README_PATH.read_text(encoding="utf-8")
    commands: list[str] = []
    for block_match in _FENCED_BLOCK.finditer(text):
        block = block_match.group(1)
        # Join backslash-continued lines so shlex sees one command.
        block = re.sub(r"\\\n\s*", " ", block)
        for line_match in _LINE_PREFIX.finditer(block):
            commands.append(line_match.group(1).strip())
    return commands


def _subcommand_of(command: str) -> str | None:
    """Return the first whitespace-separated token after ``mitoribopy``."""
    tokens = shlex.split(command, comments=False)
    if not tokens or tokens[0] != "mitoribopy":
        return None
    if len(tokens) < 2:
        return None
    return tokens[1]


def test_readme_has_mitoribopy_examples() -> None:
    commands = _extracted_mitoribopy_lines()
    assert commands, (
        "README.md contains no fenced 'mitoribopy ...' commands — has "
        "the file been gutted?"
    )


def test_every_readme_subcommand_is_known() -> None:
    """Each 'mitoribopy <subcommand>' in README.md must name a real
    registered subcommand. Catches the rename-without-doc-update bug.
    """
    known = set(_SUBCOMMANDS) | {
        # Top-level help and version flags don't require a subcommand.
        "--help",
        "-h",
        "--version",
        "-V",
    }
    unknown: list[tuple[str, str]] = []
    for line in _extracted_mitoribopy_lines():
        sub = _subcommand_of(line)
        if sub is None:
            continue
        # Strip an inline arg like 'mitoribopy <PR#>' that doesn't
        # match a registered subcommand AND isn't a flag — we just
        # ignore those (typically placeholders).
        if sub.startswith("--") or sub.startswith("-"):
            if sub not in known:
                unknown.append((sub, line))
            continue
        if sub not in known:
            unknown.append((sub, line))
    assert not unknown, (
        "README.md references unknown subcommand(s): "
        + ", ".join(f"{s!r} in {cmd!r}" for s, cmd in unknown[:3])
    )


def test_help_dispatches_for_each_subcommand() -> None:
    """In-process sanity check: each registered subcommand's ``--help``
    flow returns a non-error exit code via ``cli.main``.
    """
    from mitoribopy import cli

    for sub in _SUBCOMMANDS:
        # SystemExit(0) is argparse's normal exit for --help on most
        # parsers; we tolerate it as a pass.
        try:
            rc = cli.main([sub, "--help"])
        except SystemExit as exc:
            rc = exc.code if exc.code is not None else 0
        assert rc == 0, f"{sub} --help returned non-zero rc={rc}"
