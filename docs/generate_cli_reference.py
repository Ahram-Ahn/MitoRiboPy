"""Regenerate ``docs/reference/cli.md`` from the live argparse parsers.

Per the publication-readiness audit (report §5.2):

  Add a small script:

      python docs/generate_cli_reference.py

  It should write docs/reference/cli.md by importing each subcommand's
  parser and rendering ``parser.format_help()`` into a single Markdown
  document. CI can then diff the regenerated file against the tracked
  one to detect README ↔ --help drift.

Why import-rather-than-subprocess
---------------------------------
Running ``mitoribopy <subcommand> --help`` in a subprocess would also
work, but importing keeps the regeneration step self-contained: no
``mitoribopy`` install required (just ``PYTHONPATH=src`` from the
repo root), and the script can be run from any CI image that has the
declared dev dependencies.

Usage
-----
::

    python docs/generate_cli_reference.py             # write the file
    python docs/generate_cli_reference.py --check     # exit 1 if stale

The ``--check`` mode is what CI should call.
"""

from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path
from typing import Callable

# Make ``src/mitoribopy`` importable when this script is run directly
# from a checkout (``python docs/generate_cli_reference.py``) without
# requiring an editable install.
_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT / "src"))

from mitoribopy import __version__  # noqa: E402
from mitoribopy.cli import align as align_cli  # noqa: E402
from mitoribopy.cli import all_ as all_cli  # noqa: E402
from mitoribopy.cli import benchmark as benchmark_cli  # noqa: E402
from mitoribopy.cli import migrate_config as migrate_config_cli  # noqa: E402
from mitoribopy.cli import periodicity as periodicity_cli  # noqa: E402
from mitoribopy.cli import rnaseq as rnaseq_cli  # noqa: E402
from mitoribopy.cli import summarize as summarize_cli  # noqa: E402
from mitoribopy.cli import validate_config as validate_config_cli  # noqa: E402
from mitoribopy.cli import validate_figures as validate_figures_cli  # noqa: E402
from mitoribopy.cli import validate_reference as validate_reference_cli  # noqa: E402
from mitoribopy.config import DEFAULT_CONFIG  # noqa: E402
from mitoribopy.pipeline import runner as pipeline_runner  # noqa: E402


# (subcommand-name, help-summary, factory) tuples.
_ParserFactory = Callable[[], argparse.ArgumentParser]


def _rpf_parser_factory() -> argparse.ArgumentParser:
    """The rpf parser is built once per call from runtime defaults so a
    test fixture that mutates DEFAULT_CONFIG cannot leak into the
    generated reference."""
    return pipeline_runner.build_parser(dict(DEFAULT_CONFIG))


SUBCOMMANDS: list[tuple[str, str, _ParserFactory]] = [
    ("align", align_cli.ALIGN_SUBCOMMAND_HELP, align_cli.build_parser),
    ("rpf", "Ribo-seq analysis from BED/BAM inputs.", _rpf_parser_factory),
    ("rnaseq", rnaseq_cli.RNASEQ_SUBCOMMAND_HELP, rnaseq_cli.build_parser),
    ("all", all_cli.ALL_SUBCOMMAND_HELP, all_cli.build_parser),
    (
        "periodicity",
        periodicity_cli.PERIODICITY_SUBCOMMAND_HELP,
        periodicity_cli.build_parser,
    ),
    (
        "migrate-config",
        migrate_config_cli.MIGRATE_CONFIG_HELP,
        migrate_config_cli.build_parser,
    ),
    (
        "validate-config",
        validate_config_cli.VALIDATE_CONFIG_HELP,
        validate_config_cli.build_parser,
    ),
    (
        "validate-reference",
        validate_reference_cli.VALIDATE_REFERENCE_HELP,
        validate_reference_cli.build_parser,
    ),
    (
        "validate-figures",
        validate_figures_cli.VALIDATE_FIGURES_HELP,
        validate_figures_cli.build_parser,
    ),
    (
        "summarize",
        summarize_cli.SUMMARIZE_SUBCOMMAND_HELP,
        summarize_cli.build_parser,
    ),
    (
        "benchmark",
        benchmark_cli.BENCHMARK_SUBCOMMAND_HELP,
        benchmark_cli.build_parser,
    ),
]


HEADER = """\
# CLI reference

> **GENERATED FILE — do not hand-edit.** Regenerate with
> `python docs/generate_cli_reference.py`. CI runs the same script in
> `--check` mode to fail when this file drifts from the live argparse
> parsers in `src/mitoribopy/`.

This document is the canonical machine reference for every flag in
every `mitoribopy` subcommand. For prose, examples, and decision
trees see [the README](../../README.md) and the tutorials under
[`docs/tutorials/`](../tutorials/).

Generated against MitoRiboPy v{version}.

## Subcommand summary

| Subcommand | What it does |
|---|---|
"""


def _format_subcommand_block(name: str, factory: _ParserFactory) -> str:
    """Render one ``## mitoribopy <name>`` section as Markdown."""
    parser = factory()
    buf = io.StringIO()
    parser.print_help(file=buf)
    help_text = buf.getvalue().rstrip() + "\n"
    return (
        f"## `mitoribopy {name}`\n\n"
        "```text\n"
        f"{help_text}"
        "```\n\n"
    )


def render() -> str:
    """Return the full Markdown document."""
    out: list[str] = [HEADER.format(version=__version__)]
    for name, summary, _ in SUBCOMMANDS:
        # Markdown table escaping: collapse pipes inside the summary so
        # they don't break the table.
        clean = summary.replace("|", "\\|").replace("\n", " ").strip()
        out.append(f"| [`mitoribopy {name}`](#mitoribopy-{name.replace('-', '')}) | {clean} |\n")
    out.append("\n---\n\n")
    for name, _, factory in SUBCOMMANDS:
        out.append(_format_subcommand_block(name, factory))
    return "".join(out)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        default=False,
        help=(
            "Do not write the file. Compute the rendered reference, "
            "compare against the tracked file, and exit 1 if they "
            "differ. Suitable for CI."
        ),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=_REPO_ROOT / "docs" / "reference" / "cli.md",
        help="Output path. Default: docs/reference/cli.md.",
    )
    args = parser.parse_args(argv)

    rendered = render()

    if args.check:
        existing = args.out.read_text(encoding="utf-8") if args.out.is_file() else ""
        if rendered != existing:
            sys.stderr.write(
                f"[generate_cli_reference] {args.out} is stale. Run "
                "`python docs/generate_cli_reference.py` and commit "
                "the result.\n"
            )
            return 1
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(rendered, encoding="utf-8")
    sys.stderr.write(f"[generate_cli_reference] wrote {args.out}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
