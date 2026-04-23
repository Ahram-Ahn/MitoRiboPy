"""MitoRiboPy command-line entrypoint with subcommand dispatch.

Top-level synopsis::

    mitoribopy [--version] [--help]
    mitoribopy align   [options]   # Phase 3: FASTQ -> BAM + BED + read counts
    mitoribopy rpf     [options]   # Ribo-seq analysis pipeline from BED / BAM
    mitoribopy rnaseq  [options]   # Phase 5: DE-table + rpf -> TE / dTE tables + plots
    mitoribopy all     [options]   # Phase 6: align + rpf + (optional) rnaseq

Backward compatibility (v0.3.x only, removed in v0.4.0)
-------------------------------------------------------
Invoking ``mitoribopy`` without a subcommand and passing flags directly
(``mitoribopy -s h -f ref.fa ...``) is treated as ``mitoribopy rpf ...``
with a one-time deprecation warning printed to ``stderr``.
"""

from __future__ import annotations

import sys
from typing import Callable, Iterable

from .. import __version__
# Re-export ``run_pipeline_cli`` for monkeypatching in tests that existed
# before the subcommand refactor.
from ..pipeline.runner import run_pipeline_cli  # noqa: F401 (re-export)

from . import align as _align
from . import all_ as _all
from . import rnaseq as _rnaseq
from . import rpf as _rpf

__all__ = ["main", "run_pipeline_cli"]


_SUBCOMMANDS: dict[str, Callable[[Iterable[str]], int]] = {
    "align": _align.run,
    "rpf": _rpf.run,
    "rnaseq": _rnaseq.run,
    "all": _all.run,
}


_SUBCOMMAND_SUMMARIES: list[tuple[str, str]] = [
    ("align", _align.ALIGN_SUBCOMMAND_HELP),
    ("rpf", _rpf.RPF_SUBCOMMAND_HELP),
    ("rnaseq", _rnaseq.RNASEQ_SUBCOMMAND_HELP),
    ("all", _all.ALL_SUBCOMMAND_HELP),
]


_TOP_HELP = """\
usage: mitoribopy [--version] <subcommand> [options]

MitoRiboPy - mitochondrial ribosome profiling analysis pipeline.

Subcommands:
{rows}

Top-level options:
  --version, -V  Print MitoRiboPy version and exit.
  --help,    -h  Print this help message.

Run 'mitoribopy <subcommand> --help' for per-subcommand options.
"""


def _normalize_args(argv: Iterable[str]) -> list[str]:
    """Strip legacy 'run' / '--' separators at the front of *argv*.

    Preserved verbatim from the pre-subcommand CLI so that existing tests
    and user scripts invoking ``mitoribopy run -- <args>`` continue to work.
    """
    args = list(argv)
    if args and args[0] == "run":
        args = args[1:]
    if args and args[0] == "--":
        args = args[1:]
    return args


def _print_top_help() -> None:
    rows = "\n".join(f"  {name:<7}{summary}" for name, summary in _SUBCOMMAND_SUMMARIES)
    print(_TOP_HELP.format(rows=rows))


def _emit_legacy_fallback_warning(first_flag: str) -> None:
    sys.stderr.write(
        "[mitoribopy] DEPRECATION: invoking 'mitoribopy "
        f"{first_flag} ...' without a subcommand is deprecated. "
        "Use 'mitoribopy rpf ...' instead. This fallback will be removed "
        "in v0.4.0.\n"
    )


def main(argv: Iterable[str] | None = None) -> int:
    """Dispatch the top-level ``mitoribopy`` entry point."""
    raw_args = sys.argv[1:] if argv is None else list(argv)
    args = _normalize_args(raw_args)

    if args in (["--version"], ["-V"]):
        print(f"MitoRiboPy {__version__}")
        return 0

    if not args or args[0] in ("-h", "--help"):
        _print_top_help()
        return 0

    first, rest = args[0], args[1:]

    if first in _SUBCOMMANDS:
        return _SUBCOMMANDS[first](rest)

    if first.startswith("-"):
        _emit_legacy_fallback_warning(first)
        return _SUBCOMMANDS["rpf"](args)

    sys.stderr.write(
        f"mitoribopy: error: unknown subcommand '{first}'. "
        "Known subcommands: align, rpf, rnaseq, all. "
        "Try 'mitoribopy --help'.\n"
    )
    return 2
