"""MitoRiboPy command-line entrypoint with subcommand dispatch.

Top-level synopsis::

    mitoribopy [--version] [--help]
    mitoribopy align              [options]   # FASTQ -> BAM + BED + read counts
    mitoribopy rpf                [options]   # Ribo-seq analysis from BED / BAM
    mitoribopy rnaseq             [options]   # DE table + rpf -> TE / dTE
    mitoribopy all                [options]   # align + rpf + (optional) rnaseq
    mitoribopy migrate-config     <path>      # rewrite legacy YAML keys
    mitoribopy validate-config    <path>      # schema + cross-stage checks
    mitoribopy validate-reference [options]   # FASTA / annotation sanity checks
    mitoribopy validate-figures   <run_dir>   # figure_qc.tsv producer
    mitoribopy summarize          <run_dir>   # SUMMARY.md / outputs_index.tsv
    mitoribopy benchmark          [options]   # wall / RSS benchmark harness

A subcommand is required. The pre-v0.6.0 fallback that re-routed
``mitoribopy -s h -f ref.fa ...`` to ``mitoribopy rpf ...`` was
removed in v0.6.0 (publication freeze). Use the explicit form::

    mitoribopy rpf -s h.sapiens -f ref.fa ...
"""

from __future__ import annotations

import sys
from typing import Callable, Iterable

from .. import __version__

from . import align as _align
from . import all_ as _all
from . import benchmark as _benchmark
from . import migrate_config as _migrate_config
from . import rnaseq as _rnaseq
from . import rpf as _rpf
from . import summarize as _summarize
from . import validate_config as _validate_config
from . import validate_figures as _validate_figures
from . import validate_reference as _validate_reference

__all__ = ["main", "run_pipeline_cli"]


def __getattr__(name: str):
    """Lazy attribute access for the back-compat ``run_pipeline_cli``
    re-export.

    Eagerly importing ``mitoribopy.pipeline.runner`` at module-init
    time creates a cli<->pipeline circular when callers import
    ``mitoribopy.pipeline.context`` (or any other pipeline submodule)
    before they ever reach the ``cli`` package. Deferring the import
    until the attribute is actually accessed keeps the convenience
    while letting pipeline-only callers (tests, scripts) avoid the
    cli load.
    """
    if name == "run_pipeline_cli":
        from ..pipeline.runner import run_pipeline_cli as _impl
        return _impl
    raise AttributeError(f"module 'mitoribopy.cli' has no attribute {name!r}")


_SUBCOMMANDS: dict[str, Callable[[Iterable[str]], int]] = {
    "align": _align.run,
    "rpf": _rpf.run,
    "rnaseq": _rnaseq.run,
    "all": _all.run,
    "migrate-config": _migrate_config.run,
    "validate-config": _validate_config.run,
    "validate-reference": _validate_reference.run,
    "validate-figures": _validate_figures.run,
    "summarize": _summarize.run,
    "benchmark": _benchmark.run,
}


_SUBCOMMAND_SUMMARIES: list[tuple[str, str]] = [
    ("align", _align.ALIGN_SUBCOMMAND_HELP),
    ("rpf", _rpf.RPF_SUBCOMMAND_HELP),
    ("rnaseq", _rnaseq.RNASEQ_SUBCOMMAND_HELP),
    ("all", _all.ALL_SUBCOMMAND_HELP),
    ("migrate-config", _migrate_config.MIGRATE_CONFIG_HELP),
    ("validate-config", _validate_config.VALIDATE_CONFIG_HELP),
    ("validate-reference", _validate_reference.VALIDATE_REFERENCE_HELP),
    ("validate-figures", _validate_figures.VALIDATE_FIGURES_HELP),
    ("summarize", _summarize.SUMMARIZE_SUBCOMMAND_HELP),
    ("benchmark", _benchmark.BENCHMARK_SUBCOMMAND_HELP),
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
    name_w = max(len(name) for name, _ in _SUBCOMMAND_SUMMARIES) + 2
    rows = "\n".join(
        f"  {name:<{name_w}}{summary}" for name, summary in _SUBCOMMAND_SUMMARIES
    )
    print(_TOP_HELP.format(rows=rows))


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
        sys.stderr.write(
            f"mitoribopy: error: missing subcommand. Got '{first}' as the "
            "first token, but a subcommand is required as of v0.6.0. "
            f"Did you mean 'mitoribopy rpf {first} ...'? "
            "Run 'mitoribopy --help' for the full list.\n"
        )
        return 2

    sys.stderr.write(
        f"mitoribopy: error: unknown subcommand '{first}'. "
        "Known subcommands: " + ", ".join(_SUBCOMMANDS) + ". "
        "Try 'mitoribopy --help'.\n"
    )
    return 2
