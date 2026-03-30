"""MitoRiboPy command-line entrypoint.

Phase I behavior:
- forwards arguments to the existing legacy pipeline parser in ``main.py``
- preserves current command behavior while package refactor is in progress

Examples
--------
Run with legacy arguments:
    mitoribopy -s h -f human.fa --directory input_bed --align stop

Optional compatibility subcommand form:
    mitoribopy run -s h -f human.fa --directory input_bed --align stop
"""

from __future__ import annotations

import sys
from typing import Iterable

from . import __version__
from .pipeline.runner import run_legacy_pipeline


def _normalize_args(argv: Iterable[str]) -> list[str]:
    args = list(argv)
    if args and args[0] == "run":
        args = args[1:]
    if args and args[0] == "--":
        args = args[1:]
    return args


def main(argv: Iterable[str] | None = None) -> int:
    args = _normalize_args(sys.argv[1:] if argv is None else argv)
    if args in (["--version"], ["-V"]):
        print(f"MitoRiboPy {__version__}")
        return 0
    try:
        return run_legacy_pipeline(args)
    except RuntimeError as exc:
        print(f"[mitoribopy] {exc}", file=sys.stderr)
        return 2
