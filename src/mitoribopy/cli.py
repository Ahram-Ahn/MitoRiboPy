"""MitoRiboPy command-line entrypoint."""

from __future__ import annotations

import sys
from typing import Iterable

from . import __version__
from .pipeline.runner import run_pipeline_cli


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
        return run_pipeline_cli(args)
    except RuntimeError as exc:
        print(f"[mitoribopy] {exc}", file=sys.stderr)
        return 2
