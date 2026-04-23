"""``mitoribopy all`` subcommand (Phase 6 stub).

Phase 2 registers the argument surface and ``--help`` text. The end-to-end
orchestrator that wires ``align`` -> ``rpf`` -> (optional) ``rnaseq``
with a shared config and a ``run_manifest.json`` is implemented in
Phase 6.
"""

from __future__ import annotations

import argparse
import sys
from typing import Iterable

from . import common


ALL_SUBCOMMAND_HELP = (
    "End-to-end orchestrator: align + rpf, optionally + rnaseq when a "
    "DE table is provided. Writes a run_manifest.json with tool versions, "
    "parameters, and input/output hashes. (Phase 6; not yet implemented.)"
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mitoribopy all",
        description=ALL_SUBCOMMAND_HELP,
    )
    common.add_common_arguments(parser)
    planned = parser.add_argument_group("Planned in Phase 6 (not yet active)")
    planned.add_argument(
        "--resume",
        action="store_true",
        default=False,
        help="Skip steps whose output files already exist (Phase 6).",
    )
    planned.add_argument(
        "--manifest",
        default="run_manifest.json",
        metavar="PATH",
        help="Path to the run manifest written at the end of the run (Phase 6).",
    )
    planned.add_argument(
        "--de-table",
        default=None,
        metavar="PATH",
        help="Optional DE table; when given, rnaseq integration is run last (Phase 6).",
    )
    planned.add_argument(
        "--output",
        default=None,
        metavar="DIR",
        help="Shared output directory for align + rpf + rnaseq results (Phase 6).",
    )
    return parser


def run(argv: Iterable[str]) -> int:
    """Entry point for ``mitoribopy all <args>``."""
    parser = build_parser()
    args = parser.parse_args(list(argv))
    common.apply_common_arguments(args)

    if args.dry_run:
        return common.emit_dry_run(
            "all",
            [
                "resolve shared config for align + rpf + rnaseq",
                "run align (FASTQ -> BAM + BED + read counts)",
                "run rpf (BED -> offsets, translation profile, coverage plots)",
                "if --de-table given: run rnaseq (TE, ΔTE, plots) with "
                "reference-consistency gate",
                "write run_manifest.json with tool versions and input/output hashes",
            ],
        )

    common.warn_stub_subcommand("all", "Phase 6")
    print(
        "[mitoribopy all] ERROR: all is not yet implemented. Planned for "
        "Phase 6. Re-run with --dry-run to inspect the planned steps.",
        file=sys.stderr,
    )
    return 2
