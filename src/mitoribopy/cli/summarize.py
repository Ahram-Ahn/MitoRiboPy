"""``mitoribopy summarize <run-root>`` (P1.8).

Reads the manifest of a finished run and (re-)emits ``SUMMARY.md`` and
``summary_qc.tsv``. The same path is also auto-invoked by
``mitoribopy all`` after the per-stage work completes so SUMMARY.md is
always produced; the standalone command is for regenerating summaries
on archival runs.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Iterable

from . import common


SUMMARIZE_SUBCOMMAND_HELP = (
    "Regenerate SUMMARY.md and summary_qc.tsv from a finished MitoRiboPy "
    "run by reading the run_manifest.json and per-stage TSV outputs. "
    "Useful for re-rendering summaries on archival runs without "
    "re-executing any pipeline stage."
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mitoribopy summarize",
        description=SUMMARIZE_SUBCOMMAND_HELP,
        formatter_class=common.MitoRiboPyHelpFormatter,
    )
    parser.add_argument(
        "run_root",
        metavar="RUN_DIR",
        help="The output directory of a finished `mitoribopy all` run.",
    )
    parser.add_argument(
        "--manifest",
        default="run_manifest.json",
        metavar="NAME",
        help="Manifest filename within RUN_DIR.",
    )
    return parser


def run(argv: Iterable[str]) -> int:
    args = build_parser().parse_args(list(argv))

    run_root = Path(args.run_root)
    if not run_root.is_dir():
        print(
            f"[mitoribopy summarize] ERROR: {run_root} is not a directory.",
            file=sys.stderr,
        )
        return 2

    manifest_path = run_root / args.manifest
    if not manifest_path.is_file():
        print(
            f"[mitoribopy summarize] ERROR: manifest not found at "
            f"{manifest_path}. Was --manifest set correctly?",
            file=sys.stderr,
        )
        return 2

    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(
            f"[mitoribopy summarize] ERROR: failed to parse manifest "
            f"{manifest_path}: {exc}",
            file=sys.stderr,
        )
        return 2

    from ..summary import build_summary_qc, render_summary_md, write_summary_qc

    qc_rows = build_summary_qc(run_root, manifest)
    qc_path = write_summary_qc(qc_rows, run_root / "summary_qc.tsv")

    md = render_summary_md(
        run_root=run_root,
        manifest=manifest,
        summary_qc_rows=qc_rows,
        warnings=list(manifest.get("warnings") or []),
    )
    md_path = run_root / "SUMMARY.md"
    md_path.write_text(md, encoding="utf-8")

    sys.stderr.write(
        f"[mitoribopy summarize] wrote {qc_path} and {md_path}.\n"
    )
    return 0
