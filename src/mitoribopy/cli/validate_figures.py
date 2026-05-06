"""``mitoribopy validate-figures`` — mechanical QC for plot outputs.

Walks a finished run root, scores every PNG / SVG plot against image-file
checks plus any checks available from a ``.metadata.json`` sidecar
(point counts, label overlap, label outside axes, min font size,
source-data link), and writes ``figure_qc.tsv``. Each plot gets a
``status`` of ``pass``, ``warn``, or ``fail``; failures are
mirrored into ``warnings.tsv`` so they show up in the manifest's
warnings array as well.

Exit codes:

* ``0`` — every plot passed.
* ``1`` — at least one ``warn`` row, no fails (``--strict`` upgrades).
* ``2`` — at least one ``fail`` row.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

from . import common


VALIDATE_FIGURES_HELP = (
    "Mechanically validate every plot under a finished MitoRiboPy run "
    "root: check label / legend / stat-box overlap, label clipping, "
    "point counts vs source TSV, SVG text editability, and PNG dpi. "
    "Writes <RUN_DIR>/figure_qc.tsv. Exit "
    "0 / 1 / 2 (all pass / warn-only / fail; --strict upgrades warn → fail)."
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mitoribopy validate-figures",
        description=VALIDATE_FIGURES_HELP,
        formatter_class=common.MitoRiboPyHelpFormatter,
    )
    parser.add_argument(
        "run_root",
        metavar="RUN_DIR",
        help="The output directory of a finished `mitoribopy all` run.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        default=False,
        help=(
            "Treat warn-level conditions (low PNG dpi, non-editable SVG "
            "text) as fail. Hard contract violations (point-count "
            "mismatch, label-overlap > 0) are always fail regardless."
        ),
    )
    parser.add_argument(
        "--require-png-dpi",
        type=int,
        default=300,
        metavar="DPI",
        help="Minimum acceptable PNG dpi (default 300).",
    )
    parser.add_argument(
        "--out",
        default=None,
        metavar="PATH",
        help=(
            "Override the figure_qc.tsv output path. Defaults to "
            "<RUN_DIR>/figure_qc.tsv."
        ),
    )
    return parser


def run(argv: Iterable[str]) -> int:
    args = build_parser().parse_args(list(argv))

    run_root = Path(args.run_root)
    if not run_root.is_dir():
        sys.stderr.write(
            f"[mitoribopy validate-figures] ERROR: {run_root} is not a "
            "directory.\n"
        )
        return 2

    # Local imports keep matplotlib off the dispatch path until needed.
    from ..io import warnings_log
    from ..plotting.figure_validator import (
        build_figure_qc_rows,
        validate_figures,
        write_figure_qc,
    )

    records = validate_figures(
        run_root,
        strict=args.strict,
        require_png_dpi=args.require_png_dpi,
    )
    rows = build_figure_qc_rows(records, run_root=run_root)

    target = Path(args.out) if args.out else run_root / "figure_qc.tsv"
    figure_qc_path = write_figure_qc(rows, target)

    n_total = len(records)
    n_warn = sum(1 for r in records if r.status == "warn")
    n_fail = sum(1 for r in records if r.status == "fail")
    sys.stderr.write(
        f"[mitoribopy validate-figures] {figure_qc_path}: "
        f"{n_total} plot(s), {n_warn} warn, {n_fail} fail.\n"
    )

    # Mirror every warn / fail into warnings.tsv (and the in-memory
    # warnings_log so a follow-up `mitoribopy summarize` picks them up).
    for r in records:
        if r.status == "pass":
            continue
        warnings_log.record(
            stage="VALIDATE_FIGURES",
            sample_id=None,
            severity="warn" if r.status == "warn" else "error",
            code="FIGURE_QC_FAIL" if r.status == "fail" else "FIGURE_QC_WARN",
            message=(
                f"{r.plot_path.name}: {'; '.join(r.warnings) or 'no detail'}"
            ),
            suggested_action=(
                "Re-render the affected plot(s) and re-run "
                "validate-figures."
            ),
        )
    # Always (re-)flush the run-root warnings.tsv so the new rows land
    # next to the manifest. Non-fatal if the path is read-only.
    try:
        warnings_log.flush_tsv(run_root / "warnings.tsv")
    except OSError:
        pass

    if n_fail:
        return 2
    if n_warn:
        return 1
    return 0
