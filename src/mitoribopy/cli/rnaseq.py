"""``mitoribopy rnaseq`` subcommand (Phase 5 stub).

Phase 2 registers the argument surface and ``--help`` text for this
subcommand. The DE-table + Ribo-seq integration (TE, ΔTE, plots, with an
SHA256 reference-consistency gate) is implemented in Phase 5 under
:mod:`mitoribopy.rnaseq`.
"""

from __future__ import annotations

import argparse
import sys
from typing import Iterable

from . import common


RNASEQ_SUBCOMMAND_HELP = (
    "Integrate a pre-computed differential-expression (DESeq2 / Xtail / "
    "Anota2Seq) table with Ribo-seq outputs to produce TE and ΔTE tables "
    "and plots. (Phase 5; not yet implemented.)"
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mitoribopy rnaseq",
        description=RNASEQ_SUBCOMMAND_HELP,
    )
    common.add_common_arguments(parser)
    planned = parser.add_argument_group("Planned in Phase 5 (not yet active)")
    planned.add_argument(
        "--de-table",
        default=None,
        metavar="PATH",
        help="DESeq2/Xtail/Anota2Seq results table (CSV or TSV; Phase 5).",
    )
    planned.add_argument(
        "--de-format",
        choices=["auto", "deseq2", "xtail", "anota2seq"],
        default="auto",
        help="DE table format; 'auto' detects from column headers (Phase 5).",
    )
    planned.add_argument(
        "--gene-id-convention",
        choices=["ensembl", "refseq", "hgnc", "mt_prefixed", "bare"],
        default=None,
        help=(
            "Gene identifier convention used in the DE table. "
            "REQUIRED in Phase 5 (no default)."
        ),
    )
    planned.add_argument(
        "--ribo-dir",
        default=None,
        metavar="DIR",
        help="Directory produced by a prior 'mitoribopy rpf' run (Phase 5).",
    )
    planned.add_argument(
        "--reference-gtf",
        default=None,
        metavar="PATH",
        help=(
            "GTF used by BOTH Ribo-seq and RNA-seq; hash-verified against "
            "the rpf run's recorded reference (Phase 5)."
        ),
    )
    planned.add_argument(
        "--reference-checksum",
        default=None,
        metavar="SHA256",
        help=(
            "SHA-256 of the shared reference. Used instead of --reference-gtf "
            "when the file is not available locally (Phase 5)."
        ),
    )
    planned.add_argument(
        "--condition-labels",
        nargs="+",
        default=None,
        metavar="LABEL",
        help="Optional condition labels for plotting (Phase 5).",
    )
    planned.add_argument(
        "--output",
        default=None,
        metavar="DIR",
        help="Output directory for TE/ΔTE tables and plots (Phase 5).",
    )
    return parser


def run(argv: Iterable[str]) -> int:
    """Entry point for ``mitoribopy rnaseq <args>``."""
    parser = build_parser()
    args = parser.parse_args(list(argv))
    common.apply_common_arguments(args)

    if args.dry_run:
        return common.emit_dry_run(
            "rnaseq",
            [
                "load DE table via --de-table; detect format via --de-format",
                "hash-verify --reference-gtf / --reference-checksum against "
                "the Ribo-seq run manifest (HARD FAIL on mismatch)",
                "match the 13 mt-mRNA rows to Ribo-seq gene counts; warn if fewer",
                "compute TE per gene per sample = RPF_count / mRNA_count",
                "compute ΔTE between conditions using DE log2FC + Ribo-seq log2FC",
                "emit TE and ΔTE tables and scatter / volcano-style plots",
            ],
        )

    common.warn_stub_subcommand("rnaseq", "Phase 5")
    print(
        "[mitoribopy rnaseq] ERROR: rnaseq is not yet implemented. "
        "Planned for Phase 5. Re-run with --dry-run to inspect the planned steps.",
        file=sys.stderr,
    )
    return 2
