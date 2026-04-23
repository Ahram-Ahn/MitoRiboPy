"""``mitoribopy align`` subcommand (Phase 3 stub).

Phase 2 registers the argument surface and ``--help`` text for this
subcommand so users see it in the top-level ``mitoribopy --help`` output.
The actual FASTQ -> BED preprocessing pipeline is implemented in Phase 3
under :mod:`mitoribopy.align`.

Invoking this subcommand currently emits a warning, prints the planned
steps, and exits with status 2.
"""

from __future__ import annotations

import argparse
from typing import Iterable

from . import common


ALIGN_SUBCOMMAND_HELP = (
    "Preprocess FASTQ inputs through trimming, contaminant subtraction, "
    "mt-genome alignment, MAPQ filtering, UMI deduplication, and produce "
    "BAM + BED + read-counts outputs. (Phase 3; not yet implemented.)"
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mitoribopy align",
        description=ALIGN_SUBCOMMAND_HELP,
    )
    common.add_common_arguments(parser)
    # Phase 3 flag surface is designed here so users and CI scripts can
    # depend on stable flag names before the implementation lands.
    planned = parser.add_argument_group("Planned in Phase 3 (not yet active)")
    planned.add_argument(
        "--fastq-dir",
        default=None,
        metavar="DIR",
        help="Directory of input FASTQ files (Phase 3).",
    )
    planned.add_argument(
        "--adapter",
        default="TGGAATTCTCGGGTGCCAAGG",
        metavar="SEQ",
        help=(
            "3' adapter sequence for cutadapt "
            "(default: Illumina TruSeq small-RNA; Phase 3)."
        ),
    )
    planned.add_argument(
        "--umi-length",
        type=int,
        default=0,
        metavar="N",
        help="UMI length in nt; 0 disables UMI handling (Phase 3).",
    )
    planned.add_argument(
        "--min-length",
        type=int,
        default=15,
        metavar="NT",
        help="Minimum read length kept after trimming (mt-RPF default 15; Phase 3).",
    )
    planned.add_argument(
        "--max-length",
        type=int,
        default=45,
        metavar="NT",
        help="Maximum read length kept after trimming (mt-RPF default 45; Phase 3).",
    )
    planned.add_argument(
        "--contam-index",
        default=None,
        metavar="BT2_PREFIX",
        help="bowtie2 index prefix of rRNA/tRNA contaminants (Phase 3).",
    )
    planned.add_argument(
        "--mt-index",
        default=None,
        metavar="BT2_PREFIX",
        help="bowtie2 index prefix of the mt-genome or transcriptome (Phase 3).",
    )
    planned.add_argument(
        "--mapq",
        type=int,
        default=10,
        metavar="Q",
        help="MAPQ threshold for samtools view (default 10; Phase 3).",
    )
    planned.add_argument(
        "--output",
        default=None,
        metavar="DIR",
        help="Output directory for BAM/BED/read-counts (Phase 3).",
    )
    return parser


def run(argv: Iterable[str]) -> int:
    """Entry point for ``mitoribopy align <args>``."""
    parser = build_parser()
    args = parser.parse_args(list(argv))
    common.apply_common_arguments(args)

    if args.dry_run:
        return common.emit_dry_run(
            "align",
            [
                "check external tools (cutadapt, bowtie2, samtools, umi_tools) on PATH",
                "cutadapt: trim --adapter, filter --min-length/--max-length, "
                "extract --umi-length UMIs into read name",
                "bowtie2: subtract --contam-index contaminants",
                "bowtie2: align to --mt-index (end-to-end --very-sensitive -L 18)",
                "samtools view -q --mapq for MAPQ filtering",
                "umi_tools dedup (only if --umi-length > 0)",
                "generate per-sample read_counts table",
                "convert BAM -> BED for downstream rpf consumption",
            ],
        )

    common.warn_stub_subcommand("align", "Phase 3")
    print(
        "[mitoribopy align] ERROR: align is not yet implemented. "
        "Planned for Phase 3. Re-run with --dry-run to inspect the planned steps.",
        file=__import__("sys").stderr,
    )
    return 2
