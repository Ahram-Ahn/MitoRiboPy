"""``mitoribopy align`` subcommand - FASTQ -> BAM + BED + read counts.

End-to-end orchestrator for Phase 3. For every input FASTQ:

1. cutadapt trim (kit-aware; optional UMI extraction)
2. bowtie2 contaminant subtraction (user-supplied index; fails loudly if missing)
3. bowtie2 mt-transcriptome alignment (Path A; end-to-end --very-sensitive -L 18)
4. MAPQ filter (pysam; default 10 to suppress NUMT cross-talk)
5. deduplication per Approach D (umi_tools / skip / mark-duplicates)
6. BAM -> BED6 (strand-aware)

Writes ``read_counts.tsv`` and ``run_settings.json`` at the run root for
provenance. The ``--dry-run`` path prints the resolved settings and
planned actions without touching any files.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Iterable

from .. import __version__
from ..align import (
    adapter_detect,
    align as align_step,
    bam_utils,
    contam as contam_step,
    dedup as dedup_step,
    read_counts as counts_step,
    tool_check,
    trim as trim_step,
)
from ..align._types import KIT_PRESETS, SampleCounts
from ..console import log_info, log_warning
from . import common


ALIGN_SUBCOMMAND_HELP = (
    "Preprocess FASTQ inputs: cutadapt trim + bowtie2 contaminant subtraction + "
    "bowtie2 mt-transcriptome alignment + MAPQ filter + dedup + BAM->BED6. "
    "Produces drop-in inputs for 'mitoribopy rpf'."
)


_FASTQ_GLOB_PATTERNS = ("*.fq.gz", "*.fastq.gz", "*.fq", "*.fastq")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mitoribopy align",
        description=ALIGN_SUBCOMMAND_HELP,
    )
    common.add_common_arguments(parser)

    inputs = parser.add_argument_group("Inputs")
    inputs.add_argument(
        "--fastq-dir",
        default=None,
        metavar="DIR",
        help=(
            "Directory containing input FASTQ files (*.fq, *.fq.gz, "
            "*.fastq, *.fastq.gz). Pass --fastq-dir OR --fastq (or both)."
        ),
    )
    inputs.add_argument(
        "--fastq",
        action="append",
        default=None,
        metavar="PATH",
        help="Individual FASTQ input; repeatable.",
    )
    inputs.add_argument(
        "--contam-index",
        default=None,
        required=False,
        metavar="BT2_PREFIX",
        help=(
            "bowtie2 index prefix of the contaminant panel (rRNA, tRNA, "
            "any nuclear spike-ins to subtract). Required for non-dry-run "
            "invocations. Build with 'bowtie2-build contaminants.fa "
            "<prefix>'."
        ),
    )
    inputs.add_argument(
        "--mt-index",
        default=None,
        required=False,
        metavar="BT2_PREFIX",
        help=(
            "bowtie2 index prefix of the mt-transcriptome (one FASTA "
            "record per mt-mRNA, header matching annotation "
            "sequence_name). Required for non-dry-run invocations."
        ),
    )
    inputs.add_argument(
        "--output",
        default=None,
        required=False,
        metavar="DIR",
        help="Output directory for BAM/BED/read_counts (required).",
    )

    library = parser.add_argument_group("Library prep")
    library.add_argument(
        "--kit-preset",
        choices=sorted(KIT_PRESETS),
        default="custom",
        help=(
            "Library-prep kit. 'custom' (default) forces --adapter to be "
            "supplied explicitly to prevent silently wrong adapter "
            "defaults from shifting the mt-RPF length distribution."
        ),
    )
    library.add_argument(
        "--adapter",
        default=None,
        metavar="SEQ",
        help=(
            "3' adapter sequence. Overrides the kit preset's adapter. "
            "Required when --kit-preset custom."
        ),
    )
    library.add_argument(
        "--umi-length",
        type=int,
        default=None,
        metavar="N",
        help="UMI length in nt. Overrides the kit preset's default.",
    )
    library.add_argument(
        "--umi-position",
        choices=["5p", "3p"],
        default=None,
        help="UMI position within the insert (overrides kit preset).",
    )
    library.add_argument(
        "--adapter-detection",
        choices=["auto", "off", "strict"],
        default="auto",
        metavar="MODE",
        help=(
            "Sanity-check the supplied --kit-preset against the data. "
            "'auto' (default): scan the head of the first FASTQ; pick "
            "the matching preset when --kit-preset is left at 'custom', "
            "warn on mismatch otherwise. 'strict': always scan and "
            "HARD-FAIL on any disagreement (best for batch / CI runs). "
            "'off': trust --kit-preset / --adapter without inspection."
        ),
    )
    library.add_argument(
        "--library-strandedness",
        choices=["forward", "reverse", "unstranded"],
        default="forward",
        help=(
            "Library strandedness. 'forward' (default, NEBNext/TruSeq "
            "small-RNA kits) enforces --norc at alignment time; 'reverse' "
            "enforces --nofw; 'unstranded' leaves bowtie2 permissive."
        ),
    )
    library.add_argument(
        "--min-length",
        type=int,
        default=15,
        metavar="NT",
        help="Minimum read length kept after trimming (mt-RPF default 15).",
    )
    library.add_argument(
        "--max-length",
        type=int,
        default=45,
        metavar="NT",
        help="Maximum read length kept after trimming (mt-RPF default 45).",
    )
    library.add_argument(
        "--quality",
        type=int,
        default=20,
        metavar="Q",
        help="cutadapt -q Phred+33 3' quality trim threshold.",
    )

    alignment = parser.add_argument_group("Alignment")
    alignment.add_argument(
        "--mapq",
        type=int,
        default=10,
        metavar="Q",
        help=(
            "MAPQ threshold for the post-alignment filter (default 10). "
            "The main reason for this filter is NUMT cross-talk "
            "suppression."
        ),
    )
    alignment.add_argument(
        "--seed",
        type=int,
        default=42,
        metavar="N",
        help="bowtie2 --seed value (deterministic output).",
    )

    dedup = parser.add_argument_group("Deduplication")
    dedup.add_argument(
        "--dedup-strategy",
        choices=["auto", "umi-tools", "skip", "mark-duplicates"],
        default="auto",
        help=(
            "'auto' (default) -> umi-tools when UMIs are present, else "
            "skip. 'mark-duplicates' is coordinate-only and destroys "
            "codon-occupancy signal on mt-Ribo-seq; gated behind the "
            f"{dedup_step.CONFIRM_MARK_DUPLICATES_FLAG} confirmation flag."
        ),
    )
    dedup.add_argument(
        "--umi-dedup-method",
        choices=["unique", "percentile", "cluster", "adjacency", "directional"],
        default="unique",
        help=(
            "umi_tools --method. 'unique' (default) collapses only on "
            "exact coord+UMI match; other methods may over-collapse in "
            "low-complexity mt regions."
        ),
    )
    dedup.add_argument(
        dedup_step.CONFIRM_MARK_DUPLICATES_FLAG,
        dest="confirm_mark_duplicates",
        action="store_true",
        default=False,
        help=(
            "Required confirmation to opt into --dedup-strategy "
            "mark-duplicates. Only pass this if you have an independent "
            "reason (carrier nuclear RNA, >22 PCR cycles, very low input)."
        ),
    )

    return parser


# ---------------------------------------------------------------------------
# FASTQ enumeration + sample naming
# ---------------------------------------------------------------------------


def _enumerate_fastqs(
    fastq_dir: str | None, extra_fastqs: list[str] | None
) -> list[Path]:
    paths: list[Path] = []
    if fastq_dir:
        root = Path(fastq_dir)
        if not root.is_dir():
            raise FileNotFoundError(f"--fastq-dir is not a directory: {root}")
        for pattern in _FASTQ_GLOB_PATTERNS:
            paths.extend(sorted(root.glob(pattern)))
    if extra_fastqs:
        paths.extend(Path(p) for p in extra_fastqs)

    # De-duplicate while preserving order.
    seen: set[str] = set()
    unique: list[Path] = []
    for path in paths:
        resolved = str(path.resolve())
        if resolved in seen:
            continue
        seen.add(resolved)
        unique.append(path)
    return unique


def _sample_name(fastq_path: Path) -> str:
    """Strip .fq / .fq.gz / .fastq / .fastq.gz from the filename."""
    name = fastq_path.name
    for suffix in (".fq.gz", ".fastq.gz", ".fq", ".fastq"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return fastq_path.stem


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


def _planned_actions(
    resolved_kit,
    resolved_dedup: str,
    strandedness: str,
    mapq: int,
    min_length: int,
    max_length: int,
    n_inputs: int,
) -> list[str]:
    return [
        f"validate external tools on PATH (cutadapt, bowtie2; plus "
        f"umi_tools / picard as needed for dedup='{resolved_dedup}')",
        (
            f"cutadapt trim: kit={resolved_kit.kit}, "
            f"adapter={resolved_kit.adapter}, umi_length="
            f"{resolved_kit.umi_length} ({resolved_kit.umi_position}), "
            f"min={min_length} max={max_length} nt"
        ),
        "bowtie2 contam subtract (--contam-index); unaligned reads pass "
        f"through. strand={strandedness}",
        "bowtie2 mt-transcriptome align (--mt-index) --end-to-end "
        f"--very-sensitive -L 18; strand={strandedness}",
        f"MAPQ filter at q >= {mapq} (pysam; NUMT suppression)",
        f"dedup strategy: {resolved_dedup}",
        "BAM -> BED6 (strand-aware)",
        f"write read_counts.tsv and run_settings.json for {n_inputs} sample(s)",
    ]


def _write_run_settings(
    output_dir: Path,
    *,
    args: argparse.Namespace,
    resolved_kit,
    resolved_dedup: str,
    tool_versions: dict[str, str | None],
    inputs: list[Path],
) -> Path:
    settings = {
        "subcommand": "align",
        "mitoribopy_version": __version__,
        "kit_preset": resolved_kit.kit,
        "adapter": resolved_kit.adapter,
        "adapter_detection_mode": args.adapter_detection,
        "umi_length": resolved_kit.umi_length,
        "umi_position": resolved_kit.umi_position,
        "library_strandedness": args.library_strandedness,
        "min_length": args.min_length,
        "max_length": args.max_length,
        "quality": args.quality,
        "mapq_threshold": args.mapq,
        "seed": args.seed,
        "dedup_strategy": resolved_dedup,
        "umi_dedup_method": args.umi_dedup_method,
        "confirm_mark_duplicates": args.confirm_mark_duplicates,
        "contam_index": args.contam_index,
        "mt_index": args.mt_index,
        "threads": getattr(args, "threads", None),
        "log_level": getattr(args, "log_level", "INFO"),
        "tool_versions": tool_versions,
        "inputs": [str(p) for p in inputs],
    }
    path = output_dir / "run_settings.json"
    path.write_text(json.dumps(settings, indent=2, sort_keys=True), encoding="utf-8")
    return path


def _process_one_sample(
    *,
    fastq_in: Path,
    sample: str,
    output_dir: Path,
    args: argparse.Namespace,
    resolved_kit,
    resolved_dedup: str,
) -> SampleCounts:
    trimmed_dir = output_dir / "trimmed"
    contam_dir = output_dir / "contam_filtered"
    aligned_dir = output_dir / "aligned"
    deduped_dir = output_dir / "deduped"
    bed_dir = output_dir / "bed"

    trimmed_fq = trimmed_dir / f"{sample}.trimmed.fq.gz"
    cutadapt_log = trimmed_dir / f"{sample}.cutadapt.json"
    nocontam_fq = contam_dir / f"{sample}.nocontam.fq.gz"
    aligned_bam = aligned_dir / f"{sample}.bam"
    mapq_bam = aligned_dir / f"{sample}.mapq.bam"
    dedup_bam = deduped_dir / f"{sample}.dedup.bam"
    dedup_log = deduped_dir / f"{sample}.dedup.log"
    bed_out = bed_dir / f"{sample}.bed"

    threads = getattr(args, "threads", None) or 1

    cutadapt = trim_step.run_cutadapt(
        fastq_in=fastq_in,
        fastq_out=trimmed_fq,
        resolved=resolved_kit,
        min_length=args.min_length,
        max_length=args.max_length,
        quality=args.quality,
        threads=threads,
        log_json=cutadapt_log,
    )

    contam = contam_step.subtract_contaminants(
        fastq_in=trimmed_fq,
        contam_index=Path(args.contam_index),
        fastq_out_unaligned=nocontam_fq,
        strandedness=args.library_strandedness,
        threads=threads,
        seed=args.seed,
    )

    alignment = align_step.align_mt(
        fastq_in=nocontam_fq,
        mt_index=Path(args.mt_index),
        bam_out=aligned_bam,
        strandedness=args.library_strandedness,
        threads=threads,
        seed=args.seed,
    )

    mapq_kept = bam_utils.filter_bam_mapq(
        bam_in=aligned_bam,
        bam_out=mapq_bam,
        mapq_threshold=args.mapq,
    )

    dedup = dedup_step.run_dedup(
        strategy=args.dedup_strategy,
        umi_length=resolved_kit.umi_length,
        confirm_mark_duplicates=args.confirm_mark_duplicates,
        bam_in=mapq_bam,
        bam_out=dedup_bam,
        log_path=dedup_log,
        umi_tools_method=args.umi_dedup_method,
    )

    bam_utils.bam_to_bed6(bam_in=dedup.bam_path, bed_out=bed_out)

    return counts_step.assemble_sample_counts(
        sample=sample,
        cutadapt=cutadapt,
        contam=contam,
        align=alignment,
        mapq_count=mapq_kept,
        dedup=dedup,
    )


def _apply_adapter_detection(
    args: argparse.Namespace,
    first_fastq: Path,
) -> int | None:
    """Scan ``first_fastq`` and apply ``--adapter-detection`` policy to ``args``.

    May mutate ``args.kit_preset`` when the user left it at ``custom`` and
    detection succeeded. Returns ``None`` to continue, or a non-zero exit
    code to abort (strict-mode hard fail).
    """
    if not first_fastq.exists():
        log_warning(
            "ADAPTER",
            f"Cannot scan {first_fastq} (file not found); skipping detection.",
        )
        return None

    try:
        detection = adapter_detect.detect_adapter(first_fastq)
    except OSError as exc:
        log_warning(
            "ADAPTER",
            f"Could not scan {first_fastq.name}: {exc}; skipping detection.",
        )
        return None

    user_specified_preset = args.kit_preset != "custom" or args.adapter is not None
    detected = detection.preset_name
    rates_str = adapter_detect.format_per_kit_rates(detection.per_kit_rates)

    if args.adapter_detection == "strict":
        if detected is None:
            print(
                "[mitoribopy align] ERROR: --adapter-detection strict could not "
                f"identify a known adapter in {first_fastq.name} "
                f"({detection.n_reads_scanned} reads scanned). "
                f"Per-kit match rates: {rates_str}. "
                "Re-run with --adapter-detection off to bypass.",
                file=sys.stderr,
            )
            return 2
        if user_specified_preset and detected != args.kit_preset:
            print(
                "[mitoribopy align] ERROR: --adapter-detection strict: "
                f"data looks like '{detected}' "
                f"({detection.match_rate * 100:.1f}% match) but you supplied "
                f"--kit-preset {args.kit_preset}. "
                f"Per-kit match rates: {rates_str}. "
                f"Re-run with --kit-preset {detected} "
                "(or --adapter-detection off).",
                file=sys.stderr,
            )
            return 2
        if not user_specified_preset:
            args.kit_preset = detected
            log_info(
                "ADAPTER",
                f"Auto-selected --kit-preset {detected} "
                f"({detection.match_rate * 100:.1f}% match, "
                f"{detection.n_reads_scanned} reads scanned).",
            )
        return None

    # auto mode (off mode never enters this function)
    if detected is None:
        if not user_specified_preset:
            log_warning(
                "ADAPTER",
                f"No known adapter detected in {first_fastq.name} "
                "and --kit-preset is 'custom' with no --adapter. "
                f"Per-kit match rates: {rates_str}. "
                "Pass --kit-preset / --adapter explicitly or use "
                "--adapter-detection off.",
            )
        return None

    if not user_specified_preset:
        args.kit_preset = detected
        log_info(
            "ADAPTER",
            f"Auto-selected --kit-preset {detected} "
            f"({detection.match_rate * 100:.1f}% match, "
            f"{detection.n_reads_scanned} reads scanned).",
        )
    elif detected != args.kit_preset:
        log_warning(
            "ADAPTER",
            f"--kit-preset {args.kit_preset} disagrees with adapter scan "
            f"(detected '{detected}' at {detection.match_rate * 100:.1f}%). "
            f"Per-kit match rates: {rates_str}. "
            "Re-run with --adapter-detection strict to fail-fast on this.",
        )
    return None


def run(argv: Iterable[str]) -> int:
    """Entry point for ``mitoribopy align <args>``."""
    parser = build_parser()
    args = parser.parse_args(list(argv))
    common.apply_common_arguments(args)

    # Enumerate FASTQs first so adapter detection can sanity-check
    # --kit-preset against the actual data before any expensive resolution
    # or subprocess work.
    inputs: list[Path]
    try:
        inputs = _enumerate_fastqs(args.fastq_dir, args.fastq)
    except FileNotFoundError as exc:
        print(f"[mitoribopy align] ERROR: {exc}", file=sys.stderr)
        return 2

    # Adapter sanity check. Skipped in dry-run (no real data inspection)
    # and when no FASTQ paths exist yet (validation will fail loudly later).
    if args.adapter_detection != "off" and not args.dry_run and inputs:
        bail = _apply_adapter_detection(args, inputs[0])
        if bail is not None:
            return bail

    # Resolve configuration up front so invalid combinations fail before
    # any subprocess runs.
    try:
        resolved_kit = trim_step.resolve_kit_settings(
            args.kit_preset,
            adapter=args.adapter,
            umi_length=args.umi_length,
            umi_position=args.umi_position,
        )
    except (KeyError, ValueError) as exc:
        print(f"[mitoribopy align] ERROR: {exc}", file=sys.stderr)
        return 2

    try:
        resolved_dedup = dedup_step.resolve_dedup_strategy(
            args.dedup_strategy,
            umi_length=resolved_kit.umi_length,
            confirm_mark_duplicates=args.confirm_mark_duplicates,
        )
    except ValueError as exc:
        print(f"[mitoribopy align] ERROR: {exc}", file=sys.stderr)
        return 2

    if args.dry_run:
        for line in _planned_actions(
            resolved_kit,
            resolved_dedup,
            strandedness=args.library_strandedness,
            mapq=args.mapq,
            min_length=args.min_length,
            max_length=args.max_length,
            n_inputs=len(inputs),
        ):
            pass
        return common.emit_dry_run(
            "align",
            _planned_actions(
                resolved_kit,
                resolved_dedup,
                strandedness=args.library_strandedness,
                mapq=args.mapq,
                min_length=args.min_length,
                max_length=args.max_length,
                n_inputs=len(inputs),
            ),
        )

    # Validate non-dry-run required inputs.
    missing: list[str] = []
    if not args.output:
        missing.append("--output")
    if not args.contam_index:
        missing.append("--contam-index")
    if not args.mt_index:
        missing.append("--mt-index")
    if not inputs:
        missing.append("--fastq-dir / --fastq")
    if missing:
        print(
            "[mitoribopy align] ERROR: missing required argument(s): "
            + ", ".join(missing),
            file=sys.stderr,
        )
        return 2

    # Tool check.
    required = ["cutadapt", "bowtie2", "bowtie2-build"]
    optional: list[str] = []
    if resolved_dedup == "umi-tools":
        required.append("umi_tools")
    if resolved_dedup == "mark-duplicates":
        required.append("picard")
    for tool in ("fastqc",):
        optional.append(tool)
    try:
        tools = tool_check.ensure_tools_available(
            required=required, optional=optional
        )
    except tool_check.ToolNotFoundError as exc:
        print(f"[mitoribopy align] ERROR: {exc}", file=sys.stderr)
        return 2

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    tool_versions = {
        name: (info.version if info is not None else None)
        for name, info in tools.items()
    }
    _write_run_settings(
        output_dir,
        args=args,
        resolved_kit=resolved_kit,
        resolved_dedup=resolved_dedup,
        tool_versions=tool_versions,
        inputs=inputs,
    )

    rows: list[SampleCounts] = []
    for fastq_in in inputs:
        sample = _sample_name(fastq_in)
        counts = _process_one_sample(
            fastq_in=fastq_in,
            sample=sample,
            output_dir=output_dir,
            args=args,
            resolved_kit=resolved_kit,
            resolved_dedup=resolved_dedup,
        )
        rows.append(counts)

    # Deterministic sort by sample name.
    rows.sort(key=lambda row: row.sample)
    counts_step.write_read_counts_table(rows, output_dir / "read_counts.tsv")
    return 0
