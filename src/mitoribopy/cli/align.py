"""``mitoribopy align`` subcommand - FASTQ -> BAM + BED + read counts.

End-to-end per-sample orchestrator. For every input FASTQ:

1. cutadapt trim (kit-aware; optional UMI extraction)
2. bowtie2 contaminant subtraction (user-supplied index; fails loudly if missing)
3. bowtie2 mt-transcriptome alignment (end-to-end --very-sensitive -L 18)
4. MAPQ filter (pysam; default 10 to suppress NUMT cross-talk)
5. deduplication (umi_coordinate / skip)
6. BAM -> BED6 (strand-aware)

Writes ``read_counts.tsv`` and ``run_settings.json`` at the run root for
provenance. The ``--dry-run`` path prints the resolved settings and
planned actions without touching any files.
"""

from __future__ import annotations

import argparse
import json
import sys
import threading
from concurrent.futures import FIRST_EXCEPTION, Future, ThreadPoolExecutor, wait
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
from ..align._types import (
    ResolvedKit,
    SampleCounts,
    SampleOverride,
)
from ..align.sample_resolve import (
    SampleResolution,
    SampleResolutionError,
    read_sample_overrides_tsv,
    required_dedup_tools,
    resolution_summary_lines,
    resolve_sample_resolutions,
    write_kit_resolution_tsv,
)
from ..console import (
    configure_file_logging,
    log_error,
    log_info,
    log_progress,
    log_warning,
)
from ..progress import (
    SampleCounter,
    StageTimings,
    Stopwatch,
    format_duration,
    render_summary_lines,
    stage_timer,
)
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
        formatter_class=common.MitoRiboPyHelpFormatter,
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
        "--adapter",
        default=None,
        metavar="SEQ",
        help=(
            "3' adapter sequence. By default the pipeline auto-detects "
            "the adapter from the head of each FASTQ; pass --adapter "
            "<SEQ> when detection cannot identify your library or when "
            "you want to pin a specific sequence. Mutually exclusive "
            "with --pretrimmed."
        ),
    )
    library.add_argument(
        "--pretrimmed",
        action="store_true",
        default=False,
        help=(
            "Declare that the input FASTQ has already been adapter-"
            "trimmed (e.g. SRA-deposited data). cutadapt skips the -a "
            "flag and only enforces length and quality filtering. "
            "Auto-detection also infers this when no known adapter "
            "signature is present; pass this flag to assert it "
            "explicitly. Mutually exclusive with --adapter."
        ),
    )
    library.add_argument(
        "--umi-length",
        type=int,
        default=None,
        metavar="N",
        help=(
            "UMI length in nt. Overrides the kit preset's default. For "
            "--umi-position=both this MUST equal --umi-length-5p + "
            "--umi-length-3p (it is the canonical concatenated QNAME "
            "UMI length umi_tools dedups on)."
        ),
    )
    library.add_argument(
        "--umi-position",
        choices=["5p", "3p", "both"],
        default=None,
        help=(
            "UMI position within the insert (overrides kit preset). "
            "'5p' / '3p' are single-end UMIs; 'both' is a dual-end UMI "
            "library (e.g. xGen Duplex, Twist) — supply --umi-length-5p "
            "and --umi-length-3p in that mode."
        ),
    )
    library.add_argument(
        "--umi-length-5p",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Per-end 5' UMI length in nt. Used only when "
            "--umi-position=both; ignored otherwise."
        ),
    )
    library.add_argument(
        "--umi-length-3p",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Per-end 3' UMI length in nt. Used only when "
            "--umi-position=both; ignored otherwise."
        ),
    )
    library.add_argument(
        "--sample-overrides",
        default=None,
        metavar="TSV",
        help=(
            "Path to a TSV with per-sample overrides for "
            "adapter / pretrimmed / umi_length / umi_position / "
            "dedup_strategy. Required header columns: 'sample' plus at "
            "least one of the override columns. The 'sample' value must "
            "match the FASTQ basename with the .fq[.gz] / .fastq[.gz] "
            "suffix removed. Empty cells (or NA / None / null) fall "
            "through to the global CLI default for that field, so a "
            "single sample can override only its UMI without restating "
            "the rest. Useful for mixed-UMI / mixed-adapter batches."
        ),
    )
    library.add_argument(
        "--keep-intermediates",
        action="store_true",
        default=False,
        help=(
            "Keep the per-step intermediate files (trimmed FASTQ, "
            "contam-filtered FASTQ, pre-MAPQ BAM). By default these are "
            "deleted as soon as the next step has consumed them, since "
            "they are large, regenerable, and not needed by any "
            "downstream stage. Pass this flag when debugging a sample "
            "or comparing per-step intermediate counts."
        ),
    )
    library.add_argument(
        "--tmpdir",
        default=None,
        help=(
            "Optional override for the directory used for per-step "
            "scratch files (trimmed FASTQ, contam-filtered FASTQ, "
            "intermediate BAMs). Defaults to a subdirectory of "
            "--output. Set this to a fast local SSD when running on a "
            "cluster with slow shared storage, or to a pre-mounted "
            "tmpfs to avoid hitting disk altogether for short runs."
        ),
    )
    library.add_argument(
        "--allow-count-invariant-warning",
        action="store_true",
        default=False,
        help=(
            "DEVELOPER / DEBUG ONLY. Demote read_counts.tsv invariant "
            "violations from errors to warnings. The default is to fail "
            "the run on any violation, since a real violation indicates "
            "a bug somewhere upstream. NEVER use for a publication run "
            "(`mitoribopy all --strict` will reject this flag too)."
        ),
    )
    library.add_argument(
        "--strict-publication-mode",
        action="store_true",
        default=False,
        help=(
            "Strict Mode for align. Reject runs that rely on inferred-"
            "rather-than-declared metadata: inferred-pretrimmed kits and "
            "ambiguous adapter detection (confidence margin < 0.10). "
            "Use this to make metadata/config problems fail loudly."
        ),
    )
    library.add_argument(
        "--resume",
        action="store_true",
        default=False,
        help=(
            "Skip samples that have already completed in a previous "
            "invocation against this --output directory. Each completed "
            "sample writes a small JSON file under "
            "<output>/.sample_done/; on resume, samples whose JSON is "
            "present and parses are reloaded instead of re-run. Use this "
            "after a crash or kill mid-batch to avoid redoing the "
            "samples that already finished. The orchestrator "
            "('mitoribopy all --resume') sets this automatically when "
            "the align stage's read_counts.tsv is missing."
        ),
    )
    library.add_argument(
        "--adapter-detection",
        choices=["auto", "off", "strict"],
        default="auto",
        metavar="MODE",
        help=(
            "Per-sample adapter detection policy. 'auto' (default) scans "
            "every input FASTQ and picks the matching adapter family per "
            "sample; samples whose scan fails fall back to --adapter "
            "when supplied, or to 'pretrimmed' when no fallback is set "
            "and the data looks already-trimmed. 'strict' scans and "
            "HARD-FAILS on any sample whose detected adapter conflicts "
            "with an explicit --adapter or where no adapter can be "
            "identified. 'off' skips the scan and trusts --adapter / "
            "--pretrimmed for every sample (one of those is required)."
        ),
    )
    library.add_argument(
        "--adapter-detect-reads",
        type=int,
        default=5000,
        metavar="N",
        help=(
            "Number of FASTQ reads to scan per sample during adapter "
            "auto-detection. Increase for noisy libraries where the first "
            "5000 reads have unusual adapter distributions; decrease for a "
            "faster pre-flight pass on cleaner data."
        ),
    )
    library.add_argument(
        "--adapter-detect-min-rate",
        type=float,
        default=0.30,
        metavar="FRAC",
        help=(
            "Minimum fraction of scanned reads that must contain an adapter "
            "prefix for the kit to be considered detected. Lower for "
            "sparsely-adapted libraries (e.g. 0.10); raise for stricter "
            "calls."
        ),
    )
    library.add_argument(
        "--adapter-detect-min-len",
        type=int,
        default=12,
        metavar="N",
        help=(
            "Adapter prefix length used as the search needle (nt). Default "
            "12. Lower (e.g. 8) tolerates noisy adapter regions; raise "
            "(e.g. 16) for stricter matches."
        ),
    )
    library.add_argument(
        "--adapter-detect-pretrimmed-threshold",
        type=float,
        default=0.05,
        metavar="FRAC",
        help=(
            "When EVERY kit's match rate is at or below this value, the "
            "FASTQ is classified as already adapter-trimmed and resolved "
            "to the 'pretrimmed' kit (cutadapt skips the -a flag). Default "
            "0.05 (5%%)."
        ),
    )
    library.add_argument(
        "--no-pretrimmed-inference",
        dest="allow_pretrimmed_inference",
        action="store_false",
        default=True,
        help=(
            "Disable the auto-fallback to 'pretrimmed' when adapter "
            "detection finds no known kit. With this flag, detection "
            "failure with no --adapter / --pretrimmed fallback raises "
            "an error instead."
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
        # Public choices describe the statistical operation
        # (`umi_coordinate`) rather than the implementation tool.
        # `umi-tools` / `umi_tools` continue to parse as deprecated
        # aliases; they normalise to `umi_coordinate` at parse time.
        choices=[
            "auto",
            "umi_coordinate",
            "umi-tools",
            "umi_tools",
            "skip",
        ],
        default="auto",
        type=lambda v: (
            "umi_coordinate"
            if str(v).lower() in {"umi-tools", "umi_tools"}
            else str(v)
        ),
        help=(
            "Canonical: auto | umi_coordinate | skip. 'auto' (default) "
            "-> umi_coordinate when UMIs are present, else skip. "
            "umi_coordinate collapses reads on (coordinate, UMI); the "
            "implementation calls into umi_tools but the statistical "
            "operation is coordinate+UMI dedup. The legacy aliases "
            "'umi-tools' / 'umi_tools' are accepted and rewritten to "
            "'umi_coordinate' (canonical_config.yaml records the "
            "canonical name). Coordinate-only mark-duplicates is not "
            "supported: it destroys codon-occupancy signal on low-"
            "complexity mt-Ribo-seq libraries."
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

    execution = parser.add_argument_group("Execution")
    execution.add_argument(
        "--max-parallel-samples",
        type=str,
        default=None,
        metavar="N|auto",
        help=(
            "Number of samples to process concurrently. Accepts an integer "
            "or 'auto' (default). 'auto' picks min(n_samples, "
            "threads/min_per_sample, memory_gb/est_per_sample) so a "
            "modern CPU is not left idle on a multi-sample run. With "
            "--threads T, each worker uses max(1, T // N) tool threads so "
            "total CPU use stays around T. Pass --max-parallel-samples 1 "
            "(or --single-sample-mode) for legacy serial behaviour. "
            "Per-sample work (cutadapt + bowtie2 + dedup + BAM->BED) is "
            "embarrassingly parallel; the joint 'mitoribopy rpf' stage is "
            "unaffected."
        ),
    )
    execution.add_argument(
        "--parallel-samples",
        dest="max_parallel_samples",
        type=str,
        default=None,
        metavar="N|auto",
        help=argparse.SUPPRESS,
    )
    execution.add_argument(
        "--single-sample-mode",
        action="store_true",
        default=False,
        help=(
            "Force serial execution (alias for --max-parallel-samples 1). "
            "Use when you want one sample's logs interleaved cleanly or "
            "when memory pressure rules out concurrency."
        ),
    )
    execution.add_argument(
        "--memory-gb",
        type=str,
        default=None,
        metavar="GB|auto",
        help=(
            "Total memory budget (in GiB) the auto scheduler may use. "
            "Accepts a float or 'auto' (no memory cap). Used only when "
            "--max-parallel-samples auto."
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
    resolutions: list[SampleResolution],
    *,
    strandedness: str,
    mapq: int,
    min_length: int,
    max_length: int,
) -> list[str]:
    dedup_strategies = sorted({res.dedup_strategy for res in resolutions}) or ["skip"]
    actions: list[str] = [
        "validate external tools on PATH (cutadapt, bowtie2; plus "
        f"umi_tools as needed for dedup={'/'.join(dedup_strategies)})",
        f"resolve per-sample kit + dedup for {len(resolutions)} sample(s):",
    ]
    actions.extend("  " + line for line in resolution_summary_lines(resolutions))
    actions.extend(
        [
            f"cutadapt trim per sample (min={min_length} max={max_length} nt)",
            "bowtie2 contam subtract (--contam-index); unaligned reads pass "
            f"through. strand={strandedness}",
            "bowtie2 mt-transcriptome align (--mt-index) --end-to-end "
            f"--very-sensitive -L 18; strand={strandedness}",
            f"MAPQ filter at q >= {mapq} (pysam; NUMT suppression)",
            "dedup per sample (auto-resolved: umi-tools when UMIs present, else skip)",
            "BAM -> BED6 (strand-aware)",
            f"write kit_resolution.tsv, read_counts.tsv, run_settings.json "
            f"for {len(resolutions)} sample(s)",
        ]
    )
    return actions


def _per_sample_settings(resolutions: list[SampleResolution]) -> list[dict]:
    """JSON-friendly per-sample settings block for run_settings.json."""
    rows: list[dict] = []
    for resolution in resolutions:
        rows.append(
            {
                "sample": resolution.sample,
                "fastq": str(resolution.fastq),
                "applied_kit": resolution.kit.kit,
                "adapter": resolution.kit.adapter,
                "umi_length": resolution.kit.umi_length,
                "umi_position": resolution.kit.umi_position,
                "dedup_strategy": resolution.dedup_strategy,
                "detected_kit": resolution.detected_kit,
                "detection_match_rate": resolution.detection_match_rate,
                "detection_ambiguous": resolution.detection_ambiguous,
                "source": resolution.source,
            }
        )
    return rows


def _format_execution_table(
    resolutions: list[SampleResolution],
) -> list[str]:
    """Render a compact per-sample execution table for the pre-run banner.

    P5.6: lets the user verify the kit / adapter / UMI decisions and
    expected output BED paths before any expensive work runs. Columns
    mirror :data:`mitoribopy.align.sample_resolve._KIT_RESOLUTION_COLUMNS`
    but are formatted for stderr (truncated long fields, fixed-width
    columns).
    """
    if not resolutions:
        return ["execution table: no samples resolved."]

    cols = (
        "sample",
        "input_fastq",
        "applied_kit",
        "adapter",
        "umi",
        "dedup",
        "strand",
        "expected_output_bed",
    )
    header = (
        f"{cols[0]:<14} {cols[1]:<30} {cols[2]:<22} "
        f"{cols[3]:<22} {cols[4]:<8} {cols[5]:<10} "
        f"{cols[6]:<11} {cols[7]:<22}"
    )
    lines = ["execution table (kit_resolution.tsv has the full layout):"]
    lines.append(header)
    lines.append("-" * len(header))
    for r in resolutions:
        adapter = (r.kit.adapter or "").ljust(22)[:22]
        fastq = str(r.fastq.name)[:30].ljust(30)
        umi = f"{r.kit.umi_length}{r.kit.umi_position}".ljust(8)
        bed_name = f"{r.sample}.mt.bed".ljust(22)[:22]
        lines.append(
            f"{r.sample:<14} {fastq} {r.kit.kit:<22} "
            f"{adapter} {umi} {r.dedup_strategy:<10} "
            f"{'(see cfg)':<11} {bed_name}"
        )
    return lines


def _resolution_base_source(resolution: SampleResolution) -> str:
    """Return the resolver source without the per-sample override prefix."""
    return resolution.source.split(":", 1)[-1]


def _warn_adapter_detection_summary(
    resolutions: list[SampleResolution],
) -> None:
    """Emit explicit adapter-detection warnings for each scanned sample."""
    for resolution in resolutions:
        source = _resolution_base_source(resolution)
        if resolution.detected_kit:
            adapter = resolution.kit.adapter or ""
            detail = (
                f"{resolution.sample}: adapter detected: "
                f"{resolution.detected_kit} "
                f"({resolution.detection_match_rate * 100:.1f}% match)"
            )
            if adapter:
                detail += f"; resolved adapter={adapter}"
            if source == "user_adapter":
                detail += "; using the user-supplied adapter sequence"
            log_warning("ADAPTER", detail + ".")
            continue

        if source in {"inferred_pretrimmed", "user_adapter"}:
            if source == "user_adapter":
                suffix = " Using the user-supplied adapter sequence."
            else:
                suffix = (
                    " Input may already be adapter-filtered/trimmed; "
                    "if this is raw data with an unlisted adapter, pass "
                    "--adapter <SEQ> or set a per-sample adapter."
                )
            log_warning(
                "ADAPTER",
                f"{resolution.sample}: no known adapter was detected."
                + suffix,
            )


def _strict_publication_mode_errors(
    resolutions: list[SampleResolution],
) -> list[str]:
    """Return the list of align Strict Mode rejection reasons.

    Empty list means the resolutions passed the strict metadata checks.
    The CLI exits 2 when this list is non-empty.

    Rules:

    * No sample picked an inferred ``pretrimmed`` kit. The user must
      either pass ``--pretrimmed`` explicitly or supply ``--adapter``.
    * No sample has ``detection_ambiguous`` true OR a confidence
      margin below 0.10 — both signal a coin-flip between kits.
    * No sample landed on a UMI-bearing kit by inference; the
      ``umi_source`` must be ``declared`` or ``none`` (i.e. the user
      either declared the UMI or there isn't one).
    """
    errors: list[str] = []
    for r in resolutions:
        if r.source.startswith("inferred_pretrimmed"):
            errors.append(
                f"Strict Mode: sample {r.sample!r} picked "
                "kit=pretrimmed by inference (no known adapter signature). "
                "Declare --pretrimmed explicitly or supply --adapter."
            )
        if r.detection_ambiguous or (
            r.detection_match_rate > 0.0
            and getattr(r, "confidence_margin", 0.0) < 0.10
        ):
            errors.append(
                f"Strict Mode: sample {r.sample!r} adapter "
                f"detection is ambiguous "
                f"(margin={getattr(r, 'confidence_margin', 0.0):.2f}, "
                f"second_best={getattr(r, 'second_best_kit', None)!r}). "
                "Pin --adapter explicitly to remove the ambiguity."
            )
        umi_src = getattr(r, "umi_source", "preset_default")
        if r.kit.umi_length > 0 and umi_src not in ("declared", "none"):
            errors.append(
                f"Strict Mode: sample {r.sample!r} resolved "
                f"to a UMI-bearing kit (umi_length={r.kit.umi_length}) "
                f"with umi_source={umi_src!r}. Declare umi_length and "
                "umi_position in the sample sheet."
            )
    return errors


def _legacy_global_dedup(resolutions: list[SampleResolution]) -> str:
    """Pick a representative dedup label for legacy single-value consumers.

    Per-sample dedup means a single run can mix ``umi-tools`` and
    ``skip``. For backward compatibility some tests still read a
    single ``dedup_strategy`` field; we report the union as ``mixed``
    when it isn't uniform so an external reader notices.
    """
    strategies = {res.dedup_strategy for res in resolutions}
    if not strategies:
        return "skip"
    if len(strategies) == 1:
        return next(iter(strategies))
    return "mixed"


def _write_run_settings(
    output_dir: Path,
    *,
    args: argparse.Namespace,
    resolutions: list[SampleResolution],
    tool_versions: dict[str, str | None],
    inputs: list[Path],
) -> Path:
    settings = {
        "subcommand": "align",
        "mitoribopy_version": __version__,
        "adapter_default": args.adapter,
        "pretrimmed_default": bool(getattr(args, "pretrimmed", False)),
        "adapter_detection_mode": args.adapter_detection,
        "library_strandedness": args.library_strandedness,
        "min_length": args.min_length,
        "max_length": args.max_length,
        "quality": args.quality,
        "mapq_threshold": args.mapq,
        "seed": args.seed,
        "dedup_strategy": _legacy_global_dedup(resolutions),
        "umi_dedup_method": args.umi_dedup_method,
        "contam_index": args.contam_index,
        "mt_index": args.mt_index,
        "threads": getattr(args, "threads", None),
        "log_level": getattr(args, "log_level", "INFO"),
        "tool_versions": tool_versions,
        "inputs": [str(p) for p in inputs],
        "per_sample": _per_sample_settings(resolutions),
    }
    path = output_dir / "run_settings.json"
    path.write_text(json.dumps(settings, indent=2, sort_keys=True), encoding="utf-8")
    return path


def _per_worker_threads(threads: int | None, max_parallel: int) -> int:
    """Divide the global --threads budget across parallel sample workers.

    With --threads T and --max-parallel-samples N, each worker's external
    tools (cutadapt / bowtie2 / umi_tools) get ``max(1, T // N)`` threads
    so total CPU use stays around T regardless of N. Floors at 1 so a
    pathological T < N does not produce 0-thread tool invocations.
    """
    t = int(threads) if threads is not None else 1
    n = max(1, int(max_parallel))
    return max(1, t // n)


def _process_one_sample(
    *,
    output_dir: Path,
    args: argparse.Namespace,
    resolution: SampleResolution,
    threads: int,
    timings: StageTimings | None = None,
) -> SampleCounts:
    sample = resolution.sample
    fastq_in = resolution.fastq
    resolved_kit = resolution.kit
    resolved_dedup = resolution.dedup_strategy
    keep_intermediates = bool(getattr(args, "keep_intermediates", False))

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

    threads = max(1, int(threads))
    sample_sw = Stopwatch()
    sample_sw.__enter__()

    with stage_timer(timings, sample, "trim") as sw:
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
    _log_stage_done(
        sample,
        "trim",
        sw.seconds,
        f"kit={resolved_kit.kit}, kept "
        f"{_short_count(cutadapt.reads_passing_filters)}/"
        f"{_short_count(cutadapt.input_reads)} reads",
    )

    with stage_timer(timings, sample, "contam-filter") as sw:
        contam = contam_step.subtract_contaminants(
            fastq_in=trimmed_fq,
            contam_index=Path(args.contam_index),
            fastq_out_unaligned=nocontam_fq,
            strandedness=args.library_strandedness,
            threads=threads,
            seed=args.seed,
        )
    _log_stage_done(
        sample,
        "contam-filter",
        sw.seconds,
        f"{_short_count(contam.unaligned_reads)} passed, "
        f"{_short_count(contam.aligned_to_contam)} removed",
    )
    # The trimmed FASTQ is no longer read after contam-filter consumes
    # it. Drop it now so a 1000-sample run does not fill the disk with
    # ~50% of the original input size in regenerable .fq.gz copies.
    _maybe_unlink(trimmed_fq, keep=keep_intermediates)

    with stage_timer(timings, sample, "mt-align") as sw:
        alignment = align_step.align_mt(
            fastq_in=nocontam_fq,
            mt_index=Path(args.mt_index),
            bam_out=aligned_bam,
            strandedness=args.library_strandedness,
            threads=threads,
            seed=args.seed,
        )
    _log_stage_done(
        sample,
        "mt-align",
        sw.seconds,
        f"aligned {_short_count(alignment.aligned)}/"
        f"{_short_count(alignment.total_reads)} to mt-tx",
    )
    _maybe_unlink(nocontam_fq, keep=keep_intermediates)

    with stage_timer(timings, sample, "mapq-filter") as sw:
        mapq_kept = bam_utils.filter_bam_mapq(
            bam_in=aligned_bam,
            bam_out=mapq_bam,
            mapq_threshold=args.mapq,
        )
    _log_stage_done(
        sample,
        "mapq-filter",
        sw.seconds,
        f"MAPQ>={args.mapq}, kept {_short_count(mapq_kept)} reads",
    )
    # The pre-MAPQ BAM is identical to mapq_bam minus the filtered reads;
    # nothing downstream reads it again.
    _maybe_unlink(aligned_bam, keep=keep_intermediates)
    _maybe_unlink(aligned_bam.with_suffix(".bam.bai"), keep=keep_intermediates)

    with stage_timer(timings, sample, "dedup") as sw:
        if resolved_dedup == "skip":
            # Avoid the hardlink/copy of skip_dedup. bam_to_bed6 reads BAM
            # content, not a path, so feeding mapq.bam directly is identical
            # and saves the duplicated dedup.bam + dedup.bam.bai pair on
            # every no-UMI sample.
            dedup = dedup_step.skip_dedup_in_place(mapq_bam)
        else:
            dedup = dedup_step.run_dedup(
                strategy=resolved_dedup,
                umi_length=resolved_kit.umi_length,
                bam_in=mapq_bam,
                bam_out=dedup_bam,
                log_path=dedup_log,
                umi_tools_method=args.umi_dedup_method,
            )
    if resolved_dedup == "skip":
        dedup_detail = f"skip; {_short_count(dedup.output_reads)} reads through"
    else:
        dedup_detail = (
            f"{resolved_dedup}; kept {_short_count(dedup.output_reads)}/"
            f"{_short_count(dedup.input_reads)} reads"
        )
    _log_stage_done(sample, "dedup", sw.seconds, dedup_detail)

    with stage_timer(timings, sample, "bam-to-bed") as sw:
        bam_utils.bam_to_bed6(bam_in=dedup.bam_path, bed_out=bed_out)
    _log_stage_done(sample, "bam-to-bed", sw.seconds, f"wrote {bed_out.name}")

    sample_sw.__exit__(None, None, None)
    log_info(
        "ALIGN",
        f"{sample}: ✓ done in {format_duration(sample_sw.seconds)}",
    )

    return counts_step.assemble_sample_counts(
        sample=sample,
        cutadapt=cutadapt,
        contam=contam,
        align=alignment,
        mapq_count=mapq_kept,
        dedup=dedup,
    )


def _sample_done_path(sample_done_dir: Path, sample: str) -> Path:
    """Return the per-sample done-marker JSON path."""
    return Path(sample_done_dir) / f"{sample}.json"


def _write_sample_done(sample_done_dir: Path, counts: SampleCounts) -> Path:
    """Persist a per-sample completion marker so --resume can skip it later.

    The JSON body is the full :class:`SampleCounts` dict so we can rebuild
    the read_counts.tsv row without re-running the sample. Written
    atomically (write to .tmp, then rename) so a kill mid-write cannot
    leave a half-flushed marker that would silently feed a wrong row
    into the aggregated table.
    """
    sample_done_dir = Path(sample_done_dir)
    sample_done_dir.mkdir(parents=True, exist_ok=True)
    target = _sample_done_path(sample_done_dir, counts.sample)
    tmp = target.with_suffix(target.suffix + ".tmp")
    tmp.write_text(
        json.dumps(asdict(counts), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    tmp.replace(target)
    return target


def _load_sample_done(
    sample_done_dir: Path, sample: str
) -> SampleCounts | None:
    """Load a per-sample done marker; return None if absent or malformed."""
    path = _sample_done_path(sample_done_dir, sample)
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    try:
        return SampleCounts(**payload)
    except TypeError:
        # Schema drift: an old marker from a previous version that no
        # longer fits SampleCounts. Treat as missing so the sample
        # re-runs cleanly.
        return None


def _maybe_unlink(path: Path, *, keep: bool) -> None:
    """Best-effort delete an intermediate file unless the user asked to keep them."""
    if keep:
        return
    try:
        Path(path).unlink()
    except FileNotFoundError:
        pass
    except OSError:
        # Cross-fs hardlinks or read-only mounts can fail; not worth
        # crashing the run over a leftover ~MB-sized file.
        pass


def _log_sample_stage(sample: str, stage: str, detail: str) -> None:
    """Emit a concise, consistent per-sample status line.

    Retained for callers that have not been migrated to the timed
    :func:`_log_stage_done` form. New per-sample stages should prefer
    :func:`_log_stage_done` so the duration is captured uniformly.
    """
    log_info("ALIGN", f"{sample}: {stage} - {detail}")


def _log_stage_done(
    sample: str, stage: str, seconds: float, detail: str
) -> None:
    """Emit a single compact per-stage line with elapsed wall time.

    Format::

        [ALIGN] sampleA: trim          12.3s — kit=auto, kept 941k/1.0M reads

    A fixed 11-char stage column and a 7-char right-aligned duration
    keep multi-sample logs scannable at a glance.
    """
    duration = format_duration(seconds)
    log_info(
        "ALIGN",
        f"{sample}: {stage:<11} {duration:>7} — {detail}",
    )


_COUNT_SUFFIXES = (("G", 1_000_000_000), ("M", 1_000_000), ("k", 1_000))


def _short_count(n: int | float | None) -> str:
    """Render an integer-ish count as a compact string.

    ``941_382 -> '941k'``, ``1_000_000 -> '1.0M'``, ``523 -> '523'``.
    Returns ``'-'`` for ``None`` / negatives so a missing field never
    shows up as a confusing zero in a stage log.
    """
    if n is None:
        return "-"
    try:
        v = float(n)
    except (TypeError, ValueError):
        return str(n)
    if v < 0:
        return "-"
    if v < 1000:
        return f"{int(v)}"
    for suffix, divisor in _COUNT_SUFFIXES:
        if v >= divisor:
            scaled = v / divisor
            if scaled >= 100:
                return f"{int(round(scaled))}{suffix}"
            if scaled >= 10:
                return f"{scaled:.1f}{suffix}".replace(".0", "")
            return f"{scaled:.1f}{suffix}"
    return f"{int(v)}"


def _resolve_per_sample(
    *,
    args: argparse.Namespace,
    inputs: list[Path],
    detector=None,
) -> list[SampleResolution]:
    """Run the per-sample resolver for the supplied inputs.

    Lifts ``--adapter-detection`` from ``args`` and translates it into
    the resolver's policy. The detection tuning flags
    (``--adapter-detect-reads``, ``--adapter-detect-min-rate``,
    ``--adapter-detect-min-len``, ``--adapter-detect-pretrimmed-threshold``)
    are folded into a closure-bound detector so tests can still inject
    their own detector callable while production runs honour every CLI
    knob.
    """
    samples = [(_sample_name(fq), fq) for fq in inputs]
    if detector is None:
        # Capture the user's detection knobs in a thin closure so the
        # underlying detect_adapter sees them per call without forcing
        # every resolver-call site to forward four extra kwargs.
        detect_n_reads = getattr(args, "adapter_detect_reads", 5000)
        detect_min_rate = getattr(args, "adapter_detect_min_rate", 0.30)
        detect_min_len = getattr(args, "adapter_detect_min_len", 12)
        detect_pretrimmed_threshold = getattr(
            args, "adapter_detect_pretrimmed_threshold", 0.05
        )

        def detector(fastq_path):  # noqa: E306 — local helper
            return adapter_detect.detect_adapter(
                fastq_path,
                n_reads=detect_n_reads,
                min_match_rate=detect_min_rate,
                min_match_len=detect_min_len,
                pretrimmed_threshold=detect_pretrimmed_threshold,
            )

    overrides_path = getattr(args, "sample_overrides", None)
    sample_overrides: dict[str, SampleOverride] | None = None
    if overrides_path:
        sample_overrides = read_sample_overrides_tsv(Path(overrides_path))

    return resolve_sample_resolutions(
        samples,
        adapter=args.adapter,
        pretrimmed=bool(getattr(args, "pretrimmed", False)),
        umi_length=args.umi_length,
        umi_position=args.umi_position,
        umi_length_5p=getattr(args, "umi_length_5p", None),
        umi_length_3p=getattr(args, "umi_length_3p", None),
        dedup_strategy=args.dedup_strategy,
        adapter_detection_mode=args.adapter_detection,
        allow_pretrimmed_inference=getattr(
            args, "allow_pretrimmed_inference", True
        ),
        detector=detector,
        sample_overrides=sample_overrides,
    )


def run(argv: Iterable[str]) -> int:
    """Entry point for ``mitoribopy align <args>``."""
    argv_list = list(argv)
    parser = build_parser()

    # Pre-parse just --config so we can fold any YAML / JSON / TOML
    # values into the parser's defaults BEFORE the real parse pass.
    # CLI flags still win on conflict because they appear in argv and
    # parse_args treats them as explicitly-set values that override the
    # default. Same pattern as `mitoribopy rpf`.
    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument("--config", default=None)
    pre_args, _ = pre_parser.parse_known_args(argv_list)
    if pre_args.config:
        try:
            cfg = common.load_config_file(pre_args.config)
        except (FileNotFoundError, RuntimeError, ValueError) as exc:
            print(f"[mitoribopy align] ERROR: {exc}", file=sys.stderr)
            return 2
        known_dests = {
            action.dest for action in parser._actions
            if action.dest and action.dest != "help"
        }
        applied = {
            key: value for key, value in cfg.items() if key in known_dests
        }
        unknown = sorted(set(cfg) - set(applied))
        if unknown:
            log_warning(
                "ALIGN",
                "Ignoring unknown --config keys: " + ", ".join(unknown),
            )
        if applied:
            parser.set_defaults(**applied)

    args = parser.parse_args(argv_list)
    common.apply_common_arguments(args)

    # Enumerate FASTQs first so the per-sample resolver can scan each
    # one and fail-fast before any expensive subprocess work.
    inputs: list[Path]
    try:
        inputs = _enumerate_fastqs(args.fastq_dir, args.fastq)
    except FileNotFoundError as exc:
        print(f"[mitoribopy align] ERROR: {exc}", file=sys.stderr)
        return 2

    # In --dry-run we skip per-sample scanning (we have not yet
    # validated that paths exist and reading bytes is undesirable). We
    # still produce a representative dry-run plan from the global args.
    resolutions: list[SampleResolution] = []
    if not args.dry_run and inputs:
        try:
            resolutions = _resolve_per_sample(args=args, inputs=inputs)
        except SampleResolutionError as exc:
            print(f"[mitoribopy align] ERROR: {exc}", file=sys.stderr)
            return 2
        _warn_adapter_detection_summary(resolutions)
        # Surface ambiguous detections (NEB family) so the user reviews
        # whether their kit really has UMIs.
        for resolution in resolutions:
            if resolution.detection_ambiguous:
                log_warning(
                    "ADAPTER",
                    f"{resolution.sample}: adapter scan was ambiguous between "
                    "kits sharing this adapter; confirm whether your library "
                    "has UMIs (resolved kit: "
                    f"{resolution.kit.kit}, umi_length="
                    f"{resolution.kit.umi_length}).",
                )

    if args.dry_run:
        # Build a synthetic resolution per planned input so the dry-run
        # plan reflects what the per-sample resolver WOULD do. We do
        # not actually scan FASTQs in dry-run; the per-sample plan
        # surfaces the resolver intent without touching disk.
        pretrimmed_flag = bool(getattr(args, "pretrimmed", False))
        if pretrimmed_flag:
            placeholder_kit_global = trim_step.resolve_kit_settings(
                "pretrimmed",
                adapter=None,
                umi_length=args.umi_length,
                umi_position=args.umi_position,
                umi_length_5p=getattr(args, "umi_length_5p", None),
                umi_length_3p=getattr(args, "umi_length_3p", None),
            )
            global_source = "dry_run_pretrimmed"
        elif args.adapter is not None:
            placeholder_kit_global = ResolvedKit(
                kit="custom",
                adapter=args.adapter,
                umi_length=int(args.umi_length or 0),
                umi_position=args.umi_position or "5p",  # type: ignore[arg-type]
                umi_length_5p=int(getattr(args, "umi_length_5p", 0) or 0),
                umi_length_3p=int(getattr(args, "umi_length_3p", 0) or 0),
            )
            global_source = "dry_run_user_adapter"
        else:
            placeholder_kit_global = ResolvedKit(
                kit="auto-detect-at-runtime",
                adapter="(detect at runtime)",
                umi_length=0,
                umi_position="5p",
            )
            global_source = "dry_run_auto"

        if inputs:
            planned: list[SampleResolution] = []
            for fastq in inputs:
                sample = _sample_name(fastq)
                if pretrimmed_flag or args.adapter is not None:
                    dedup_label = dedup_step.resolve_dedup_strategy(
                        args.dedup_strategy,
                        umi_length=placeholder_kit_global.umi_length,
                    )
                else:
                    dedup_label = "auto"
                planned.append(
                    SampleResolution(
                        sample=sample,
                        fastq=fastq,
                        kit=placeholder_kit_global,
                        dedup_strategy=dedup_label,  # type: ignore[arg-type]
                        detected_kit=None,
                        detection_match_rate=0.0,
                        detection_ambiguous=False,
                        source=global_source,
                    )
                )
            try:
                actions = _planned_actions(
                    planned,
                    strandedness=args.library_strandedness,
                    mapq=args.mapq,
                    min_length=args.min_length,
                    max_length=args.max_length,
                )
            except (KeyError, ValueError) as exc:
                print(f"[mitoribopy align] ERROR: {exc}", file=sys.stderr)
                return 2
        else:
            # No inputs supplied; still validate the dedup combo so
            # the user gets fast feedback on bad flag combinations.
            try:
                dedup_label = dedup_step.resolve_dedup_strategy(
                    args.dedup_strategy,
                    umi_length=placeholder_kit_global.umi_length,
                )
            except (KeyError, ValueError) as exc:
                print(f"[mitoribopy align] ERROR: {exc}", file=sys.stderr)
                return 2
            placeholder = SampleResolution(
                sample="(no inputs)",
                fastq=Path(""),
                kit=placeholder_kit_global,
                dedup_strategy=dedup_label,  # type: ignore[arg-type]
                detected_kit=None,
                detection_match_rate=0.0,
                detection_ambiguous=False,
                source="dry_run_no_inputs",
            )
            actions = _planned_actions(
                [placeholder],
                strandedness=args.library_strandedness,
                mapq=args.mapq,
                min_length=args.min_length,
                max_length=args.max_length,
            )
            actions.append(
                "no input FASTQs supplied; per-sample resolution would run "
                f"with adapter={args.adapter!r}, pretrimmed={pretrimmed_flag}, "
                f"adapter_detection={args.adapter_detection}."
            )
        return common.emit_dry_run("align", actions)

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

    # Tool check uses the union of every sample's resolved dedup
    # strategy so we only require umi_tools when at least one sample
    # actually needs it.
    required = ["cutadapt", "bowtie2", "bowtie2-build"]
    required.extend(sorted(required_dedup_tools(resolutions)))
    optional: list[str] = ["fastqc"]
    try:
        tools = tool_check.ensure_tools_available(
            required=required, optional=optional
        )
    except tool_check.ToolNotFoundError as exc:
        print(f"[mitoribopy align] ERROR: {exc}", file=sys.stderr)
        return 2

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = configure_file_logging(output_dir / "mitoribopy.log")
    log_info("ALIGN", f"Log file: {log_path}")

    tool_versions = {
        name: (info.version if info is not None else None)
        for name, info in tools.items()
    }

    # Persist the per-sample resolution table for provenance + UX.
    write_kit_resolution_tsv(resolutions, output_dir / "kit_resolution.tsv")
    log_info(
        "ALIGN",
        f"Wrote kit_resolution.tsv with per-sample kit + dedup decisions for "
        f"{len(resolutions)} sample(s).",
    )

    # P5.6: pre-run per-sample execution table. The kit_resolution.tsv
    # we just wrote covers everything the user needs to verify, but
    # showing a compact table on stderr before any expensive work
    # starts catches a wrong adapter / missing UMI declaration in
    # seconds rather than at minutes-into-the-run.
    for line in _format_execution_table(resolutions):
        log_info("ALIGN", line)

    # P5.6: align Strict Mode gate.
    if getattr(args, "strict_publication_mode", False):
        strict_errors = _strict_publication_mode_errors(resolutions)
        # Refactor-4 (report §3.2.E): a publication run must never demote
        # the read-count invariant check to a warning, even by accident.
        if getattr(args, "allow_count_invariant_warning", False):
            strict_errors.append(
                "--allow-count-invariant-warning is incompatible with "
                "Strict Mode (it is a developer / debug escape hatch)."
            )
        if strict_errors:
            for line in strict_errors:
                log_error("ALIGN", line)
            return 2

    log_info(
        "ALIGN",
        f"Starting align pipeline for {len(inputs)} sample(s) "
        f"(strand='{args.library_strandedness}'; per-sample resolution active).",
    )
    _write_run_settings(
        output_dir,
        args=args,
        resolutions=resolutions,
        tool_versions=tool_versions,
        inputs=inputs,
    )

    sample_done_dir = output_dir / ".sample_done"
    resume_active = bool(getattr(args, "resume", False))
    if resume_active:
        sample_done_dir.mkdir(parents=True, exist_ok=True)

    from ..pipeline.resource_plan import plan_parallelism, write_resource_plan

    requested_threads = getattr(args, "threads", None) or "auto"
    requested_parallel_raw = getattr(args, "max_parallel_samples", None)
    if getattr(args, "single_sample_mode", False):
        requested_parallel: int | str = 1
    elif requested_parallel_raw is None:
        requested_parallel = "auto"
    else:
        requested_parallel = requested_parallel_raw
    requested_memory = getattr(args, "memory_gb", None)

    resource_plan = plan_parallelism(
        n_samples=len(resolutions) or 1,
        requested_threads=requested_threads,
        requested_parallel=requested_parallel,
        memory_gb=requested_memory,
    )
    plan_path = write_resource_plan(resource_plan, output_dir)
    log_info(
        "ALIGN",
        f"Resource plan: {resource_plan.parallel_samples} parallel sample(s) x "
        f"{resource_plan.per_sample_threads} thread(s) per worker "
        f"(threads={resource_plan.total_threads}, "
        f"reason={resource_plan.reason!r}); written to {plan_path}.",
    )

    max_parallel = resource_plan.parallel_samples
    per_worker_threads = resource_plan.per_sample_threads

    # Split resolutions into resume-cached (instantly available) and
    # pending (need real work). Cached samples short-circuit the executor
    # so they don't occupy a worker slot for what is just a JSON read.
    rows: list[SampleCounts] = []
    pending: list[SampleResolution] = []
    for resolution in resolutions:
        cached = (
            _load_sample_done(sample_done_dir, resolution.sample)
            if resume_active
            else None
        )
        if cached is not None:
            log_info(
                "ALIGN",
                f"{resolution.sample}: resumed from "
                f".sample_done/{resolution.sample}.json (already complete).",
            )
            rows.append(cached)
        else:
            pending.append(resolution)

    effective_parallel = min(max_parallel, max(1, len(pending)))
    log_info(
        "ALIGN",
        f"Concurrency: {effective_parallel} parallel sample(s) x "
        f"{per_worker_threads} thread(s) per tool "
        f"(--threads={requested_threads}, --max-parallel-samples={max_parallel}, "
        f"{len(pending)} sample(s) to run, {len(rows)} resumed).",
    )

    completed_counter = {"done": len(rows)}
    counter_lock = threading.Lock()
    total = len(resolutions)
    timings = StageTimings()
    sample_counter: SampleCounter | None = None
    wall_sw = Stopwatch()
    wall_sw.__enter__()

    def _run_one(resolution: SampleResolution) -> SampleCounts:
        counts = _process_one_sample(
            output_dir=output_dir,
            args=args,
            resolution=resolution,
            threads=per_worker_threads,
            timings=timings,
        )
        _write_sample_done(sample_done_dir, counts)
        with counter_lock:
            completed_counter["done"] += 1
            done = completed_counter["done"]
        log_info(
            "ALIGN",
            f"{resolution.sample}: complete ({done}/{total}) "
            f"(post_trim={counts.post_trim}, "
            f"mt_aligned={counts.mt_aligned}, "
            f"post_dedup={counts.mt_aligned_after_dedup}).",
        )
        if sample_counter is not None:
            sample_counter.advance(resolution.sample)
        return counts

    # Optional tqdm "samples done" counter for the parallel path; the
    # serial path keeps its existing per-sample log_progress banner so
    # tests / non-TTY runs see the same output as before.
    counter_disabled = effective_parallel <= 1 or len(pending) == 0
    sample_counter = SampleCounter(
        total=len(pending),
        desc="ALIGN samples",
        disable=counter_disabled,
    )
    sample_counter.__enter__()
    try:
        if effective_parallel <= 1:
            # Serial path: preserves deterministic log ordering and
            # avoids spinning a thread pool for the common single-
            # sample case.
            for index, resolution in enumerate(pending, start=1):
                log_progress(
                    "ALIGN",
                    len(rows) + index,
                    total,
                    f"Processing sample {len(rows) + index}/{total} "
                    f"({resolution.sample}).",
                )
                rows.append(_run_one(resolution))
        else:
            # Parallel path. Submit all pending samples to a thread
            # pool; iterate in completion order. On first exception,
            # cancel not-yet-started futures and re-raise -- preserves
            # the existing fail-fast semantics. Already-running
            # futures finish naturally because we cannot interrupt
            # subprocess.run() cleanly.
            log_info(
                "ALIGN",
                f"Submitting {len(pending)} sample(s) to a "
                f"ThreadPoolExecutor(max_workers={effective_parallel}).",
            )
            for resolution in pending:
                log_info("ALIGN", f"{resolution.sample}: queued.")
            with ThreadPoolExecutor(max_workers=effective_parallel) as executor:
                future_to_sample: dict[Future[SampleCounts], str] = {}
                for resolution in pending:
                    future = executor.submit(_run_one, resolution)
                    future_to_sample[future] = resolution.sample
                try:
                    done_set, _ = wait(
                        future_to_sample.keys(), return_when=FIRST_EXCEPTION
                    )
                    # If a future raised, surface the first exception
                    # now; other in-flight futures will be left to
                    # finish via the executor context manager's
                    # shutdown(wait=True).
                    for future in done_set:
                        exc = future.exception()
                        if exc is not None:
                            sample = future_to_sample[future]
                            log_warning(
                                "ALIGN",
                                f"{sample}: failed with "
                                f"{type(exc).__name__}: {exc}; "
                                "cancelling pending samples.",
                            )
                            for pending_future in future_to_sample:
                                if not pending_future.done():
                                    pending_future.cancel()
                            raise exc
                    # No exception path -- collect the rest of the
                    # results.
                    for future in future_to_sample:
                        rows.append(future.result())
                except BaseException:
                    # Best-effort cancel + propagate. The 'with' block
                    # waits for already-running workers to finish
                    # (subprocess.run cannot be interrupted mid-call).
                    for pending_future in future_to_sample:
                        if not pending_future.done():
                            pending_future.cancel()
                    raise
    finally:
        sample_counter.__exit__(None, None, None)

    # Deterministic sort by sample name.
    rows.sort(key=lambda row: row.sample)
    # P5.6: enforce read_counts invariants. Default behaviour aborts
    # the run on any violation; --allow-count-invariant-warning demotes
    # to warnings.tsv entries. The previous behaviour was to silently
    # write the bad numbers into the TSV.
    try:
        counts_step.enforce_count_invariants(
            rows,
            allow_warning=getattr(args, "allow_count_invariant_warning", False),
        )
    except counts_step.CountInvariantError as exc:
        log_error("ALIGN", f"{exc}")
        return 2
    counts_step.write_read_counts_table(rows, output_dir / "read_counts.tsv")

    # P5.7 (publication review): emit umi_qc.tsv so reviewers can
    # confirm that UMI extraction + coordinate+UMI dedup behaved as
    # advertised. We pair every SampleCounts row with the resolved kit
    # so the table records the umi_length / umi_position / dedup_method
    # actually used, not just the requested values.
    umi_qc_rows: list[dedup_step.UmiQCRow] = []
    resolutions_by_sample = {r.sample: r for r in resolutions}
    requested_method = str(getattr(args, "umi_dedup_method", "unique") or "unique")
    for row in rows:
        resolution = resolutions_by_sample.get(row.sample)
        umi_length = (
            int(resolution.kit.umi_length)
            if resolution is not None and resolution.kit is not None
            else 0
        )
        umi_position = (
            str(resolution.kit.umi_position)
            if resolution is not None
            and resolution.kit is not None
            and getattr(resolution.kit, "umi_position", None) is not None
            else None
        )
        dedup_strategy = (
            resolution.dedup_strategy
            if resolution is not None and resolution.dedup_strategy
            else "auto"
        )
        effective_strategy = dedup_step.resolve_dedup_strategy(
            dedup_strategy, umi_length=umi_length,
        )
        method = requested_method if effective_strategy == "umi-tools" else "skip"
        umi_length_5p = (
            int(getattr(resolution.kit, "umi_length_5p", 0) or 0)
            if resolution is not None and resolution.kit is not None
            else 0
        )
        umi_length_3p = (
            int(getattr(resolution.kit, "umi_length_3p", 0) or 0)
            if resolution is not None and resolution.kit is not None
            else 0
        )
        umi_qc_rows.append(
            dedup_step.build_umi_qc_row(
                sample_id=row.sample,
                umi_length=umi_length,
                umi_position=umi_position,
                dedup_strategy=str(effective_strategy),
                dedup_method=method,
                pre_count=int(row.mt_aligned_after_mapq),
                post_count=int(row.mt_aligned_after_dedup),
                umi_length_5p=umi_length_5p,
                umi_length_3p=umi_length_3p,
            )
        )
    dedup_step.write_umi_qc_tsv(umi_qc_rows, output_dir / "umi_qc.tsv")

    wall_sw.__exit__(None, None, None)
    log_info(
        "ALIGN",
        f"Finished align pipeline for {len(rows)} sample(s) in "
        f"{format_duration(wall_sw.seconds)}. Wrote read_counts.tsv, "
        "kit_resolution.tsv, and run_settings.json to "
        f"{output_dir}.",
    )

    # Per-stage timing summary table (skipped for resume-only runs
    # where every sample was loaded from cache and no stage timings
    # were recorded).
    if timings.per_sample_view():
        for line in render_summary_lines(
            timings,
            wall_seconds=wall_sw.seconds,
            samples_total=len(timings.per_sample_view()),
        ):
            log_info("ALIGN", line)
    return 0
