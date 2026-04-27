"""``mitoribopy align`` subcommand - FASTQ -> BAM + BED + read counts.

End-to-end orchestrator for Phase 3. For every input FASTQ:

1. cutadapt trim (kit-aware; optional UMI extraction)
2. bowtie2 contaminant subtraction (user-supplied index; fails loudly if missing)
3. bowtie2 mt-transcriptome alignment (Path A; end-to-end --very-sensitive -L 18)
4. MAPQ filter (pysam; default 10 to suppress NUMT cross-talk)
5. deduplication (umi_tools / skip)
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
    sample_resolve,
    tool_check,
    trim as trim_step,
)
from ..align._types import (
    KIT_PRESET_ALIASES,
    KIT_PRESETS,
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
from ..console import configure_file_logging, log_info, log_progress, log_warning
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
        "--kit-preset",
        choices=sorted(KIT_PRESETS) + sorted(KIT_PRESET_ALIASES),
        default="auto",
        metavar="PRESET",
        help=(
            "Library-prep adapter family. 'auto' (default) runs per-sample "
            "adapter detection and resolves the kit independently for each "
            "input FASTQ; mixed-kit and mixed-UMI runs are supported. "
            "Canonical presets organized by adapter family:\n"
            "  illumina_smallrna     Illumina TruSeq Small RNA adapter, no UMI\n"
            "  illumina_truseq       Illumina TruSeq Read 1 adapter, no UMI "
            "(NEBNext Small RNA, TruSeq Stranded Total, SMARTer Pico v3, "
            "SEQuoia Express, …)\n"
            "  illumina_truseq_umi   Same adapter + 8 nt 5' UMI (NEBNext Ultra "
            "II UMI, SEQuoia Complete UMI, …)\n"
            "  qiaseq_mirna          QIAseq miRNA adapter + 12 nt 3' UMI\n"
            "  pretrimmed            FASTQs already adapter-trimmed; cutadapt "
            "skips the -a flag\n"
            "  custom                requires --adapter <SEQ>\n"
            "Legacy vendor names (truseq_smallrna, nebnext_smallrna, "
            "nebnext_ultra_umi, …) remain accepted as aliases."
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
        "--sample-overrides",
        default=None,
        metavar="TSV",
        help=(
            "Path to a TSV with per-sample overrides for "
            "kit_preset / adapter / umi_length / umi_position / "
            "dedup_strategy. Required header columns: 'sample' plus at "
            "least one of the override columns. The 'sample' value must "
            "match the FASTQ basename with the .fq[.gz] / .fastq[.gz] "
            "suffix removed. Empty cells (or NA / None / null) fall "
            "through to the global CLI default for that field, so a "
            "single sample can override only its UMI without restating "
            "the rest of the kit. Useful for mixed-UMI batches where "
            "each sample's UMI length / position differs."
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
            "every input FASTQ and picks the matching preset per sample; "
            "samples whose scan fails fall back to the user's --kit-preset / "
            "--adapter when supplied, or to the 'pretrimmed' kit when no "
            "fallback is set and the data looks already-trimmed. 'strict' "
            "scans and HARD-FAILS on any sample whose scan disagrees with "
            "an explicit --kit-preset or where no preset can be identified. "
            "'off' skips the scan and trusts --kit-preset / --adapter for "
            "every sample (requires an explicit preset)."
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
            "failure with no --kit-preset fallback raises an error "
            "instead (the v0.4.0 behaviour before this change)."
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
        choices=["auto", "umi-tools", "skip"],
        default="auto",
        help=(
            "'auto' (default) -> umi-tools when UMIs are present, else "
            "skip. The legacy 'mark-duplicates' option was removed in "
            "v0.4.5: coordinate-only dedup destroys codon-occupancy "
            "signal on low-complexity mt-Ribo-seq libraries."
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
        type=int,
        default=1,
        metavar="N",
        help=(
            "Number of samples to process concurrently. With --threads T, "
            "each worker uses max(1, T // N) tool threads so total CPU "
            "use stays around T. Default 1 (serial; backward-compatible). "
            "Per-sample work (cutadapt + bowtie2 + dedup + BAM->BED) is "
            "embarrassingly parallel; the joint 'mitoribopy rpf' stage is "
            "unaffected."
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


def _legacy_global_dedup(resolutions: list[SampleResolution]) -> str:
    """Pick a representative dedup label for legacy single-value consumers.

    The per-sample switch in v0.4.0 means a single run can mix
    ``umi-tools`` and ``skip``. For backward compatibility some tests
    still read a single ``dedup_strategy`` field; we report the union as
    ``mixed`` when it isn't uniform so an external reader notices.
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
        "kit_preset_default": args.kit_preset,
        "adapter_default": args.adapter,
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

    _log_sample_stage(
        sample,
        "trim",
        f"kit={resolved_kit.kit}, adapter={resolved_kit.adapter}, "
        f"umi_length={resolved_kit.umi_length}",
    )
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
    _log_sample_stage(
        sample,
        "trim",
        f"kept {cutadapt.reads_passing_filters}/{cutadapt.input_reads} reads after cutadapt",
    )

    _log_sample_stage(sample, "contam-filter", f"subtracting contaminants with index {args.contam_index}")
    contam = contam_step.subtract_contaminants(
        fastq_in=trimmed_fq,
        contam_index=Path(args.contam_index),
        fastq_out_unaligned=nocontam_fq,
        strandedness=args.library_strandedness,
        threads=threads,
        seed=args.seed,
    )
    _log_sample_stage(
        sample,
        "contam-filter",
        f"{contam.unaligned_reads} reads passed through; {contam.aligned_to_contam} removed as contaminants",
    )
    # The trimmed FASTQ is no longer read after contam-filter consumes
    # it. Drop it now so a 1000-sample run does not fill the disk with
    # ~50% of the original input size in regenerable .fq.gz copies.
    _maybe_unlink(trimmed_fq, keep=keep_intermediates)

    _log_sample_stage(sample, "mt-align", f"aligning against {args.mt_index}")
    alignment = align_step.align_mt(
        fastq_in=nocontam_fq,
        mt_index=Path(args.mt_index),
        bam_out=aligned_bam,
        strandedness=args.library_strandedness,
        threads=threads,
        seed=args.seed,
    )
    _log_sample_stage(
        sample,
        "mt-align",
        f"aligned {alignment.aligned}/{alignment.total_reads} reads to the mt-transcriptome",
    )
    _maybe_unlink(nocontam_fq, keep=keep_intermediates)

    _log_sample_stage(sample, "mapq-filter", f"keeping reads with MAPQ >= {args.mapq}")
    mapq_kept = bam_utils.filter_bam_mapq(
        bam_in=aligned_bam,
        bam_out=mapq_bam,
        mapq_threshold=args.mapq,
    )
    _log_sample_stage(sample, "mapq-filter", f"kept {mapq_kept} reads")
    # The pre-MAPQ BAM is identical to mapq_bam minus the filtered reads;
    # nothing downstream reads it again.
    _maybe_unlink(aligned_bam, keep=keep_intermediates)
    _maybe_unlink(aligned_bam.with_suffix(".bam.bai"), keep=keep_intermediates)

    _log_sample_stage(sample, "dedup", f"strategy={resolved_dedup}")
    if resolved_dedup == "skip":
        # Avoid the hardlink/copy of skip_dedup. bam_to_bed6 reads BAM
        # content, not a path, so feeding mapq.bam directly is identical
        # and saves the duplicated dedup.bam + dedup.bam.bai pair on
        # every no-UMI sample.
        dedup = dedup_step.skip_dedup_in_place(mapq_bam)
        _log_sample_stage(
            sample,
            "dedup",
            f"strategy=skip; passing {dedup.output_reads} reads through unchanged",
        )
    else:
        dedup = dedup_step.run_dedup(
            strategy=resolved_dedup,
            umi_length=resolved_kit.umi_length,
            bam_in=mapq_bam,
            bam_out=dedup_bam,
            log_path=dedup_log,
            umi_tools_method=args.umi_dedup_method,
        )
        _log_sample_stage(
            sample,
            "dedup",
            f"kept {dedup.output_reads}/{dedup.input_reads} reads after deduplication",
        )

    _log_sample_stage(sample, "bam-to-bed", f"writing {bed_out.name}")
    bam_utils.bam_to_bed6(bam_in=dedup.bam_path, bed_out=bed_out)
    _log_sample_stage(sample, "bam-to-bed", "finished BED6 export")

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
    """Emit a concise, consistent per-sample status line."""
    log_info("ALIGN", f"{sample}: {stage} - {detail}")


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
        kit_preset=args.kit_preset,
        adapter=args.adapter,
        umi_length=args.umi_length,
        umi_position=args.umi_position,
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
        # plan reflects what the per-sample resolver WOULD do. We trust
        # the user-supplied kit_preset; when 'auto' we cannot scan, so
        # report 'pending: auto-detect at run time'.
        if inputs:
            planned: list[SampleResolution] = []
            for fastq in inputs:
                sample = _sample_name(fastq)
                if args.kit_preset == "auto":
                    # Use a synthetic placeholder so the plan still
                    # surfaces the per-sample loop intent.
                    placeholder_kit = ResolvedKit(
                        kit="auto-detect-at-runtime",
                        adapter="(detect at runtime)",
                        umi_length=0,
                        umi_position="5p",
                    )
                    dedup_label = "auto"
                    source_label = "dry_run_auto"
                else:
                    placeholder_kit = trim_step.resolve_kit_settings(
                        args.kit_preset,
                        adapter=args.adapter,
                        umi_length=args.umi_length,
                        umi_position=args.umi_position,
                    )
                    dedup_label = dedup_step.resolve_dedup_strategy(
                        args.dedup_strategy,
                        umi_length=placeholder_kit.umi_length,
                    )
                    source_label = "dry_run_explicit"
                planned.append(
                    SampleResolution(
                        sample=sample,
                        fastq=fastq,
                        kit=placeholder_kit,
                        dedup_strategy=dedup_label,  # type: ignore[arg-type]
                        detected_kit=None,
                        detection_match_rate=0.0,
                        detection_ambiguous=False,
                        source=source_label,
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
            # No inputs supplied; still validate the kit/dedup combo so
            # the user gets fast feedback on bad flag combinations.
            try:
                if args.kit_preset == "auto":
                    placeholder_kit = ResolvedKit(
                        kit="auto",
                        adapter="(per-sample auto-detect)",
                        umi_length=0,
                        umi_position="5p",
                    )
                else:
                    placeholder_kit = trim_step.resolve_kit_settings(
                        args.kit_preset,
                        adapter=args.adapter,
                        umi_length=args.umi_length,
                        umi_position=args.umi_position,
                    )
                dedup_label = dedup_step.resolve_dedup_strategy(
                    args.dedup_strategy,
                    umi_length=placeholder_kit.umi_length,
                )
            except (KeyError, ValueError) as exc:
                print(f"[mitoribopy align] ERROR: {exc}", file=sys.stderr)
                return 2
            placeholder = SampleResolution(
                sample="(no inputs)",
                fastq=Path(""),
                kit=placeholder_kit,
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
                f"no input FASTQs supplied; per-sample resolution would run "
                f"with kit_preset={args.kit_preset}, "
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

    requested_threads = getattr(args, "threads", None) or 1
    max_parallel = max(1, int(getattr(args, "max_parallel_samples", 1) or 1))
    per_worker_threads = _per_worker_threads(requested_threads, max_parallel)

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

    def _run_one(resolution: SampleResolution) -> SampleCounts:
        counts = _process_one_sample(
            output_dir=output_dir,
            args=args,
            resolution=resolution,
            threads=per_worker_threads,
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
        return counts

    if effective_parallel <= 1:
        # Serial path: preserves deterministic log ordering and avoids
        # spinning a thread pool for the common single-sample case.
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
        # Parallel path. Submit all pending samples to a thread pool;
        # iterate in completion order. On first exception, cancel
        # not-yet-started futures and re-raise -- preserves the existing
        # fail-fast semantics. Already-running futures finish naturally
        # because we cannot interrupt subprocess.run() cleanly.
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
                # If a future raised, surface the first exception now;
                # other in-flight futures will be left to finish via the
                # executor context manager's shutdown(wait=True).
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
                        # Cancel any future that hasn't started yet.
                        for pending_future in future_to_sample:
                            if not pending_future.done():
                                pending_future.cancel()
                        raise exc
                # No exception path -- collect the rest of the results.
                for future in future_to_sample:
                    rows.append(future.result())
            except BaseException:
                # Best-effort cancel + propagate. The 'with' block waits
                # for already-running workers to finish (subprocess.run
                # cannot be interrupted mid-call).
                for pending_future in future_to_sample:
                    if not pending_future.done():
                        pending_future.cancel()
                raise

    # Deterministic sort by sample name.
    rows.sort(key=lambda row: row.sample)
    counts_step.write_read_counts_table(rows, output_dir / "read_counts.tsv")
    log_info(
        "ALIGN",
        f"Finished align pipeline for {len(rows)} sample(s). Wrote "
        "read_counts.tsv, kit_resolution.tsv, and run_settings.json to "
        f"{output_dir}.",
    )
    return 0
