"""``mitoribopy rnaseq`` subcommand — RNA-seq + Ribo-seq → TE / ΔTE + plots.

Two ways to drive it:

* **Default (from-FASTQ).** Pass raw RNA-seq + Ribo-seq FASTQs and a
  transcriptome FASTA. The subcommand auto-detects SE vs PE from the
  filename mate token, reuses the existing ``mitoribopy.align`` adapter
  / kit / UMI machinery, runs cutadapt + bowtie2 per sample, counts
  reads per transcript, runs pyDESeq2 on the RNA side, then falls
  through into the TE / ΔTE / plot path. Requires the ``[fastq]``
  optional-dependency extra (``pip install 'mitoribopy[fastq]'``).

* **Alternative (pre-computed DE).** Pass an existing DE table (DESeq2
  / Xtail / Anota2Seq) plus a prior ``mitoribopy rpf`` run via
  ``--de-table`` + ``--ribo-dir``. A SHA256 reference-consistency gate
  rejects mismatched references before any math runs. Use this when
  you already ran DE externally on the full transcriptome — the
  recommended path for publication-grade DE statistics.

The two paths are mutually exclusive (passing both ``--rna-fastq`` and
``--de-table`` exits with code 2) and produce the same ``te.tsv``,
``delta_te.tsv``, and ``plots/`` output shape; the from-FASTQ path
additionally writes the intermediate counts matrices and a generated
``de_table.tsv``.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Iterable

from .. import __version__
from ..console import log_warning
from ..rnaseq import (
    compute_delta_te,
    compute_reference_checksum,
    compute_te,
    detect_de_format,
    load_de_table,
    load_ribo_counts,
    match_mt_mrnas,
    verify_reference_consistency,
    ReferenceMismatchError,
)
from ..rnaseq._types import DTeRow, TeRow
from ..rnaseq.plots import (
    plot_de_volcano,
    plot_delta_te_volcano,
    plot_ma,
    plot_mrna_vs_rpf_scatter,
    plot_pca,
    plot_te_bar_grouped,
    plot_te_compare_scatter,
    plot_te_heatmap,
    plot_te_log2fc_bar,
)
from . import common


RNASEQ_SUBCOMMAND_HELP = (
    "Translation efficiency (TE / delta-TE) from paired RNA-seq + Ribo-seq. "
    "Default flow: pass --rna-fastq + --ribo-fastq + --reference-fasta and "
    "the subcommand runs trimming, bowtie2 alignment, per-transcript counting, "
    "and pyDESeq2 itself before emitting te.tsv, delta_te.tsv, and plots. "
    "Alternative: pass --de-table from a prior external DESeq2 / Xtail / "
    "Anota2Seq run together with --ribo-dir; this path is mutually exclusive "
    "with --rna-fastq and enforces a SHA256 reference-consistency gate."
)


# Explicit mode names. de_table is the publication-grade flow (external
# full-transcriptome DE table); from_fastq is the in-tree exploratory
# flow (mt-only pyDESeq2); rna_only emits an mt-mRNA count matrix
# from RNA FASTQs without computing TE (useful when RPF data is not
# available yet). 'none' means the rnaseq stage is configured but no
# inputs were supplied — typically used by `mitoribopy all` to record
# "stage section present, intentionally inert".
RNASEQ_MODES: tuple[str, ...] = ("de_table", "from_fastq", "rna_only", "none")


def _normalize_mode(value: str | None) -> str | None:
    """Accept hyphenated forms (`from-fastq`) and normalize to underscore form."""
    if value is None:
        return None
    return value.replace("-", "_").lower()


def _resolve_rnaseq_mode(args) -> tuple[str | None, str | None]:
    """Resolve the rnaseq mode from --rnaseq-mode + supplied inputs.

    Returns ``(mode, error_message)``. If ``error_message`` is set, the
    caller should print it and exit 2. The mode strings are members of
    :data:`RNASEQ_MODES`.

    Inference rule when --rnaseq-mode is unset:

    * ``de_table`` if ``--de-table`` is set
    * ``from_fastq`` if ``--rna-fastq`` or ``--sample-sheet`` is set
    * ``none`` otherwise

    When --rnaseq-mode IS set, the supplied inputs must match the chosen
    mode; mismatches return a clear error.
    """
    requested = _normalize_mode(getattr(args, "rnaseq_mode", None))
    has_de_table = bool(args.de_table)
    has_fastq_inputs = bool(args.rna_fastq) or bool(args.sample_sheet)
    inferred = (
        "de_table" if has_de_table
        else "from_fastq" if has_fastq_inputs
        else "none"
    )

    if requested is None:
        # Ambiguous: user supplied inputs for BOTH flows without picking
        # one. The two paths are mutually exclusive — fail loudly so
        # users see exactly which flag pair conflicts.
        if has_de_table and has_fastq_inputs:
            return None, (
                "rnaseq inputs are ambiguous: --de-table (publication "
                "flow) and --rna-fastq / --sample-sheet (exploratory "
                "flow) are mutually exclusive. Pick one, or set "
                "--rnaseq-mode={de_table|from_fastq} to disambiguate."
            )
        return inferred, None

    if requested not in RNASEQ_MODES:
        return None, (
            f"--rnaseq-mode={requested!r} is not a valid mode; "
            f"choose one of {list(RNASEQ_MODES)}."
        )

    # Reject mismatched mode + inputs.
    if requested == "de_table" and has_fastq_inputs:
        return None, (
            "--rnaseq-mode=de_table conflicts with --rna-fastq / "
            "--sample-sheet inputs (the de_table flow consumes a "
            "pre-computed external DE table, not raw FASTQs). Drop the "
            "FASTQ inputs or switch to --rnaseq-mode=from_fastq."
        )
    if requested == "from_fastq" and has_de_table:
        return None, (
            "--rnaseq-mode=from_fastq conflicts with --de-table (the "
            "from-FASTQ flow runs pyDESeq2 itself; --de-table belongs "
            "to the de_table flow). Drop --de-table or switch to "
            "--rnaseq-mode=de_table."
        )
    if requested == "de_table" and not has_de_table:
        return None, (
            "--rnaseq-mode=de_table requires --de-table <path>."
        )
    if requested == "from_fastq" and not has_fastq_inputs:
        return None, (
            "--rnaseq-mode=from_fastq requires either --rna-fastq or "
            "--sample-sheet."
        )

    return requested, None


def _emit_exploratory_mode_banner() -> None:
    """One-line stderr note when from_fastq mode is active.

    Distinct from the pseudo-replicate banner, which fires only when
    n=1 conditions are split. This banner fires for *every* from_fastq
    run because the in-tree mt-only pyDESeq2 path is exploratory by
    construction; publication-grade DE should run externally on the
    full transcriptome and feed back via --de-table.
    """
    sys.stderr.write(
        "[mitoribopy rnaseq] NOTE: rnaseq mode = from_fastq (exploratory). "
        "This runs pyDESeq2 on the mt-mRNA subset only; for "
        "publication-grade DE, run DESeq2 / Xtail / Anota2Seq externally "
        "on the full transcriptome and re-run with mode=de_table + "
        "--de-table <path>.\n"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mitoribopy rnaseq",
        description=RNASEQ_SUBCOMMAND_HELP,
        formatter_class=common.MitoRiboPyHelpFormatter,
    )
    common.add_common_arguments(parser)

    # ----- Mode (explicit; recommended) ------------------------------------
    parser.add_argument(
        "--rnaseq-mode",
        dest="rnaseq_mode",
        choices=list(RNASEQ_MODES),
        default=None,
        metavar="MODE",
        help=(
            "Explicit rnaseq stage mode. Recommended over relying on "
            "input-presence inference. Modes:\n"
            "  de_table   PUBLICATION: consume external DE table + prior "
            "rpf run (--de-table + --ribo-dir).\n"
            "  from_fastq EXPLORATORY: run cutadapt + bowtie2 + pyDESeq2 "
            "on the mt-mRNA subset (--rna-fastq or --sample-sheet).\n"
            "  none       stage configured but inert (no inputs).\n"
            "When omitted, the mode is inferred from supplied inputs."
        ),
    )

    # ----- Default flow: from raw FASTQ ------------------------------------
    fastq = parser.add_argument_group("Inputs (from raw FASTQ — default flow)")
    fastq.add_argument(
        "--sample-sheet",
        default=None,
        metavar="PATH",
        help=(
            "Unified sample sheet TSV (one row per FASTQ; columns "
            "documented in mitoribopy.sample_sheet). When provided, "
            "--rna-fastq, --ribo-fastq, and --condition-map are derived "
            "from the sheet and passing any of those flags alongside "
            "--sample-sheet is an error. Pairing between RNA-seq and "
            "Ribo-seq is by sample_id; index-based pairing is no longer "
            "supported."
        ),
    )
    fastq.add_argument(
        "--rna-fastq",
        nargs="+",
        default=None,
        metavar="PATH",
        help=(
            "RNA-seq FASTQ files or directories. The default driver of "
            "the rnaseq subcommand: pass these and the rest of the "
            "pipeline (trim → bowtie2 → count → pyDESeq2 → TE / ΔTE / "
            "plots) runs end-to-end. Mutually exclusive with --de-table "
            "and --sample-sheet."
        ),
    )
    fastq.add_argument(
        "--ribo-fastq",
        nargs="+",
        default=None,
        metavar="PATH",
        help=(
            "Ribo-seq FASTQ files or directories. When omitted the run "
            "short-circuits after writing de_table.tsv "
            "(mode='from-fastq-rna-only' in run_settings.json)."
        ),
    )
    fastq.add_argument(
        "--reference-fasta",
        default=None,
        metavar="PATH",
        help=(
            "Transcriptome FASTA used to build (or reuse) a content-"
            "addressed bowtie2 index. Required for the default flow. "
            "The SHA256 is recorded under from_fastq.reference_checksum "
            "in run_settings.json."
        ),
    )
    fastq.add_argument(
        "--bowtie2-index",
        default=None,
        metavar="PREFIX",
        help=(
            "Pre-built bowtie2 index prefix. When set, bowtie2-build is "
            "skipped; --reference-fasta is still required for hashing."
        ),
    )
    fastq.add_argument(
        "--workdir",
        default=None,
        metavar="DIR",
        help="Scratch directory for trim / index / BAM artefacts. "
             "Defaults to <output>/work.",
    )
    fastq.add_argument(
        "--align-threads",
        type=int,
        default=4,
        metavar="N",
        help=(
            "Threads passed to cutadapt and bowtie2. Separate from "
            "--threads (which caps BLAS / pyDESeq2 thread pools)."
        ),
    )
    fastq.add_argument(
        "--recount-ribo-fastq",
        dest="recount_ribo_fastq",
        action="store_true",
        default=False,
        help=(
            "When 'mitoribopy all' has just produced rpf_counts.tsv from "
            "the rpf stage, the rnaseq from-FASTQ flow defaults to "
            "REUSING those counts instead of re-aligning the Ribo FASTQs "
            "(rpf already did the work; re-running cutadapt + bowtie2 + "
            "counting is duplicate effort and a source of subtle "
            "inconsistencies between the two stages' counts). Pass "
            "--recount-ribo-fastq to force a second pass over the Ribo "
            "FASTQs from inside the rnaseq stage."
        ),
    )
    fastq.add_argument(
        "--upstream-rpf-counts",
        dest="upstream_rpf_counts",
        default=None,
        metavar="PATH",
        help=argparse.SUPPRESS,  # internal; set by `mitoribopy all` orchestrator
    )
    fastq.add_argument(
        "--upstream-rna-counts",
        dest="upstream_rna_counts",
        default=None,
        metavar="PATH",
        help=argparse.SUPPRESS,  # internal; set by `mitoribopy all` orchestrator
    )
    fastq.add_argument(
        "--align-only",
        dest="align_only",
        action="store_true",
        default=False,
        help=(
            "Run only the FASTQ-trim + bowtie2 + counts phase and exit "
            "before DESeq2 / TE / plots. Used by `mitoribopy all` to "
            "kick off RNA-seq alignment in parallel with the Ribo-seq "
            "align stage so total wall time is dominated by the slower "
            "of the two assays instead of the sum. Writes "
            "`rna_counts.tsv` (and `rpf_counts.tsv` if Ribo FASTQs are "
            "supplied) to the output dir; subsequent rnaseq invocations "
            "with --upstream-rna-counts pick up the prebuilt matrix."
        ),
    )
    fastq.add_argument(
        "--allow-pseudo-replicates-for-demo-not-publication",
        dest="allow_pseudo_replicates",
        action="store_true",
        default=False,
        help=(
            "Opt INTO the auto-pseudo-replicate fallback for conditions "
            "with only 1 sample. Without this flag, an n=1 design is a "
            "hard error (exit 2) — the safe default for publication-grade "
            "DE. With this flag, each n=1 condition is split into rep1 / "
            "rep2 by FASTQ record parity so pyDESeq2 has n>=2 to fit "
            "dispersion on, BUT the resulting p-values and padj are NOT "
            "biologically defensible (the two halves are mechanical "
            "subsamples of one library). Use ONLY for demos / tutorials. "
            "The run_settings.json records pseudo_replicate_mode=true "
            "and an EXPLORATORY.md sidecar is written to the output dir."
        ),
    )
    # Deprecated alias kept for back-compat. Pseudo-replicate mode is
    # now opt-IN by default, so this flag is a no-op accepted only to
    # keep older config files / scripts working.
    fastq.add_argument(
        "--no-auto-pseudo-replicate",
        dest="no_auto_pseudo_replicate",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS,
    )

    # ----- Required in both flows ------------------------------------------
    gene = parser.add_argument_group("Gene identifiers")
    gene.add_argument(
        "--gene-id-convention",
        choices=["ensembl", "refseq", "hgnc", "mt_prefixed", "bare"],
        required=False,
        help=(
            "Gene identifier convention used in the FASTA / DE table. "
            "REQUIRED (no default): mismatched conventions silently "
            "produce zero-match runs."
        ),
    )
    gene.add_argument(
        "--organism",
        choices=["h.sapiens", "s.cerevisiae", "h", "y", "human", "yeast"],
        default="h.sapiens",
        help=(
            "Organism for the mt-mRNA registry. Use 'h.sapiens' for "
            "human, 's.cerevisiae' for budding yeast. Short / spelled-"
            "out forms ('h', 'y', 'human', 'yeast') are accepted as "
            "synonyms."
        ),
    )

    conditions = parser.add_argument_group(
        "Conditions (required in default flow; optional in --de-table flow)"
    )
    conditions.add_argument(
        "--condition-map",
        default=None,
        metavar="PATH",
        help=(
            "TSV with columns 'sample' and 'condition'. Required by the "
            "default flow (drives the pyDESeq2 contrast). In the "
            "--de-table flow it is optional and enables a replicate-"
            "based Ribo log2FC for delta-TE; without it the delta-TE "
            "rows carry only the mRNA log2FC and a "
            "'single_replicate_no_statistics' note."
        ),
    )
    conditions.add_argument(
        "--condition-a",
        default=None,
        metavar="NAME",
        help=(
            "Reference condition (denominator of the WT-vs-X contrast). "
            "Used as the baseline in DE / TE comparison plots. "
            "``--base-sample`` is the preferred spelling (matches the "
            "rpf config's ``base_sample`` key); both flags resolve to "
            "the same internal field."
        ),
    )
    conditions.add_argument(
        "--condition-b",
        default=None,
        metavar="NAME",
        help=(
            "Comparison condition (numerator of the contrast). "
            "``--compare-sample`` is the preferred alias."
        ),
    )
    conditions.add_argument(
        "--base-sample",
        default=None,
        metavar="NAME",
        dest="base_sample",
        help=(
            "Reference condition for the contrast and for all per-gene "
            "comparison plots (mRNA / RPF DE volcanoes, TE compare "
            "scatter, TE log2FC bar). Alias for ``--condition-a``; "
            "named to mirror the rpf config's ``base_sample`` key. "
            "When both ``--base-sample`` and ``--condition-a`` are "
            "provided they must agree."
        ),
    )
    conditions.add_argument(
        "--compare-sample",
        default=None,
        metavar="NAME",
        dest="compare_sample",
        help=(
            "Comparison condition (alias for ``--condition-b``). When "
            "both ``--compare-sample`` and ``--condition-b`` are "
            "provided they must agree."
        ),
    )

    out = parser.add_argument_group("Output")
    out.add_argument(
        "--output",
        required=False,
        metavar="DIR",
        help="Output directory for te.tsv, delta_te.tsv, and plots.",
    )

    # ----- Alternative flow: bring your own DE table -----------------------
    de = parser.add_argument_group(
        "Alternative inputs (bring your own DE table; mutually exclusive with --rna-fastq)"
    )
    de.add_argument(
        "--de-table",
        required=False,
        metavar="PATH",
        help=(
            "Pre-computed DESeq2 / Xtail / Anota2Seq results table "
            "(CSV or TSV). Use this when you already ran DE externally "
            "on the full transcriptome — the recommended path for "
            "publication-grade DE statistics. Requires --ribo-dir / "
            "--ribo-counts and one of --reference-gtf / "
            "--reference-checksum. Mutually exclusive with --rna-fastq."
        ),
    )
    de.add_argument(
        "--de-format",
        choices=["auto", "deseq2", "xtail", "anota2seq", "custom"],
        default="auto",
        help="DE table format. 'auto' detects from column headers.",
    )
    de.add_argument("--de-gene-col", default=None, metavar="NAME")
    de.add_argument("--de-log2fc-col", default=None, metavar="NAME")
    de.add_argument("--de-padj-col", default=None, metavar="NAME")
    de.add_argument("--de-basemean-col", default=None, metavar="NAME")

    ribo = parser.add_argument_group(
        "Ribo-seq inputs (--de-table flow)"
    )
    ribo.add_argument(
        "--ribo-dir",
        required=False,
        metavar="DIR",
        help=(
            "Directory produced by a prior 'mitoribopy rpf' run; expected "
            "to contain rpf_counts.tsv and run_settings.json / "
            "run_manifest.json with a recorded reference_checksum."
        ),
    )
    ribo.add_argument(
        "--ribo-counts",
        required=False,
        metavar="PATH",
        help=(
            "Explicit path to rpf_counts.tsv. Defaults to "
            "<ribo-dir>/rpf_counts.tsv when --ribo-dir is set."
        ),
    )

    ref = parser.add_argument_group(
        "Reference-consistency gate (--de-table flow; exactly one)"
    )
    ref.add_argument(
        "--reference-gtf",
        default=None,
        metavar="PATH",
        help=(
            "Reference GTF / FASTA used by RNA-seq; we hash this and "
            "verify it matches the hash recorded in the rpf run's "
            "manifest. EXACTLY ONE of --reference-gtf / "
            "--reference-checksum must be provided in the --de-table flow."
        ),
    )
    ref.add_argument(
        "--reference-checksum",
        default=None,
        metavar="SHA256",
        help=(
            "Precomputed SHA-256 of the shared reference (use instead of "
            "--reference-gtf when the file is not on this host)."
        ),
    )

    return parser


# ---------- helpers ---------------------------------------------------------


def _resolve_de_column_map(args) -> "object | None":
    """Return a user-provided DeColumnMap when --de-format custom, else None."""
    from ..rnaseq._types import DeColumnMap

    if args.de_format != "custom":
        return None
    if not args.de_gene_col:
        raise ValueError("--de-format custom requires --de-gene-col.")
    return DeColumnMap(
        gene_id=args.de_gene_col,
        log2fc=args.de_log2fc_col,
        padj=args.de_padj_col,
        basemean=args.de_basemean_col,
    )


def _load_condition_map(path: Path | None) -> dict[str, str]:
    if not path:
        return {}
    path = Path(path)
    mapping: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample = (row.get("sample") or "").strip()
            condition = (row.get("condition") or "").strip()
            if sample and condition:
                mapping[sample] = condition
    return mapping


_TE_TSV_COLUMNS: tuple[str, ...] = (
    "sample_id",
    "condition",
    "assay",
    "gene",
    "rpf_count",
    "rna_abundance",
    "te",
    "log2_te",
    "note",
)


_DELTA_TE_TSV_COLUMNS: tuple[str, ...] = (
    "gene",
    "base_condition",
    "compare_condition",
    "mrna_log2fc",
    "rpf_log2fc",
    "delta_te_log2",
    "padj_mrna",
    "padj_rpf",
    "padj_delta_te",
    "method",
    "note",
)


def _fmt_optional(value, fmt: str = ".6g") -> str:
    if value is None:
        return ""
    return format(value, fmt)


def _write_te_table(rows: Iterable[TeRow], path: Path) -> None:
    """Write ``te.tsv`` (schema 2.0, §5 spec)."""
    from ..io.schema_versions import schema_header_line

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(schema_header_line("te.tsv"))
        handle.write("\t".join(_TE_TSV_COLUMNS) + "\n")
        for row in rows:
            cells = {
                "sample_id": row.sample_id,
                "condition": row.condition or "",
                "assay": row.assay,
                "gene": row.gene,
                "rpf_count": str(row.rpf_count),
                "rna_abundance": _fmt_optional(row.rna_abundance),
                "te": _fmt_optional(row.te),
                "log2_te": _fmt_optional(row.log2_te),
                "note": row.note,
            }
            handle.write(
                "\t".join(cells[c] for c in _TE_TSV_COLUMNS) + "\n"
            )


def _write_delta_te_table(rows: Iterable[DTeRow], path: Path) -> None:
    """Write ``delta_te.tsv`` (schema 2.0, §5 spec)."""
    from ..io.schema_versions import schema_header_line

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(schema_header_line("delta_te.tsv"))
        handle.write("\t".join(_DELTA_TE_TSV_COLUMNS) + "\n")
        for row in rows:
            cells = {
                "gene": row.gene,
                "base_condition": row.base_condition or "",
                "compare_condition": row.compare_condition or "",
                "mrna_log2fc": _fmt_optional(row.mrna_log2fc),
                "rpf_log2fc": _fmt_optional(row.rpf_log2fc),
                "delta_te_log2": _fmt_optional(row.delta_te_log2),
                "padj_mrna": _fmt_optional(row.padj_mrna),
                "padj_rpf": _fmt_optional(row.padj_rpf),
                "padj_delta_te": _fmt_optional(row.padj_delta_te),
                "method": row.method,
                "note": row.note,
            }
            handle.write(
                "\t".join(cells[c] for c in _DELTA_TE_TSV_COLUMNS) + "\n"
            )


# ---------- from-FASTQ orchestrator -----------------------------------------


def _conditions_with_one_sample(samples: list, condition_map: dict) -> list[str]:
    """Return sorted condition names that map to exactly one sample.

    Samples whose ``sample`` field is missing from ``condition_map`` are
    ignored: the rest of the pipeline already errors on them downstream.
    """
    counts: dict[str, int] = {}
    for s in samples:
        cond = condition_map.get(s.sample)
        if cond is None:
            continue
        counts[cond] = counts.get(cond, 0) + 1
    return sorted(c for c, n in counts.items() if n == 1)


def _emit_pseudo_replicate_banner(
    conditions: list[str], *, when: str
) -> None:
    bar = "=" * 72
    head = (
        "PSEUDO-REPLICATE MODE ACTIVE — EXPLORATORY ONLY"
        if when == "start"
        else "PSEUDO-REPLICATE RUN COMPLETE — EXPLORATORY ONLY"
    )
    body = (
        f"Conditions split by FASTQ-record parity: {', '.join(conditions)}.\n"
        "padj / p-values in plots and tables come from mechanical halves "
        "of one library and are NOT biologically defensible. Use these "
        "outputs for tutorials, smoke tests, or QC only — never for "
        "publication or for a 'significant gene' claim. To produce "
        "publication-grade DE, supply biological replicates or run DE "
        "externally and pass --de-table."
    )
    sys.stderr.write(f"\n{bar}\n{head}\n{bar}\n{body}\n{bar}\n\n")


def _write_exploratory_sidecar(
    output_dir: Path, conditions: list[str]
) -> None:
    """Drop EXPLORATORY.md next to te.tsv / delta_te.tsv."""
    sidecar = Path(output_dir) / "EXPLORATORY.md"
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    sidecar.write_text(
        "# EXPLORATORY RUN — pseudo-replicate mode active\n\n"
        f"Conditions split by FASTQ-record parity: {', '.join(conditions)}\n\n"
        "This run was launched with `--allow-pseudo-replicates-for-"
        "demo-not-publication`. The condition(s) above had only one "
        "sample, so the FASTQ for each was streamed-split in two by "
        "record parity to give pyDESeq2 the n>=2 it needs to fit "
        "dispersion. The two 'replicates' are therefore mechanical "
        "halves of the same library, NOT biological replicates.\n\n"
        "## Outputs to NOT cite\n\n"
        "- p-values, padj, and 'significant' markers in `de_table.tsv`, "
        "`rpf_de_table.tsv`, `delta_te.tsv`, and the volcano plots.\n"
        "- Dispersion estimates in pyDESeq2 fits.\n\n"
        "## Outputs that remain meaningful\n\n"
        "- Read counts, alignment statistics, and the FASTA / FASTQ "
        "provenance recorded in `run_settings.json`.\n"
        "- Descriptive log2 fold-change point estimates (treat as "
        "exploratory; do not gate publication claims on them).\n\n"
        "To produce defensible statistics, supply biological "
        "replicates or run DE externally and re-run `mitoribopy "
        "rnaseq` with `--de-table`.\n",
        encoding="utf-8",
    )


def _auto_split_singletons(
    samples: list,
    condition_map: dict,
    workdir: Path,
    label: str,
) -> tuple[list, dict]:
    """Auto-split any sample whose condition has only one member.

    Returns ``(updated_samples, updated_condition_map)``. Both inputs
    are left intact; the caller swaps in the returned objects. Emits a
    stderr WARNING per sample that gets split so the behaviour is
    never silent.
    """
    from ..rnaseq.split_replicates import split_sample_into_pseudo_replicates

    by_condition: dict = {}
    for s in samples:
        cond = condition_map.get(s.sample)
        if cond is None:
            continue
        by_condition.setdefault(cond, []).append(s)

    new_samples: list = []
    new_map = dict(condition_map)
    for s in samples:
        cond = condition_map.get(s.sample)
        if cond is not None and len(by_condition.get(cond, [])) == 1:
            rep1, rep2 = split_sample_into_pseudo_replicates(s, workdir)
            new_samples.extend([rep1, rep2])
            new_map[rep1.sample] = cond
            new_map[rep2.sample] = cond
            sys.stderr.write(
                f"[mitoribopy rnaseq] WARNING: {label} condition "
                f"{cond!r} has only 1 sample ({s.sample!r}); auto-"
                f"splitting reads by record parity into pseudo-"
                f"replicates {rep1.sample!r} / {rep2.sample!r}. These "
                "are mechanical halves of the same library, NOT "
                "biological replicates — DESeq2 dispersion estimates "
                "will be artificially low and the resulting padj is "
                "NOT publication-grade.\n"
            )
        else:
            new_samples.append(s)
    return new_samples, new_map


def _provenance_entry(result) -> dict:
    rk = result.resolved_kit
    return {
        "sample": result.sample,
        "paired": result.paired,
        "total_reads": result.total_reads,
        "aligned_reads": result.aligned_reads,
        "kit": rk.kit,
        "adapter": rk.adapter,
        "umi_length": rk.umi_length,
        "umi_position": rk.umi_position,
    }


def _synth_rna_results_from_upstream_counts(
    counts_path: Path,
) -> "list":
    """Synthesise per-sample RNA `SampleAlignmentResult`s from rna counts.

    Mirrors :func:`_synth_ribo_results_from_upstream_counts` but for the
    RNA assay. Used by ``mitoribopy all`` when it has already produced
    ``rna_counts.tsv`` in a prior parallel align-only invocation: the
    rnaseq DE stage skips bowtie2 entirely and consumes the matrix
    directly. The matrix is wide (genes × samples) and we pivot to the
    per-sample dict shape ``SampleAlignmentResult.counts`` expects.
    """
    import pandas as pd

    from ..align._types import ResolvedKit
    from ..rnaseq.alignment import SampleAlignmentResult

    counts_path = Path(counts_path)
    # comment="#" skips the leading `# schema_version: X.Y` line that
    # every MitoRiboPy TSV writer prepends; without it pandas treats
    # that comment as the first data row, the column layout collapses
    # to one column, and the downstream DESeq2 fit raises "No 'rna'
    # samples available for DESeq2."
    df = pd.read_csv(counts_path, sep="\t", comment="#")
    if df.empty:
        return []
    # Conventional layout: first column is gene id, remaining columns
    # are samples. We tolerate either 'gene' or 'gene_id' as the index.
    gene_col = df.columns[0]
    sample_cols = [c for c in df.columns[1:] if c.lower() != "gene_id"]

    sentinel_kit = ResolvedKit(
        kit="upstream-rna-counts",
        adapter=None,
        umi_length=0,
        umi_position="5p",
    )
    sentinel_bam = Path(str(counts_path))

    results: list = []
    for sample in sample_cols:
        per_sample: dict[str, int] = {
            str(g): int(v)
            for g, v in zip(df[gene_col].astype(str), df[sample].fillna(0))
        }
        total = sum(per_sample.values())
        results.append(
            SampleAlignmentResult(
                sample=str(sample),
                bam_path=sentinel_bam,
                counts=per_sample,
                paired=False,
                total_reads=total,
                aligned_reads=total,
                resolved_kit=sentinel_kit,
            )
        )
    return results


def _synth_ribo_results_from_upstream_counts(
    counts_path: Path,
) -> "list":
    """Synthesise per-sample Ribo `SampleAlignmentResult`s from rpf counts.

    When `mitoribopy all` has already produced `rpf_counts.tsv`, the
    rnaseq stage can build a counts matrix without re-aligning the
    Ribo FASTQs. The synthesised results carry only what the
    downstream `build_sample_sheet` + `run_deseq2(assay='ribo')` and
    `_provenance_entry` callers actually read: `sample`, `counts`,
    `paired`, `total_reads`, `aligned_reads`, `resolved_kit`. The BAM
    path is a sentinel (no BAM exists for these samples in the rnaseq
    workdir).
    """
    from ..align._types import ResolvedKit
    from ..rnaseq.alignment import SampleAlignmentResult
    from ..rnaseq.counts import load_ribo_counts

    # load_ribo_counts returns {gene: {sample: count}}; pivot to
    # {sample: {gene: count}} so we can build per-sample counts dicts.
    by_gene = load_ribo_counts(counts_path)
    by_sample: dict[str, dict[str, int]] = {}
    for gene, sample_counts in by_gene.items():
        for sample, n in sample_counts.items():
            by_sample.setdefault(sample, {})[gene] = int(n)

    sentinel_kit = ResolvedKit(
        kit="upstream-rpf-counts",
        adapter=None,
        umi_length=0,
        umi_position="5p",
    )
    sentinel_bam = Path(str(counts_path))

    results: list = []
    for sample in sorted(by_sample):
        per_sample = by_sample[sample]
        total = sum(per_sample.values())
        results.append(
            SampleAlignmentResult(
                sample=sample,
                bam_path=sentinel_bam,
                counts=per_sample,
                paired=False,
                total_reads=total,
                aligned_reads=total,
                resolved_kit=sentinel_kit,
            )
        )
    return results


def _run_from_fastq(args, output_dir: Path) -> int:
    """Trim / align / count / DESeq2 from raw FASTQs.

    Mutates ``args`` in place so the existing pre-computed-DE path can
    fall through unchanged: sets ``args.de_table``, ``args.de_format``,
    ``args.ribo_counts``, plus a few private (underscore-prefixed)
    attributes used to skip the reference-consistency hard-gate and
    enrich ``run_settings.json``.

    Returns ``0`` and continues only when Ribo-seq FASTQs are also
    provided; for an RNA-only run it short-circuits after writing
    ``de_table.tsv`` + ``run_settings.json`` and returns ``-1`` to
    signal "stop here".

    Reuse path (P0.2): when ``args.upstream_rpf_counts`` is set
    (typically by ``mitoribopy all`` after the rpf stage just wrote
    ``rpf_counts.tsv``) AND ``args.recount_ribo_fastq`` is False, the
    Ribo side is built from the upstream counts file instead of
    re-running cutadapt + bowtie2 + counting on the Ribo FASTQs.
    """
    from ..rnaseq.alignment import (
        SampleAlignmentResult,
        align_sample,
        build_bowtie2_index,
        write_counts_matrix,
        write_long_counts,
    )
    from ..rnaseq.de_analysis import (
        build_sample_sheet,
        deseq2_to_de_table,
        run_deseq2,
        write_de_table_tsv,
    )
    from ..rnaseq.fastq_pairing import detect_samples, enumerate_fastqs

    output_dir.mkdir(parents=True, exist_ok=True)
    workdir = Path(args.workdir) if args.workdir else output_dir / "work"
    workdir.mkdir(parents=True, exist_ok=True)

    upstream_rna_counts = getattr(args, "upstream_rna_counts", None)
    align_only = bool(getattr(args, "align_only", False))
    reuse_upstream_rna = bool(upstream_rna_counts) and not align_only

    if reuse_upstream_rna:
        # The orchestrator (`mitoribopy all`) handed us the prebuilt
        # rna_counts.tsv from a parallel align-only invocation. Skip
        # FASTQ enumeration entirely and synthesise the result objects.
        rna_paths = []
        rna_samples = []
        sys.stderr.write(
            "[mitoribopy rnaseq] reusing upstream rna counts "
            f"({upstream_rna_counts}); skipping RNA-FASTQ re-alignment.\n"
        )
        # The align-only worker writes condition_map.augmented.tsv next
        # to rna_counts.tsv when it had to split pseudo-replicates.
        # Prefer that augmented map over the user-supplied one so the
        # DE step sees the rep1/rep2 entries the count matrix is keyed
        # on.
        augmented_map = output_dir / "condition_map.augmented.tsv"
        if augmented_map.is_file() and (
            not args.condition_map or Path(args.condition_map) != augmented_map
        ):
            sys.stderr.write(
                "[mitoribopy rnaseq] using augmented condition map "
                f"{augmented_map} (written by upstream align-only step).\n"
            )
            args.condition_map = str(augmented_map)
    else:
        rna_paths = enumerate_fastqs(args.rna_fastq)
        rna_samples = detect_samples(rna_paths)

    upstream_rpf_counts = getattr(args, "upstream_rpf_counts", None)
    recount = bool(getattr(args, "recount_ribo_fastq", False))
    reuse_upstream_rpf = bool(upstream_rpf_counts) and not recount

    if reuse_upstream_rpf:
        # The orchestrator (`mitoribopy all`) handed us the rpf stage's
        # rpf_counts.tsv. Skip Ribo FASTQ enumeration entirely; downstream
        # we will build a synthetic ribo_results list from the counts
        # file. The pseudo-replicate logic still runs only on rna_samples.
        ribo_samples = []
        sys.stderr.write(
            "[mitoribopy rnaseq] reusing upstream rpf counts "
            f"({upstream_rpf_counts}); skipping Ribo-FASTQ re-alignment. "
            "Pass --recount-ribo-fastq to force a second pass.\n"
        )
    elif args.ribo_fastq:
        ribo_paths = enumerate_fastqs(args.ribo_fastq)
        ribo_samples = detect_samples(ribo_paths)
    else:
        ribo_samples = []

    condition_map = _load_condition_map(
        Path(args.condition_map) if args.condition_map else None
    )

    # Pseudo-replicate gate.
    #
    # Pseudo-replicates auto-split a condition with n=1 into two halves
    # by FASTQ record parity, which makes pyDESeq2 happy but produces
    # p-values that look real and are not — the two halves are
    # mechanical subsamples of one library. The default is opt-IN:
    # detect singletons up-front, fail with a clear message unless the
    # user has passed --allow-pseudo-replicates-for-demo-not-publication.
    if getattr(args, "no_auto_pseudo_replicate", False):
        sys.stderr.write(
            "[mitoribopy rnaseq] DEPRECATED: --no-auto-pseudo-replicate "
            "is a no-op (pseudo-replicate mode is opt-in via "
            "--allow-pseudo-replicates-for-demo-not-publication). "
            "Drop the flag from your scripts; it will be removed in a "
            "future release.\n"
        )

    allow_pseudo = bool(getattr(args, "allow_pseudo_replicates", False))
    rna_singletons = _conditions_with_one_sample(rna_samples, condition_map)
    ribo_singletons = (
        _conditions_with_one_sample(ribo_samples, condition_map)
        if ribo_samples else []
    )
    pseudo_conditions = sorted(set(rna_singletons) | set(ribo_singletons))

    if pseudo_conditions and not allow_pseudo:
        sys.stderr.write(
            "[mitoribopy rnaseq] ERROR: the following condition(s) have "
            "only 1 sample and pyDESeq2 cannot fit dispersion on n=1 "
            f"designs: {', '.join(repr(c) for c in pseudo_conditions)}.\n"
            "  Resolve by ONE of:\n"
            "    1. supply biological replicates (recommended);\n"
            "    2. run external DE and use the --de-table flow;\n"
            "    3. pass --allow-pseudo-replicates-for-demo-not-publication\n"
            "       for a non-publication exploratory run (FASTQ-record\n"
            "       parity halves; padj/p-values are NOT biologically\n"
            "       defensible).\n"
        )
        return 2

    original_map = dict(condition_map)
    pseudo_replicate_mode = bool(pseudo_conditions and allow_pseudo)
    if pseudo_replicate_mode:
        _emit_pseudo_replicate_banner(pseudo_conditions, when="start")
        rna_samples, condition_map = _auto_split_singletons(
            rna_samples, condition_map, workdir / "rna_split", "RNA-seq"
        )
        if ribo_samples:
            ribo_samples, condition_map = _auto_split_singletons(
                ribo_samples, condition_map, workdir / "ribo_split", "Ribo-seq"
            )
        _write_exploratory_sidecar(output_dir, pseudo_conditions)

    # Persist the (possibly augmented) condition map so the downstream
    # delta-TE step picks up the new rep1 / rep2 entries when it re-
    # reads --condition-map from disk.
    if condition_map != original_map:
        augmented = output_dir / "condition_map.augmented.tsv"
        with augmented.open("w", encoding="utf-8") as handle:
            handle.write("sample\tcondition\n")
            for s, c in sorted(condition_map.items()):
                handle.write(f"{s}\t{c}\n")
        args.condition_map = str(augmented)

    if args.bowtie2_index:
        bt2_index = Path(args.bowtie2_index)
    else:
        bt2_index = build_bowtie2_index(
            Path(args.reference_fasta), workdir / "bt2_cache"
        )

    rna_results: list[SampleAlignmentResult] = []
    if reuse_upstream_rna:
        rna_results = _synth_rna_results_from_upstream_counts(
            Path(upstream_rna_counts)
        )
    else:
        for s in rna_samples:
            rna_results.append(
                align_sample(
                    s,
                    bt2_index=bt2_index,
                    workdir=workdir / "rna",
                    threads=args.align_threads,
                )
            )

    ribo_results: list[SampleAlignmentResult] = []
    if reuse_upstream_rpf:
        # Build a synthetic ribo_results from rpf_counts.tsv. The
        # second pyDESeq2 fit (assay='ribo') still runs, just on the
        # rpf-stage counts instead of a fresh re-alignment.
        ribo_results = _synth_ribo_results_from_upstream_counts(
            Path(upstream_rpf_counts)
        )
    else:
        for s in ribo_samples:
            ribo_results.append(
                align_sample(
                    s,
                    bt2_index=bt2_index,
                    workdir=workdir / "ribo",
                    threads=args.align_threads,
                )
            )

    if rna_results and not reuse_upstream_rna:
        write_counts_matrix(rna_results, output_dir / "rna_counts.tsv")
    if ribo_results and not reuse_upstream_rpf:
        # In the normal (non-reuse) path we materialise rpf_counts.tsv
        # under the rnaseq output dir so downstream consumers find it
        # at a stable location. In the reuse path the upstream file is
        # already on disk under <run_root>/rpf/rpf_counts.tsv and we
        # point args.ribo_counts at it directly below.
        write_counts_matrix(ribo_results, output_dir / "rpf_counts_matrix.tsv")
        write_long_counts(ribo_results, output_dir / "rpf_counts.tsv")

    if align_only:
        # The `mitoribopy all` orchestrator launched us in parallel with
        # the Ribo align stage solely to produce the counts matrices it
        # needs for the downstream DE/TE step. Exit cleanly before
        # touching pyDESeq2 / pseudo-replicate logic / plots.
        # Return -1 (the caller's "fully done; short-circuit" sentinel)
        # so the reference-consistency gate further down does NOT fire —
        # the rpf side hasn't written its output yet, so the gate would
        # spuriously reject this worker.
        sys.stderr.write(
            "[mitoribopy rnaseq] --align-only complete; wrote "
            f"{output_dir / 'rna_counts.tsv'}"
            + (
                f" + {output_dir / 'rpf_counts.tsv'}"
                if (ribo_results and not reuse_upstream_rpf)
                else ""
            )
            + ". Skipping DE/TE/plots.\n"
        )
        return -1

    samples, counts_df, metadata_df = build_sample_sheet(
        rna_results, ribo_results, condition_map
    )
    plot_pca(
        counts_df, metadata_df, output_dir / "plots" / "sample_pca.png"
    )
    results_df = run_deseq2(
        counts_df,
        metadata_df,
        contrast_factor="condition",
        contrast_a=args.condition_a,
        contrast_b=args.condition_b,
        assay="rna",
    )
    de_table = deseq2_to_de_table(results_df)
    de_table_path = output_dir / "de_table.tsv"
    write_de_table_tsv(de_table, de_table_path)

    # Run pyDESeq2 a second time on the Ribo-seq subset so we get a
    # statistically grounded RPF log2FC + padj. Only attempted when at
    # least one Ribo sample landed in the metadata; otherwise skipped
    # silently. ``run_deseq2`` raises ValueError if the subset has
    # fewer than two condition levels — caught here so the run continues
    # without an RPF DE table (the fallback plots still work).
    rpf_de_table_path: Path | None = None
    if ribo_results:
        try:
            rpf_results_df = run_deseq2(
                counts_df,
                metadata_df,
                contrast_factor="condition",
                contrast_a=args.condition_a,
                contrast_b=args.condition_b,
                assay="ribo",
            )
        except ValueError as exc:
            sys.stderr.write(
                "[mitoribopy rnaseq] WARNING: Skipping Ribo-seq DESeq2 "
                f"({exc}); RPF volcano plot will be omitted.\n"
            )
        else:
            rpf_de_table = deseq2_to_de_table(rpf_results_df)
            rpf_de_table_path = output_dir / "rpf_de_table.tsv"
            write_de_table_tsv(rpf_de_table, rpf_de_table_path)

    reference_checksum = compute_reference_checksum(Path(args.reference_fasta))

    args.de_table = str(de_table_path)
    args.de_format = "deseq2"
    if reuse_upstream_rpf:
        args.ribo_counts = str(upstream_rpf_counts)
    elif ribo_results:
        args.ribo_counts = str(output_dir / "rpf_counts.tsv")
    else:
        args.ribo_counts = None
    args._from_fastq = True
    args._reference_checksum = reference_checksum
    args._rpf_de_table_path = (
        str(rpf_de_table_path) if rpf_de_table_path is not None else None
    )
    args._pseudo_replicate_mode = pseudo_replicate_mode
    args._pseudo_replicate_conditions = list(pseudo_conditions) if pseudo_replicate_mode else []
    if reuse_upstream_rpf:
        provenance_mode = "from-fastq-rna-with-upstream-rpf-counts"
    elif ribo_results:
        provenance_mode = "from-fastq"
    else:
        provenance_mode = "from-fastq-rna-only"
    args._fastq_provenance = {
        "mode": provenance_mode,
        "reference_fasta": str(args.reference_fasta),
        "bowtie2_index": str(bt2_index),
        "upstream_rpf_counts": (
            str(upstream_rpf_counts) if reuse_upstream_rpf else None
        ),
        "rna_samples": [_provenance_entry(r) for r in rna_results],
        "ribo_samples": [_provenance_entry(r) for r in ribo_results],
        "pseudo_replicate_mode": pseudo_replicate_mode,
        "pseudo_replicate_conditions": (
            list(pseudo_conditions) if pseudo_replicate_mode else []
        ),
    }

    if not ribo_results:
        # RNA-only short-circuit: write a slim run_settings.json and stop.
        slim_settings = {
            "subcommand": "rnaseq",
            "mitoribopy_version": __version__,
            "rnaseq_mode": getattr(args, "_rnaseq_mode", "from_fastq"),
            "from_fastq": args._fastq_provenance,
            "reference_checksum": reference_checksum,
            "de_table": args.de_table,
            "de_format_resolved": "deseq2",
            "pseudo_replicate_mode": pseudo_replicate_mode,
            "pseudo_replicate_conditions": (
                list(pseudo_conditions) if pseudo_replicate_mode else []
            ),
            "sample_sheet": getattr(args, "_sample_sheet_path", None),
        }
        (output_dir / "run_settings.json").write_text(
            json.dumps(slim_settings, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        if pseudo_replicate_mode:
            _emit_pseudo_replicate_banner(pseudo_conditions, when="end")
        return -1

    return 0


# ---------- sample sheet -> per-flag input derivation -----------------------


def _apply_sample_sheet(args) -> int:
    """Derive --rna-fastq / --ribo-fastq / --condition-map from a sheet.

    Mutates ``args`` in place. Returns 0 on success or 2 on user-facing
    errors (mixing with mutually-exclusive flags, sheet validation
    failure, unknown fields). The caller is expected to propagate the
    return code.
    """
    from ..sample_sheet import (
        SampleSheetError,
        check_sheet_conflicts,
        format_sheet_conflict_error,
        load_sample_sheet,
    )

    # Single source of truth: every per-flag input the sheet shadows
    # gets reported via the shared formatter so the wording matches
    # what `mitoribopy all` emits for the equivalent YAML conflict.
    args_as_dict = {
        "rna_fastq": args.rna_fastq,
        "ribo_fastq": args.ribo_fastq,
        "condition_map": args.condition_map,
        "de_table": args.de_table,
        "ribo_dir": getattr(args, "ribo_dir", None),
        "ribo_counts": getattr(args, "ribo_counts", None),
    }
    conflict_keys = check_sheet_conflicts(
        args_as_dict,
        conflict_keys=tuple(args_as_dict.keys()),
    )
    if conflict_keys:
        flag_names = [f"--{k.replace('_', '-')}" for k in conflict_keys]
        sys.stderr.write(
            format_sheet_conflict_error(
                "mitoribopy rnaseq",
                flag_names,
                sheet_label="--sample-sheet",
            )
            + "\n"
        )
        return 2

    try:
        sheet = load_sample_sheet(args.sample_sheet)
    except SampleSheetError as exc:
        sys.stderr.write(f"[mitoribopy rnaseq] ERROR: {exc}\n")
        return 2

    rna_rows = sheet.by_assay("rna")
    ribo_rows = sheet.by_assay("ribo")

    if not rna_rows:
        sys.stderr.write(
            "[mitoribopy rnaseq] ERROR: --sample-sheet contains no "
            "active rows with assay='rna'. The from-FASTQ flow needs at "
            "least one RNA-seq sample.\n"
        )
        return 2

    args.rna_fastq = [str(r.fastq_1) for r in rna_rows]
    args.ribo_fastq = (
        [str(r.fastq_1) for r in ribo_rows] if ribo_rows else None
    )

    # Materialise a condition-map TSV next to the sheet so the
    # downstream condition_map loader (which re-reads from disk after
    # the pseudo-replicate splitter possibly augments it) finds a
    # canonical file. This keeps the on-disk artefact discoverable for
    # users debugging a run.
    cmap_path = Path(args.sample_sheet).with_suffix(".condition_map.tsv")
    with cmap_path.open("w", encoding="utf-8") as handle:
        handle.write("sample\tcondition\n")
        for row in sheet.active():
            handle.write(f"{row.sample_id}\t{row.condition}\n")
    args.condition_map = str(cmap_path)
    args._sample_sheet_path = str(sheet.path)
    return 0


# ---------- orchestrator ----------------------------------------------------


def run(argv: Iterable[str]) -> int:
    argv_list = list(argv)
    parser = build_parser()

    # Pre-parse --config so we can fold YAML / JSON / TOML values into
    # the parser's defaults BEFORE the real parse pass. CLI flags still
    # win on conflict because they appear in argv_list as explicitly-
    # set values that override the new defaults. Same pattern as
    # `mitoribopy align` / `mitoribopy rpf`.
    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument("--config", default=None)
    pre_args, _ = pre_parser.parse_known_args(argv_list)
    if pre_args.config:
        try:
            cfg = common.load_config_file(pre_args.config)
        except (FileNotFoundError, RuntimeError, ValueError) as exc:
            print(f"[mitoribopy rnaseq] ERROR: {exc}", file=sys.stderr)
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
                "RNASEQ",
                "Ignoring unknown --config keys: " + ", ".join(unknown),
            )
        if applied:
            parser.set_defaults(**applied)

    args = parser.parse_args(argv_list)
    common.apply_common_arguments(args)

    # --sample-sheet is a higher-level input that derives --rna-fastq /
    # --ribo-fastq / --condition-map for the from-FASTQ flow. Reject
    # mixing it with the per-flag alternatives so users do not get a
    # surprising shadow when both define the same field.
    if args.sample_sheet is not None:
        rc = _apply_sample_sheet(args)
        if rc != 0:
            return rc

    # --base-sample / --compare-sample are aliases for --condition-a /
    # --condition-b. Reconcile here so the rest of the run uses a single
    # canonical pair (args.condition_a / args.condition_b). When both an
    # alias and the legacy flag are passed they must agree.
    alias_pairs = (
        ("base_sample", "condition_a"),
        ("compare_sample", "condition_b"),
    )
    for alias_attr, canonical_attr in alias_pairs:
        alias_value = getattr(args, alias_attr, None)
        canonical_value = getattr(args, canonical_attr, None)
        if alias_value and canonical_value and alias_value != canonical_value:
            print(
                f"[mitoribopy rnaseq] ERROR: --{alias_attr.replace('_', '-')}"
                f"={alias_value!r} conflicts with "
                f"--{canonical_attr.replace('_', '-')}={canonical_value!r}; "
                "pick one.",
                file=sys.stderr,
            )
            return 2
        if alias_value and not canonical_value:
            setattr(args, canonical_attr, alias_value)
        elif canonical_value and not alias_value:
            setattr(args, alias_attr, canonical_value)

    # --- Default-flow vs --de-table flow dispatch -----------------------
    # The default flow (from raw FASTQ) is selected unless the user
    # passes --de-table, in which case we run the alternative pre-
    # computed-DE flow. Passing both is a hard error so users do not
    # accidentally silently shadow one with the other. The explicit
    # --rnaseq-mode flag (when set) overrides inference and is the
    # documented public contract; inference remains for back-compat.
    mode, mode_error = _resolve_rnaseq_mode(args)
    if mode_error is not None:
        print(f"[mitoribopy rnaseq] ERROR: {mode_error}", file=sys.stderr)
        return 2
    args._rnaseq_mode = mode
    use_de_table = mode == "de_table"
    if mode == "from_fastq":
        _emit_exploratory_mode_banner()
    if args.dry_run:
        base = args.condition_a or "<--base-sample>"
        compare = args.condition_b or "<--compare-sample>"
        if use_de_table:
            actions = [
                "resolve --de-format (auto-detect from headers) and load DE table",
                "hash-verify --reference-gtf / --reference-checksum against "
                "the rpf run's recorded reference_checksum (HARD FAIL on mismatch)",
                "load rpf_counts.tsv from --ribo-dir or --ribo-counts",
                f"match DE gene_ids against mt-mRNA registry ({args.organism}) "
                f"using --gene-id-convention {args.gene_id_convention}",
                "compute TE per sample per gene = (RPF + 0.5) / (mRNA + 0.5)",
                "compute delta-TE = log2(RPF_fc) - log2(mRNA_fc) with "
                "replicate-based Ribo log2FC when --condition-map is given",
                f"emit te.tsv, delta_te.tsv, and 6 always-on plots "
                f"(mrna_vs_rpf, delta_te_volcano, ma, de_volcano_mrna "
                f"({base} vs {compare}), te_bar_by_condition, "
                f"te_heatmap); plus te_compare_scatter and "
                "te_log2fc_bar when --condition-map is provided. "
                "Each PNG is co-emitted as an SVG sidecar at the "
                "same stem (.png + .svg) for vector / Illustrator use",
            ]
        else:
            actions = [
                "auto-detect SE / PE for every --rna-fastq + --ribo-fastq input",
                "auto-detect adapter + UMI per sample (PE+UMI is NotImplementedError)",
                "build (or reuse) bowtie2 index from --reference-fasta "
                "(content-addressed cache: <workdir>/bt2_cache/transcriptome_<sha12>)",
                "trim + align every sample (cutadapt + bowtie2); count primary "
                "mapped reads per transcript",
                "fail fast on any condition with n=1 (publication-safe "
                "default); pass --allow-pseudo-replicates-for-demo-not-"
                "publication to opt into FASTQ-record-parity splitting "
                "for demos / tutorials only",
                f"build sample sheet, run pyDESeq2 contrast on the RNA "
                f"subset ({base} vs {compare})",
                "serialise mRNA pyDESeq2 result to de_table.tsv (DESeq2 schema)",
                f"run pyDESeq2 a second time on the Ribo subset ({base} "
                f"vs {compare}) and write rpf_de_table.tsv (skipped if "
                "the Ribo subset has fewer than two condition levels)",
                f"match DE gene_ids against mt-mRNA registry ({args.organism}) "
                f"using --gene-id-convention {args.gene_id_convention}",
                "compute TE + delta-TE; emit te.tsv, delta_te.tsv, and "
                "10 plots (sample_pca, mrna_vs_rpf, delta_te_volcano, "
                "ma, de_volcano_mrna, de_volcano_rpf, te_bar_by_condition, "
                "te_heatmap, te_compare_scatter, te_log2fc_bar). Each "
                "PNG is co-emitted as an SVG sidecar at the same stem "
                "(.png + .svg) for vector / Illustrator use; PNGs render "
                "at 300 dpi with the Okabe-Ito colour-blind-safe palette",
                "also write the wide / long counts matrices "
                "(rna_counts.tsv, rpf_counts_matrix.tsv, rpf_counts.tsv) "
                "and an rpf_de_table.tsv from a second pyDESeq2 fit on "
                "the Ribo subset",
                "record FASTA SHA256 + per-sample provenance in run_settings.json",
            ]
        return common.emit_dry_run("rnaseq", actions)

    from_fastq = args.rna_fastq is not None
    if from_fastq and use_de_table:
        print(
            "[mitoribopy rnaseq] ERROR: --rna-fastq and --de-table are "
            "mutually exclusive. Pick one: --rna-fastq for the default "
            "from-FASTQ flow (cutadapt + bowtie2 + pyDESeq2 inside the "
            "subcommand), or --de-table to bring your own pre-computed "
            "DE table from an external DESeq2 / Xtail / Anota2Seq run.",
            file=sys.stderr,
        )
        return 2

    missing: list[str] = []
    if use_de_table:
        # Alternative flow: DE table + prior rpf run.
        if not args.gene_id_convention:
            missing.append("--gene-id-convention")
        if not args.ribo_dir and not args.ribo_counts:
            missing.append("--ribo-dir or --ribo-counts")
        if not args.output:
            missing.append("--output")
        if not (args.reference_gtf or args.reference_checksum):
            missing.append("--reference-gtf or --reference-checksum")
    else:
        # Default flow: from raw FASTQ.
        if not args.rna_fastq:
            missing.append("--rna-fastq")
        if not args.reference_fasta:
            missing.append("--reference-fasta")
        if not args.gene_id_convention:
            missing.append("--gene-id-convention")
        if not args.output:
            missing.append("--output")
        if not args.condition_map:
            missing.append("--condition-map")
        # condition_a / condition_b are populated from --base-sample /
        # --compare-sample by the alias-reconciliation step above, so
        # we only need to flag the canonical pair here. Naming the
        # alias too keeps the error actionable for users who read
        # docs starting from --base-sample.
        if not args.condition_a:
            missing.append("--condition-a (or --base-sample)")
        if not args.condition_b:
            missing.append("--condition-b (or --compare-sample)")
    if missing:
        flow = "--de-table flow" if use_de_table else "default flow (from raw FASTQ)"
        print(
            f"[mitoribopy rnaseq] ERROR: missing required argument(s) for "
            f"the {flow}: " + ", ".join(missing),
            file=sys.stderr,
        )
        if not use_de_table and not args.rna_fastq:
            print(
                "[mitoribopy rnaseq] HINT: pass --rna-fastq + --reference-fasta "
                "to drive the default flow, or --de-table + --ribo-dir + "
                "--reference-gtf to bring your own pre-computed DE table.",
                file=sys.stderr,
            )
        return 2

    output_dir = Path(args.output)

    if from_fastq:
        try:
            short_circuit = _run_from_fastq(args, output_dir)
        except (FileNotFoundError, ValueError, NotImplementedError, RuntimeError) as exc:
            print(f"[mitoribopy rnaseq] ERROR: {exc}", file=sys.stderr)
            return 2
        if short_circuit == -1:
            return 0

    ribo_dir = Path(args.ribo_dir) if args.ribo_dir else None
    ribo_counts_path = (
        Path(args.ribo_counts)
        if args.ribo_counts
        else (ribo_dir / "rpf_counts.tsv" if ribo_dir else None)
    )

    # --- Reference-consistency gate (HARD FAIL) -------------------------
    if not getattr(args, "_from_fastq", False):
        if ribo_dir is None:
            print(
                "[mitoribopy rnaseq] ERROR: --ribo-dir is required for the "
                "reference-consistency gate; pass it alongside --ribo-counts if "
                "your counts file lives elsewhere.",
                file=sys.stderr,
            )
            return 2

        try:
            verified_checksum = verify_reference_consistency(
                ribo_dir=ribo_dir,
                reference_path=Path(args.reference_gtf)
                if args.reference_gtf
                else None,
                reference_checksum=args.reference_checksum,
            )
        except (ReferenceMismatchError, FileNotFoundError, ValueError) as exc:
            print(f"[mitoribopy rnaseq] ERROR: {exc}", file=sys.stderr)
            return 2
    else:
        verified_checksum = args._reference_checksum

    # --- DE table -------------------------------------------------------
    try:
        column_map = _resolve_de_column_map(args)
    except ValueError as exc:
        print(f"[mitoribopy rnaseq] ERROR: {exc}", file=sys.stderr)
        return 2

    try:
        de_table = load_de_table(Path(args.de_table), column_map=column_map)
    except (FileNotFoundError, ValueError) as exc:
        print(f"[mitoribopy rnaseq] ERROR: {exc}", file=sys.stderr)
        return 2

    # --- Gene ID match --------------------------------------------------
    de_gene_ids = [row["gene_id"] for row in de_table.rows]
    organism = "h" if args.organism.startswith("h") else "y"
    match = match_mt_mrnas(de_gene_ids, args.gene_id_convention, organism=organism)
    if match["missing"]:
        sys.stderr.write(
            "[mitoribopy rnaseq] WARNING: "
            f"{len(match['missing'])} expected mt-mRNA(s) not found in DE "
            f"table under convention '{args.gene_id_convention}': "
            + ", ".join(match["missing"])
            + "\n"
        )

    # --- Ribo counts ----------------------------------------------------
    try:
        ribo_counts = load_ribo_counts(ribo_counts_path)
    except (FileNotFoundError, ValueError) as exc:
        print(f"[mitoribopy rnaseq] ERROR: {exc}", file=sys.stderr)
        return 2

    # --- TE + delta-TE --------------------------------------------------
    mrna_abundances = {
        row["gene_id"]: row["basemean"]
        for row in de_table.rows
        if row["basemean"] is not None
    }
    condition_map = _load_condition_map(
        Path(args.condition_map) if args.condition_map else None
    )
    # When the from_fastq orchestrator built `de_table.tsv` itself we
    # want the TE/dTE rows to be tagged ``exploratory_mt_only``, NOT
    # ``publication_grade``. The orchestrator sets ``_from_fastq`` /
    # ``_pseudo_replicate_mode`` on args before falling through to
    # this de_table-style writer.
    is_from_fastq = bool(getattr(args, "_from_fastq", False))
    pseudo_rep = bool(getattr(args, "_pseudo_replicate_mode", False))
    te_mode = "from_fastq" if is_from_fastq else "de_table"
    te_rows = compute_te(
        ribo_counts,
        mrna_abundances,
        condition_map=condition_map or None,
        assay="ribo",
        mode=te_mode,
        pseudo_replicate=pseudo_rep,
    )

    dte_rows = compute_delta_te(
        ribo_counts,
        de_table,
        condition_map=condition_map or None,
        base_condition=args.condition_a,
        compare_condition=args.condition_b,
        mode=te_mode,
        pseudo_replicate=pseudo_rep,
    )

    # --- Write outputs --------------------------------------------------
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    _write_te_table(te_rows, output_dir / "te.tsv")
    _write_delta_te_table(dte_rows, output_dir / "delta_te.tsv")

    contrast_label = (
        f"{args.condition_a} vs {args.condition_b}"
        if args.condition_a and args.condition_b
        else None
    )

    plot_mrna_vs_rpf_scatter(
        dte_rows,
        output_dir / "plots" / "mrna_vs_rpf.png",
        title=(
            f"log2FC RPF vs mRNA — {contrast_label}"
            if contrast_label
            else "log2FC RPF vs mRNA"
        ),
    )
    plot_delta_te_volcano(
        dte_rows,
        output_dir / "plots" / "delta_te_volcano.png",
        title=(
            f"Delta-TE (log2) volcano — {contrast_label}"
            if contrast_label
            else "Delta-TE (log2) volcano"
        ),
    )
    plot_ma(
        de_table,
        output_dir / "plots" / "ma.png",
        title=(
            f"MA plot (mRNA DE) — {contrast_label}"
            if contrast_label
            else "MA plot (mRNA DE)"
        ),
    )
    plot_de_volcano(
        de_table,
        output_dir / "plots" / "de_volcano_mrna.png",
        contrast_label=contrast_label,
        title=(
            f"mRNA DE volcano — {contrast_label}"
            if contrast_label
            else "mRNA DE volcano"
        ),
    )

    rpf_de_table_path = getattr(args, "_rpf_de_table_path", None)
    if rpf_de_table_path:
        try:
            rpf_de_table = load_de_table(Path(rpf_de_table_path))
        except (FileNotFoundError, ValueError) as exc:
            sys.stderr.write(
                "[mitoribopy rnaseq] WARNING: Could not load Ribo-seq "
                f"DE table from {rpf_de_table_path!r} ({exc}); skipping "
                "RPF volcano plot.\n"
            )
        else:
            plot_de_volcano(
                rpf_de_table,
                output_dir / "plots" / "de_volcano_rpf.png",
                contrast_label=contrast_label,
                title=(
                    f"Ribo-seq DE volcano — {contrast_label}"
                    if contrast_label
                    else "Ribo-seq DE volcano"
                ),
            )

    plot_te_bar_grouped(
        te_rows, condition_map, output_dir / "plots" / "te_bar_by_condition.png"
    )
    plot_te_heatmap(
        te_rows, condition_map, output_dir / "plots" / "te_heatmap.png"
    )

    if condition_map and args.condition_a and args.condition_b:
        plot_te_compare_scatter(
            te_rows,
            condition_map,
            args.condition_a,
            args.condition_b,
            output_dir / "plots" / "te_compare_scatter.png",
        )
        plot_te_log2fc_bar(
            te_rows,
            condition_map,
            args.condition_a,
            args.condition_b,
            output_dir / "plots" / "te_log2fc_bar.png",
        )

    settings = {
        "subcommand": "rnaseq",
        "mitoribopy_version": __version__,
        "rnaseq_mode": getattr(args, "_rnaseq_mode", None),
        "de_table": str(args.de_table),
        "rpf_de_table": getattr(args, "_rpf_de_table_path", None),
        "de_format_resolved": de_table.format,
        "de_column_map": asdict(de_table.column_map),
        "gene_id_convention": args.gene_id_convention,
        "organism": organism,
        "ribo_dir": str(ribo_dir) if ribo_dir is not None else None,
        "ribo_counts": str(ribo_counts_path),
        "reference_checksum": verified_checksum,
        "mt_mrna_match": match,
        "condition_map_path": args.condition_map,
        "condition_a": args.condition_a,
        "condition_b": args.condition_b,
        "base_sample": args.condition_a,
        "compare_sample": args.condition_b,
        "te_row_count": len(te_rows),
        "delta_te_row_count": len(dte_rows),
    }
    fastq_provenance = getattr(args, "_fastq_provenance", None)
    if fastq_provenance is not None:
        settings["from_fastq"] = fastq_provenance
    pseudo_mode = bool(getattr(args, "_pseudo_replicate_mode", False))
    pseudo_conditions = list(getattr(args, "_pseudo_replicate_conditions", []) or [])
    settings["pseudo_replicate_mode"] = pseudo_mode
    settings["pseudo_replicate_conditions"] = pseudo_conditions
    settings["sample_sheet"] = getattr(args, "_sample_sheet_path", None)
    (output_dir / "run_settings.json").write_text(
        json.dumps(settings, indent=2, sort_keys=True), encoding="utf-8"
    )
    if pseudo_mode:
        _emit_pseudo_replicate_banner(pseudo_conditions, when="end")
    return 0
