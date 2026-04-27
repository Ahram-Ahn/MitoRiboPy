"""``mitoribopy rnaseq`` subcommand - DE + rpf -> TE / delta-TE + plots.

Phase 5 of the v0.3.0 refactor. Consumes a pre-computed differential
expression table (DESeq2 / Xtail / Anota2Seq) + a prior ``mitoribopy
rpf`` run and emits ``te.tsv``, ``delta_te.tsv``, and diagnostic plots.
A SHA256 reference-consistency gate rejects mismatched references
before any math runs.
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
from ..rnaseq.plots import plot_delta_te_volcano, plot_mrna_vs_rpf_scatter
from . import common


RNASEQ_SUBCOMMAND_HELP = (
    "Integrate a pre-computed differential-expression (DESeq2 / Xtail / "
    "Anota2Seq) table with a prior 'mitoribopy rpf' run and emit TE and "
    "delta-TE tables + plots. Enforces a SHA256 reference-consistency "
    "gate so Ribo-seq and RNA-seq must share the same transcript set."
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mitoribopy rnaseq",
        description=RNASEQ_SUBCOMMAND_HELP,
        formatter_class=common.MitoRiboPyHelpFormatter,
    )
    common.add_common_arguments(parser)

    de = parser.add_argument_group("DE table")
    de.add_argument(
        "--de-table",
        required=False,
        metavar="PATH",
        help="DESeq2 / Xtail / Anota2Seq results table (CSV or TSV).",
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

    gene = parser.add_argument_group("Gene identifiers")
    gene.add_argument(
        "--gene-id-convention",
        choices=["ensembl", "refseq", "hgnc", "mt_prefixed", "bare"],
        required=False,
        help=(
            "Gene identifier convention used by the DE table. REQUIRED "
            "(no default): mismatched conventions between DE and Ribo-seq "
            "sides silently produce zero-match runs."
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

    ribo = parser.add_argument_group("Ribo-seq inputs")
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

    ref = parser.add_argument_group("Reference-consistency gate")
    ref.add_argument(
        "--reference-gtf",
        default=None,
        metavar="PATH",
        help=(
            "Reference GTF / FASTA used by RNA-seq; we hash this and "
            "verify it matches the hash recorded in the rpf run's "
            "manifest. EXACTLY ONE of --reference-gtf / "
            "--reference-checksum must be provided."
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

    conditions = parser.add_argument_group("Conditions (optional)")
    conditions.add_argument(
        "--condition-map",
        default=None,
        metavar="PATH",
        help=(
            "TSV with columns 'sample' and 'condition' assigning Ribo-seq "
            "samples to conditions. Required to emit an internal Ribo-seq "
            "log2FC (used for delta-TE). Without this the delta-TE rows "
            "carry only the mRNA log2FC and a 'single_replicate_"
            "no_statistics' note."
        ),
    )
    conditions.add_argument("--condition-a", default=None, metavar="NAME")
    conditions.add_argument("--condition-b", default=None, metavar="NAME")

    out = parser.add_argument_group("Output")
    out.add_argument(
        "--output",
        required=False,
        metavar="DIR",
        help="Output directory for te.tsv, delta_te.tsv, and plots.",
    )

    fastq = parser.add_argument_group("From-FASTQ mode")
    fastq.add_argument(
        "--rna-fastq",
        nargs="+",
        default=None,
        metavar="PATH",
        help=(
            "RNA-seq FASTQ files or directories. Activates from-FASTQ "
            "mode: the subcommand will trim, align, count, and run "
            "pyDESeq2 itself. Mutually exclusive with --de-table."
        ),
    )
    fastq.add_argument(
        "--ribo-fastq",
        nargs="+",
        default=None,
        metavar="PATH",
        help=(
            "Ribo-seq FASTQ files or directories. Optional in from-FASTQ "
            "mode; when omitted the run short-circuits after writing "
            "de_table.tsv."
        ),
    )
    fastq.add_argument(
        "--reference-fasta",
        default=None,
        metavar="PATH",
        help="Transcriptome FASTA used to build (or re-use) a bowtie2 index.",
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
        help="Threads passed to cutadapt and bowtie2 in from-FASTQ mode.",
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


def _write_te_table(rows: Iterable[TeRow], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("sample\tgene\trpf_count\tmrna_abundance\tte\n")
        for row in rows:
            handle.write(
                f"{row.sample}\t{row.gene}\t{row.rpf_count}\t"
                f"{row.mrna_abundance:.6g}\t{row.te:.6g}\n"
            )


def _write_delta_te_table(rows: Iterable[DTeRow], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("gene\tmrna_log2fc\trpf_log2fc\tdelta_te_log2\tpadj\tnote\n")
        for row in rows:
            handle.write(
                f"{row.gene}\t"
                f"{'' if row.mrna_log2fc is None else f'{row.mrna_log2fc:.6g}'}\t"
                f"{'' if row.rpf_log2fc is None else f'{row.rpf_log2fc:.6g}'}\t"
                f"{'' if row.delta_te_log2 is None else f'{row.delta_te_log2:.6g}'}\t"
                f"{'' if row.padj is None else f'{row.padj:.6g}'}\t"
                f"{row.note}\n"
            )


# ---------- from-FASTQ orchestrator -----------------------------------------


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

    rna_paths = enumerate_fastqs(args.rna_fastq)
    rna_samples = detect_samples(rna_paths)
    if args.ribo_fastq:
        ribo_paths = enumerate_fastqs(args.ribo_fastq)
        ribo_samples = detect_samples(ribo_paths)
    else:
        ribo_samples = []

    if args.bowtie2_index:
        bt2_index = Path(args.bowtie2_index)
    else:
        bt2_index = build_bowtie2_index(
            Path(args.reference_fasta), workdir / "bt2_cache"
        )

    rna_results: list[SampleAlignmentResult] = []
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
    for s in ribo_samples:
        ribo_results.append(
            align_sample(
                s,
                bt2_index=bt2_index,
                workdir=workdir / "ribo",
                threads=args.align_threads,
            )
        )

    if rna_results:
        write_counts_matrix(rna_results, output_dir / "rna_counts.tsv")
    if ribo_results:
        write_counts_matrix(ribo_results, output_dir / "rpf_counts_matrix.tsv")
        write_long_counts(ribo_results, output_dir / "rpf_counts.tsv")

    condition_map = _load_condition_map(
        Path(args.condition_map) if args.condition_map else None
    )

    samples, counts_df, metadata_df = build_sample_sheet(
        rna_results, ribo_results, condition_map
    )
    results_df = run_deseq2(
        counts_df,
        metadata_df,
        contrast_factor="condition",
        contrast_a=args.condition_a,
        contrast_b=args.condition_b,
    )
    de_table = deseq2_to_de_table(results_df)
    de_table_path = output_dir / "de_table.tsv"
    write_de_table_tsv(de_table, de_table_path)

    reference_checksum = compute_reference_checksum(Path(args.reference_fasta))

    args.de_table = str(de_table_path)
    args.de_format = "deseq2"
    args.ribo_counts = (
        str(output_dir / "rpf_counts.tsv") if ribo_results else None
    )
    args._from_fastq = True
    args._reference_checksum = reference_checksum
    args._fastq_provenance = {
        "mode": "from-fastq" if ribo_results else "from-fastq-rna-only",
        "reference_fasta": str(args.reference_fasta),
        "bowtie2_index": str(bt2_index),
        "rna_samples": [_provenance_entry(r) for r in rna_results],
        "ribo_samples": [_provenance_entry(r) for r in ribo_results],
    }

    if not ribo_results:
        # RNA-only short-circuit: write a slim run_settings.json and stop.
        slim_settings = {
            "subcommand": "rnaseq",
            "mitoribopy_version": __version__,
            "from_fastq": args._fastq_provenance,
            "reference_checksum": reference_checksum,
            "de_table": args.de_table,
            "de_format_resolved": "deseq2",
        }
        (output_dir / "run_settings.json").write_text(
            json.dumps(slim_settings, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        return -1

    return 0


# ---------- orchestrator ----------------------------------------------------


def run(argv: Iterable[str]) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv))
    common.apply_common_arguments(args)

    # --- Resolve required flags up front --------------------------------
    if args.dry_run:
        return common.emit_dry_run(
            "rnaseq",
            [
                "resolve --de-format (auto-detect from headers) and load DE table",
                "hash-verify --reference-gtf / --reference-checksum against "
                "the rpf run's recorded reference_checksum (HARD FAIL on mismatch)",
                "load rpf_counts.tsv from --ribo-dir or --ribo-counts",
                f"match DE gene_ids against mt-mRNA registry ({args.organism}) "
                f"using --gene-id-convention {args.gene_id_convention}",
                "compute TE per sample per gene = (RPF + 0.5) / (mRNA + 0.5)",
                "compute delta-TE = log2(RPF_fc) - log2(mRNA_fc) with "
                "replicate-based Ribo log2FC when --condition-map is given",
                "emit te.tsv, delta_te.tsv, and plots (scatter + volcano)",
            ],
        )

    from_fastq = args.rna_fastq is not None
    if from_fastq and args.de_table is not None:
        print(
            "[mitoribopy rnaseq] ERROR: --de-table cannot be combined with "
            "--rna-fastq; from-FASTQ mode generates the DE table itself.",
            file=sys.stderr,
        )
        return 2

    missing: list[str] = []
    if from_fastq:
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
        if not args.condition_a:
            missing.append("--condition-a")
        if not args.condition_b:
            missing.append("--condition-b")
    else:
        if not args.gene_id_convention:
            missing.append("--gene-id-convention")
        if not args.de_table:
            missing.append("--de-table")
        if not args.ribo_dir and not args.ribo_counts:
            missing.append("--ribo-dir or --ribo-counts")
        if not args.output:
            missing.append("--output")
        if not (args.reference_gtf or args.reference_checksum):
            missing.append("--reference-gtf or --reference-checksum")
    if missing:
        print(
            "[mitoribopy rnaseq] ERROR: missing required argument(s): "
            + ", ".join(missing),
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
    te_rows = compute_te(ribo_counts, mrna_abundances)

    condition_map = _load_condition_map(
        Path(args.condition_map) if args.condition_map else None
    )
    dte_rows = compute_delta_te(
        ribo_counts,
        de_table,
        condition_map=condition_map or None,
        condition_a=args.condition_a,
        condition_b=args.condition_b,
    )

    # --- Write outputs --------------------------------------------------
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    _write_te_table(te_rows, output_dir / "te.tsv")
    _write_delta_te_table(dte_rows, output_dir / "delta_te.tsv")
    plot_mrna_vs_rpf_scatter(dte_rows, output_dir / "plots" / "mrna_vs_rpf.png")
    plot_delta_te_volcano(dte_rows, output_dir / "plots" / "delta_te_volcano.png")

    settings = {
        "subcommand": "rnaseq",
        "mitoribopy_version": __version__,
        "de_table": str(args.de_table),
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
        "te_row_count": len(te_rows),
        "delta_te_row_count": len(dte_rows),
    }
    fastq_provenance = getattr(args, "_fastq_provenance", None)
    if fastq_provenance is not None:
        settings["from_fastq"] = fastq_provenance
    (output_dir / "run_settings.json").write_text(
        json.dumps(settings, indent=2, sort_keys=True), encoding="utf-8"
    )
    return 0
