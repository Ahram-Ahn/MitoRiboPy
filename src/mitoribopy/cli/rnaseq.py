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
        choices=["h", "y", "human", "yeast"],
        default="h",
        help="Organism for the mt-mRNA registry (default human).",
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

    missing = []
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

    ribo_dir = Path(args.ribo_dir) if args.ribo_dir else None
    ribo_counts_path = (
        Path(args.ribo_counts)
        if args.ribo_counts
        else (ribo_dir / "rpf_counts.tsv" if ribo_dir else None)
    )

    # --- Reference-consistency gate (HARD FAIL) -------------------------
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
        "ribo_dir": str(ribo_dir),
        "ribo_counts": str(ribo_counts_path),
        "reference_checksum": verified_checksum,
        "mt_mrna_match": match,
        "condition_map_path": args.condition_map,
        "condition_a": args.condition_a,
        "condition_b": args.condition_b,
        "te_row_count": len(te_rows),
        "delta_te_row_count": len(dte_rows),
    }
    (output_dir / "run_settings.json").write_text(
        json.dumps(settings, indent=2, sort_keys=True), encoding="utf-8"
    )
    return 0
