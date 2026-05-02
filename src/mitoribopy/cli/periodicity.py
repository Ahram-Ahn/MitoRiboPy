"""``mitoribopy periodicity`` — standalone metagene Fourier QC.

Re-runs the Fourier-spectrum periodicity bundle (per-(sample, length,
gene_set, region) metagene amplitude curve + period-3 spectral ratio +
three-figure plot bundle) from a pre-assigned site table, without
invoking the upstream ``rpf`` pipeline. Use this when you want to re-
score periodicity from a saved site_table — e.g., to try a different
site assignment (P-site vs A-site) or a different window width.

Required input columns (TSV / CSV):
    sample, gene, transcript_id, read_length, site_type,
    site_pos, cds_start, cds_end

Optional columns:
    count   (defaults to 1 per row when missing)

The command does NOT recalculate offsets — it scores already-assigned
sites. For an end-to-end run, use ``mitoribopy rpf`` or ``mitoribopy
all`` instead.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd

from ..analysis.fourier_spectrum import (
    DEFAULT_DROP_CODONS_AFTER_START,
    DEFAULT_DROP_CODONS_BEFORE_STOP,
    DEFAULT_MIN_MEAN_COVERAGE,
    DEFAULT_MIN_TOTAL_COUNTS,
    DEFAULT_WINDOW_NT,
)
from ..analysis.periodicity_qc import run_periodicity_qc_bundle
from . import common


PERIODICITY_SUBCOMMAND_HELP = (
    "Quantify 3-nt periodicity by running the metagene Fourier "
    "analysis on a pre-assigned site table."
)


_REQUIRED_COLUMNS = (
    "sample",
    "gene",
    "transcript_id",
    "read_length",
    "site_type",
    "site_pos",
    "cds_start",
    "cds_end",
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mitoribopy periodicity",
        description=PERIODICITY_SUBCOMMAND_HELP,
        formatter_class=common.MitoRiboPyHelpFormatter,
        epilog=(
            "Frame definition: (site_pos - cds_start) mod 3. The Fourier "
            "analysis is anchored at the start codon (orf_start) and the "
            "stop codon (orf_stop). site_pos must be transcript-oriented "
            "(forward-strand-relative) and 0-based."
        ),
    )
    parser.add_argument(
        "--site-table",
        required=True,
        metavar="PATH",
        help=(
            "Per-read site table; required columns: "
            + ", ".join(_REQUIRED_COLUMNS)
            + ". `count` is optional (defaults to 1)."
        ),
    )
    parser.add_argument(
        "--output",
        required=True,
        metavar="DIR",
        help="Directory for the Fourier QC bundle.",
    )
    parser.add_argument(
        "--site",
        choices=("p", "a"),
        default="p",
        help="Ribosomal site to score. Default: p.",
    )
    parser.add_argument(
        "--fourier-window-nt",
        type=int,
        default=DEFAULT_WINDOW_NT,
        metavar="N",
        help=(
            "Window size (nt) per region. Default: 99 (33 codons). "
            "Must be a multiple of 3 for clean period-3 bin alignment."
        ),
    )
    parser.add_argument(
        "--drop-codons-after-start",
        type=int,
        default=DEFAULT_DROP_CODONS_AFTER_START,
        metavar="N",
        help=(
            "Codons after the AUG to skip in the orf_start window. "
            "Default: 5 (skip the initiation peak)."
        ),
    )
    parser.add_argument(
        "--drop-codons-before-stop",
        type=int,
        default=DEFAULT_DROP_CODONS_BEFORE_STOP,
        metavar="N",
        help=(
            "Codons before the stop codon to skip in the orf_stop "
            "window. Default: 1 (skip the termination peak)."
        ),
    )
    parser.add_argument(
        "--min-mean-coverage",
        type=float,
        default=DEFAULT_MIN_MEAN_COVERAGE,
        metavar="X",
        help=(
            "Skip per-gene windows whose mean coverage is below X. "
            "Default: 0.1."
        ),
    )
    parser.add_argument(
        "--min-total-counts",
        type=int,
        default=DEFAULT_MIN_TOTAL_COUNTS,
        metavar="N",
        help=(
            "Skip per-gene windows whose total site count is below N. "
            "Default: 30."
        ),
    )
    parser.add_argument(
        "--no-plots",
        dest="render_plots",
        action="store_false",
        default=True,
        help="Skip the per-(sample, read_length) figures; TSVs still written.",
    )
    return parser


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _load_site_table(path: str | Path) -> pd.DataFrame:
    """Load TSV or CSV based on suffix."""
    p = Path(path)
    sep = "," if p.suffix.lower() == ".csv" else "\t"
    return pd.read_csv(p, sep=sep, comment="#")


def _validate(df: pd.DataFrame) -> list[str]:
    missing = [c for c in _REQUIRED_COLUMNS if c not in df.columns]
    errs: list[str] = []
    if missing:
        errs.append(
            "site-table is missing required column(s): " + ", ".join(missing)
        )
        return errs
    for col in ("site_pos", "cds_start", "cds_end", "read_length"):
        if not pd.api.types.is_numeric_dtype(df[col]):
            try:
                pd.to_numeric(df[col])
            except Exception:
                errs.append(f"column {col!r} must be numeric")
    if "count" in df.columns and (df["count"] < 0).any():
        errs.append("column 'count' contains negative values")
    return errs


def _to_bed_with_psite(df: pd.DataFrame, *, site_filter: str) -> pd.DataFrame:
    """Project the site table into the BED-shaped layout the Fourier
    extractor consumes (sample_name, chrom, P_site, read_length).
    """
    sub = df[df["site_type"].astype(str).str.lower() == site_filter].copy()
    if sub.empty:
        return pd.DataFrame()
    sub["sample_name"] = sub["sample"].astype(str)
    sub["chrom"] = sub["transcript_id"].astype(str)
    sub["P_site"] = pd.to_numeric(sub["site_pos"]).astype(int)
    sub["read_length"] = pd.to_numeric(sub["read_length"]).astype(int)
    return sub[["sample_name", "chrom", "P_site", "read_length"]]


def _annotation_from_site_table(df: pd.DataFrame) -> pd.DataFrame:
    """Synthesise the annotation from the cds_start / cds_end columns
    when the user passes only a site table.

    Convention: ``cds_end`` is the inclusive last nt of the stop codon
    (so cds_end + 1 = transcript end when there is no 3' UTR). The first
    nt of the stop trinucleotide is therefore ``cds_end - 2``.
    """
    ann = (
        df.groupby(["transcript_id", "gene"], dropna=False)
        .agg(
            start_codon=("cds_start", "min"),
            cds_end_max=("cds_end", "max"),
        )
        .reset_index()
    )
    ann["start_codon"] = ann["start_codon"].astype(int)
    cds_end = ann["cds_end_max"].astype(int)
    ann["l_cds"] = (cds_end - ann["start_codon"]).clip(lower=0)
    ann["stop_codon"] = (cds_end - 2).clip(lower=0)
    ann["l_tr"] = cds_end + 1
    ann["transcript"] = ann["transcript_id"].astype(str)
    ann["sequence_name"] = ann["transcript_id"].astype(str)
    ann["sequence_aliases"] = ""
    return ann[[
        "transcript", "sequence_name", "sequence_aliases",
        "start_codon", "l_cds", "stop_codon", "l_tr",
    ]]


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def run(argv: Iterable[str]) -> int:
    args = build_parser().parse_args(list(argv))

    try:
        df = _load_site_table(args.site_table)
    except (FileNotFoundError, OSError, ValueError) as exc:
        print(f"[mitoribopy periodicity] ERROR: {exc}", file=sys.stderr)
        return 2

    errs = _validate(df)
    if errs:
        for line in errs:
            print(f"[mitoribopy periodicity] ERROR: {line}", file=sys.stderr)
        return 2

    site_filter = args.site.lower()
    bed = _to_bed_with_psite(df, site_filter=site_filter)
    annotation = _annotation_from_site_table(df)
    samples = sorted(bed["sample_name"].astype(str).unique()) if not bed.empty else []

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    result = run_periodicity_qc_bundle(
        bed_with_psite=bed if not bed.empty else None,
        annotation_df=annotation if not annotation.empty else None,
        samples=samples,
        output_dir=output_dir,
        site_type=site_filter,  # type: ignore[arg-type]
        window_nt=int(args.fourier_window_nt),
        drop_codons_after_start=int(args.drop_codons_after_start),
        drop_codons_before_stop=int(args.drop_codons_before_stop),
        min_mean_coverage=float(args.min_mean_coverage),
        min_total_counts=int(args.min_total_counts),
        render_plots=bool(args.render_plots),
    )

    print(
        "[mitoribopy periodicity] wrote: "
        f"{result['paths']['fourier_spectrum_combined']}, "
        f"{result['paths']['fourier_period3_score_combined']}",
        file=sys.stderr,
    )
    return 0
