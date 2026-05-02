"""``mitoribopy periodicity`` — standalone 3-nt periodicity QC.

Re-runs the publication-grade periodicity bundle (frame_counts_by_sample_length,
qc_summary, gene_periodicity, qc_summary.md) from a pre-assigned site
table, without invoking the upstream ``rpf`` pipeline. Use this when you
want to re-score periodicity with different thresholds, or to score an
externally produced site table that the rpf stage already saved.

Required input columns (TSV / CSV):
    sample, gene, transcript_id, read_length, site_type,
    site_pos, cds_start, cds_end

Optional columns:
    count   (defaults to 1 per row when missing)

The command does NOT recalculate offsets — it validates already-assigned
sites. For an end-to-end run, use ``mitoribopy rpf`` or
``mitoribopy all`` instead.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from ..analysis.periodicity_qc import (
    QC_THRESHOLDS_DEFAULT,
    build_frame_counts_by_gene,
    build_frame_counts_by_sample_length,
    build_gene_periodicity,
    build_qc_summary,
    calculate_entropy_bias,
    calculate_frame_enrichment,
    write_qc_summary_markdown,
)
from . import common


PERIODICITY_SUBCOMMAND_HELP = (
    "Quantify and summarise 3-nt periodicity from an already-assigned "
    "P-site / A-site coordinate table (no offset recomputation)."
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
            "Frame definition: (site_pos - cds_start) mod 3. Frame 0 is "
            "the annotated coding frame. site_pos must be transcript-"
            "oriented (forward-strand-relative) and 0-based."
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
        help="Directory for the periodicity QC bundle.",
    )
    parser.add_argument(
        "--site",
        choices=("p", "a"),
        default="p",
        help="Ribosomal site to score. Default: p.",
    )
    parser.add_argument(
        "--expected-frame",
        type=int,
        default=0,
        help="Expected CDS frame after site assignment. Default: 0.",
    )
    parser.add_argument(
        "--min-reads-per-length",
        type=int,
        default=int(QC_THRESHOLDS_DEFAULT["min_reads_per_length"]),
        help="Minimum CDS sites required to score a read length.",
    )
    parser.add_argument(
        "--min-reads-per-gene",
        type=int,
        default=int(QC_THRESHOLDS_DEFAULT["min_reads_per_gene"]),
        help="Minimum CDS sites required to score a gene.",
    )
    parser.add_argument(
        "--good-frame-fraction",
        type=float,
        default=float(QC_THRESHOLDS_DEFAULT["good_frame_fraction"]),
        help="Lower bound on expected_frame_fraction for qc_call=good.",
    )
    parser.add_argument(
        "--warn-frame-fraction",
        type=float,
        default=float(QC_THRESHOLDS_DEFAULT["warn_frame_fraction"]),
        help="Lower bound on expected_frame_fraction for qc_call=warn.",
    )
    parser.add_argument(
        "--exclude-start-codons",
        type=int,
        default=6,
        help="Codons after CDS start to exclude (initiation pause). Default: 6.",
    )
    parser.add_argument(
        "--exclude-stop-codons",
        type=int,
        default=3,
        help="Codons before CDS stop to exclude (termination pause). Default: 3.",
    )
    parser.add_argument(
        "--include-overlaps",
        action="store_true",
        default=False,
        help=(
            "Keep rows with is_overlap=true; default masks them when the "
            "column is present."
        ),
    )
    parser.add_argument(
        "--phase-score",
        action="store_true",
        default=False,
        help="Add a ribotricer-style gene-level phase_score column.",
    )
    parser.add_argument(
        "--fourier-spectrum",
        dest="fourier_spectrum",
        action="store_true",
        default=True,
        help=(
            "Compute the Wakigawa-style amplitude spectrum per "
            "(sample, read_length, gene, region) and emit "
            "fourier_spectrum.tsv, fourier_period3_score.tsv, and "
            "fourier_period3_summary.tsv. ENABLED BY DEFAULT — pass "
            "--no-fourier-spectrum to skip. Replaces the legacy "
            "fft_period3_power.tsv output (single-scalar power ratio)."
        ),
    )
    parser.add_argument(
        "--no-fourier-spectrum",
        dest="fourier_spectrum",
        action="store_false",
        help="Disable the Fourier-spectrum bundle.",
    )
    parser.add_argument(
        "--fourier-window-nt",
        type=int,
        default=100,
        metavar="N",
        help=(
            "Window size (nt) on each side of the canonical stop "
            "codon for the Fourier spectrum. Wakigawa default: 100."
        ),
    )
    parser.add_argument(
        "--fourier-period-min",
        type=float,
        default=2.0,
        metavar="P",
        help=(
            "Minimum period (nt) retained in the spectrum table and "
            "plot. Default: 2."
        ),
    )
    parser.add_argument(
        "--fourier-period-max",
        type=float,
        default=10.0,
        metavar="P",
        help=(
            "Maximum period (nt) retained in the spectrum table and "
            "plot. Default: 10."
        ),
    )
    parser.add_argument(
        "--fourier-no-plots",
        dest="fourier_plots",
        action="store_false",
        default=True,
        help=(
            "Skip the per-(sample, read_length) two-panel overlay "
            "plots; the TSVs are still written."
        ),
    )
    parser.add_argument(
        "--metagene-nt",
        type=int,
        default=90,
        help=(
            "Number of nucleotides upstream/downstream of start/stop "
            "codons for metagene plots. Default: 90 (per spec)."
        ),
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


def _to_frame_by_length(
    df: pd.DataFrame,
    *,
    site_filter: str,
    expected_frame: int,
    exclude_start_codons: int,
    exclude_stop_codons: int,
    include_overlaps: bool,
) -> pd.DataFrame:
    """Aggregate the site table into the per-(sample, read_length) layout
    that :func:`build_frame_counts_by_sample_length` expects.
    """
    sub = df[df["site_type"].astype(str).str.lower() == site_filter].copy()
    if "count" not in sub.columns:
        sub["count"] = 1
    if "is_overlap" in sub.columns and not include_overlaps:
        sub = sub[~sub["is_overlap"].astype(bool)]
    for col in ("site_pos", "cds_start", "cds_end", "read_length", "count"):
        sub[col] = pd.to_numeric(sub[col])
    sub["count"] = sub["count"].astype(int)

    in_cds = (
        (sub["site_pos"] >= sub["cds_start"] + 3 * exclude_start_codons)
        & (sub["site_pos"] < sub["cds_end"] - 3 * exclude_stop_codons)
    )
    cds = sub.loc[in_cds].copy()
    cds["frame"] = ((cds["site_pos"] - cds["cds_start"]) % 3).astype(int)

    rows: list[dict] = []
    if cds.empty:
        return pd.DataFrame(
            columns=[
                "sample_id", "read_length", "n_reads_total", "n_reads_cds",
                "frame0_fraction", "frame1_fraction", "frame2_fraction",
                "dominant_frame", "frame0_dominance", "periodicity_score",
                "frame_entropy", "include_for_downstream", "exclusion_reason",
            ]
        )

    n_total_by = (
        sub.groupby(["sample", "read_length"])["count"].sum().to_dict()
    )
    grouped = cds.groupby(["sample", "read_length", "frame"], observed=True)
    counts = grouped["count"].sum().reset_index()
    pivot = (
        counts.pivot_table(
            index=["sample", "read_length"],
            columns="frame",
            values="count",
            fill_value=0,
        )
        .reindex(columns=[0, 1, 2], fill_value=0)
        .reset_index()
    )

    for _, r in pivot.iterrows():
        sample = str(r["sample"])
        read_length = int(r["read_length"])
        n0, n1, n2 = int(r[0]), int(r[1]), int(r[2])
        n_cds = n0 + n1 + n2
        n_total = int(n_total_by.get((sample, read_length), n_cds))
        if n_cds == 0:
            f0 = f1 = f2 = 0.0
        else:
            f0, f1, f2 = n0 / n_cds, n1 / n_cds, n2 / n_cds
        fractions = [f0, f1, f2]
        dominant_frame = int(np.argmax(fractions))
        dominance = fractions[expected_frame] - max(
            f for i, f in enumerate(fractions) if i != expected_frame
        )
        rows.append({
            "sample_id": sample,
            "read_length": read_length,
            "n_reads_total": n_total,
            "n_reads_cds": n_cds,
            "frame0_fraction": f0,
            "frame1_fraction": f1,
            "frame2_fraction": f2,
            "dominant_frame": dominant_frame,
            "frame0_dominance": float(dominance),
            "periodicity_score": float(dominance),
            "frame_entropy": calculate_entropy_bias((f0, f1, f2)),
            "include_for_downstream": True,
            "exclusion_reason": "none",
        })
    return pd.DataFrame(rows)


def _to_bed_with_psite(df: pd.DataFrame, *, site_filter: str) -> pd.DataFrame:
    """Project the site table into the layout that
    :func:`build_gene_periodicity` consumes (sample_name, chrom, P_site,
    read_length).
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
    """Synthesise the annotation table from the cds_start / cds_end
    columns when the user passes only the site table.

    The Fourier-spectrum module also needs ``stop_codon`` (first nt of
    the stop codon) and ``l_tr`` (transcript length); both are derived
    from cds_end here. We assume the site table's ``cds_end`` is the
    inclusive last nt of the stop codon (cds_end + 1 = transcript end
    when no 3' UTR is annotated). For a synthetic site table that
    carries no 3' UTR, ``stop_codon = cds_end - 2`` (the first of the
    stop trinucleotide) and ``l_tr = cds_end + 1`` so the Fourier
    window check skips genes whose 3' UTR doesn't fit.
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
    # Stop codon = first nt of the 3-nt stop trinucleotide. For a CDS
    # whose end inclusively names the last stop nt, the first stop nt
    # is two positions earlier.
    ann["stop_codon"] = (cds_end - 2).clip(lower=0)
    ann["l_tr"] = cds_end + 1
    ann["transcript"] = ann["transcript_id"].astype(str)
    ann["sequence_name"] = ann["transcript_id"].astype(str)
    return ann[[
        "transcript", "sequence_name",
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
    thresholds = {
        "good_frame_fraction": float(args.good_frame_fraction),
        "warn_frame_fraction": float(args.warn_frame_fraction),
        "min_reads_per_length": int(args.min_reads_per_length),
        "min_reads_per_gene": int(args.min_reads_per_gene),
    }

    frame_by_length = _to_frame_by_length(
        df,
        site_filter=site_filter,
        expected_frame=int(args.expected_frame),
        exclude_start_codons=int(args.exclude_start_codons),
        exclude_stop_codons=int(args.exclude_stop_codons),
        include_overlaps=bool(args.include_overlaps),
    )

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    by_length = build_frame_counts_by_sample_length(
        frame_by_length,
        expected_frame=int(args.expected_frame),
        thresholds=thresholds,
    )
    by_length_path = output_dir / "frame_counts_by_sample_length.tsv"
    by_length.to_csv(by_length_path, sep="\t", index=False, na_rep="")

    qc_summary = build_qc_summary(
        by_length, site_type=site_filter, thresholds=thresholds,
    )
    qc_summary_path = output_dir / "qc_summary.tsv"
    qc_summary.to_csv(qc_summary_path, sep="\t", index=False, na_rep="")
    qc_summary_md = output_dir / "qc_summary.md"
    write_qc_summary_markdown(qc_summary, qc_summary_md, site_type=site_filter)

    gene_path = output_dir / "gene_periodicity.tsv"
    frame_counts_by_gene_path = output_dir / "frame_counts_by_gene.tsv"
    fourier_spectrum_path = output_dir / "fourier_spectrum.tsv"
    fourier_score_path = output_dir / "fourier_period3_score.tsv"
    fourier_summary_path = output_dir / "fourier_period3_summary.tsv"
    bed = _to_bed_with_psite(df, site_filter=site_filter)
    annotation = _annotation_from_site_table(df)
    written_extras: list[Path] = []
    if not bed.empty:
        samples = sorted(bed["sample_name"].astype(str).unique())
        gene_table = build_gene_periodicity(
            bed,
            annotation,
            samples=samples,
            expected_frame=int(args.expected_frame),
            thresholds=thresholds,
            compute_phase_score=bool(args.phase_score),
            exclude_start_codons=int(args.exclude_start_codons),
            exclude_stop_codons=int(args.exclude_stop_codons),
            annotate_overlap=True,
        )
        gene_table.to_csv(gene_path, sep="\t", index=False, na_rep="")
        # Spec output: per-(sample, gene, read_length) frame counts.
        fc_gene = build_frame_counts_by_gene(
            bed,
            annotation,
            samples=samples,
            expected_frame=int(args.expected_frame),
            thresholds=thresholds,
            exclude_start_codons=int(args.exclude_start_codons),
            exclude_stop_codons=int(args.exclude_stop_codons),
            annotate_overlap=True,
            site_type=site_filter,
        )
        fc_gene.to_csv(frame_counts_by_gene_path, sep="\t", index=False, na_rep="")
        written_extras.append(frame_counts_by_gene_path)
        # Wakigawa-style Fourier amplitude spectrum + derived scalars.
        # ENABLED BY DEFAULT (--no-fourier-spectrum to disable).
        if bool(getattr(args, "fourier_spectrum", True)):
            from ..analysis.fourier_spectrum import (
                build_fourier_spectrum_table,
                build_period3_score_table,
                build_period3_summary_table,
            )

            spectrum_table = build_fourier_spectrum_table(
                bed,
                annotation,
                samples=samples,
                window_nt=int(getattr(args, "fourier_window_nt", 100)),
                period_range=(
                    float(getattr(args, "fourier_period_min", 2.0)),
                    float(getattr(args, "fourier_period_max", 10.0)),
                ),
                site=site_filter,  # type: ignore[arg-type]
            )
            spectrum_table.to_csv(
                fourier_spectrum_path, sep="\t", index=False, na_rep="",
            )
            written_extras.append(fourier_spectrum_path)
            score_table = build_period3_score_table(spectrum_table)
            score_table.to_csv(
                fourier_score_path, sep="\t", index=False, na_rep="",
            )
            written_extras.append(fourier_score_path)
            summary_table = build_period3_summary_table(score_table)
            summary_table.to_csv(
                fourier_summary_path, sep="\t", index=False, na_rep="",
            )
            written_extras.append(fourier_summary_path)
            if (
                bool(getattr(args, "fourier_plots", True))
                and not spectrum_table.empty
            ):
                from ..plotting.fourier_spectrum_plots import (
                    render_fourier_spectrum_panels,
                )

                render_fourier_spectrum_panels(
                    spectrum_table,
                    output_dir=output_dir / "fourier_spectrum",
                )

    metadata = {
        "site": site_filter,
        "expected_frame": int(args.expected_frame),
        "frame_formula": "(site_pos - cds_start) % 3",
        "thresholds": thresholds,
        "exclude_start_codons": int(args.exclude_start_codons),
        "exclude_stop_codons": int(args.exclude_stop_codons),
        "include_overlaps": bool(args.include_overlaps),
        "phase_score": bool(args.phase_score),
        "fourier_spectrum": bool(getattr(args, "fourier_spectrum", True)),
        "fourier_window_nt": int(getattr(args, "fourier_window_nt", 100)),
        "fourier_period_range": [
            float(getattr(args, "fourier_period_min", 2.0)),
            float(getattr(args, "fourier_period_max", 10.0)),
        ],
        "fourier_plots": bool(getattr(args, "fourier_plots", True)),
        "metagene_nt": int(args.metagene_nt),
        "input_site_table": str(Path(args.site_table).resolve()),
    }
    (output_dir / "periodicity.metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    extras_str = (", " + ", ".join(str(p) for p in written_extras)) if written_extras else ""
    print(
        "[mitoribopy periodicity] wrote: "
        f"{by_length_path}, {qc_summary_path}, {qc_summary_md}, "
        f"{gene_path}{extras_str}",
        file=sys.stderr,
    )
    return 0
