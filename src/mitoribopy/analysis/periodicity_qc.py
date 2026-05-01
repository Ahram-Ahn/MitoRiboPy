"""Publication-grade 3-nt periodicity QC bundle (v0.6.2).

This module is layered on top of :mod:`mitoribopy.analysis.periodicity`
and provides the additional artefacts requested in the periodicity
implementation specification:

* ``frame_counts_by_sample_length.tsv`` — frame fractions, expected /
  dominant frame, frame enrichment, entropy bias, and a soft QC call
  per (sample, read length).
* ``gene_periodicity.tsv`` — same metrics resolved per gene with an
  optional gene-level ribotricer-style phase score.
* ``qc_summary.tsv`` / ``qc_summary.md`` — one row / one short markdown
  block per sample with the overall QC verdict.

The module is intentionally pure (DataFrame in, DataFrame out) so the
CLI and the test suite can drive it without disk I/O. The
:func:`run_periodicity_qc_bundle` orchestrator is the only function
that touches the filesystem.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


__all__ = [
    "GeneRow",
    "QCSummaryRow",
    "QC_THRESHOLDS_DEFAULT",
    "build_frame_counts_by_sample_length",
    "build_gene_periodicity",
    "build_qc_summary",
    "calculate_frame_enrichment",
    "calculate_entropy_bias",
    "calculate_phase_score",
    "calculate_fft_period3_ratio",
    "run_periodicity_qc_bundle",
    "write_qc_summary_markdown",
]


# Default thresholds match the values requested in the periodicity
# implementation spec. They map to the QC labels:
#
#   good     >= good_frame_fraction
#   warn     >= warn_frame_fraction
#   poor     <  warn_frame_fraction (and not low_depth)
#   low_depth   total assigned sites < min_reads_per_length
#
QC_THRESHOLDS_DEFAULT: dict[str, float] = {
    "good_frame_fraction": 0.60,
    "warn_frame_fraction": 0.50,
    "min_reads_per_length": 1000,
    "min_reads_per_gene": 50,
    "min_phase_score_good": 0.50,
}


_EPS = 1e-12


@dataclass(frozen=True)
class GeneRow:
    """One row of ``gene_periodicity.tsv``."""

    sample_id: str
    gene: str
    transcript_id: str
    read_length: int  # -1 means "all lengths pooled"
    n_frame0: int
    n_frame1: int
    n_frame2: int
    n_total: int
    frame0_fraction: float
    frame1_fraction: float
    frame2_fraction: float
    expected_frame_fraction: float
    expected_frame_enrichment: float
    dominant_frame: int
    dominant_frame_fraction: float
    entropy_bias: float
    phase_score: float
    qc_call: str


@dataclass(frozen=True)
class QCSummaryRow:
    """One row of ``qc_summary.tsv``."""

    sample_id: str
    site_type: str
    n_total_sites: int
    best_read_length: int
    best_read_length_expected_frame_fraction: float
    best_read_length_dominant_frame: int
    best_read_length_dominant_fraction: float
    global_expected_frame_fraction: float
    global_entropy_bias: float
    n_good_lengths: int
    n_warn_lengths: int
    n_poor_lengths: int
    n_low_depth_lengths: int
    overall_qc_call: str
    notes: str


# ---------------------------------------------------------------------------
# Pure metric helpers
# ---------------------------------------------------------------------------


def calculate_entropy_bias(fractions: tuple[float, float, float]) -> float:
    """Return ``1 - H_3(p)`` where H_3 is the base-3 entropy of ``fractions``.

    A value of 0 means uniform across frames; 1 means all reads in one
    frame. Values outside [0, 1] are clipped to that range to stay
    interpretable across small floating-point drift.
    """
    log3 = np.log(3.0)
    h = 0.0
    for p in fractions:
        if p > 0:
            h -= p * (np.log(p) / log3)
    return float(max(0.0, min(1.0, 1.0 - h)))


def calculate_frame_enrichment(
    fractions: tuple[float, float, float], expected_frame: int = 0
) -> float:
    """Return frame enrichment relative to the mean of the other two frames.

    For ``expected_frame=0`` this is ``f0 / mean(f1, f2)``; generalised
    for other expected frames. Returns ``inf`` when both other frames
    are zero (every read in the expected frame).
    """
    f = list(fractions)
    expected = f.pop(expected_frame)
    other_mean = (sum(f) / len(f)) if f else 0.0
    if other_mean <= _EPS:
        return float("inf") if expected > 0 else 0.0
    return float(expected / other_mean)


_PHASE_VECTORS: dict[int, np.ndarray] = {
    0: np.array([1.0, 0.0]),
    1: np.array([-0.5, np.sqrt(3.0) / 2]),
    2: np.array([-0.5, -np.sqrt(3.0) / 2]),
}


def calculate_phase_score(
    frame_counts_per_codon: Iterable[tuple[int, int, int]],
) -> float:
    """Return the ribotricer-style consistency of frame dominance.

    Each codon contributes a unit vector pointing toward its dominant
    frame; the score is the magnitude of the average over all codons.
    1.0 means every codon dominates the same frame; 0.0 means no
    consistent direction.
    """
    vectors: list[np.ndarray] = []
    for counts in frame_counts_per_codon:
        n0, n1, n2 = float(counts[0]), float(counts[1]), float(counts[2])
        total = n0 + n1 + n2
        if total <= 0:
            continue
        v = (
            n0 * _PHASE_VECTORS[0]
            + n1 * _PHASE_VECTORS[1]
            + n2 * _PHASE_VECTORS[2]
        )
        norm = np.linalg.norm(v)
        if norm > 0:
            vectors.append(v / norm)
    if not vectors:
        return float("nan")
    return float(np.linalg.norm(np.mean(np.stack(vectors), axis=0)))


def calculate_fft_period3_ratio(coverage) -> float:
    """Return the period-3 FFT power ratio of a 1-D coverage array.

    Returns NaN for arrays with fewer than 30 nt or zero variance.
    Compares the squared magnitude at frequency 1/3 against the median
    of the rest of the spectrum (excluding DC and the period-3 bin).
    """
    x = np.asarray(coverage, dtype=float)
    if x.size < 30:
        return float("nan")
    x = x - np.nanmean(x)
    if not np.isfinite(x).any() or float(np.nanvar(x)) <= _EPS:
        return float("nan")
    power = np.abs(np.fft.rfft(np.nan_to_num(x))) ** 2
    freqs = np.fft.rfftfreq(x.size, d=1)
    if freqs.size < 3:
        return float("nan")
    k3 = int(np.argmin(np.abs(freqs - 1.0 / 3.0)))
    background = np.delete(power, [0, k3])
    if background.size == 0:
        return float("nan")
    bg = float(np.nanmedian(background))
    return float(power[k3] / max(bg, _EPS))


# ---------------------------------------------------------------------------
# Frame-by-(sample, length) and gene-level tables
# ---------------------------------------------------------------------------


def _qc_call_from_fraction(
    n_total: int,
    expected_fraction: float,
    *,
    min_reads: int,
    good_threshold: float,
    warn_threshold: float,
) -> str:
    if n_total < min_reads:
        return "low_depth"
    if expected_fraction >= good_threshold:
        return "good"
    if expected_fraction >= warn_threshold:
        return "warn"
    return "poor"


def build_frame_counts_by_sample_length(
    frame_by_length: pd.DataFrame,
    *,
    expected_frame: int = 0,
    thresholds: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Promote a per-(sample, length) frame table to the publication schema.

    Adds the ``expected_frame_*`` columns, ``entropy_bias``, and a
    soft ``qc_call`` derived from the configured thresholds.
    """
    th = {**QC_THRESHOLDS_DEFAULT, **(thresholds or {})}
    if frame_by_length.empty:
        return frame_by_length.assign(
            expected_frame=expected_frame,
            expected_frame_fraction=0.0,
            expected_frame_enrichment=0.0,
            entropy_bias=0.0,
            qc_call="low_depth",
        )

    df = frame_by_length.copy()
    fractions = df[["frame0_fraction", "frame1_fraction", "frame2_fraction"]]
    expected_col = f"frame{expected_frame}_fraction"
    df["expected_frame"] = expected_frame
    df["expected_frame_fraction"] = fractions[expected_col].astype(float)
    df["expected_frame_enrichment"] = fractions.apply(
        lambda row: calculate_frame_enrichment(
            (
                float(row["frame0_fraction"]),
                float(row["frame1_fraction"]),
                float(row["frame2_fraction"]),
            ),
            expected_frame=expected_frame,
        ),
        axis=1,
    )
    df["entropy_bias"] = fractions.apply(
        lambda row: calculate_entropy_bias(
            (
                float(row["frame0_fraction"]),
                float(row["frame1_fraction"]),
                float(row["frame2_fraction"]),
            )
        ),
        axis=1,
    )
    df["qc_call"] = df.apply(
        lambda row: _qc_call_from_fraction(
            int(row["n_reads_cds"]),
            float(row["expected_frame_fraction"]),
            min_reads=int(th["min_reads_per_length"]),
            good_threshold=float(th["good_frame_fraction"]),
            warn_threshold=float(th["warn_frame_fraction"]),
        ),
        axis=1,
    )
    return df


def build_gene_periodicity(
    bed_with_psite_and_gene: pd.DataFrame,
    annotation_df: pd.DataFrame,
    *,
    samples: Iterable[str],
    expected_frame: int = 0,
    thresholds: dict[str, float] | None = None,
    compute_phase_score: bool = False,
) -> pd.DataFrame:
    """Per-(sample, gene) frame fractions + soft QC call.

    Expects the BED table from
    :func:`mitoribopy.analysis.periodicity.compute_p_site_positions` —
    rows must already carry ``P_site``, ``sample_name``, ``read_length``,
    and ``chrom``.
    """
    th = {**QC_THRESHOLDS_DEFAULT, **(thresholds or {})}
    if bed_with_psite_and_gene.empty or "P_site" not in bed_with_psite_and_gene.columns:
        return pd.DataFrame(
            columns=[
                "sample_id", "gene", "transcript_id", "read_length",
                "n_frame0", "n_frame1", "n_frame2", "n_total",
                "frame0_fraction", "frame1_fraction", "frame2_fraction",
                "expected_frame", "expected_frame_fraction",
                "expected_frame_enrichment", "dominant_frame",
                "dominant_frame_fraction", "entropy_bias",
                "phase_score", "qc_call",
            ],
        )

    ann = annotation_df.set_index("transcript")[["start_codon", "l_cds"]]
    chrom_to_tx = annotation_df.set_index("sequence_name")["transcript"].to_dict()

    df = bed_with_psite_and_gene.copy()
    df["transcript"] = df["chrom"].map(chrom_to_tx).fillna(df["chrom"])
    df = df[df["transcript"].isin(ann.index)]
    if df.empty:
        return pd.DataFrame()

    starts = df["transcript"].map(ann["start_codon"]).astype(int)
    cds_lens = df["transcript"].map(ann["l_cds"]).astype(int)
    in_cds_mask = (df["P_site"] >= starts) & (df["P_site"] < starts + cds_lens)
    df = df.loc[in_cds_mask].copy()
    if df.empty:
        return pd.DataFrame()
    df["frame"] = (
        (df["P_site"].astype(int) - starts.loc[df.index]) % 3
    ).astype(int)
    df["codon_index"] = (
        (df["P_site"].astype(int) - starts.loc[df.index]) // 3
    ).astype(int)

    rows: list[dict] = []
    samples_set = set(samples)
    df = df[df["sample_name"].astype(str).isin(samples_set)]
    if df.empty:
        return pd.DataFrame()

    for (sample, transcript), group in df.groupby(["sample_name", "transcript"]):
        counts = group["frame"].value_counts().reindex([0, 1, 2], fill_value=0)
        n_total = int(counts.sum())
        if n_total <= 0:
            continue
        f0, f1, f2 = (float(counts[i] / n_total) for i in (0, 1, 2))
        dominant_frame = int(counts.idxmax())
        dominant_fraction = float(counts.max() / n_total)
        expected_fraction = float(counts[expected_frame] / n_total)
        enrichment = calculate_frame_enrichment(
            (f0, f1, f2), expected_frame=expected_frame
        )
        entropy = calculate_entropy_bias((f0, f1, f2))

        phase = float("nan")
        if compute_phase_score and n_total >= max(int(th["min_reads_per_gene"]), 1):
            codon_counts = (
                group.groupby("codon_index")["frame"]
                .value_counts()
                .unstack(fill_value=0)
                .reindex(columns=[0, 1, 2], fill_value=0)
                .astype(int)
                .values
                .tolist()
            )
            phase = calculate_phase_score(codon_counts)

        qc_call = _qc_call_from_fraction(
            n_total,
            expected_fraction,
            min_reads=int(th["min_reads_per_gene"]),
            good_threshold=float(th["good_frame_fraction"]),
            warn_threshold=float(th["warn_frame_fraction"]),
        )
        rows.append({
            "sample_id": str(sample),
            "gene": str(transcript),
            "transcript_id": str(transcript),
            "read_length": -1,
            "n_frame0": int(counts[0]),
            "n_frame1": int(counts[1]),
            "n_frame2": int(counts[2]),
            "n_total": n_total,
            "frame0_fraction": f0,
            "frame1_fraction": f1,
            "frame2_fraction": f2,
            "expected_frame": int(expected_frame),
            "expected_frame_fraction": expected_fraction,
            "expected_frame_enrichment": enrichment,
            "dominant_frame": dominant_frame,
            "dominant_frame_fraction": dominant_fraction,
            "entropy_bias": entropy,
            "phase_score": phase,
            "qc_call": qc_call,
        })
    return pd.DataFrame(rows)


def build_qc_summary(
    frame_by_sample_length: pd.DataFrame,
    *,
    site_type: str = "p",
    thresholds: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Collapse the per-length table to one row per sample."""
    th = {**QC_THRESHOLDS_DEFAULT, **(thresholds or {})}
    if frame_by_sample_length.empty:
        return pd.DataFrame(
            columns=[
                "sample_id", "site_type", "n_total_sites",
                "best_read_length",
                "best_read_length_expected_frame_fraction",
                "best_read_length_dominant_frame",
                "best_read_length_dominant_fraction",
                "global_expected_frame_fraction",
                "global_entropy_bias",
                "n_good_lengths", "n_warn_lengths",
                "n_poor_lengths", "n_low_depth_lengths",
                "overall_qc_call", "notes",
            ],
        )

    rows: list[dict] = []
    for sample, group in frame_by_sample_length.groupby("sample_id"):
        n_total_sites = int(group["n_reads_cds"].astype(int).sum())
        if n_total_sites <= 0:
            rows.append(_empty_qc_row(str(sample), site_type))
            continue
        best_idx = group["expected_frame_fraction"].astype(float).idxmax()
        best_row = group.loc[best_idx]
        # Global frame fractions: depth-weighted across lengths.
        weights = group["n_reads_cds"].astype(float)
        weights_sum = float(weights.sum())
        if weights_sum <= 0:
            rows.append(_empty_qc_row(str(sample), site_type))
            continue
        f0 = float((group["frame0_fraction"].astype(float) * weights).sum() / weights_sum)
        f1 = float((group["frame1_fraction"].astype(float) * weights).sum() / weights_sum)
        f2 = float((group["frame2_fraction"].astype(float) * weights).sum() / weights_sum)
        global_expected = (f0, f1, f2)[0]
        global_entropy = calculate_entropy_bias((f0, f1, f2))
        calls = group["qc_call"].astype(str).value_counts().to_dict()
        n_good = int(calls.get("good", 0))
        n_warn = int(calls.get("warn", 0))
        n_poor = int(calls.get("poor", 0))
        n_low_depth = int(calls.get("low_depth", 0))
        if n_good >= 1 and n_total_sites >= int(th["min_reads_per_length"]):
            overall = "good"
            note = "at least one read length cleared the good threshold"
        elif n_good + n_warn >= 1:
            overall = "warn"
            note = "no length passed good but at least one is in warn band"
        elif n_low_depth == len(group):
            overall = "low_depth"
            note = "every read length is below min_reads_per_length"
        else:
            overall = "poor"
            note = "no length cleared the warn threshold"
        rows.append({
            "sample_id": str(sample),
            "site_type": site_type,
            "n_total_sites": n_total_sites,
            "best_read_length": int(best_row["read_length"]),
            "best_read_length_expected_frame_fraction": float(
                best_row["expected_frame_fraction"]
            ),
            "best_read_length_dominant_frame": int(best_row["dominant_frame"]),
            "best_read_length_dominant_fraction": float(
                best_row["dominant_frame"] == best_row.get("expected_frame", 0)
            ) and float(best_row["expected_frame_fraction"])
                or float(max(
                    best_row["frame0_fraction"],
                    best_row["frame1_fraction"],
                    best_row["frame2_fraction"],
                )),
            "global_expected_frame_fraction": float(global_expected),
            "global_entropy_bias": float(global_entropy),
            "n_good_lengths": n_good,
            "n_warn_lengths": n_warn,
            "n_poor_lengths": n_poor,
            "n_low_depth_lengths": n_low_depth,
            "overall_qc_call": overall,
            "notes": note,
        })
    return pd.DataFrame(rows)


def _empty_qc_row(sample_id: str, site_type: str) -> dict:
    return {
        "sample_id": sample_id,
        "site_type": site_type,
        "n_total_sites": 0,
        "best_read_length": -1,
        "best_read_length_expected_frame_fraction": 0.0,
        "best_read_length_dominant_frame": -1,
        "best_read_length_dominant_fraction": 0.0,
        "global_expected_frame_fraction": 0.0,
        "global_entropy_bias": 0.0,
        "n_good_lengths": 0,
        "n_warn_lengths": 0,
        "n_poor_lengths": 0,
        "n_low_depth_lengths": 0,
        "overall_qc_call": "low_depth",
        "notes": "no assigned sites",
    }


# ---------------------------------------------------------------------------
# Markdown summary writer + bundle orchestrator
# ---------------------------------------------------------------------------


def write_qc_summary_markdown(
    qc_summary: pd.DataFrame,
    path: Path,
    *,
    site_type: str = "p",
) -> Path:
    """Write a per-sample markdown summary of the periodicity QC."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = ["# 3-nt periodicity QC summary", ""]
    if qc_summary.empty:
        lines.append("(no samples had assignable P-sites)")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path
    for _, row in qc_summary.iterrows():
        lines.extend([
            f"## {row['sample_id']}",
            "",
            f"- Site: **{site_type.upper()}-site**",
            f"- Overall call: **{row['overall_qc_call']}** "
            f"({row['notes']})",
            f"- Best read length: **{int(row['best_read_length'])} nt** "
            f"(expected-frame fraction "
            f"{float(row['best_read_length_expected_frame_fraction']):.3f}, "
            f"dominant frame "
            f"{int(row['best_read_length_dominant_frame'])})",
            f"- Global expected-frame fraction: "
            f"{float(row['global_expected_frame_fraction']):.3f}",
            f"- Global entropy bias: "
            f"{float(row['global_entropy_bias']):.3f}",
            f"- Read-length QC counts: good={int(row['n_good_lengths'])}, "
            f"warn={int(row['n_warn_lengths'])}, "
            f"poor={int(row['n_poor_lengths'])}, "
            f"low_depth={int(row['n_low_depth_lengths'])}",
            "",
        ])
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    return path


def run_periodicity_qc_bundle(
    *,
    frame_by_length: pd.DataFrame,
    bed_with_psite: pd.DataFrame | None,
    annotation_df: pd.DataFrame | None,
    samples: Iterable[str],
    output_dir: Path,
    site_type: str = "p",
    expected_frame: int = 0,
    thresholds: dict[str, float] | None = None,
    compute_phase_score: bool = False,
) -> dict:
    """Write the full periodicity QC bundle and return the in-memory tables.

    The orchestrator is intentionally tolerant of missing inputs: when
    ``bed_with_psite`` / ``annotation_df`` are not provided the
    gene-level table is skipped but the per-(sample, length) and QC
    summary still write.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    by_length = build_frame_counts_by_sample_length(
        frame_by_length,
        expected_frame=expected_frame,
        thresholds=thresholds,
    )
    by_length_path = output_dir / "frame_counts_by_sample_length.tsv"
    by_length.to_csv(by_length_path, sep="\t", index=False, na_rep="")

    qc_summary = build_qc_summary(
        by_length, site_type=site_type, thresholds=thresholds,
    )
    qc_summary_path = output_dir / "qc_summary.tsv"
    qc_summary.to_csv(qc_summary_path, sep="\t", index=False, na_rep="")
    qc_summary_md = output_dir / "qc_summary.md"
    write_qc_summary_markdown(qc_summary, qc_summary_md, site_type=site_type)

    gene_table = pd.DataFrame()
    gene_path = output_dir / "gene_periodicity.tsv"
    if bed_with_psite is not None and annotation_df is not None:
        gene_table = build_gene_periodicity(
            bed_with_psite,
            annotation_df,
            samples=samples,
            expected_frame=expected_frame,
            thresholds=thresholds,
            compute_phase_score=compute_phase_score,
        )
        gene_table.to_csv(gene_path, sep="\t", index=False, na_rep="")
    elif not gene_path.exists():
        # Header-only sentinel so outputs_index can advertise the path
        # consistently across runs.
        gene_path.write_text(
            "sample_id\tgene\ttranscript_id\tread_length\tn_frame0\t"
            "n_frame1\tn_frame2\tn_total\tframe0_fraction\t"
            "frame1_fraction\tframe2_fraction\texpected_frame\t"
            "expected_frame_fraction\texpected_frame_enrichment\t"
            "dominant_frame\tdominant_frame_fraction\tentropy_bias\t"
            "phase_score\tqc_call\n",
            encoding="utf-8",
        )

    return {
        "frame_counts_by_sample_length": by_length,
        "qc_summary": qc_summary,
        "gene_periodicity": gene_table,
        "paths": {
            "frame_counts_by_sample_length": by_length_path,
            "qc_summary": qc_summary_path,
            "qc_summary_md": qc_summary_md,
            "gene_periodicity": gene_path,
        },
    }
