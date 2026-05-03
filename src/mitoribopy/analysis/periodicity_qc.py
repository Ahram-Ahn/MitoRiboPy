"""Fourier periodicity orchestrator.

Replaces the legacy v0.6.x QC bundle (frame_counts_*, gene_periodicity,
qc_summary, phase_score, frame heatmap plots, etc.). The current QC
contract is the metagene Fourier spectrum + period-3 spectral ratio
documented in :mod:`mitoribopy.analysis.fourier_spectrum`.

Outputs written by :func:`run_periodicity_qc_bundle`:

  fourier_spectrum_combined.tsv         — per-(sample, length, gene_set,
                                          region) metagene amplitude
                                          curve over period 2-10 nt.
  fourier_period3_score_combined.tsv    — per-(sample, length, gene_set,
                                          region) headline scalars
                                          (amp_at_3nt, spectral_ratio_3nt,
                                          snr_call).
  fourier_spectrum/<sample>/*.{png,svg} — three two-panel figures per
                                          (sample, length): combined,
                                          ATP86, ND4L4.
  periodicity.metadata.json             — every knob and gene-set
                                          definition that produced the
                                          tables.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .fourier_spectrum import (
    DEFAULT_BOOTSTRAP_N,
    DEFAULT_CI_ALPHA,
    DEFAULT_DROP_CODONS_AFTER_START,
    DEFAULT_DROP_CODONS_BEFORE_STOP,
    DEFAULT_MIN_MEAN_COVERAGE,
    DEFAULT_MIN_TOTAL_COUNTS,
    DEFAULT_PERIOD_GRID,
    DEFAULT_PERMUTATIONS_N,
    DEFAULT_RANDOM_SEED,
    DEFAULT_WINDOW_NT,
    GENE_SETS,
    REGIONS,
    Site,
    build_fourier_period3_score_combined_table,
    build_fourier_spectrum_combined_table,
    extract_per_gene_normalized_tracks,
)


__all__ = [
    "run_periodicity_qc_bundle",
]


def run_periodicity_qc_bundle(
    *,
    bed_with_psite: pd.DataFrame | None,
    annotation_df: pd.DataFrame | None,
    samples: Iterable[str],
    output_dir: Path,
    site_type: Site = "p",
    window_nt: int = DEFAULT_WINDOW_NT,
    drop_codons_after_start: int = DEFAULT_DROP_CODONS_AFTER_START,
    drop_codons_before_stop: int = DEFAULT_DROP_CODONS_BEFORE_STOP,
    min_mean_coverage: float = DEFAULT_MIN_MEAN_COVERAGE,
    min_total_counts: int = DEFAULT_MIN_TOTAL_COUNTS,
    render_plots: bool = True,
    n_bootstrap: int = DEFAULT_BOOTSTRAP_N,
    n_permutations: int = DEFAULT_PERMUTATIONS_N,
    ci_alpha: float = DEFAULT_CI_ALPHA,
    random_seed: int = DEFAULT_RANDOM_SEED,
    compute_stats: bool = True,
) -> dict:
    """Write the metagene Fourier bundle and return in-memory tables.

    Tolerant of missing inputs: when ``bed_with_psite`` or
    ``annotation_df`` is ``None`` the function writes empty (header-
    only) TSVs so downstream outputs_index can advertise the paths
    consistently.

    The statistical hardening knobs (``n_bootstrap``, ``n_permutations``,
    ``ci_alpha``, ``random_seed``, ``compute_stats``) are forwarded to
    :func:`build_fourier_period3_score_combined_table`. Set
    ``compute_stats=False`` for a fast smoke run that produces only the
    point-estimate columns.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    spectrum_path = output_dir / "fourier_spectrum_combined.tsv"
    score_path = output_dir / "fourier_period3_score_combined.tsv"
    plots_dir = output_dir / "fourier_spectrum"

    spectrum_table = pd.DataFrame()
    score_table = pd.DataFrame()

    if bed_with_psite is not None and annotation_df is not None:
        tracks = extract_per_gene_normalized_tracks(
            bed_with_psite,
            annotation_df,
            samples=samples,
            window_nt=int(window_nt),
            drop_codons_after_start=int(drop_codons_after_start),
            drop_codons_before_stop=int(drop_codons_before_stop),
            min_mean_coverage=float(min_mean_coverage),
            min_total_counts=int(min_total_counts),
            site=str(site_type).lower() or "a",  # type: ignore[arg-type]
        )
        spectrum_table = build_fourier_spectrum_combined_table(
            tracks, periods=DEFAULT_PERIOD_GRID,
        )
        score_table = build_fourier_period3_score_combined_table(
            tracks,
            periods=DEFAULT_PERIOD_GRID,
            n_bootstrap=int(n_bootstrap),
            n_permutations=int(n_permutations),
            ci_alpha=float(ci_alpha),
            random_seed=int(random_seed),
            compute_stats=bool(compute_stats),
        )

    spectrum_table.to_csv(spectrum_path, sep="\t", index=False, na_rep="")
    score_table.to_csv(score_path, sep="\t", index=False, na_rep="")

    if render_plots and not spectrum_table.empty:
        from ..plotting.fourier_spectrum_plots import (
            render_fourier_spectrum_panels,
        )
        render_fourier_spectrum_panels(
            spectrum_table,
            score_table,
            output_dir=plots_dir,
            source_data_relpath=spectrum_path.name,
        )

    metadata = {
        "fourier_window_nt": int(window_nt),
        "drop_codons_after_start": int(drop_codons_after_start),
        "drop_codons_before_stop": int(drop_codons_before_stop),
        "min_mean_coverage": float(min_mean_coverage),
        "min_total_counts": int(min_total_counts),
        "site_type": str(site_type),
        "regions": list(REGIONS),
        "gene_sets": list(GENE_SETS),
        "period_grid_count": int(np.size(DEFAULT_PERIOD_GRID)),
        "period_grid_first": float(DEFAULT_PERIOD_GRID[0]),
        "period_grid_last": float(DEFAULT_PERIOD_GRID[-1]),
        "method": "metagene_dft",
        # Statistical hardening (v0.9.0+): record exactly what produced
        # the CI / p columns so a downstream reviewer can reproduce
        # them by name.
        "compute_stats": bool(compute_stats),
        "n_bootstrap": int(n_bootstrap),
        "n_permutations": int(n_permutations),
        "ci_alpha": float(ci_alpha),
        "ci_method": "percentile_over_genes" if compute_stats else "disabled",
        "null_method": "circular_shift_per_gene" if compute_stats else "disabled",
        "random_seed": int(random_seed),
    }
    (output_dir / "periodicity.metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    return {
        "fourier_spectrum_combined": spectrum_table,
        "fourier_period3_score_combined": score_table,
        "paths": {
            "fourier_spectrum_combined": spectrum_path,
            "fourier_period3_score_combined": score_path,
            "fourier_spectrum_plots_dir": plots_dir,
        },
    }
