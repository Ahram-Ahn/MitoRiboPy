"""Export structure-density tables derived from footprint-density CSV files."""

from __future__ import annotations

import os
from math import log2

import numpy as np
import pandas as pd

from ..console import log_info, log_warning


def run_structure_density_export(
    translation_profile_dir: str,
    structure_density_norm_perc: float = 0.98,
    output_dir: str = "analysis_results/structure_density",
    site_column: str = "A_site",
) -> None:
    """Export scaled structure-density values from footprint-density CSV files."""
    log_info(
        "DENSITY",
        "Searching for footprint-density CSV files in "
        f"'{translation_profile_dir}' using column '{site_column}'.",
    )
    os.makedirs(output_dir, exist_ok=True)

    footprint_csvs: list[str] = []
    for root, _, files in os.walk(translation_profile_dir):
        for filename in files:
            if filename.endswith("_footprint_density.csv"):
                footprint_csvs.append(os.path.join(root, filename))

    if not footprint_csvs:
        log_warning("DENSITY", "No footprint-density CSV files found; skipping export.")
        return

    log_info(
        "DENSITY",
        f"Found {len(footprint_csvs)} footprint-density file(s); "
        f"normalization percentile={structure_density_norm_perc}.",
    )

    for csv_file in sorted(footprint_csvs):
        base_name = os.path.basename(csv_file)
        transcript = base_name.replace("_footprint_density.csv", "")
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(csv_file)))
        sample_output_dir = os.path.join(output_dir, sample_name)
        os.makedirs(sample_output_dir, exist_ok=True)

        try:
            df = pd.read_csv(csv_file)
        except Exception as exc:
            log_warning("DENSITY", f"Error reading {csv_file}: {exc}")
            continue

        if site_column not in df.columns:
            log_warning(
                "DENSITY",
                f"Column '{site_column}' not found in {csv_file}; skipping.",
            )
            continue

        coverage = df[site_column].to_numpy()
        log2_density = np.array(
            [0.0 if value <= 0 else log2(value + 1) for value in coverage],
            dtype=float,
        )
        df["log2_density"] = log2_density

        nonzero = log2_density[log2_density > 0]
        if len(nonzero) > 0:
            cap_value = np.percentile(nonzero, structure_density_norm_perc * 100)
        else:
            cap_value = 1.0
        cap_value = max(float(cap_value), 1e-9)

        scaled_density = np.clip(log2_density, None, cap_value) / cap_value
        df["scaled_density"] = scaled_density

        scaled_density_txt = os.path.join(sample_output_dir, f"{transcript}_scaled_density.txt")
        with open(scaled_density_txt, "w", encoding="utf-8") as output_handle:
            for value in scaled_density:
                output_handle.write(f"{value:.4f}\n")

        structure_density_csv = os.path.join(
            sample_output_dir,
            f"{transcript}_structure_density.csv",
        )
        df.to_csv(structure_density_csv, index=False)

        log_info(
            "DENSITY",
            f"{sample_name}/{transcript} -> wrote {scaled_density_txt} and {structure_density_csv}",
        )

    log_info("DENSITY", f"All structure-density outputs saved to {output_dir}")
