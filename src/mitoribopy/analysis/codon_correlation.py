"""Codon-level correlation analysis between a base sample and comparison samples."""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import linregress

from ..console import iter_with_progress, log_info, log_warning
from ..plotting.style import apply_publication_style

apply_publication_style()


def run_codon_correlation(
    translation_profile_dir: str,
    samples: list[str],
    base_sample: str,
    column: str = "CoverageDivFreq",
    output_dir: str | None = None,
    mask_method: str = "percentile",
    mask_percentile: float = 0.99,
    mask_threshold: float | None = None,
) -> None:
    """
    Run codon correlation analysis comparing the base sample to every other sample.
    For each pair, for stop codons (AA=="*"), the value from the A-site analysis is used.
    Two versions ("all" and "masked") are produced based on thresholds.
    CSV files and SVG plots are saved for each version.
    """
    log_info("COR", f"Starting correlation vs base='{base_sample}', column='{column}'.")
    if output_dir is None:
        output_dir = os.path.join(translation_profile_dir, "codon_correlation")
    os.makedirs(output_dir, exist_ok=True)

    # Helper function: override stop codon value with A-site value
    def override_stop_values(p_csv: str, a_csv: str) -> pd.DataFrame:
        df = pd.read_csv(p_csv)
        if os.path.isfile(a_csv):
            a_df = pd.read_csv(a_csv)
            # For rows with AA=="*", update value
            for idx, row in df.iterrows():
                if row["AA"] == "*":
                    match = a_df[(a_df["Codon"] == row["Codon"]) &
                                 (a_df["AA"] == row["AA"]) &
                                 (a_df["Category"] == row["Category"])]
                    if not match.empty:
                        df.at[idx, column] = match.iloc[0]["CoverageDivFreq"]
        return df

    # Read base sample p-site CSV and override stop codons
    base_p_csv = os.path.join(translation_profile_dir, base_sample, "codon_usage", "codon_usage_total.csv")
    base_a_csv = os.path.join(translation_profile_dir, base_sample, "codon_usage", "a_site_codon_usage_total.csv")
    if not os.path.isfile(base_p_csv):
        log_warning("COR", f"Base sample CSV not found => {base_p_csv}")
        return
    base_df = override_stop_values(base_p_csv, base_a_csv)
    if column not in base_df.columns:
        log_warning("COR", f"Column '{column}' not found in base sample CSV => {base_p_csv}")
        return
    base_df = base_df.rename(columns={column: "base_val"})
    base_df = base_df[["Codon", "AA", "Category", "base_val"]]

    corr_records = []

    for sample_name in iter_with_progress(
        list(samples),
        component="COR",
        noun="sample",
        labeler=str,
    ):
        if sample_name == base_sample:
            continue
        sample_p_csv = os.path.join(translation_profile_dir, sample_name, "codon_usage", "codon_usage_total.csv")
        sample_a_csv = os.path.join(translation_profile_dir, sample_name, "codon_usage", "a_site_codon_usage_total.csv")
        if not os.path.isfile(sample_p_csv):
            log_warning("COR", f"Missing codon-usage file => {sample_p_csv}; skipping.")
            continue
        s_df = override_stop_values(sample_p_csv, sample_a_csv)
        if column not in s_df.columns:
            log_warning("COR", f"Column '{column}' not found in {sample_p_csv}; skipping.")
            continue
        s_df = s_df.rename(columns={column: "sample_val"})
        s_df = s_df[["Codon", "AA", "Category", "sample_val"]]

        # Merge based on Codon, AA, and Category
        merged = pd.merge(base_df, s_df, on=["Codon", "AA", "Category"], how="inner")
        if merged.empty:
            log_warning("COR", f"No overlapping codons found for sample {sample_name}; skipping.")
            continue

        threshold = mask_threshold
        if threshold is None and mask_method == "percentile":
            both_vals = pd.concat([merged["base_val"], merged["sample_val"]], axis=0)
            threshold = float(both_vals.quantile(mask_percentile))

        versions = ["all", "masked"] if threshold is not None else ["all"]
        for version in versions:
            if version == "masked":
                merged_current = merged[(merged["base_val"] <= threshold) & (merged["sample_val"] <= threshold)]
                if merged_current.empty:
                    log_warning(
                        "COR",
                        f"No data remain after masking for sample {sample_name} ({version}).",
                    )
                    continue
            else:
                merged_current = merged.copy()
            merged_current = merged_current.copy()

            out_csv = os.path.join(output_dir, f"{base_sample}_vs_{sample_name}_{version}.csv")
            merged_current.to_csv(out_csv, index=False)
            log_info("COR", f"Wrote correlation CSV => {out_csv}")

            # Regression using CoverageDivFreq values
            slope, intercept, r_value, _, _ = linregress(merged_current["base_val"], merged_current["sample_val"])
            merged_current["predicted"] = slope * merged_current["base_val"] + intercept
            merged_current["residual"] = abs(merged_current["sample_val"] - merged_current["predicted"])
            top10 = merged_current.nlargest(10, "residual")

            plt.figure(figsize=(8, 8))
            sns.scatterplot(x="base_val", y="sample_val", hue="Category", data=merged_current, s=100)
            x_vals = merged_current["base_val"]
            plt.plot(x_vals, slope * x_vals + intercept, color="red", linestyle="--", label=f"Regression (r={r_value:.3f})")
            for _, row in top10.iterrows():
                plt.text(row["base_val"] + 0.001, row["sample_val"] + 0.001,
                         f"{row['Codon']} - {row['AA']}", fontsize=9, fontweight="bold")
            plt.title(f"{base_sample} vs. {sample_name} ({version})\ncolumn={column}", fontsize=14, fontweight="bold")
            plt.xlabel(f"{base_sample} ({column})", fontsize=12, fontweight="bold")
            plt.ylabel(f"{sample_name} ({column})", fontsize=12, fontweight="bold")
            plt.legend(title="Category", bbox_to_anchor=(1.05, 1), loc="upper left")
            plt.tight_layout()

            out_svg = os.path.join(output_dir, f"{base_sample}_vs_{sample_name}_{version}.svg")
            plt.savefig(out_svg)
            plt.close()
            log_info("COR", f"Plot saved => {out_svg}")

            corr_records.append({
                "Base_Sample": base_sample,
                "Other_Sample": sample_name,
                "Version": version,
                "Pearson_r": r_value,
                "NumPoints": len(merged_current)
            })

    cor_df = pd.DataFrame(corr_records)
    out_summary = os.path.join(output_dir, f"codon_correlation_summary_{base_sample}.csv")
    cor_df.to_csv(out_summary, index=False)
    log_info("COR", f"Correlation summary saved => {out_summary}")
    log_info("COR", "Done correlation for all samples.")

# Example usage:
# run_codon_correlation(translation_profile_dir="path/to/analysis_results",
#                       samples=["sample1", "sample2", "sample3"],
#                       base_sample="sample1",
#                       column="CoverageDivFreq")
