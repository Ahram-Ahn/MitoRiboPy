"""Codon-level correlation analysis between a base sample and comparison samples."""

from __future__ import annotations

import math
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import linregress

from ..console import iter_with_progress, log_info, log_warning
from ..plotting.style import apply_publication_style

apply_publication_style()


# ---------------------------------------------------------------------------
# Publication-quality scatter helper (no external deps)
# ---------------------------------------------------------------------------

# Okabe-Ito colour-blind safe palette, mapped onto the codon Category
# levels we expect to see (start / stop / standard amino-acid groups).
# Falls through to a Set2 fallback for unknown categories so a future
# annotation file with new categories does not crash the plot.
_CATEGORY_PALETTE: dict[str, str] = {
    "start": "#E69F00",
    "stop": "#D55E00",
    "polar": "#0072B2",
    "non-polar": "#009E73",
    "acidic": "#CC79A7",
    "basic": "#56B4E9",
    "aromatic": "#F0E442",
    "other": "#999999",
}


def _categorical_colour(category: str, fallback_cycle) -> str:
    """Return a deterministic colour for *category*."""
    key = str(category).strip().lower()
    if key in _CATEGORY_PALETTE:
        return _CATEGORY_PALETTE[key]
    # Fall back to a deterministic seaborn Set2 colour for novel
    # categories so the legend remains stable across runs.
    return fallback_cycle[hash(key) % len(fallback_cycle)]


def _place_residual_labels(
    ax,
    labels_df: pd.DataFrame,
    *,
    x_col: str,
    y_col: str,
    label_col: str,
    axis_min: float,
    axis_max: float,
) -> None:
    """Greedy label placement that avoids overlap without adjustText.

    For each label point we evaluate four cardinal anchor positions
    (NE / NW / SE / SW) at a fixed nudge distance from the point, score
    each by its distance to every already-placed label and to every
    point cluster, and keep the best candidate. A leader line connects
    the label to its point so the reader can match them up unambiguously.
    """
    if labels_df.empty:
        return

    span = max(axis_max - axis_min, 1e-9)
    nudge = span * 0.04
    placed_xy: list[tuple[float, float]] = []
    point_xy = list(zip(labels_df[x_col].astype(float), labels_df[y_col].astype(float)))
    bbox = dict(
        boxstyle="round,pad=0.25",
        facecolor="white",
        edgecolor="0.5",
        linewidth=0.6,
        alpha=0.85,
    )

    for x, y, label in zip(
        labels_df[x_col].astype(float),
        labels_df[y_col].astype(float),
        labels_df[label_col],
    ):
        candidates = [
            (x + nudge, y + nudge, "left", "bottom"),
            (x - nudge, y + nudge, "right", "bottom"),
            (x + nudge, y - nudge, "left", "top"),
            (x - nudge, y - nudge, "right", "top"),
        ]
        best = None
        best_score = -math.inf
        for cx, cy, ha, va in candidates:
            if not (axis_min <= cx <= axis_max and axis_min <= cy <= axis_max):
                continue
            # Repulsion score: distance to nearest other label + distance
            # to nearest other point. Higher is better.
            label_d = (
                min(math.hypot(cx - lx, cy - ly) for lx, ly in placed_xy)
                if placed_xy
                else span
            )
            point_d = min(
                math.hypot(cx - px, cy - py)
                for px, py in point_xy
                if (px, py) != (x, y)
            )
            score = min(label_d, point_d)
            if score > best_score:
                best = (cx, cy, ha, va)
                best_score = score
        if best is None:
            best = (x + nudge, y + nudge, "left", "bottom")

        cx, cy, ha, va = best
        placed_xy.append((cx, cy))
        ax.annotate(
            str(label),
            xy=(x, y),
            xytext=(cx, cy),
            ha=ha,
            va=va,
            fontsize=8,
            fontweight="bold",
            bbox=bbox,
            arrowprops=dict(
                arrowstyle="-",
                color="0.4",
                lw=0.6,
                shrinkA=0,
                shrinkB=2,
            ),
        )


def run_codon_correlation(
    translation_profile_dir: str,
    samples: list[str],
    base_sample: str,
    column: str = "CoverageDivFreq",
    output_dir: str | None = None,
    mask_method: str = "percentile",
    mask_percentile: float = 0.99,
    mask_threshold: float | None = None,
    site: str = "p",
) -> None:
    """Run codon correlation analysis comparing the base sample to every
    other sample for the requested site.

    ``site`` selects which codon-usage table to read:
      * ``"p"``: read ``p_site_codon_usage_total.csv``; for stop-codon
        rows (``AA=="*"``), substitute the A-site value because the
        P-site CSV applies stop-codon masking.
      * ``"a"``: read ``a_site_codon_usage_total.csv`` directly (no
        stop-codon override; A-site usage never masks the stop).

    Two versions (``all`` and ``masked``) are produced based on
    thresholds. CSV + SVG + PNG outputs are saved for each.
    """
    site = str(site).lower()
    if site not in {"p", "a"}:
        site = "p"
    site_label = "P-site" if site == "p" else "A-site"
    log_info(
        "COR",
        f"Starting {site_label} correlation vs base='{base_sample}', column='{column}'.",
    )
    if output_dir is None:
        output_dir = os.path.join(translation_profile_dir, "codon_correlation")
    os.makedirs(output_dir, exist_ok=True)

    primary_basename = (
        "p_site_codon_usage_total.csv" if site == "p" else "a_site_codon_usage_total.csv"
    )
    a_site_basename = "a_site_codon_usage_total.csv"

    def override_stop_values(primary_csv: str, a_csv: str) -> pd.DataFrame:
        """Replace stop-codon rows with A-site values when masking the P-site.

        For ``site=="a"`` no override is needed (A-site usage already
        contains the stop-codon row unmasked); we still pass through
        :func:`pd.read_csv` for a consistent return type.
        """
        df = pd.read_csv(primary_csv)
        if site == "p" and os.path.isfile(a_csv):
            a_df = pd.read_csv(a_csv)
            for idx, row in df.iterrows():
                if row["AA"] == "*":
                    match = a_df[
                        (a_df["Codon"] == row["Codon"])
                        & (a_df["AA"] == row["AA"])
                        & (a_df["Category"] == row["Category"])
                    ]
                    if not match.empty:
                        df.at[idx, column] = match.iloc[0]["CoverageDivFreq"]
        return df

    base_p_csv = os.path.join(
        translation_profile_dir, base_sample, "codon_usage", primary_basename
    )
    base_a_csv = os.path.join(
        translation_profile_dir, base_sample, "codon_usage", a_site_basename
    )
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
        sample_p_csv = os.path.join(
            translation_profile_dir, sample_name, "codon_usage", primary_basename
        )
        sample_a_csv = os.path.join(
            translation_profile_dir, sample_name, "codon_usage", a_site_basename
        )
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

            # Regression on the CoverageDivFreq values.
            slope, intercept, r_value, p_value, _ = linregress(
                merged_current["base_val"], merged_current["sample_val"]
            )
            merged_current["predicted"] = (
                slope * merged_current["base_val"] + intercept
            )
            merged_current["residual"] = abs(
                merged_current["sample_val"] - merged_current["predicted"]
            )
            # Label the 10 codons farthest from the regression line.
            top10 = merged_current.nlargest(10, "residual").copy()
            top10["label"] = top10["Codon"].astype(str) + " (" + top10["AA"].astype(str) + ")"

            # ----- Publication-quality scatter --------------------------
            #
            # Design choices:
            # * square aspect ratio with shared axis limits + identity
            #   line so deviation from y=x is visually obvious;
            # * alpha + small marker edge to mute overplotting at the
            #   low-coverage corner;
            # * Okabe-Ito colour-blind safe palette per Category;
            # * leader-line labels for the 10 highest-residual codons,
            #   placed greedily to avoid overlap (no adjustText dep);
            # * compact stat box with r, R^2, p, slope, intercept, N;
            # * SVG (vector, for figures) + 300 dpi PNG (for slides)
            #   side by side.
            fallback_palette = list(sns.color_palette("Set2", 8).as_hex())
            categories = list(merged_current["Category"].astype(str).unique())
            colour_map = {
                cat: _categorical_colour(cat, fallback_palette) for cat in categories
            }

            fig, ax = plt.subplots(figsize=(7.5, 7.5))
            for category, group in merged_current.groupby("Category", sort=False):
                ax.scatter(
                    group["base_val"],
                    group["sample_val"],
                    s=42,
                    alpha=0.75,
                    color=colour_map[str(category)],
                    edgecolor="white",
                    linewidth=0.6,
                    label=str(category),
                )

            # Square axes anchored on the data range so y=x is informative.
            data_min = float(
                min(merged_current["base_val"].min(), merged_current["sample_val"].min())
            )
            data_max = float(
                max(merged_current["base_val"].max(), merged_current["sample_val"].max())
            )
            pad = max((data_max - data_min) * 0.05, 1e-6)
            axis_min = data_min - pad
            axis_max = data_max + pad
            ax.set_xlim(axis_min, axis_max)
            ax.set_ylim(axis_min, axis_max)
            ax.set_aspect("equal", adjustable="box")

            identity = np.linspace(axis_min, axis_max, 100)
            ax.plot(
                identity,
                identity,
                color="0.6",
                linestyle=":",
                linewidth=1.0,
                label="y = x",
            )
            x_vals = np.linspace(axis_min, axis_max, 100)
            ax.plot(
                x_vals,
                slope * x_vals + intercept,
                color="#D55E00",
                linestyle="--",
                linewidth=1.4,
                label=f"OLS fit (r = {r_value:.3f})",
            )

            _place_residual_labels(
                ax,
                top10,
                x_col="base_val",
                y_col="sample_val",
                label_col="label",
                axis_min=axis_min,
                axis_max=axis_max,
            )

            r_squared = r_value * r_value
            stat_lines = [
                f"r = {r_value:.3f}",
                f"R$^2$ = {r_squared:.3f}",
                f"p = {p_value:.2e}",
                f"slope = {slope:.3f}",
                f"intercept = {intercept:.3f}",
                f"N = {len(merged_current)}",
            ]
            ax.text(
                0.03,
                0.97,
                "\n".join(stat_lines),
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=9,
                family="monospace",
                bbox=dict(
                    boxstyle="round,pad=0.4",
                    facecolor="white",
                    edgecolor="0.7",
                    linewidth=0.6,
                    alpha=0.9,
                ),
            )

            title_suffix = "all codons" if version == "all" else "outlier-masked"
            ax.set_title(
                f"Codon usage: {base_sample} vs {sample_name} ({title_suffix})",
                fontsize=12,
                fontweight="bold",
                pad=12,
            )
            ax.set_xlabel(f"{base_sample}  ({column})", fontsize=11)
            ax.set_ylabel(f"{sample_name}  ({column})", fontsize=11)
            ax.grid(True, which="major", alpha=0.25, linewidth=0.6)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.legend(
                title="Category",
                title_fontsize=10,
                fontsize=9,
                loc="lower right",
                frameon=True,
                framealpha=0.9,
                edgecolor="0.7",
            )
            fig.tight_layout()

            out_svg = os.path.join(
                output_dir, f"{base_sample}_vs_{sample_name}_{version}.svg"
            )
            out_png = os.path.join(
                output_dir, f"{base_sample}_vs_{sample_name}_{version}.png"
            )
            fig.savefig(out_svg)
            fig.savefig(out_png, dpi=300)
            plt.close(fig)
            log_info("COR", f"Plot saved => {out_svg} (+ 300 dpi PNG)")

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
