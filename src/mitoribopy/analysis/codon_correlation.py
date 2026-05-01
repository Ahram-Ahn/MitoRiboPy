"""Codon-level correlation analysis between a base sample and comparison samples.

Publication-readiness redesign (2026-05): the headline statistic is no
longer a Pearson correlation on raw count-like columns. Raw-count
correlations are dominated by sequencing depth and codon-usage
frequency in the reference, so "high correlation" mostly means "the
samples sequenced to similar depth", not "codon occupancy is similar".

The new pipeline:

* **Metric.** Default ``log2_density_rpm`` — log2 of codon density
  normalized to a million assigned P-sites, with a small additive
  pseudocount to handle zeros.
* **Regression.** Default ``theil_sen`` (median-based, robust to
  outliers; uses :func:`scipy.stats.theilslopes`).
* **Support filter.** Codons with very low raw P-site support in
  *either* sample are flagged as low-support and excluded from the
  primary-figure label set (controlled by ``support_min_raw``).
* **MA / Bland-Altman panel.** Always written alongside the scatter so
  reviewers can see codon-specific shifts that a correlation hides.
* **Label scoring.** Outlier labels are picked by
  ``|log2 fold change| * log10(1 + min_raw_support)`` instead of pure
  residual magnitude — this prevents a low-support codon from dominating
  the figure on rounding noise.

Backwards compatibility: the legacy ``column="CoverageDivFreq"`` argument
is still honoured (it now selects the raw column we transform on); the
legacy per-version CSV / scatter is still emitted; passing
``metric="raw_count"`` reproduces the old behaviour but raises the
``W_CODON_RAW_COUNT_PRIMARY`` warning.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import kendalltau, linregress, spearmanr, theilslopes

from ..console import iter_with_progress, log_info, log_warning
from ..plotting.style import apply_publication_style

apply_publication_style()


# ---------------------------------------------------------------------------
# Metric / regression policy
# ---------------------------------------------------------------------------


_VALID_METRICS = ("log2_density_rpm", "log2_rpm", "linear", "raw_count")
_VALID_REGRESSIONS = ("theil_sen", "ols", "none")
_DEFAULT_PSEUDOCOUNT = 0.5
_DEFAULT_SUPPORT_MIN_RAW = 10
W_CODON_RAW_COUNT_PRIMARY = "W_CODON_RAW_COUNT_PRIMARY"


def _coerce_pseudocount(pseudocount, fallback_series: pd.Series) -> float:
    """Resolve ``pseudocount='auto'`` to half the smallest non-zero value."""
    if isinstance(pseudocount, (int, float)) and not isinstance(pseudocount, bool):
        return float(pseudocount)
    if pseudocount in (None, "auto", "AUTO"):
        nonzero = fallback_series[fallback_series > 0]
        if nonzero.empty:
            return _DEFAULT_PSEUDOCOUNT
        return max(_DEFAULT_PSEUDOCOUNT, float(nonzero.min()) * 0.5)
    return _DEFAULT_PSEUDOCOUNT


def _transform_metric(
    base_vals: pd.Series,
    sample_vals: pd.Series,
    *,
    metric: str,
    pseudocount: float,
) -> tuple[pd.Series, pd.Series, dict]:
    """Apply the chosen metric transform and return (base_t, sample_t, meta).

    ``meta`` records what was done so it can be embedded in the
    ``codon_correlation.metadata.json`` sidecar.
    """
    metric = str(metric).lower()
    if metric == "raw_count":
        return base_vals.astype(float), sample_vals.astype(float), {
            "metric": "raw_count",
            "transform": "identity",
            "pseudocount": 0.0,
        }
    if metric == "linear":
        return base_vals.astype(float), sample_vals.astype(float), {
            "metric": "linear",
            "transform": "identity",
            "pseudocount": 0.0,
        }
    # Both log2 modes apply the same per-sample-RPM scaling (so a
    # depth-difference does not move the headline scatter).
    scale_base = max(float(base_vals.sum()) / 1e6, 1e-12)
    scale_sample = max(float(sample_vals.sum()) / 1e6, 1e-12)
    base_rpm = base_vals.astype(float) / scale_base
    sample_rpm = sample_vals.astype(float) / scale_sample
    base_t = np.log2(base_rpm + pseudocount)
    sample_t = np.log2(sample_rpm + pseudocount)
    return base_t, sample_t, {
        "metric": metric,
        "transform": "log2(value/RPM_scale + pseudocount)",
        "pseudocount": float(pseudocount),
        "rpm_scale_base": float(scale_base),
        "rpm_scale_sample": float(scale_sample),
    }


def _fit_regression(
    x: np.ndarray, y: np.ndarray, *, method: str,
) -> dict:
    """Fit the requested regression and return diagnostics.

    Always returns ``slope``, ``intercept``, ``predicted`` (over ``x``),
    ``residual`` (``y - predicted``), and a ``method`` string. For
    ``method="none"`` the line is the identity ``y = x``, which lets the
    rest of the pipeline keep working.
    """
    method = str(method).lower()
    n = int(min(len(x), len(y)))
    if n < 3 or method == "none":
        slope, intercept = 1.0, 0.0
        predicted = slope * x + intercept
        return {
            "method": "identity" if method == "none" else f"{method}_too_few_points",
            "slope": slope,
            "intercept": intercept,
            "predicted": predicted,
            "residual": y - predicted,
        }
    if method == "theil_sen":
        slope, intercept, lo_slope, hi_slope = theilslopes(y, x)
        predicted = slope * x + intercept
        return {
            "method": "theil_sen",
            "slope": float(slope),
            "intercept": float(intercept),
            "slope_ci_low": float(lo_slope),
            "slope_ci_high": float(hi_slope),
            "predicted": predicted,
            "residual": y - predicted,
        }
    # OLS fallback (back-compat with the old behaviour).
    slope, intercept, r_value, p_value, _ = linregress(x, y)
    predicted = slope * x + intercept
    return {
        "method": "ols",
        "slope": float(slope),
        "intercept": float(intercept),
        "predicted": predicted,
        "residual": y - predicted,
        "r": float(r_value),
        "p": float(p_value),
    }


def _label_score(
    log2_fc: pd.Series, support_min_raw: pd.Series,
) -> pd.Series:
    """``|log2 FC| * log10(1 + min support)`` — biases labels to high-support codons."""
    support_term = np.log10(1.0 + support_min_raw.astype(float))
    return (log2_fc.abs().astype(float) * support_term).fillna(0.0)


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
    *,
    metric: str = "log2_density_rpm",
    regression: str = "theil_sen",
    pseudocount: float | str = "auto",
    support_min_raw: int = _DEFAULT_SUPPORT_MIN_RAW,
    label_top_n: int = 10,
    raw_panel: str = "qc_only",
) -> None:
    """Run codon correlation analysis comparing the base sample to every
    other sample for the requested site.

    Parameters
    ----------
    metric
        Primary plotting / scoring transform. ``log2_density_rpm`` (the
        default) is the publication-recommended choice; ``raw_count``
        reproduces the old behaviour but emits a structured warning.
    regression
        Regression line drawn through the scatter and used to compute
        residuals. ``theil_sen`` (default) is robust to outliers;
        ``ols`` reproduces the old linear regression; ``none`` skips
        the regression and uses an identity line.
    pseudocount
        Additive pseudocount before the log2. ``"auto"`` resolves to
        ``max(0.5, 0.5 * min_nonzero_value)`` per (base, sample) pair.
    support_min_raw
        Minimum raw value required in BOTH samples for a codon to be
        considered "primary"; lower-support codons remain in the data
        table but are excluded from labels and the residual ranking.
    label_top_n
        Number of codons to label on the scatter, ranked by
        ``|log2_fold_change| * log10(1 + min_raw_support)``.
    raw_panel
        ``"qc_only"`` (default): write the raw-count scatter under
        ``raw_count_qc/`` so it is clearly a QC artefact, not a
        publication figure. ``"off"`` skips the raw panel entirely.

    ``site`` selects which codon-usage table to read:
      * ``"p"``: read ``p_site_codon_usage_total.csv``; for stop-codon
        rows (``AA=="*"``), substitute the A-site value because the
        P-site CSV applies stop-codon masking.
      * ``"a"``: read ``a_site_codon_usage_total.csv`` directly (no
        stop-codon override; A-site usage never masks the stop).

    Two versions (``all`` and ``masked``) are produced based on
    thresholds. The function writes:

    * ``codon_correlation_metrics.tsv`` (one row per sample/codon) — the
      authoritative tabular output, suitable for re-plotting.
    * ``{base}_vs_{sample}_{version}.csv`` — legacy per-pair CSVs.
    * ``{base}_vs_{sample}_{version}.svg/.png`` — three-panel publication
      figure (log2 density scatter, MA plot, residual plot).
    * ``codon_correlation.metadata.json`` — sidecar describing the
      metric, regression method, pseudocount, and warnings.
    """
    if metric not in _VALID_METRICS:
        raise ValueError(
            f"metric={metric!r} not in {_VALID_METRICS}"
        )
    if regression not in _VALID_REGRESSIONS:
        raise ValueError(
            f"regression={regression!r} not in {_VALID_REGRESSIONS}"
        )
    raw_metric = metric == "raw_count"
    if raw_metric:
        log_warning(
            "COR",
            f"[{W_CODON_RAW_COUNT_PRIMARY}] metric='raw_count' selected "
            "for primary codon-correlation figure. Raw-count codon "
            "correlation is depth- and abundance-dominated; use "
            "metric='log2_density_rpm' (default) or 'log2_rpm' for "
            "publication figures.",
        )
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
    metrics_records: list[pd.DataFrame] = []

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

            # Resolve pseudocount (auto -> 0.5 * smallest non-zero, floored at 0.5).
            pc = _coerce_pseudocount(
                pseudocount,
                pd.concat([merged_current["base_val"], merged_current["sample_val"]]),
            )
            base_t, sample_t, transform_meta = _transform_metric(
                merged_current["base_val"],
                merged_current["sample_val"],
                metric=metric,
                pseudocount=pc,
            )
            merged_current["base_metric"] = base_t.to_numpy()
            merged_current["sample_metric"] = sample_t.to_numpy()
            merged_current["log2_fold_change"] = (
                merged_current["sample_metric"] - merged_current["base_metric"]
                if metric.startswith("log2")
                else np.log2(
                    (merged_current["sample_val"].astype(float) + pc)
                    / (merged_current["base_val"].astype(float) + pc)
                )
            )
            merged_current["mean_log2_density"] = 0.5 * (
                merged_current["base_metric"] + merged_current["sample_metric"]
            )
            merged_current["support_min_raw"] = (
                merged_current[["base_val", "sample_val"]].min(axis=1).astype(float)
            )
            merged_current["support_total_raw"] = (
                merged_current["base_val"].astype(float)
                + merged_current["sample_val"].astype(float)
            )
            merged_current["include_primary"] = (
                merged_current["support_min_raw"] >= float(support_min_raw)
            )
            merged_current["exclusion_reason"] = np.where(
                merged_current["include_primary"], "none", "low_support"
            )

            primary = merged_current[merged_current["include_primary"]]
            x = primary["base_metric"].to_numpy(dtype=float)
            y = primary["sample_metric"].to_numpy(dtype=float)

            fit = _fit_regression(x, y, method=regression)
            slope = float(fit["slope"])
            intercept = float(fit["intercept"])
            # Predict over the full table (including low-support codons)
            # so the residual column covers every row.
            full_pred = slope * merged_current["base_metric"].to_numpy() + intercept
            merged_current["predicted"] = full_pred
            merged_current["robust_residual"] = (
                merged_current["sample_metric"].to_numpy() - full_pred
            )
            # Legacy "residual" column kept as the absolute version for
            # backwards compatibility with existing downstream parsers.
            merged_current["residual"] = np.abs(merged_current["robust_residual"])

            # Robust label policy: pick by |log2 FC| weighted by support
            # so a low-support codon doesn't dominate the figure.
            merged_current["label_score"] = _label_score(
                merged_current["log2_fold_change"],
                merged_current["support_min_raw"],
            )
            label_pool = merged_current[merged_current["include_primary"]].copy()
            top10 = label_pool.nlargest(int(label_top_n), "label_score").copy()
            top10["label"] = (
                top10["Codon"].astype(str) + " (" + top10["AA"].astype(str) + ")"
            )

            # Headline correlations (r/spearman/kendall) on the
            # transformed, primary-support subset.
            if len(x) >= 3:
                r_value = float(np.corrcoef(x, y)[0, 1])
                try:
                    sp_r, _ = spearmanr(x, y)
                    spearman_r = float(sp_r)
                except Exception:
                    spearman_r = float("nan")
                try:
                    kt, _ = kendalltau(x, y)
                    kendall_tau = float(kt)
                except Exception:
                    kendall_tau = float("nan")
            else:
                r_value = spearman_r = kendall_tau = float("nan")
            r_squared = r_value * r_value if not math.isnan(r_value) else float("nan")

            out_csv = os.path.join(output_dir, f"{base_sample}_vs_{sample_name}_{version}.csv")
            merged_current.to_csv(out_csv, index=False)
            log_info("COR", f"Wrote correlation CSV => {out_csv}")

            # Append to the long-format metrics table written once at
            # the end of the run.
            metrics_row = merged_current.copy()
            metrics_row.insert(0, "version", version)
            metrics_row.insert(0, "compare_sample", sample_name)
            metrics_row.insert(0, "base_sample", base_sample)
            metrics_row.insert(0, "site", site_label)
            metrics_records.append(metrics_row)

            # ----- Three-panel publication figure -----------------------
            #
            # Panel A: log2-density scatter with identity + robust fit
            #          and support-aware labels.
            # Panel B: MA / Bland-Altman plot — surfaces codon-specific
            #          shifts a correlation hides.
            # Panel C: residual plot — distribution of (sample - fit) so
            #          a scale or biased shift is obvious.
            # When metric='raw_count', panel A is rendered in the original
            # raw scale and the figure is written under raw_count_qc/ to
            # mark it as a QC artefact, not a publication figure.
            fallback_palette = list(sns.color_palette("Set2", 8).as_hex())
            categories = list(merged_current["Category"].astype(str).unique())
            colour_map = {
                cat: _categorical_colour(cat, fallback_palette) for cat in categories
            }

            primary_mask = merged_current["include_primary"].to_numpy()

            fig, (ax_scatter, ax_ma, ax_resid) = plt.subplots(
                1, 3, figsize=(18, 6)
            )

            x_plot = merged_current["base_metric"].to_numpy()
            y_plot = merged_current["sample_metric"].to_numpy()

            # Panel A — scatter
            for category, group in merged_current.groupby("Category", sort=False):
                ax_scatter.scatter(
                    group["base_metric"],
                    group["sample_metric"],
                    s=42,
                    alpha=np.where(group["include_primary"], 0.85, 0.30),
                    color=colour_map[str(category)],
                    edgecolor="white",
                    linewidth=0.6,
                    label=str(category),
                )
            data_min = float(min(x_plot.min(), y_plot.min()))
            data_max = float(max(x_plot.max(), y_plot.max()))
            pad = max((data_max - data_min) * 0.05, 1e-6)
            axis_min = data_min - pad
            axis_max = data_max + pad
            ax_scatter.set_xlim(axis_min, axis_max)
            ax_scatter.set_ylim(axis_min, axis_max)
            ax_scatter.set_aspect("equal", adjustable="box")
            line_x = np.linspace(axis_min, axis_max, 100)
            ax_scatter.plot(line_x, line_x, color="0.6", linestyle=":", linewidth=1.0, label="y = x")
            ax_scatter.plot(
                line_x,
                slope * line_x + intercept,
                color="#D55E00",
                linestyle="--",
                linewidth=1.4,
                label=f"{regression} (slope={slope:.3f})",
            )
            _place_residual_labels(
                ax_scatter, top10,
                x_col="base_metric",
                y_col="sample_metric",
                label_col="label",
                axis_min=axis_min,
                axis_max=axis_max,
            )
            stat_lines = [
                f"metric = {metric}",
                f"regression = {regression}",
                f"r (transformed) = {r_value:.3f}",
                f"Spearman rho = {spearman_r:.3f}",
                f"Kendall tau = {kendall_tau:.3f}",
                f"slope = {slope:.3f}",
                f"intercept = {intercept:.3f}",
                f"N primary = {int(primary_mask.sum())}",
                f"N total = {len(merged_current)}",
            ]
            ax_scatter.text(
                0.03, 0.97, "\n".join(stat_lines),
                transform=ax_scatter.transAxes, ha="left", va="top",
                fontsize=9, family="monospace",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                          edgecolor="0.7", linewidth=0.6, alpha=0.9),
            )
            scatter_xlabel = (
                f"{base_sample} ({metric})" if metric != "raw_count"
                else f"{base_sample} ({column})"
            )
            scatter_ylabel = (
                f"{sample_name} ({metric})" if metric != "raw_count"
                else f"{sample_name} ({column})"
            )
            ax_scatter.set_xlabel(scatter_xlabel, fontsize=11)
            ax_scatter.set_ylabel(scatter_ylabel, fontsize=11)
            ax_scatter.set_title("A. Codon density scatter", fontsize=12, fontweight="bold")
            ax_scatter.grid(True, which="major", alpha=0.25, linewidth=0.6)
            # Place the legend OUTSIDE the data axes (to the right of
            # panel A) so it cannot occlude data points or labels.
            # bbox_inches="tight" on savefig (below) keeps the figure
            # framing snug around the legend.
            ax_scatter.legend(
                title="Category", title_fontsize=10, fontsize=9,
                loc="upper left", bbox_to_anchor=(1.02, 1.0),
                borderaxespad=0.0, frameon=True, framealpha=0.9,
                edgecolor="0.7",
            )

            # Panel B — MA plot
            for category, group in merged_current.groupby("Category", sort=False):
                ax_ma.scatter(
                    group["mean_log2_density"],
                    group["log2_fold_change"],
                    s=36,
                    alpha=np.where(group["include_primary"], 0.85, 0.25),
                    color=colour_map[str(category)],
                    edgecolor="white",
                    linewidth=0.5,
                    label=str(category),
                )
            ax_ma.axhline(0.0, color="0.5", linestyle=":", linewidth=1.0)
            ax_ma.axhline(1.0, color="#D55E00", linestyle="--", linewidth=0.8, alpha=0.7)
            ax_ma.axhline(-1.0, color="#D55E00", linestyle="--", linewidth=0.8, alpha=0.7)
            ax_ma.set_xlabel("Mean log2 density (A)", fontsize=11)
            ax_ma.set_ylabel("log2 fold change (M)", fontsize=11)
            ax_ma.set_title("B. MA / Bland-Altman", fontsize=12, fontweight="bold")
            ax_ma.grid(True, which="major", alpha=0.25, linewidth=0.6)

            # Panel C — residuals
            ax_resid.scatter(
                merged_current["base_metric"],
                merged_current["robust_residual"],
                s=30,
                alpha=np.where(primary_mask, 0.85, 0.25),
                color="#0072B2",
                edgecolor="white",
                linewidth=0.4,
            )
            ax_resid.axhline(0.0, color="0.5", linestyle=":", linewidth=1.0)
            ax_resid.set_xlabel(scatter_xlabel, fontsize=11)
            ax_resid.set_ylabel("Residual (sample - fit)", fontsize=11)
            ax_resid.set_title("C. Robust residuals", fontsize=12, fontweight="bold")
            ax_resid.grid(True, which="major", alpha=0.25, linewidth=0.6)

            title_suffix = "all codons" if version == "all" else "outlier-masked"
            fig.suptitle(
                f"Codon density: {base_sample} vs {sample_name} ({title_suffix})",
                fontsize=14, fontweight="bold",
            )
            fig.tight_layout(rect=(0, 0, 1, 0.95))

            plot_dir = (
                Path(output_dir) / "raw_count_qc"
                if (raw_metric and str(raw_panel).lower() == "qc_only")
                else Path(output_dir)
            )
            plot_dir.mkdir(parents=True, exist_ok=True)
            out_svg = plot_dir / f"{base_sample}_vs_{sample_name}_{version}.svg"
            out_png = plot_dir / f"{base_sample}_vs_{sample_name}_{version}.png"
            fig.savefig(out_svg, bbox_inches="tight")
            fig.savefig(out_png, dpi=300, bbox_inches="tight")
            plt.close(fig)
            log_info("COR", f"Plot saved => {out_svg} (+ 300 dpi PNG)")

            corr_records.append({
                "Base_Sample": base_sample,
                "Other_Sample": sample_name,
                "Site": site_label,
                "Version": version,
                "Metric": metric,
                "Regression": regression,
                "Pseudocount": float(pc),
                "Pearson_r_metric": r_value,
                "Spearman_r": spearman_r,
                "Kendall_tau": kendall_tau,
                "RobustSlope": slope,
                "RobustIntercept": intercept,
                "NumCodons": len(merged_current),
                "NumCodonsPrimary": int(primary_mask.sum()),
                "PseudocountResolved": float(pc),
            })

    if metrics_records:
        metrics_df = pd.concat(metrics_records, ignore_index=True)
        # Friendly column ordering for the long-format TSV.
        ordered = [
            "site", "base_sample", "compare_sample", "version",
            "Codon", "AA", "Category",
            "base_val", "sample_val",
            "base_metric", "sample_metric",
            "log2_fold_change", "mean_log2_density",
            "support_min_raw", "support_total_raw",
            "predicted", "robust_residual", "residual",
            "label_score", "include_primary", "exclusion_reason",
        ]
        ordered = [c for c in ordered if c in metrics_df.columns]
        metrics_df = metrics_df[ordered + [c for c in metrics_df.columns if c not in ordered]]
        metrics_path = os.path.join(output_dir, "codon_correlation_metrics.tsv")
        metrics_df.to_csv(metrics_path, sep="\t", index=False)
        log_info("COR", f"Wrote codon correlation metrics => {metrics_path}")

    cor_df = pd.DataFrame(corr_records)
    out_summary = os.path.join(output_dir, f"codon_correlation_summary_{base_sample}.csv")
    cor_df.to_csv(out_summary, index=False)
    log_info("COR", f"Correlation summary saved => {out_summary}")

    metadata_payload = {
        "metric": metric,
        "regression": regression,
        "support_min_raw": int(support_min_raw),
        "label_top_n": int(label_top_n),
        "raw_panel": str(raw_panel),
        "site": site_label,
        "base_sample": base_sample,
        "compare_samples": [s for s in samples if s != base_sample],
        "primary_column": column,
        "warnings": (
            [W_CODON_RAW_COUNT_PRIMARY] if raw_metric else []
        ),
    }
    metadata_path = os.path.join(output_dir, "codon_correlation.metadata.json")
    Path(metadata_path).write_text(
        json.dumps(metadata_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    log_info("COR", "Done correlation for all samples.")

# Example usage:
# run_codon_correlation(translation_profile_dir="path/to/analysis_results",
#                       samples=["sample1", "sample2", "sample3"],
#                       base_sample="sample1",
#                       column="CoverageDivFreq")
