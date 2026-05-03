"""Tests for the v0.7.0 Theil-Sen slope CI surfacing.

The slope CI was always computed (scipy's ``theilslopes`` returns the
``lo_slope`` / ``hi_slope`` bounds by default), but pre-v0.7.0 those
values were not propagated to the user-facing surfaces — the figure
annotation showed only the point estimate, and the codon-correlation
summary CSV / per-plot sidecar omitted them. This test locks the new
behaviour in.
"""

from __future__ import annotations

import math

import numpy as np

from mitoribopy.analysis.codon_correlation import _fit_regression


def test_theil_sen_returns_slope_ci_columns() -> None:
    rng = np.random.default_rng(0)
    x = np.linspace(0, 10, 30)
    y = 2.0 * x + 1.0 + rng.normal(size=x.size, scale=0.5)
    fit = _fit_regression(x, y, method="theil_sen")
    assert fit["method"] == "theil_sen"
    # Slope point estimate hugs the true slope of 2.0 within noise.
    assert abs(fit["slope"] - 2.0) < 0.3
    # CI keys are present and bracket the point estimate.
    assert "slope_ci_low" in fit
    assert "slope_ci_high" in fit
    assert fit["slope_ci_low"] <= fit["slope"] <= fit["slope_ci_high"]
    # CI is non-degenerate for a noisy fit.
    assert fit["slope_ci_high"] > fit["slope_ci_low"]


def test_theil_sen_ci_is_finite_for_clean_signal() -> None:
    """A perfect linear signal should still produce finite (possibly
    very tight) CI bounds — not NaN, not Inf."""
    x = np.linspace(0, 5, 20)
    y = 3.0 * x + 2.0
    fit = _fit_regression(x, y, method="theil_sen")
    assert math.isfinite(fit["slope_ci_low"])
    assert math.isfinite(fit["slope_ci_high"])
    assert fit["slope_ci_low"] <= 3.0 <= fit["slope_ci_high"]


def test_ols_does_not_emit_slope_ci_columns() -> None:
    """OLS doesn't carry a CI from linregress in this code path; the
    fields must be ABSENT (so the call site can detect via .get and
    write NaN to the sidecar)."""
    x = np.linspace(0, 10, 20)
    y = 2.0 * x + 1.0
    fit = _fit_regression(x, y, method="ols")
    assert fit["method"] == "ols"
    assert "slope_ci_low" not in fit
    assert "slope_ci_high" not in fit


def test_too_few_points_returns_identity_with_no_ci() -> None:
    fit = _fit_regression(np.array([1.0]), np.array([1.0]), method="theil_sen")
    assert fit["method"].startswith("theil_sen_too_few_points") \
        or fit["method"] == "identity"
    assert "slope_ci_low" not in fit
    assert "slope_ci_high" not in fit
