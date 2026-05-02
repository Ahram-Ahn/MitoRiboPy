"""Per-plot ``.metadata.json`` sidecar coverage for the three plot families
that previously triggered ``FIGURE_QC_WARN`` because no sidecar was
written next to the rendered PNG.

The contract is tied to :mod:`mitoribopy.plotting.figure_validator`:
``mitoribopy validate-figures`` looks up ``<plot>.metadata.json`` via
:func:`metadata_sidecar_path`, so each rendered plot must produce one.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # Headless backend for the in-memory plot fixtures.

import pandas as pd

from mitoribopy.analysis import codon_correlation as cc
from mitoribopy.plotting.figure_validator import metadata_sidecar_path
from mitoribopy.plotting.translation_profile_plots import (
    plot_codon_usage_dataframe,
    plot_site_depth_profile,
)


_REQUIRED_KEYS = (
    "plot_type",
    "stage",
    "source_data",
    "n_points_expected",
    "n_points_drawn",
)


def _assert_required_keys(payload: dict) -> None:
    for key in _REQUIRED_KEYS:
        assert key in payload, f"sidecar missing required key {key!r}: {payload}"


def _read_sidecar(plot_path: Path) -> dict:
    sidecar = metadata_sidecar_path(plot_path)
    assert sidecar.is_file(), f"missing sidecar for {plot_path}"
    return json.loads(sidecar.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# 1) translation_profile/<sample>/footprint_density/<gene>_p_site_depth.png
# ---------------------------------------------------------------------------


def test_plot_site_depth_profile_writes_sidecar(tmp_path: Path) -> None:
    foot_dir = tmp_path / "rpf" / "translation_profile" / "WT" / "footprint_density"
    foot_dir.mkdir(parents=True)
    df = pd.DataFrame(
        {
            "Position": list(range(1, 11)),
            "Nucleotide": list("ACGTACGTAC"),
            "_clipped_depth": [0, 1, 2, 3, 0, 0, 4, 1, 0, 2],
        }
    )
    out_png = foot_dir / "MT-CO1_p_site_depth.png"
    plot_site_depth_profile(
        df,
        str(out_png),
        site_column="_clipped_depth",
        site_label="P-site",
        sample_name="WT",
        transcript_name="MT-CO1",
        color="forestgreen",
        source_data="MT-CO1_footprint_density.csv",
    )

    payload = _read_sidecar(out_png)
    _assert_required_keys(payload)
    assert payload["plot_type"] == "footprint_density_depth"
    assert payload["stage"] == "rpf"
    assert payload["source_data"] == "MT-CO1_footprint_density.csv"
    assert payload["n_points_expected"] == 10
    assert payload["n_points_drawn"] == 10


# ---------------------------------------------------------------------------
# 2) translation_profile/<sample>/codon_usage/p_site_codon_usage_<gene>_plot.png
# ---------------------------------------------------------------------------


def test_plot_codon_usage_writes_sidecar(tmp_path: Path) -> None:
    codon_dir = tmp_path / "rpf" / "translation_profile" / "WT" / "codon_usage"
    codon_dir.mkdir(parents=True)
    df = pd.DataFrame(
        {
            "Codon": ["AAA", "AAC", "AGG"],
            "AA": ["K", "N", "R"],
            "Category": ["basic", "polar", "basic"],
            "CoverageDivFreq": [10.0, 20.0, 30.0],
        }
    )
    codon_label_order = ["AAA-K", "AAC-N", "AGG-R", "CTG-L"]
    out_png = codon_dir / "p_site_codon_usage_MT-CO1_plot.png"
    plot_codon_usage_dataframe(
        df,
        str(out_png),
        title="MT-CO1 - P-site Codon Usage - WT",
        codon_label_order=codon_label_order,
        source_data="p_site_codon_usage_MT-CO1.csv",
    )

    payload = _read_sidecar(out_png)
    _assert_required_keys(payload)
    assert payload["plot_type"] == "codon_usage_bar"
    assert payload["stage"] == "rpf"
    assert payload["source_data"] == "p_site_codon_usage_MT-CO1.csv"
    # 4 labels in order, 3 of them appear in the data.
    assert payload["n_points_expected"] == 4
    assert payload["n_points_drawn"] == 3


# ---------------------------------------------------------------------------
# 3) codon_correlation/<base>_vs_<sample>_<version>.png — per-plot sidecar
#    must coexist with the folder-level codon_correlation.metadata.json.
# ---------------------------------------------------------------------------


def _make_codon_csv(path: Path, codons: list[tuple[str, str, str, float]]) -> None:
    df = pd.DataFrame(codons, columns=["Codon", "AA", "Category", "CoverageDivFreq"])
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def test_codon_correlation_writes_per_plot_sidecar(tmp_path: Path) -> None:
    base = tmp_path / "translation_profile"
    codons = [
        ("AAA", "K", "basic", 100.0),
        ("AAC", "N", "polar", 80.0),
        ("AGG", "R", "basic", 60.0),
        ("CTG", "L", "non-polar", 40.0),
        ("CTC", "L", "non-polar", 30.0),
        ("ATG", "M", "start", 20.0),
        ("TGA", "*", "stop", 5.0),
    ]
    perturbed = [
        ("AAA", "K", "basic", 25.0),
        ("AAC", "N", "polar", 80.0),
        ("AGG", "R", "basic", 60.0),
        ("CTG", "L", "non-polar", 40.0),
        ("CTC", "L", "non-polar", 30.0),
        ("ATG", "M", "start", 20.0),
        ("TGA", "*", "stop", 5.0),
    ]
    for sample, vals in (("WT", codons), ("KO", perturbed)):
        _make_codon_csv(
            base / sample / "codon_usage" / "p_site_codon_usage_total.csv",
            vals,
        )
        _make_codon_csv(
            base / sample / "codon_usage" / "a_site_codon_usage_total.csv",
            vals,
        )

    out_dir = tmp_path / "cor_out"
    cc.run_codon_correlation(
        translation_profile_dir=str(base),
        samples=["WT", "KO"],
        base_sample="WT",
        output_dir=str(out_dir),
        site="p",
        mask_method="none",
    )

    # Folder-level metadata must still exist (downstream tooling
    # references it) — the per-plot sidecar coexists with it.
    folder_metadata = out_dir / "codon_correlation.metadata.json"
    assert folder_metadata.is_file()

    plot_png = out_dir / "WT_vs_KO_all.png"
    assert plot_png.is_file()
    payload = _read_sidecar(plot_png)
    _assert_required_keys(payload)
    assert payload["plot_type"] == "codon_correlation_scatter"
    assert payload["stage"] == "rpf"
    assert payload["source_data"] == "WT_vs_KO_all.csv"
    # The scatter plots one point per merged codon row; the seven test
    # codons all merge cleanly so expected == drawn == 7.
    assert payload["n_points_expected"] == payload["n_points_drawn"]
    assert payload["n_points_drawn"] >= 1
