"""Runtime pipeline configuration helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


DEFAULT_CONFIG: dict[str, Any] = {
    "strain": "y",
    "directory": ".",
    "rpf": None,
    "align": "start",
    "range": 20,
    "output": "analysis_results",
    "downstream_dir": "footprint_density",
    "min_offset": 11,
    "max_offset": 20,
    "offset_pick_reference": "p_site",
    "offset_type": "5",
    "offset_site": "p",
    "codon_overlap_mode": "full",
    "plot_dir": "plots_and_csv",
    "plot_format": "png",
    "x_breaks": None,
    "line_plot_style": "combined",
    "cap_percentile": 0.999,
    "merge_density": False,
    "psite_offset": None,
    "read_counts_file": "read_counts_summary.txt",
    "read_counts_sample_col": None,
    "read_counts_reads_col": None,
    "rpm_norm_mode": "total",
    "read_counts_reference_col": None,
    "mrna_ref_patterns": ["mt_genome", "mt-mrna", "mt_mrna"],
    "varna": False,
    "varna_norm_perc": 0.99,
    "order_samples": None,
    "cor_plot": False,
    "base_sample": None,
    "cor_mask_method": "percentile",
    "cor_mask_percentile": 0.99,
    "cor_mask_threshold": None,
    "use_rna_seq": False,
    "rna_seq_dir": None,
    "rna_order": None,
    "rna_out_dir": "rna_seq_results",
    "do_rna_ribo_ratio": False,
}


def load_user_config(config_path: str | None) -> dict[str, Any]:
    """Load JSON config and keep only recognized keys."""
    if not config_path:
        return {}

    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    with path.open("r", encoding="utf-8") as config_file:
        raw = json.load(config_file)

    if not isinstance(raw, dict):
        raise ValueError("Config file must be a JSON object (key/value dictionary).")

    allowed = set(DEFAULT_CONFIG.keys())
    unknown = sorted(key for key in raw.keys() if key not in allowed)
    if unknown:
        print(f"[CONFIG] Ignoring unknown keys: {', '.join(unknown)}")

    return {key: value for key, value in raw.items() if key in allowed}


def resolve_rpf_range(strain: str, rpf_arg: list[int] | tuple[int, int] | None) -> list[int]:
    """Resolve RPF range from CLI/config or fallback defaults by strain."""
    if rpf_arg:
        start, end = int(rpf_arg[0]), int(rpf_arg[1])
        if end < start:
            raise ValueError(f"Invalid RPF range: start={start}, end={end}")
        return list(range(start, end + 1))

    if strain == "y":
        return list(range(37, 42))
    return list(range(28, 35))

