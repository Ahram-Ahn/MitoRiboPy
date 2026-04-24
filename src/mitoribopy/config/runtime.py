"""Runtime pipeline configuration helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ..console import log_warning


DEFAULT_CONFIG: dict[str, Any] = {
    "strain": "y",
    "directory": ".",
    "rpf": None,
    "annotation_file": None,
    "codon_tables_file": None,
    "codon_table_name": None,
    "start_codons": None,
    "atp8_atp6_baseline": "ATP6",
    "nd4l_nd4_baseline": "ND4",
    "align": "start",
    "range": 20,
    "output": "analysis_results",
    "downstream_dir": "footprint_density",
    "min_offset": 11,
    "max_offset": 20,
    "min_5_offset": None,
    "max_5_offset": None,
    "min_3_offset": None,
    "max_3_offset": None,
    "offset_mask_nt": 5,
    "offset_pick_reference": "p_site",
    "offset_type": "5",
    "offset_site": "p",
    "offset_mode": "per_sample",
    "analysis_sites": "both",
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
    "unfiltered_read_length_range": [15, 50],
    "rpm_norm_mode": "total",
    "read_counts_reference_col": None,
    "mrna_ref_patterns": ["mt_genome", "mt-mrna", "mt_mrna"],
    "structure_density": False,
    "structure_density_norm_perc": 0.99,
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
    "bam_mapq": 10,
    "footprint_class": "monosome",
}


# Biological RPF windows per footprint class.
# monosome: a single ribosome footprint, ~28-34 nt in metazoan mt-Ribo-seq
#           (mt-monosomes are short compared to cytoplasmic ~28-32 nt; yeast
#           mt is longer at 37-41 nt because of the extended mS37 tail).
# disome:   two stacked ribosomes with a shared protected fragment, ~60-90 nt
#           (collided-ribosome studies, e.g. eIF5A-depletion / queueing).
FOOTPRINT_CLASS_DEFAULTS: dict[str, dict[str, tuple[int, int]]] = {
    "monosome": {
        "unfiltered_read_length_range": (15, 50),
        "rpf_h": (28, 34),   # vertebrate mt-mono footprint
        "rpf_y": (37, 41),   # yeast mt-mono footprint
    },
    "disome": {
        "unfiltered_read_length_range": (40, 110),
        "rpf_h": (60, 90),   # vertebrate mt-disome footprint (loose)
        "rpf_y": (65, 95),   # yeast mt-disome footprint
    },
    "custom": {
        # "custom" implies the user will supply --rpf and
        # --unfiltered_read_length_range explicitly; no biological default.
    },
}


def load_user_config(config_path: str | None) -> dict[str, Any]:
    """Load a JSON, YAML, or TOML config file and keep only recognized keys.

    The file format is auto-detected from the path suffix:

    * ``.json`` (or unknown extension) - stdlib :mod:`json`
    * ``.yaml`` / ``.yml``              - :mod:`yaml` (PyYAML, required dep)
    * ``.toml``                         - stdlib :mod:`tomllib` (Py 3.11+)
                                          or the optional :mod:`tomli` fallback

    Unknown keys are reported as a ``CONFIG`` warning and ignored.
    """
    if not config_path:
        return {}

    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    suffix = path.suffix.lower()
    text = path.read_text(encoding="utf-8")

    if suffix in {".yaml", ".yml"}:
        try:
            import yaml
        except ImportError as exc:  # pragma: no cover - PyYAML is a core dep
            raise RuntimeError(
                "PyYAML is required to read YAML config files. "
                "Install with: pip install PyYAML"
            ) from exc
        raw = yaml.safe_load(text) or {}
    elif suffix == ".toml":
        try:
            import tomllib  # Python 3.11+
        except ImportError:  # pragma: no cover - Python 3.10 fallback
            try:
                import tomli as tomllib  # type: ignore[no-redef]
            except ImportError as exc:
                raise RuntimeError(
                    "Reading TOML config files on Python < 3.11 requires the "
                    "'tomli' package. Install with: pip install tomli"
                ) from exc
        raw = tomllib.loads(text)
    else:
        raw = json.loads(text) if text.strip() else {}

    if not isinstance(raw, dict):
        raise ValueError(
            "Config file must parse to a mapping/object; got "
            f"{type(raw).__name__}."
        )

    allowed = set(DEFAULT_CONFIG.keys())
    unknown = sorted(key for key in raw.keys() if key not in allowed)
    if unknown:
        log_warning("CONFIG", f"Ignoring unknown keys: {', '.join(unknown)}")

    return {key: value for key, value in raw.items() if key in allowed}


def resolve_rpf_range(
    strain: str,
    rpf_arg: list[int] | tuple[int, int] | None,
    *,
    footprint_class: str = "monosome",
) -> list[int]:
    """Resolve RPF range from CLI/config or fallback defaults.

    Explicit ``rpf_arg`` always wins. Otherwise the default comes from
    ``FOOTPRINT_CLASS_DEFAULTS[footprint_class]`` keyed by strain:

    * ``strain='h'`` or ``'vm'`` -> ``rpf_h``
    * ``strain='y'`` or ``'ym'`` -> ``rpf_y``

    For ``footprint_class='custom'`` or any strain without a built-in
    (including ``'custom'`` strain), ``rpf_arg`` is required; we raise
    a clear ``ValueError`` otherwise so the caller can surface it as a
    user-facing CLI error.
    """
    if rpf_arg:
        start, end = int(rpf_arg[0]), int(rpf_arg[1])
        if end < start:
            raise ValueError(f"Invalid RPF range: start={start}, end={end}")
        return list(range(start, end + 1))

    class_defaults = FOOTPRINT_CLASS_DEFAULTS.get(footprint_class, {})
    if strain in {"h", "vm"}:
        window = class_defaults.get("rpf_h")
    elif strain in {"y", "ym"}:
        window = class_defaults.get("rpf_y")
    else:
        window = None

    if window is None:
        raise ValueError(
            f"Cannot default RPF range for strain={strain!r}, "
            f"footprint_class={footprint_class!r}. Pass -rpf MIN MAX explicitly."
        )
    return list(range(window[0], window[1] + 1))
