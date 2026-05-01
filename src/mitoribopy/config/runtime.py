"""Runtime pipeline configuration helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ..console import log_warning


DEFAULT_CONFIG: dict[str, Any] = {
    "strain": "h.sapiens",
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
    "rpf_min_count_frac": 0.20,
    "min_5_offset": None,
    "max_5_offset": None,
    "min_3_offset": None,
    "max_3_offset": None,
    "offset_mask_nt": 5,
    # 'p_site' (default): pick the offset in canonical P-site space then
    # convert into the --offset_site space for reporting. 'reported_site'
    # (formerly 'selected_site'): pick directly in --offset_site space.
    "offset_pick_reference": "p_site",
    "offset_type": "5",
    "offset_site": "p",
    "offset_mode": "per_sample",
    "analysis_sites": "both",
    "codon_overlap_mode": "full",
    "plot_dir": "offset_diagnostics",
    "plot_format": "png",
    "x_breaks": None,
    "line_plot_style": "combined",
    "cap_percentile": 0.999,
    # When True, codon-density values are smoothed with a +/-1 nt
    # window around the codon centre. Old name 'merge_density' is
    # accepted via _DEPRECATED_CONFIG_KEY_ALIASES in pipeline/runner.py.
    "codon_density_window": False,
    "psite_offset": None,
    "read_counts_file": "read_counts_summary.txt",
    "read_counts_sample_col": None,
    "read_counts_reads_col": None,
    "unfiltered_read_length_range": [15, 50],
    "rpm_norm_mode": "total",
    "read_counts_reference_col": None,
    # Substring patterns identifying mt-mRNA rows in the read-count
    # table when rpm_norm_mode == 'mt_mrna'. Old name 'mrna_ref_patterns'
    # is accepted via _DEPRECATED_CONFIG_KEY_ALIASES in pipeline/runner.py.
    "mt_mrna_substring_patterns": ["mt_genome", "mt-mrna", "mt_mrna"],
    "structure_density": False,
    "structure_density_norm_perc": 0.99,
    "order_samples": None,
    "cor_plot": False,
    "base_sample": None,
    "cor_mask_method": "percentile",
    "cor_mask_percentile": 0.99,
    "cor_mask_threshold": None,
    "cor_metric": "log2_density_rpm",
    "cor_regression": "theil_sen",
    "cor_support_min_raw": 10,
    "cor_label_top_n": 10,
    "cor_pseudocount": "auto",
    "cor_raw_panel": "qc_only",
    "read_coverage_raw": True,
    "read_coverage_rpm": True,
    "igv_export": False,
    "bam_mapq": 10,
    "footprint_class": "monosome",
}


# Biological RPF windows per footprint class.
#
# short:    truncated RNase products (~16-24 nt). These are partial
#           protections produced when the RNase digest cuts inside the
#           ribosome footprint. They sit just below the canonical
#           monosome window and are useful for mapping context-dependent
#           pausing and for QC of digest aggressiveness. Same range
#           across species (RNase truncation products do not strongly
#           vary by mt-genome lineage).
# monosome: a single ribosome footprint, ~28-34 nt in metazoan mt-Ribo-seq
#           (mt-monosomes are short compared to cytoplasmic ~28-32 nt; yeast
#           mt is longer at 37-41 nt because of the extended mS37 tail).
# disome:   two stacked ribosomes with a shared protected fragment.
#           Vertebrate mt: 50-70 nt; yeast mt: 60-90 nt.
#           Used for collided-ribosome studies (eIF5A-depletion,
#           queueing, ribosome-stalling).
# custom:   the user will supply --rpf and --unfiltered_read_length_range
#           explicitly; no biological default.
FOOTPRINT_CLASS_DEFAULTS: dict[str, dict[str, tuple[int, int]]] = {
    "short": {
        "unfiltered_read_length_range": (10, 30),
        "rpf_h.sapiens": (16, 24),
        "rpf_s.cerevisiae": (16, 24),
    },
    "monosome": {
        "unfiltered_read_length_range": (15, 50),
        "rpf_h.sapiens": (28, 34),
        "rpf_s.cerevisiae": (37, 41),
    },
    "disome": {
        "unfiltered_read_length_range": (40, 100),
        "rpf_h.sapiens": (50, 70),
        "rpf_s.cerevisiae": (60, 90),
    },
    "custom": {},
}


DEPRECATED_CONFIG_KEY_ALIASES: dict[str, str] = {
    # legacy YAML key            -> canonical key
    "merge_density": "codon_density_window",
    "mrna_ref_patterns": "mt_mrna_substring_patterns",
}


def load_user_config(config_path: str | None) -> dict[str, Any]:
    """Load a JSON, YAML, or TOML config file and keep only recognized keys.

    The file format is auto-detected from the path suffix:

    * ``.json`` (or unknown extension) - stdlib :mod:`json`
    * ``.yaml`` / ``.yml``              - :mod:`yaml` (PyYAML, required dep)
    * ``.toml``                         - stdlib :mod:`tomllib` (Py 3.11+)
                                          or the optional :mod:`tomli` fallback

    Legacy keys listed in :data:`DEPRECATED_CONFIG_KEY_ALIASES` are
    silently rewritten to their canonical names before the unknown-key
    check (the runner-level deprecation warning is emitted after this
    pass so the user sees one line per renamed key, not one line per
    pipeline invocation).

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

    # Rewrite deprecated keys before the unknown-key check so legacy
    # configs continue to load cleanly.
    canonicalized: dict[str, Any] = {}
    for key, value in raw.items():
        canonical = DEPRECATED_CONFIG_KEY_ALIASES.get(key, key)
        if canonical in canonicalized:
            # Both the legacy and canonical key were set; the canonical
            # value already loaded wins.
            continue
        canonicalized[canonical] = value

    allowed = set(DEFAULT_CONFIG.keys())
    unknown = sorted(key for key in canonicalized.keys() if key not in allowed)
    if unknown:
        log_warning("CONFIG", f"Ignoring unknown keys: {', '.join(unknown)}")

    return {key: value for key, value in canonicalized.items() if key in allowed}


def resolve_rpf_range(
    strain: str,
    rpf_arg: list[int] | tuple[int, int] | None,
    *,
    footprint_class: str = "monosome",
) -> list[int]:
    """Resolve the RPF length window from the CLI / config or a default.

    Explicit ``rpf_arg`` (a 2-tuple [min, max]) always wins. Otherwise
    the default comes from
    ``FOOTPRINT_CLASS_DEFAULTS[footprint_class]['rpf_<strain>']``. The
    only built-in defaults are for ``h.sapiens`` and ``s.cerevisiae``;
    any other strain (including ``'custom'``) MUST supply
    ``-rpf MIN MAX`` explicitly. Unknown strain or unknown footprint
    class raises ``ValueError`` so the caller can surface a clear CLI
    message.
    """
    # Local import keeps the package layering clean — config/ should
    # not pull in data/ at module load time.
    from ..data.reference_data import canonical_strain

    if rpf_arg:
        start, end = int(rpf_arg[0]), int(rpf_arg[1])
        if end < start:
            raise ValueError(f"Invalid RPF range: start={start}, end={end}")
        return list(range(start, end + 1))

    class_defaults = FOOTPRINT_CLASS_DEFAULTS.get(footprint_class, {})
    canonical = canonical_strain(strain)
    window = class_defaults.get(f"rpf_{canonical}")

    if window is None:
        raise ValueError(
            f"Cannot default the RPF range for strain={strain!r}, "
            f"footprint_class={footprint_class!r}. Either choose a "
            "footprint_class with a built-in default for this strain "
            "(short / monosome / disome work for h.sapiens and "
            "s.cerevisiae) or pass -rpf MIN MAX explicitly."
        )
    return list(range(window[0], window[1] + 1))
