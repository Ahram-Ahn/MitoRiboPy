"""Analysis modules for MitoRiboPy."""

from importlib import import_module


_EXPORTS = {
    "build_per_sample_summaries": ("offset_enrichment", "build_per_sample_summaries"),
    "run_codon_correlation": ("codon_correlation", "run_codon_correlation"),
    "run_translation_profile_analysis": (
        "translation_profile_analysis",
        "run_translation_profile_analysis",
    ),
    "compute_offsets": ("offset_enrichment", "compute_offsets"),
    "create_csv_for_offset_enrichment": (
        "offset_enrichment",
        "create_csv_for_offset_enrichment",
    ),
    "plot_offset_enrichment": ("offset_enrichment", "plot_offset_enrichment"),
    "determine_p_site_offsets": ("offset_selection", "determine_p_site_offsets"),
}


def __getattr__(name: str):
    """Resolve back-compat analysis re-exports without eager plot imports."""
    try:
        module_name, attr_name = _EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(
            f"module 'mitoribopy.analysis' has no attribute {name!r}"
        ) from exc
    module = import_module(f"{__name__}.{module_name}")
    value = getattr(module, attr_name)
    globals()[name] = value
    return value


__all__ = [
    "build_per_sample_summaries",
    "run_codon_correlation",
    "run_translation_profile_analysis",
    "compute_offsets",
    "create_csv_for_offset_enrichment",
    "plot_offset_enrichment",
    "determine_p_site_offsets",
]
