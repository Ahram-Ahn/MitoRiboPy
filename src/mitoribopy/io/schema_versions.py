"""Schema version registry for MitoRiboPy output TSVs.

Each major output table is associated with a semantic version. The
version is BOTH:

* recorded in the ``output_schemas`` map of ``run_manifest.json`` so a
  downstream script can fail fast if it expects a different layout, and
* prepended as a comment line ``# schema_version: X.Y`` to the TSV
  itself so consumers that grep the file (or open it without manifest
  context) can detect drift without having to parse the manifest.

Bump rules:

* MINOR (1.0 -> 1.1): added a new column at the END.
* MAJOR (1.x -> 2.0): renamed / removed / reordered an existing column.

Standard pandas / csv readers ignore ``#``-prefixed lines when
``comment="#"`` is passed; the TSVs remain backwards-compatible with
readers that do not pass that argument because the column header is
the FIRST non-comment line and stdlib csv stops at EOF.
"""

from __future__ import annotations


__all__ = [
    "OUTPUT_SCHEMA_VERSIONS",
    "schema_header_line",
]


OUTPUT_SCHEMA_VERSIONS: dict[str, str] = {
    "read_counts.tsv":    "1.0",
    "kit_resolution.tsv": "1.2",  # + best_adapter, best/second_best_match_rate, confidence_margin
    "umi_qc.tsv":         "1.0",  # per-sample UMI / dedup audit
    "rpf_counts.tsv":     "1.0",
    "rpf_counts_matrix.tsv": "1.0",
    "rna_counts.tsv":     "1.0",
    "te.tsv":             "2.0",  # + condition, assay, log2_te, note
    "delta_te.tsv":       "2.0",  # + base/compare_condition, padj_*, method, note
    "summary_qc.tsv":     "1.0",
    "warnings.tsv":       "2.0",  # column rename to publication-readiness spec
    "outputs_index.tsv":  "1.0",
    "figure_qc.tsv":      "1.0",  # per-plot mechanical QC (overlap / counts / SVG)
    # Wakigawa metagene Fourier periodicity bundle (v0.8.0).
    "fourier_spectrum_combined.tsv":      "1.0",
    "fourier_period3_score_combined.tsv": "1.0",
}


def schema_header_line(name: str) -> str:
    """Return the ``# schema_version: X.Y\\n`` line for the given file name.

    Raises ``KeyError`` for an unregistered name so adding a new TSV
    output requires updating :data:`OUTPUT_SCHEMA_VERSIONS` first.
    """
    return f"# schema_version: {OUTPUT_SCHEMA_VERSIONS[name]}\n"
