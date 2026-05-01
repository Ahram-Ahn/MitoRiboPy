"""Generate the top-level ``outputs_index.tsv`` for a MitoRiboPy run.

Section 8 of the assessment asks every run root to advertise every
TSV / JSON / plot it produced through one machine-readable index file
so a reviewer (or a downstream script) does not have to hunt through
``align/``, ``rpf/``, and ``rnaseq/`` to find the canonical path of
each output. The columns are fixed:

    output_type, stage, path, description, recommended_for, schema_version

Where:

* ``output_type``     — short kebab-case key (``read_counts``, ``rpf_counts``,
                        ``te_table``, ``delta_te_table``, ``offset_diagnostics``, ...)
* ``stage``           — ``align`` / ``rpf`` / ``rnaseq`` / ``all``
* ``path``            — path RELATIVE to the run root, posix-style
* ``description``     — one-sentence human description
* ``recommended_for`` — short audience hint (``downstream-scripting``,
                        ``reviewer-spot-check``, ``debugging``)
* ``schema_version``  — value from :data:`mitoribopy.io.schema_versions.OUTPUT_SCHEMA_VERSIONS`,
                        empty for non-TSV outputs (plots, JSON sidecars).

The index is *advisory*: missing entries do not fail a run. The
mapping is built by :func:`build_outputs_index_rows` from disk state,
so adding a new TSV in a stage will pick up automatically as long as
it is registered in :data:`_KNOWN_OUTPUTS` below.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .schema_versions import OUTPUT_SCHEMA_VERSIONS, schema_header_line


__all__ = [
    "OutputDescriptor",
    "build_outputs_index_rows",
    "write_outputs_index",
]


_OUTPUTS_INDEX_COLUMNS: tuple[str, ...] = (
    "output_type",
    "stage",
    "path",
    "description",
    "recommended_for",
    "schema_version",
)


@dataclass(frozen=True)
class OutputDescriptor:
    output_type: str
    stage: str
    relative_path: str
    description: str
    recommended_for: str
    schema_key: str | None = None  # key into OUTPUT_SCHEMA_VERSIONS, if any


# Registry of known outputs. Order is the order they appear in the TSV.
# Keep entries grouped by stage and roughly ordered "primary first,
# diagnostics last" so a reviewer reading the file top-to-bottom gets
# the publication-relevant tables first.
_KNOWN_OUTPUTS: tuple[OutputDescriptor, ...] = (
    # align ---------------------------------------------------------------
    OutputDescriptor(
        output_type="read_counts",
        stage="align",
        relative_path="align/read_counts.tsv",
        description=(
            "Per-sample read funnel from raw FASTQ through trim, "
            "rRNA / contaminant filter, mt-transcriptome alignment, "
            "MAPQ filter, and dedup."
        ),
        recommended_for="downstream-scripting",
        schema_key="read_counts.tsv",
    ),
    OutputDescriptor(
        output_type="kit_resolution",
        stage="align",
        relative_path="align/kit_resolution.tsv",
        description=(
            "Per-sample resolution of library kit, adapter, UMI, and "
            "dedup strategy with detection match rates."
        ),
        recommended_for="reviewer-spot-check",
        schema_key="kit_resolution.tsv",
    ),
    OutputDescriptor(
        output_type="align_run_settings",
        stage="align",
        relative_path="align/run_settings.json",
        description="Resolved align-stage configuration and tool versions.",
        recommended_for="reviewer-spot-check",
    ),
    OutputDescriptor(
        output_type="bed_dir",
        stage="align",
        relative_path="align/bed/",
        description=(
            "Per-sample BED6 files of mt-transcriptome alignments "
            "(post-dedup, post-MAPQ); RPF stage input."
        ),
        recommended_for="downstream-scripting",
    ),
    # rpf -----------------------------------------------------------------
    OutputDescriptor(
        output_type="rpf_counts",
        stage="rpf",
        relative_path="rpf/rpf_counts.tsv",
        description=(
            "RPF count matrix per (sample, gene) after RPF length "
            "filter and offset assignment."
        ),
        recommended_for="downstream-scripting",
        schema_key="rpf_counts.tsv",
    ),
    OutputDescriptor(
        output_type="rpf_counts_metadata",
        stage="rpf",
        relative_path="rpf/rpf_counts.metadata.json",
        description=(
            "Sidecar JSON for rpf_counts.tsv recording filter window "
            "and offset selection provenance."
        ),
        recommended_for="reviewer-spot-check",
    ),
    OutputDescriptor(
        output_type="offset_diagnostics",
        stage="rpf",
        relative_path="rpf/offset_diagnostics/",
        description="Per-sample P-site offset selection diagnostic plots.",
        recommended_for="reviewer-spot-check",
    ),
    OutputDescriptor(
        output_type="offset_applied",
        stage="rpf",
        relative_path="rpf/offset_applied.csv",
        description="Final per-sample, per-length offset table actually applied.",
        recommended_for="downstream-scripting",
    ),
    OutputDescriptor(
        output_type="rpf_run_settings",
        stage="rpf",
        relative_path="rpf/run_settings.json",
        description="Resolved rpf-stage configuration and tool versions.",
        recommended_for="reviewer-spot-check",
    ),
    # rnaseq --------------------------------------------------------------
    OutputDescriptor(
        output_type="te_table",
        stage="rnaseq",
        relative_path="rnaseq/te.tsv",
        description="Per-sample / per-gene translation efficiency table.",
        recommended_for="downstream-scripting",
        schema_key="te.tsv",
    ),
    OutputDescriptor(
        output_type="delta_te_table",
        stage="rnaseq",
        relative_path="rnaseq/delta_te.tsv",
        description=(
            "ΔTE (compare vs base) per gene with mRNA, RPF, and TE "
            "log2 fold-changes and adjusted p-values."
        ),
        recommended_for="downstream-scripting",
        schema_key="delta_te.tsv",
    ),
    OutputDescriptor(
        output_type="rna_counts",
        stage="rnaseq",
        relative_path="rnaseq/rna_counts.tsv",
        description="Per-sample mt-mRNA count matrix from the from-FASTQ flow.",
        recommended_for="downstream-scripting",
        schema_key="rna_counts.tsv",
    ),
    OutputDescriptor(
        output_type="exploratory_sidecar",
        stage="rnaseq",
        relative_path="rnaseq/EXPLORATORY.md",
        description=(
            "Disclaimer written when rnaseq ran in pseudo-replicate or "
            "exploratory mode; absence means publication-grade."
        ),
        recommended_for="reviewer-spot-check",
    ),
    OutputDescriptor(
        output_type="rnaseq_run_settings",
        stage="rnaseq",
        relative_path="rnaseq/run_settings.json",
        description="Resolved rnaseq-stage configuration and tool versions.",
        recommended_for="reviewer-spot-check",
    ),
    # all -----------------------------------------------------------------
    OutputDescriptor(
        output_type="run_manifest",
        stage="all",
        relative_path="run_manifest.json",
        description=(
            "Top-level run provenance: schema, versions, git commit, "
            "config canonicalisation, sample-sheet hash, tools, "
            "warnings, outputs."
        ),
        recommended_for="reviewer-spot-check",
    ),
    OutputDescriptor(
        output_type="warnings",
        stage="all",
        relative_path="warnings.tsv",
        description=(
            "Structured warnings collected across all stages "
            "(stage / sample_id / severity / code / message / suggested_action)."
        ),
        recommended_for="downstream-scripting",
        schema_key="warnings.tsv",
    ),
    OutputDescriptor(
        output_type="summary_qc",
        stage="all",
        relative_path="summary_qc.tsv",
        description="Per-sample QC summary across align + rpf + rnaseq.",
        recommended_for="downstream-scripting",
        schema_key="summary_qc.tsv",
    ),
    OutputDescriptor(
        output_type="summary_md",
        stage="all",
        relative_path="SUMMARY.md",
        description="Human-readable run summary linking to the major outputs.",
        recommended_for="reviewer-spot-check",
    ),
)


def build_outputs_index_rows(run_root: Path) -> list[dict[str, str]]:
    """Return the row dicts for the ``outputs_index.tsv`` of *run_root*.

    Outputs that don't exist on disk are skipped. The ``schema_version``
    column is populated from :data:`OUTPUT_SCHEMA_VERSIONS` when the
    descriptor names a registered TSV; otherwise it is empty.
    """
    rows: list[dict[str, str]] = []
    run_root = Path(run_root)
    for desc in _KNOWN_OUTPUTS:
        target = run_root / desc.relative_path
        if not target.exists():
            continue
        rows.append(
            {
                "output_type": desc.output_type,
                "stage": desc.stage,
                "path": desc.relative_path,
                "description": desc.description,
                "recommended_for": desc.recommended_for,
                "schema_version": (
                    OUTPUT_SCHEMA_VERSIONS.get(desc.schema_key, "")
                    if desc.schema_key
                    else ""
                ),
            }
        )
    return rows


def _scrub_cell(value: str) -> str:
    return value.replace("\t", " ").replace("\n", " ")


def write_outputs_index(run_root: Path | str) -> Path:
    """Write ``<run_root>/outputs_index.tsv`` and return the path."""
    run_root = Path(run_root)
    run_root.mkdir(parents=True, exist_ok=True)
    rows = build_outputs_index_rows(run_root)

    target = run_root / "outputs_index.tsv"
    with target.open("w", encoding="utf-8") as handle:
        handle.write(schema_header_line("outputs_index.tsv"))
        handle.write("\t".join(_OUTPUTS_INDEX_COLUMNS) + "\n")
        for row in rows:
            handle.write(
                "\t".join(_scrub_cell(row[c]) for c in _OUTPUTS_INDEX_COLUMNS)
                + "\n"
            )
    return target
