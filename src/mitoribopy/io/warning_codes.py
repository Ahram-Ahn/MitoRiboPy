"""Stable registry of every warning / error code emitted by MitoRiboPy.

The registry is the single source of truth for two consumers:

1. :mod:`mitoribopy.io.warnings_log` — every ``record(code=...)`` call
   that names a code is checked against this registry in test mode so
   a typo cannot silently ship a brand-new alias.
2. ``docs/reference/warning_codes.md`` — the user-facing table is
   regenerated from this registry by the docs generator. The Python
   side never reads the docs, so there is no risk of the table going
   stale relative to the runtime.

Adding a new code requires:

* A new :class:`WarningCode` entry in :data:`WARNING_CODES`.
* A docstring-quality ``message_template`` describing what the user
  sees and how to fix it.
* A call to :func:`mitoribopy.io.warnings_log.record` somewhere in
  the codebase (otherwise the registry entry is unreachable; the
  test suite will detect this and fail).

Severity levels:

* ``"error"`` — the run cannot finish correctly with this condition;
  ``mitoribopy all --strict`` already promotes such a record to a
  hard exit. Non-strict mode keeps running so the user still gets
  the rest of the manifest.
* ``"warn"`` — the run will produce output, but the output may not
  be publication-grade.
* ``"info"`` — informational; surfaced for audit purposes but does
  not change runtime behaviour.

Stage values match the run-manifest's stage names: ``"align"``,
``"rpf"``, ``"rnaseq"``, ``"all"`` (orchestrator-level).
"""

from __future__ import annotations

from dataclasses import dataclass


__all__ = [
    "WARNING_CODES",
    "WarningCode",
    "lookup",
]


@dataclass(frozen=True)
class WarningCode:
    """One stable entry in the warning-code registry."""

    code: str
    severity: str  # "error" | "warn" | "info"
    stage: str  # "align" | "rpf" | "rnaseq" | "all"
    sample_scoped: bool
    summary: str
    remediation: str


WARNING_CODES: tuple[WarningCode, ...] = (
    # ---- Config / sample-sheet -------------------------------------------
    WarningCode(
        code="E_CONFIG_UNKNOWN_KEY",
        severity="error",
        stage="all",
        sample_scoped=False,
        summary="Top-level YAML key not in the canonical schema.",
        remediation=(
            "Remove the key, or run `mitoribopy migrate-config` to "
            "rewrite legacy spellings into the canonical form."
        ),
    ),
    WarningCode(
        code="E_CONFIG_DEPRECATED_KEY_STRICT",
        severity="error",
        stage="all",
        sample_scoped=False,
        summary=(
            "A deprecated key was rewritten by the canonicaliser, but "
            "`--strict` was set."
        ),
        remediation=(
            "Update the YAML to use the canonical key name reported in "
            "the migrate log; re-run with `--strict`."
        ),
    ),
    WarningCode(
        code="E_SAMPLE_SHEET_DUPLICATE_ID",
        severity="error",
        stage="all",
        sample_scoped=True,
        summary="Two rows in the unified sample sheet share a sample_id.",
        remediation="Make sample_id unique across every row of the sheet.",
    ),
    WarningCode(
        code="E_SAMPLE_SHEET_UNKNOWN_COLUMN",
        severity="error",
        stage="all",
        sample_scoped=False,
        summary="Sample sheet has a column not in the documented schema.",
        remediation=(
            "Drop the column, or move it under a project-specific "
            "name that does not collide with the schema."
        ),
    ),
    # ---- Align preprocessing --------------------------------------------
    WarningCode(
        code="E_FASTQ_MATE_PAIRING_FAILED",
        severity="error",
        stage="align",
        sample_scoped=True,
        summary="Could not pair R1/R2 FASTQs from filename token.",
        remediation=(
            "Rename to `<sample>_R1.fq.gz` / `<sample>_R2.fq.gz` (or "
            "the documented mate-token aliases), or pass FASTQs as "
            "explicit pairs in the sample sheet."
        ),
    ),
    WarningCode(
        code="E_ADAPTER_AMBIGUOUS",
        severity="error",
        stage="align",
        sample_scoped=True,
        summary=(
            "Adapter auto-detection saw two candidate kits above the "
            "match threshold and could not pick one."
        ),
        remediation=(
            "Pin the 3' adapter explicitly with `--adapter <SEQ>` (CLI) "
            "or `adapter: <SEQ>` (YAML / sample sheet)."
        ),
    ),
    WarningCode(
        code="E_UMI_DECLARED_BUT_NOT_FOUND",
        severity="error",
        stage="align",
        sample_scoped=True,
        summary=(
            "Resolved kit declared a UMI but no UMI bases were "
            "extracted by cutadapt."
        ),
        remediation=(
            "Verify `umi_length` / `umi_position` against the library "
            "prep documentation; pre-trimmed reads may need "
            "`--pretrimmed` (or `pretrimmed: true` in YAML) instead."
        ),
    ),
    WarningCode(
        code="W_UMI_SHORT_COLLISION_RISK",
        severity="warn",
        stage="align",
        sample_scoped=True,
        summary="UMI shorter than 8 nt has a high random-collision rate.",
        remediation=(
            "Treat duplicate-fraction estimates as conservative; for "
            "publication-grade dedup, prefer kits with UMI ≥ 8 nt."
        ),
    ),
    WarningCode(
        code="W_UMI_HIGH_DUPLICATE_FRACTION",
        severity="warn",
        stage="align",
        sample_scoped=True,
        summary=(
            "Coordinate+UMI dedup removed > 80 % of mapped reads — "
            "library complexity is low."
        ),
        remediation=(
            "Re-check input depth / library complexity before drawing "
            "biological conclusions."
        ),
    ),
    WarningCode(
        code="W_UMI_DEDUP_SKIPPED_WITH_UMI_PRESENT",
        severity="warn",
        stage="align",
        sample_scoped=True,
        summary=(
            "UMI bases were extracted but `dedup_strategy: skip` was "
            "set explicitly."
        ),
        remediation=(
            "Switch to `dedup_strategy: umi_coordinate` (or `auto`) "
            "unless you have a specific reason to skip dedup."
        ),
    ),
    WarningCode(
        code="W_UMI_METHOD_DIRECTIONAL_ON_SHORT_UMI",
        severity="warn",
        stage="align",
        sample_scoped=True,
        summary=(
            "`umi_dedup_method: directional` was selected on a UMI "
            "shorter than 8 nt; clustering tends to over-collapse."
        ),
        remediation=(
            "Use `umi_dedup_method: unique` for short UMIs; "
            "directional clustering is recommended only for ≥ 8 nt UMIs."
        ),
    ),
    # ---- Rpf periodicity / codon correlation ----------------------------
    WarningCode(
        code="W_CODON_RAW_COUNT_PRIMARY",
        severity="warn",
        stage="rpf",
        sample_scoped=False,
        summary=(
            "Codon-correlation panel was rendered with `--cor-metric "
            "raw_count`; this is exploratory only."
        ),
        remediation=(
            "Re-run with `--cor-metric log2_density_rpm` for the "
            "publication-facing scatter."
        ),
    ),
    WarningCode(
        code="W_PERIODICITY_LOW_DEPTH",
        severity="warn",
        stage="rpf",
        sample_scoped=True,
        summary=(
            "Sample has fewer assigned P-sites than the configured "
            "min_reads_per_length threshold."
        ),
        remediation=(
            "Treat the sample's frame fractions as `low_depth`; do not "
            "use it to validate the offset choice."
        ),
    ),
    WarningCode(
        code="W_PERIODICITY_WEAK",
        severity="warn",
        stage="rpf",
        sample_scoped=True,
        summary=(
            "Expected-frame fraction is below the warn threshold "
            "(default 0.50)."
        ),
        remediation=(
            "Inspect fourier_period3_score_combined.tsv (snr_call column) "
            "and the per-(sample, read_length) Fourier overlay PNGs under "
            "rpf/qc/fourier_spectrum/ for length-specific causes; consider "
            "tightening the read-length window."
        ),
    ),
    # ---- Rnaseq / TE -----------------------------------------------------
    WarningCode(
        code="W_RNASEQ_FROM_FASTQ_EXPLORATORY",
        severity="warn",
        stage="rnaseq",
        sample_scoped=False,
        summary=(
            "rnaseq mode = from_fastq (mt-mRNA-only pyDESeq2). Output "
            "is exploratory and not publication-grade."
        ),
        remediation=(
            "For publication, run DESeq2 / Xtail / Anota2Seq on the "
            "full transcriptome externally and re-run rnaseq with "
            "`--rnaseq-mode de_table --de-table <path>`."
        ),
    ),
    WarningCode(
        code="W_RNASEQ_PSEUDO_REPLICATE",
        severity="warn",
        stage="rnaseq",
        sample_scoped=True,
        summary=(
            "n=1 condition was split into FASTQ-record-parity "
            "pseudo-replicates; padj is NOT biologically defensible."
        ),
        remediation=(
            "Use only for demo / smoke tests. For publication-grade "
            "DE add real biological replicates."
        ),
    ),
    WarningCode(
        code="E_REFERENCE_CHECKSUM_MISMATCH",
        severity="error",
        stage="rnaseq",
        sample_scoped=False,
        summary=(
            "rpf and rnaseq references disagree (different SHA256). "
            "TE comparisons across mismatched transcript sets are "
            "biologically meaningless."
        ),
        remediation=(
            "Rebuild the reference once and re-run both stages; the "
            "manifest's `reference_checksum` must match across rpf "
            "and rnaseq."
        ),
    ),
    WarningCode(
        code="E_TEMPLATE_INVALID",
        severity="error",
        stage="all",
        sample_scoped=False,
        summary=(
            "A shipped YAML / shell template failed `validate-config "
            "--strict`."
        ),
        remediation=(
            "Open an issue against MitoRiboPy itself; this is a release-"
            "engineering problem, not a user problem."
        ),
    ),
)


_BY_CODE: dict[str, WarningCode] = {entry.code: entry for entry in WARNING_CODES}


def lookup(code: str) -> WarningCode | None:
    """Return the registered :class:`WarningCode`, or ``None`` if unknown."""
    return _BY_CODE.get(code)
