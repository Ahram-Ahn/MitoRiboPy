"""Stable error codes used across MitoRiboPy CLI surfaces.

Every code is a string constant rather than an enum so it survives JSON
round-trips into ``warnings.tsv`` and ``run_manifest.json`` without
extra serialisation. The constant name and the literal value match
intentionally: a downstream script can ``grep E_OFFSET_LOW_CONFIDENCE``
across either source code or run logs and find both.

Code prefixes:

* ``E_CONFIG_*``       — YAML / CLI configuration problems detected by
                         ``validate-config`` or the per-stage parsers.
* ``E_SAMPLE_SHEET_*`` — sample-sheet schema or content errors.
* ``E_REFERENCE_*``    — reference FASTA / annotation / checksum
                         drift between stages.
* ``E_OFFSET_*``       — RPF P-site offset selection diagnostics.
* ``E_RNASEQ_*``       — rnaseq stage gating (replicates, mode mismatch).
* ``E_TOOL_*``         — external-tool resolution / version problems.
* ``E_RESUME_*``       — resume-guard hash mismatch reasons.
* ``E_OUTPUT_*``       — output-contract violations (count invariants
                         in read_counts.tsv, etc.).

These codes are part of the public CLI contract: removing or renaming
one is a breaking change. Adding a new code is additive.
"""

from __future__ import annotations


__all__ = [
    "E_CONFIG_UNKNOWN_KEY",
    "E_CONFIG_MUTUALLY_EXCLUSIVE",
    "E_CONFIG_MISSING_REQUIRED",
    "E_CONFIG_INVALID_VALUE",
    "E_SAMPLE_SHEET_MISSING_COLUMN",
    "E_SAMPLE_SHEET_DUPLICATE_ROW",
    "E_SAMPLE_SHEET_BAD_ASSAY",
    "E_SAMPLE_SHEET_MISSING_FASTQ",
    "E_REFERENCE_HASH_MISMATCH",
    "E_REFERENCE_VALIDATION_FAILED",
    "E_OFFSET_LOW_CONFIDENCE",
    "E_OFFSET_FALLBACK_USED",
    "E_RNASEQ_INSUFFICIENT_REPLICATES",
    "E_RNASEQ_GENE_ID_UNMATCHED",
    "E_RNASEQ_REFERENCE_DRIFT",
    "E_TOOL_NOT_FOUND",
    "E_TOOL_VERSION_INCOMPATIBLE",
    "E_RESUME_HASH_MISMATCH",
    "E_OUTPUT_COUNT_INVARIANT",
]


# Config-level errors -----------------------------------------------------

E_CONFIG_UNKNOWN_KEY = "E_CONFIG_UNKNOWN_KEY"
E_CONFIG_MUTUALLY_EXCLUSIVE = "E_CONFIG_MUTUALLY_EXCLUSIVE"
E_CONFIG_MISSING_REQUIRED = "E_CONFIG_MISSING_REQUIRED"
E_CONFIG_INVALID_VALUE = "E_CONFIG_INVALID_VALUE"


# Sample-sheet errors -----------------------------------------------------

E_SAMPLE_SHEET_MISSING_COLUMN = "E_SAMPLE_SHEET_MISSING_COLUMN"
E_SAMPLE_SHEET_DUPLICATE_ROW = "E_SAMPLE_SHEET_DUPLICATE_ROW"
E_SAMPLE_SHEET_BAD_ASSAY = "E_SAMPLE_SHEET_BAD_ASSAY"
E_SAMPLE_SHEET_MISSING_FASTQ = "E_SAMPLE_SHEET_MISSING_FASTQ"


# Reference / checksum errors ---------------------------------------------

E_REFERENCE_HASH_MISMATCH = "E_REFERENCE_HASH_MISMATCH"
E_REFERENCE_VALIDATION_FAILED = "E_REFERENCE_VALIDATION_FAILED"


# RPF offset diagnostics --------------------------------------------------

E_OFFSET_LOW_CONFIDENCE = "E_OFFSET_LOW_CONFIDENCE"
E_OFFSET_FALLBACK_USED = "E_OFFSET_FALLBACK_USED"


# RNA-seq stage gating ----------------------------------------------------

E_RNASEQ_INSUFFICIENT_REPLICATES = "E_RNASEQ_INSUFFICIENT_REPLICATES"
E_RNASEQ_GENE_ID_UNMATCHED = "E_RNASEQ_GENE_ID_UNMATCHED"
E_RNASEQ_REFERENCE_DRIFT = "E_RNASEQ_REFERENCE_DRIFT"


# Tool resolution ---------------------------------------------------------

E_TOOL_NOT_FOUND = "E_TOOL_NOT_FOUND"
E_TOOL_VERSION_INCOMPATIBLE = "E_TOOL_VERSION_INCOMPATIBLE"


# Resume guard ------------------------------------------------------------

E_RESUME_HASH_MISMATCH = "E_RESUME_HASH_MISMATCH"


# Output contract ---------------------------------------------------------

E_OUTPUT_COUNT_INVARIANT = "E_OUTPUT_COUNT_INVARIANT"
