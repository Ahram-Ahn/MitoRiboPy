"""Hash-validated resume for ``mitoribopy all --resume``.

The orchestrator's prior implementation skipped a stage whenever its
sentinel output file existed (``read_counts.tsv`` for align,
``rpf_counts.tsv`` for rpf, ``delta_te.tsv`` for rnaseq). That is
correct when the inputs and config have not changed. It is *silently
wrong* when the user edited the YAML, swapped a sample sheet, or
re-built the reference FASTA between runs — the stale output now
encodes one set of decisions and the new run is asked for another.

This module compares the *current* config / sample-sheet / reference
hashes against the values recorded in the prior run's
``run_manifest.json``. On mismatch, the guard refuses to skip the
affected stage(s) unless the user has passed ``--force-resume`` (also
honoured via ``MITORIBOPY_FORCE_RESUME=1`` for CI scripts that cannot
easily change argv).

The guard is intentionally conservative: it flags everything it can
detect and lets the orchestrator translate the report into a stage
re-run. False positives (e.g. cosmetic config edits that happen to
change the SHA but not the behaviour) are recoverable with
``--force-resume``; false negatives (silently using a stale output
under a new config) are not.
"""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass, field
from pathlib import Path


__all__ = [
    "FORCE_RESUME_ENV",
    "ResumeMismatch",
    "ResumeReport",
    "force_resume_requested",
    "load_prior_manifest",
    "sha256_of",
    "validate_resume",
]


FORCE_RESUME_ENV = "MITORIBOPY_FORCE_RESUME"


def force_resume_requested(*, cli_flag: bool) -> bool:
    """``True`` if the user asked to bypass the hash guard.

    The CLI flag wins; the env var is a CI-friendly fallback so
    long-lived scripts that fence on ``--resume`` can opt in without
    a code change.
    """
    if cli_flag:
        return True
    return os.environ.get(FORCE_RESUME_ENV, "").lower() in {"1", "true", "yes"}


def sha256_of(path: Path | str | None) -> str | None:
    """SHA-256 of a file, or ``None`` if the path is missing/unreadable."""
    if path is None:
        return None
    p = Path(path)
    if not p.is_file():
        return None
    digest = hashlib.sha256()
    try:
        with p.open("rb") as fh:
            for chunk in iter(lambda: fh.read(65536), b""):
                digest.update(chunk)
    except OSError:
        return None
    return digest.hexdigest()


def load_prior_manifest(run_root: Path, manifest_name: str) -> dict | None:
    """Load ``<run_root>/<manifest_name>`` if present, else ``None``.

    Returns ``None`` (not an error) when no manifest exists — that just
    means there is nothing to compare against, and the resume falls back
    to the sentinel-file check.
    """
    path = Path(run_root) / manifest_name
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


@dataclass(frozen=True)
class ResumeMismatch:
    """One field's prior-vs-current value disagreement."""

    field: str
    prior: object
    current: object

    def render(self) -> str:
        return (
            f"  - {self.field}: prior={self.prior!r}, current={self.current!r}"
        )


@dataclass(frozen=True)
class ResumeReport:
    """Aggregate result of comparing prior manifest to current inputs.

    ``mismatches`` is empty when every guarded field matches (the
    sentinel-file check decides whether to skip stages); non-empty
    when at least one field changed and the stage(s) should be
    re-run unless ``--force-resume`` is set.

    ``no_prior_manifest`` is true when ``run_manifest.json`` was
    missing in ``run_root``; the orchestrator treats that as
    "first run" and falls back to the sentinel-file check.
    """

    mismatches: tuple[ResumeMismatch, ...] = field(default_factory=tuple)
    no_prior_manifest: bool = False

    @property
    def ok(self) -> bool:
        return not self.mismatches

    def render(self) -> str:
        if self.no_prior_manifest:
            return "no prior run_manifest.json found; cannot validate resume."
        if self.ok:
            return "resume hashes match prior manifest."
        head = (
            f"resume hash mismatch ({len(self.mismatches)} field(s) "
            "differ from the prior run_manifest.json):"
        )
        body = "\n".join(m.render() for m in self.mismatches)
        return head + "\n" + body


def validate_resume(
    *,
    prior_manifest: dict | None,
    config_path: Path | str | None,
    sample_sheet_path: Path | str | None,
    reference_fasta: Path | str | None,
    mitoribopy_version: str,
    manifest_schema_version: str,
) -> ResumeReport:
    """Compare current inputs against a prior manifest's recorded hashes.

    Guarded fields (each compared independently; absent prior fields
    are skipped, so older manifests do not break resume):

    * ``config_source_sha256``     SHA of the YAML driving the run
    * ``sample_sheet_sha256``      SHA of the unified sample sheet
    * ``reference_checksum``       SHA of the reference FASTA
    * ``mitoribopy_version``       installed package version
    * ``schema_version``           manifest layout version

    The reference FASTA is hashed only when the caller passes a path;
    the orchestrator typically extracts it from ``rpf.fasta`` in the
    canonical config.
    """
    if prior_manifest is None:
        return ResumeReport(no_prior_manifest=True)

    mismatches: list[ResumeMismatch] = []

    def _check(field_name: str, prior_value: object, current_value: object) -> None:
        # Skip when the prior manifest doesn't have this field at all
        # (older manifest formats predate it).
        if prior_value is None:
            return
        if current_value != prior_value:
            mismatches.append(
                ResumeMismatch(
                    field=field_name, prior=prior_value, current=current_value
                )
            )

    _check(
        "config_source_sha256",
        prior_manifest.get("config_source_sha256"),
        sha256_of(config_path),
    )
    _check(
        "sample_sheet_sha256",
        prior_manifest.get("sample_sheet_sha256"),
        sha256_of(sample_sheet_path),
    )
    _check(
        "reference_checksum",
        prior_manifest.get("reference_checksum"),
        sha256_of(reference_fasta),
    )
    _check(
        "mitoribopy_version",
        prior_manifest.get("mitoribopy_version"),
        mitoribopy_version,
    )
    _check(
        "schema_version",
        prior_manifest.get("schema_version"),
        manifest_schema_version,
    )

    return ResumeReport(mismatches=tuple(mismatches))
