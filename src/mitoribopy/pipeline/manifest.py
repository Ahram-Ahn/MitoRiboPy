"""Manifest builder for ``mitoribopy all``.

Extracted from ``cli/all_.py`` (which had grown past 2 300 LoC) so the
manifest layout, schema version, hash helpers, and JSON-Schema
validation hook live in one place. The orchestrator imports
:func:`write_manifest` from here; the public surface
(``MANIFEST_SCHEMA_VERSION``, the per-stage status enum, the warning
code) is re-exported from ``cli/all_.py`` for backwards compatibility
with downstream consumers that import from the CLI module.

Versioning rules for ``MANIFEST_SCHEMA_VERSION`` (semver-ish):

* MAJOR — a top-level required field is removed or its semantics change.
* MINOR — a top-level field is added at the end.
* PATCH — purely cosmetic (comment / ordering changes).
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from pathlib import Path

from .. import __version__


__all__ = [
    "MANIFEST_SCHEMA_VERSION",
    "W_MANIFEST_SCHEMA_DRIFT_CODE",
    "sha256_of",
    "read_stage_settings",
    "git_commit",
    "yaml_dump",
    "build_stages_block",
    "lift_tool_versions",
    "write_manifest",
]


# Bumped whenever the run_manifest.json layout changes in a way that
# breaks downstream consumers (added fields are minor; renamed or
# removed fields are major). Read it from your own scripts to gate on a
# compatible manifest shape.
MANIFEST_SCHEMA_VERSION = "1.3.0"  # 1.3: + resource_plan (top-level execution audit, v0.6.2)


# Warning code emitted when validate_run_manifest finds drift between
# the live manifest and the bundled schema. Defined here so the
# extraction does not introduce a new import cycle through
# mitoribopy.errors.
W_MANIFEST_SCHEMA_DRIFT_CODE = "W_MANIFEST_SCHEMA_DRIFT"


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------


def sha256_of(path: Path | str | None) -> str | None:
    """Return the hex SHA256 digest of *path*, or ``None`` on any I/O error.

    ``None`` input also returns ``None`` so callers can write
    ``sha256_of(maybe_path)`` without a guard.
    """
    if path is None:
        return None
    try:
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for chunk in iter(lambda: handle.read(65536), b""):
                digest.update(chunk)
        return digest.hexdigest()
    except OSError:
        return None


def read_stage_settings(stage_dir: Path) -> dict | None:
    """Read ``<stage_dir>/run_settings.json`` if present and parseable."""
    path = Path(stage_dir) / "run_settings.json"
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def git_commit() -> str | None:
    """Best-effort current git commit — never raises.

    Returns ``None`` when not in a repo, when git is missing, or when
    the call fails for any reason. Intentionally cheap; we never want a
    manifest write to fail because git is misconfigured.
    """
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=2,
            check=False,
        )
    except (FileNotFoundError, subprocess.SubprocessError, OSError):
        return None
    if out.returncode != 0:
        return None
    sha = out.stdout.strip()
    return sha or None


def yaml_dump(payload: dict) -> str:
    """YAML dump if PyYAML is available, JSON fallback otherwise.

    The fallback is intentional: ``--print-canonical-config`` needs to
    work in lean install environments where PyYAML is not installed.
    JSON is a strict subset of YAML so the output remains valid YAML.
    """
    try:
        import yaml  # type: ignore[import-not-found]
    except ImportError:
        return json.dumps(payload, indent=2, sort_keys=True) + "\n"
    return yaml.safe_dump(
        payload, sort_keys=True, default_flow_style=False, allow_unicode=True
    )


def build_stages_block(
    stages_run: list[str],
    stages_skipped: list[str],
    runtimes: dict[str, float],
    skip_reasons: dict[str, str],
) -> dict[str, dict]:
    """Reshape parallel run/skipped lists into a ``{stage: {status, ...}}`` map.

    Status values:

    * ``completed`` — the stage executed successfully.
    * ``skipped``   — the stage was not configured or was skipped via
                      ``--skip-*`` / ``--resume``; ``reason`` carries the
                      trigger.
    * ``not_configured`` — no section for the stage in the YAML.
    """
    out: dict[str, dict] = {}
    for stage in ("align", "rpf", "rnaseq"):
        if stage in stages_run:
            out[stage] = {"status": "completed"}
            if stage in runtimes:
                out[stage]["runtime_seconds"] = round(runtimes[stage], 3)
        elif stage in stages_skipped:
            entry: dict = {"status": "skipped"}
            reason = skip_reasons.get(stage)
            if reason:
                entry["reason"] = reason
            out[stage] = entry
        else:
            out[stage] = {"status": "not_configured"}
    return out


def lift_tool_versions(
    align: dict | None, rpf: dict | None, rnaseq: dict | None,
) -> dict[str, str]:
    """Pull tool versions out of per-stage ``run_settings.json`` files.

    Stages already record (varying subsets of) their tool versions in
    their own run_settings; this just lifts whatever's there into a flat
    top-level map. Missing keys stay missing — we don't probe.
    """
    out: dict[str, str] = {}
    out["python"] = platform.python_version()
    out["mitoribopy"] = __version__
    for settings in (align, rpf, rnaseq):
        if not isinstance(settings, dict):
            continue
        tools = settings.get("tools")
        if isinstance(tools, dict):
            for k, v in tools.items():
                if v is not None and k not in out:
                    out[k] = str(v)
        # Some stages stash individual tool versions at the top level.
        for legacy_key in (
            "cutadapt_version",
            "bowtie2_version",
            "umi_tools_version",
            "pysam_version",
            "pydeseq2_version",
        ):
            v = settings.get(legacy_key)
            if v is not None:
                tool_name = legacy_key.removesuffix("_version")
                out.setdefault(tool_name, str(v))
    return out


# ---------------------------------------------------------------------------
# Manifest writer
# ---------------------------------------------------------------------------


def write_manifest(
    output_dir: Path,
    manifest_name: str,
    *,
    stages_run: list[str],
    stages_skipped: list[str],
    align_settings: dict | None,
    rpf_settings: dict | None,
    rnaseq_settings: dict | None,
    config_canonical: dict,
    config_source_path: str | None,
    sample_sheet_path: str | None,
    command_argv: list[str],
    runtimes: dict[str, float],
    skip_reasons: dict[str, str],
    total_runtime_seconds: float | None = None,
) -> Path:
    """Write the ``run_manifest.json`` for an ``all`` run.

    Layout (top-level keys):

    * ``schema_version``       — version of THIS layout (see
                                 :data:`MANIFEST_SCHEMA_VERSION`).
    * ``mitoribopy_version``   — package version that produced the run.
    * ``git_commit``           — current commit when run inside a repo,
                                 ``null`` otherwise.
    * ``command``              — the original argv joined with spaces.
    * ``config_source``        — path to the user-supplied YAML.
    * ``config_source_sha256`` — SHA256 of that file as the user wrote
                                 it (useful for "is this the same input
                                 I used last time?" diffs).
    * ``config_canonical``     — the merged + auto-wired config that
                                 actually drove the run.
    * ``sample_sheet``         — path to the unified sheet, when set.
    * ``sample_sheet_sha256``  — its SHA256.
    * ``reference_checksum``   — promoted from rpf settings (legacy
                                 field; the rnaseq subcommand reads
                                 this directly).
    * ``stages``               — ``{align: {...}, rpf: {...}, rnaseq: {...}}``
                                 with ``status`` + optional
                                 ``runtime_seconds`` / ``reason``.
    * ``align`` / ``rpf`` / ``rnaseq`` — full per-stage
                                 ``run_settings.json`` copies (kept
                                 as-is for back-compat).
    * ``tools``                — flat ``{tool: version}`` map lifted
                                 from per-stage settings.
    * ``warnings``             — structured warnings collected by
                                 ``mitoribopy.io.warnings_log`` during
                                 the run, mirrored to ``warnings.tsv``.
    * ``outputs``              — outputs_index rows (output_type,
                                 stage, path, description,
                                 recommended_for, schema_version).
                                 Same data is also written to
                                 ``outputs_index.tsv``.
    * ``runtime_seconds``      — total wall time for the orchestrator;
                                 per-stage values live under
                                 ``stages.<stage>.runtime_seconds``.
    * ``platform``             — :func:`platform.platform` string.
    * ``python_version``       — running Python version.
    * ``resource_plan``        — mirror of ``resource_plan.json``;
                                 ``null`` when the orchestrator could
                                 not write the sidecar.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Local imports keep this module free of circulars at import time.
    from ..io.outputs_index import (
        build_outputs_index_rows,
        write_outputs_index,
    )
    from ..io.schema_versions import OUTPUT_SCHEMA_VERSIONS
    from ..io.warnings_log import collected as _collected_warnings
    from ..io.warnings_log import flush_tsv as _flush_warnings_tsv

    structured_warnings = [w.as_dict() for w in _collected_warnings()]
    _flush_warnings_tsv(output_dir / "warnings.tsv")
    # Write outputs_index.tsv before the manifest so the manifest can
    # reference the same row set under "outputs".
    write_outputs_index(output_dir)
    outputs_rows = build_outputs_index_rows(output_dir)

    # Lift resource_plan.json into the manifest so a downstream script
    # does not have to re-read the sidecar file. The path is written by
    # the orchestrator before any stage runs; tolerate a missing file
    # (e.g. permissions issues) and just record null.
    resource_plan_blob: dict | None = None
    plan_path = output_dir / "resource_plan.json"
    if plan_path.is_file():
        try:
            resource_plan_blob = json.loads(plan_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            resource_plan_blob = None

    manifest: dict = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "subcommand": "all",
        "mitoribopy_version": __version__,
        "git_commit": git_commit(),
        "command": "mitoribopy all " + " ".join(command_argv),
        "config_source": config_source_path,
        "config_source_sha256": (
            sha256_of(config_source_path) if config_source_path else None
        ),
        "config_canonical": config_canonical,
        "sample_sheet": sample_sheet_path,
        "sample_sheet_sha256": (
            sha256_of(sample_sheet_path) if sample_sheet_path else None
        ),
        "stages": build_stages_block(
            stages_run, stages_skipped, runtimes, skip_reasons,
        ),
        "align": align_settings,
        "rpf": rpf_settings,
        "rnaseq": rnaseq_settings,
        "tools": lift_tool_versions(
            align_settings, rpf_settings, rnaseq_settings,
        ),
        "output_schemas": dict(OUTPUT_SCHEMA_VERSIONS),
        "warnings": structured_warnings,
        "outputs": outputs_rows,
        "runtime_seconds": (
            float(total_runtime_seconds)
            if total_runtime_seconds is not None
            else sum(runtimes.values())
        ),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "resource_plan": resource_plan_blob,
    }

    # Promote rpf's reference_checksum so future rnaseq invocations
    # reading this manifest directly can find it without drilling into
    # the rpf section.
    if rpf_settings and rpf_settings.get("reference_checksum"):
        manifest["reference_checksum"] = rpf_settings["reference_checksum"]

    path = output_dir / manifest_name
    path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8",
    )

    # Validate the manifest we just wrote against the package JSON
    # Schema. A failure here means a code change drifted away from the
    # documented reproducibility surface; we record a structured
    # warning rather than failing the run, so the user still gets the
    # rest of the pipeline outputs.
    try:
        from ..data import validate_run_manifest
        from ..io import warnings_log

        errors = validate_run_manifest(manifest)
        if errors:
            warnings_log.record(
                stage="ALL",
                code=W_MANIFEST_SCHEMA_DRIFT_CODE,
                severity="warn",
                message=(
                    "run_manifest.json failed JSON-Schema validation: "
                    + "; ".join(errors[:5])
                    + (" (+ more)" if len(errors) > 5 else "")
                ),
                suggested_action=(
                    "Compare src/mitoribopy/data/run_manifest.schema.json "
                    "against the manifest builder; update one of them so "
                    "downstream reproducibility tooling does not rely on "
                    "an undocumented field."
                ),
            )
    except (ImportError, OSError) as exc:  # pragma: no cover - defensive
        from ..console import log_warning

        log_warning(
            "ALL", f"run_manifest.json schema validation skipped: {exc}",
        )

    return path
