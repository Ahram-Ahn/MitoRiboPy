"""Process-wide structured warnings collector for MitoRiboPy.

Every warning emitted via :func:`record` is

* mirrored to stderr through :func:`mitoribopy.console.log_warning`,
* appended to the in-memory list returned by :func:`collected`,
* persisted to ``warnings.tsv`` at the run root by :func:`flush_tsv`
  (called by ``mitoribopy all`` after the per-stage work completes), and
* embedded in ``run_manifest.json`` under the ``warnings`` field
  (previously a placeholder empty list).

Records are intentionally append-only and process-global: subcommand
boundaries do NOT clear the list, so a later ``mitoribopy summarize``
or ``mitoribopy all`` invocation that re-runs the rnaseq stage can
collect warnings from every stage in a single TSV.

P5.8 schema 2.0: TSV columns are now
``stage, sample_id, severity, code, message, suggested_action`` to
match the assessment §8 spec. ``timestamp`` is retained on the
in-memory :class:`WarningRecord` (and therefore in the manifest's
``warnings`` JSON array) for diagnostics, but is no longer persisted
to the TSV.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

from .schema_versions import schema_header_line


__all__ = [
    "WarningRecord",
    "clear",
    "collected",
    "flush_tsv",
    "record",
]


_WARNING_TSV_COLUMNS: tuple[str, ...] = (
    "stage",
    "sample_id",
    "severity",
    "code",
    "message",
    "suggested_action",
)


@dataclass(frozen=True)
class WarningRecord:
    """One structured warning captured during a MitoRiboPy run.

    The dataclass keeps ``timestamp`` and the canonical ``stage`` /
    ``sample_id`` field names. Backwards-compatible properties named
    ``component`` and ``sample`` are exposed so older callers keep
    working without an immediate rewrite.
    """

    timestamp: str
    stage: str
    severity: str  # "warn" | "error" | "info"
    sample_id: str | None
    code: str | None
    message: str
    suggested_action: str | None = None

    # Back-compat aliases (read-only) ----------------------------------------
    @property
    def component(self) -> str:  # noqa: D401
        """Deprecated alias for :attr:`stage`."""
        return self.stage

    @property
    def sample(self) -> str | None:  # noqa: D401
        """Deprecated alias for :attr:`sample_id`."""
        return self.sample_id

    def as_dict(self) -> dict:
        return asdict(self)


_RECORDS: list[WarningRecord] = []


def record(
    stage: str | None = None,
    message: str = "",
    *,
    sample_id: str | None = None,
    code: str | None = None,
    severity: str = "warn",
    suggested_action: str | None = None,
    # Back-compat kwargs (deprecated; still accepted for older callers):
    component: str | None = None,
    sample: str | None = None,
) -> WarningRecord:
    """Capture a structured warning AND mirror it to the console.

    Either ``stage=`` (canonical) or ``component=`` (legacy alias) is
    accepted; the same goes for ``sample_id=`` vs ``sample=``. The
    console mirror is intentional: existing user-facing behaviour of
    ``log_warning`` is preserved, and the structured TSV / manifest
    entries are additive.
    """
    # Local import to avoid a console <-> warnings_log circular at
    # module-init time.
    from ..console import log_warning

    resolved_stage = stage if stage is not None else (component or "GENERIC")
    resolved_sample = sample_id if sample_id is not None else sample

    rec = WarningRecord(
        timestamp=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        stage=resolved_stage,
        severity=severity,
        sample_id=resolved_sample,
        code=code,
        message=message,
        suggested_action=suggested_action,
    )
    _RECORDS.append(rec)
    log_warning(resolved_stage, f"[{code or 'GENERIC'}] {message}")
    return rec


def collected() -> list[WarningRecord]:
    """Return a shallow copy of every structured warning so far."""
    return list(_RECORDS)


def clear() -> None:
    """Reset the in-memory log (used by tests)."""
    _RECORDS.clear()


def _scrub_cell(value: str) -> str:
    return value.replace("\t", " ").replace("\n", " ")


def flush_tsv(path: Path | str) -> Path:
    """Write the accumulated warnings to a TSV at *path*.

    The file is overwritten on each call. An empty log still produces
    a header-only file so downstream consumers can rely on the path
    existing after a successful run.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(schema_header_line("warnings.tsv"))
        handle.write("\t".join(_WARNING_TSV_COLUMNS) + "\n")
        for r in _RECORDS:
            row = {
                "stage": r.stage,
                "sample_id": r.sample_id or "",
                "severity": r.severity,
                "code": r.code or "",
                "message": _scrub_cell(r.message),
                "suggested_action": _scrub_cell(r.suggested_action or ""),
            }
            handle.write(
                "\t".join(row[c] for c in _WARNING_TSV_COLUMNS) + "\n"
            )
    return path
