"""Process-wide structured warnings collector for MitoRiboPy.

P1.11 / P1.12: every warning emitted via :func:`record` is

* mirrored to stderr through :func:`mitoribopy.console.log_warning`,
* appended to the in-memory list returned by :func:`collected`,
* persisted to ``warnings.tsv`` at the run root by
  :func:`flush_tsv` (called by ``mitoribopy all`` after the per-stage
  work completes), and
* embedded in ``run_manifest.json`` under the ``warnings`` field
  (previously a placeholder empty list).

Records are intentionally append-only and process-global: subcommand
boundaries do NOT clear the list, so a later ``mitoribopy summarize``
or ``mitoribopy all`` invocation that re-runs the rnaseq stage can
collect warnings from every stage in a single TSV.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from .schema_versions import schema_header_line


__all__ = [
    "WarningRecord",
    "clear",
    "collected",
    "flush_tsv",
    "record",
]


_WARNING_TSV_COLUMNS: tuple[str, ...] = (
    "timestamp",
    "component",
    "severity",
    "sample",
    "code",
    "message",
)


@dataclass(frozen=True)
class WarningRecord:
    """One structured warning captured during a MitoRiboPy run."""

    timestamp: str
    component: str
    severity: str  # "warn" | "error"
    sample: str | None
    code: str | None
    message: str

    def as_dict(self) -> dict:
        return asdict(self)


_RECORDS: list[WarningRecord] = []


def record(
    component: str,
    message: str,
    *,
    sample: str | None = None,
    code: str | None = None,
    severity: str = "warn",
) -> WarningRecord:
    """Capture a structured warning AND mirror it to the console.

    The console mirror is intentional: existing user-facing behaviour
    of ``log_warning`` is preserved, and the structured TSV /
    manifest entries are additive.
    """
    # Local import to avoid a console <-> warnings_log circular at
    # module-init time.
    from ..console import log_warning

    rec = WarningRecord(
        timestamp=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        component=component,
        severity=severity,
        sample=sample,
        code=code,
        message=message,
    )
    _RECORDS.append(rec)
    if severity == "error":
        # Error severity still routes through log_warning here; the
        # caller is responsible for an explicit log_error / exit path
        # when appropriate. This keeps the structured log lossless
        # without changing exit-code semantics.
        log_warning(component, f"[{code or 'GENERIC'}] {message}")
    else:
        log_warning(component, f"[{code or 'GENERIC'}] {message}")
    return rec


def collected() -> list[WarningRecord]:
    """Return a shallow copy of every structured warning so far."""
    return list(_RECORDS)


def clear() -> None:
    """Reset the in-memory log (used by tests)."""
    _RECORDS.clear()


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
                "timestamp": r.timestamp,
                "component": r.component,
                "severity": r.severity,
                "sample": r.sample or "",
                "code": r.code or "",
                "message": r.message.replace("\t", " ").replace("\n", " "),
            }
            handle.write(
                "\t".join(row[c] for c in _WARNING_TSV_COLUMNS) + "\n"
            )
    return path
