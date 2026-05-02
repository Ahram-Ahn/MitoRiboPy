"""Single canonicalisation entry point for parsed YAML configs.

The package historically had three different code paths that did
overlapping work on a parsed config dict:

1. :func:`mitoribopy.config.migrate.migrate` — legacy key rewrites.
2. :func:`mitoribopy.config.runtime.load_user_config` — unknown-key
   warnings (flat configs only).
3. :mod:`mitoribopy.cli.validate_config` — path / cross-section checks.

That split made it easy for a new legacy alias to land in (1) without
the validator (3) seeing it.

:func:`canonicalize_config` is the one place every CLI surface should
go through to obtain a canonical config dict plus a structured change
log. It is intentionally side-effect-free: no warnings are emitted to
the console, no files are touched. Callers that want to surface the
change log do so via the returned :class:`ConfigChange` records.

This module does NOT replace :mod:`migrate` — it composes it. Anything
that should rewrite legacy keys still lands in ``migrate.py`` so the
``mitoribopy migrate-config`` subcommand continues to work standalone.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .migrate import migrate as _migrate_legacy_keys


__all__ = [
    "ConfigChange",
    "CanonicalConfig",
    "canonicalize_config",
]


@dataclass(frozen=True)
class ConfigChange:
    """One rewrite or warning emitted by canonicalisation.

    ``severity`` is ``"info"`` for routine legacy-key rewrites and
    ``"warn"`` for things the user probably wants to know about (e.g.
    unknown top-level keys that were dropped).
    """

    severity: str
    message: str


@dataclass(frozen=True)
class CanonicalConfig:
    """The result of :func:`canonicalize_config`.

    Attributes
    ----------
    config
        The canonical config dict, suitable for handing off to the
        per-stage parsers / dataclasses.
    changes
        Ordered list of every rewrite / warning that fired during
        canonicalisation. Empty when the input was already canonical.
    """

    config: dict[str, Any]
    changes: list[ConfigChange] = field(default_factory=list)

    @property
    def has_warnings(self) -> bool:
        return any(c.severity == "warn" for c in self.changes)

    def info_messages(self) -> list[str]:
        return [c.message for c in self.changes if c.severity == "info"]

    def warning_messages(self) -> list[str]:
        return [c.message for c in self.changes if c.severity == "warn"]


# Sections we recognise at the top of a sectioned (orchestrator) config.
_KNOWN_TOP_LEVEL_SECTIONS: frozenset[str] = frozenset(
    {
        "align",
        "rpf",
        "rnaseq",
        "samples",     # unified sample sheet pointer.
        "execution",   # v0.6.2: top-level resource-plan declaration.
        "periodicity", # spec-defined periodicity QC tuning block.
        "all",         # orchestrator-level overrides.
        "shared",      # cross-stage shared knobs (bam_mapq, threads, ...).
        "resume",      # resume-guard knobs.
    }
)


def canonicalize_config(raw_config: dict[str, Any]) -> CanonicalConfig:
    """Return the canonical form of *raw_config* plus a change log.

    The pipeline is intentionally short — every step is documented:

    1. Legacy-key migration via :func:`mitoribopy.config.migrate.migrate`.
    2. Unknown top-level section detection (sectioned configs only).

    Per-stage validation (path existence, mutually-exclusive flags,
    rnaseq mode resolution) lives in :mod:`mitoribopy.cli.validate_config`
    and runs AFTER this function.
    """
    if not isinstance(raw_config, dict):
        raise TypeError(
            f"canonicalize_config expected dict, got {type(raw_config).__name__}"
        )

    canonical, log = _migrate_legacy_keys(raw_config)
    changes: list[ConfigChange] = [
        ConfigChange(severity="info", message=line) for line in log
    ]

    # Heuristic: a sectioned config has at least one of align/rpf/rnaseq
    # at the top AS A DICT. Flat per-stage templates can also use those
    # names as scalar keys (e.g. ``rpf: [29, 34]`` for the RPF length
    # window in a per-stage rpf config) — those are NOT section
    # wrappers and must NOT trigger the unknown-section check.
    is_sectioned = any(
        isinstance(canonical.get(k), dict)
        for k in ("align", "rpf", "rnaseq")
    )
    if is_sectioned:
        unknown = sorted(
            k for k in canonical.keys() if k not in _KNOWN_TOP_LEVEL_SECTIONS
        )
        for key in unknown:
            changes.append(
                ConfigChange(
                    severity="warn",
                    message=(
                        f"unknown top-level section {key!r} — recognised "
                        f"sections are {sorted(_KNOWN_TOP_LEVEL_SECTIONS)}"
                    ),
                )
            )

    return CanonicalConfig(config=canonical, changes=changes)
