"""Typed configuration models for package migration."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PhaseOneConfig:
    """Minimal config model for package identity in Phase I."""

    package_name: str = "MitoRiboPy"
    import_name: str = "mitoribopy"

