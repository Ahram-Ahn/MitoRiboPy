"""Config loader scaffolding for package migration."""

from __future__ import annotations

from .models import PhaseOneConfig


def load_phase_one_config() -> PhaseOneConfig:
    """Return package identity config used during Phase I."""
    return PhaseOneConfig()

