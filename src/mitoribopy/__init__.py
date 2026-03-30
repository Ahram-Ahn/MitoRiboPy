"""MitoRiboPy package.

Phase I package scaffold:
- keeps legacy top-level modules compatible
- introduces a package entrypoint for future refactoring
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("MitoRiboPy")
except PackageNotFoundError:  # pragma: no cover - local editable use before install
    __version__ = "0.1.0"

