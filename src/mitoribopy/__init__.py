"""MitoRiboPy package."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("MitoRiboPy")
except PackageNotFoundError:  # pragma: no cover - local editable use before install
    __version__ = "0.4.4"
