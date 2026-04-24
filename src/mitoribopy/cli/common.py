"""Shared CLI helpers for MitoRiboPy subcommands.

All MitoRiboPy subcommands (``align``, ``rpf``, ``rnaseq``, ``all``) accept
the following common flags:

* ``--config PATH``     JSON, YAML, or TOML configuration file.
* ``--dry-run``         Print planned actions and exit without running.
* ``--threads N``       Preferred thread count for external tools and
                        libraries that honor ``OMP_NUM_THREADS`` and friends.
* ``--log-level LEVEL`` Python logging level for ``mitoribopy`` console output.

This module exposes:

* :func:`add_common_arguments` to register those flags on a subparser.
* :func:`apply_common_arguments` to apply the parsed values.
* :func:`peel_common_arguments` to strip the common flags out of an argv
  list before handing the remainder to a subcommand's own parser
  (used by ``rpf`` which reuses :mod:`mitoribopy.pipeline.runner`).
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from typing import Any, Iterable

from ..console import LOGGER, log_warning

_LOG_LEVEL_CHOICES = ("DEBUG", "INFO", "WARNING", "ERROR")


class MitoRiboPyHelpFormatter(argparse.RawTextHelpFormatter):
    """Consistent help formatter with readable defaults."""

    def __init__(self, *args, **kwargs) -> None:
        kwargs.setdefault("max_help_position", 38)
        super().__init__(*args, **kwargs)

    def _get_help_string(self, action) -> str:
        help_text = action.help or ""
        if "%(default)" in help_text:
            return help_text
        if not action.option_strings or action.required:
            return help_text

        default = getattr(action, "default_display", action.default)
        if default in (None, False, argparse.SUPPRESS):
            return help_text
        if isinstance(default, (list, tuple)) and len(default) == 0:
            return help_text

        return f"{help_text} [default: {default}]"


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    """Register --config, --dry-run, --threads, and --log-level on *parser*."""
    common = parser.add_argument_group("Shared options")
    common.add_argument(
        "--config",
        default=None,
        metavar="CONFIG",
        help="Configuration file (.json, .yaml, .yml, or .toml). "
        "CLI arguments override values read from the file.",
    )
    common.add_argument(
        "--dry-run",
        action="store_true",
        default=False,
        help="Print planned actions and exit without executing.",
    )
    common.add_argument(
        "--threads",
        type=int,
        default=None,
        metavar="N",
        help="Preferred thread count for external tools and BLAS libraries.",
    )
    common.add_argument(
        "--log-level",
        choices=_LOG_LEVEL_CHOICES,
        default="INFO",
        help="Python logging level for MitoRiboPy console output.",
    )


def apply_common_arguments(args: argparse.Namespace) -> None:
    """Apply the effects of the common flags parsed from *args*.

    * ``--log-level`` sets the level on the ``mitoribopy`` logger.
    * ``--threads`` exports ``OMP_NUM_THREADS``, ``OPENBLAS_NUM_THREADS``,
      and ``MKL_NUM_THREADS`` for BLAS/OpenMP libraries and stores the
      requested count under ``MITORIBOPY_THREADS`` so subcommand code
      (for example ``align``) can forward it to external tools.
    """
    level_name = getattr(args, "log_level", "INFO") or "INFO"
    LOGGER.setLevel(getattr(logging, level_name, logging.INFO))

    threads = getattr(args, "threads", None)
    if threads is not None:
        if threads < 1:
            raise SystemExit("--threads must be a positive integer")
        os.environ["MITORIBOPY_THREADS"] = str(threads)
        for env_key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            os.environ[env_key] = str(threads)


def peel_common_arguments(argv: Iterable[str]) -> tuple[argparse.Namespace, list[str]]:
    """Peel the common flags out of *argv* without requiring --help support.

    Returns ``(ns, remaining)`` where ``ns`` has attributes ``config``,
    ``dry_run``, ``threads``, ``log_level`` and ``remaining`` is the list
    of tokens that were not consumed. Used by subcommands that delegate
    the bulk of their argument parsing to another module (``rpf`` uses
    this to hand the remainder to :mod:`mitoribopy.pipeline.runner`).
    """
    parser = argparse.ArgumentParser(add_help=False)
    add_common_arguments(parser)
    ns, remaining = parser.parse_known_args(list(argv))
    return ns, remaining


def load_config_file(path: str | None) -> dict[str, Any]:
    """Load a configuration file, auto-detected from its extension.

    Supported extensions:

    * ``.json``           - stdlib :mod:`json`
    * ``.yaml``, ``.yml`` - :mod:`yaml` (PyYAML, required)
    * ``.toml``           - stdlib :mod:`tomllib` (Python 3.11+) or
                            :mod:`tomli` (Python 3.10 fallback)

    Returns ``{}`` if *path* is falsy. Unknown extensions fall back to
    JSON for backward compatibility with v0.2.x.
    """
    if not path:
        return {}

    config_path = Path(path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    suffix = config_path.suffix.lower()
    text = config_path.read_text(encoding="utf-8")

    if suffix in {".yaml", ".yml"}:
        try:
            import yaml
        except ImportError as exc:  # pragma: no cover - PyYAML is a core dep
            raise RuntimeError(
                "PyYAML is required to read YAML config files. "
                "Install with: pip install 'mitoribopy[dev]' or pip install PyYAML."
            ) from exc
        data = yaml.safe_load(text) or {}
    elif suffix == ".toml":
        try:
            import tomllib  # Python 3.11+
        except ImportError:  # pragma: no cover - Python 3.10 fallback
            try:
                import tomli as tomllib  # type: ignore[no-redef]
            except ImportError as exc:
                raise RuntimeError(
                    "Reading TOML config files on Python < 3.11 requires the "
                    "'tomli' package. Install with: pip install tomli."
                ) from exc
        data = tomllib.loads(text)
    else:
        data = json.loads(text) if text.strip() else {}

    if not isinstance(data, dict):
        raise ValueError(
            f"Config file must parse to a mapping/object; got {type(data).__name__}."
        )
    return data


def emit_dry_run(component: str, planned_actions: list[str]) -> int:
    """Print a dry-run plan and return 0.

    Used by subcommand handlers when ``--dry-run`` is set.
    """
    print(f"[{component}] dry-run: planned actions", file=sys.stdout)
    for index, action in enumerate(planned_actions, start=1):
        print(f"  {index}. {action}", file=sys.stdout)
    return 0


def warn_stub_subcommand(subcommand: str, planned_phase: str) -> None:
    """Warn the user that a subcommand is a stub for a future phase."""
    log_warning(
        subcommand.upper(),
        f"{subcommand} is not yet implemented; planned for {planned_phase}.",
    )
