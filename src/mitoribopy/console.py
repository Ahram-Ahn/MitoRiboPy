"""Console logging and progress helpers for package CLI runs.

Logging style convention used across all MitoRiboPy subcommands:

    [<COMPONENT>] <message>

``COMPONENT`` is a short UPPERCASE tag naming the pipeline stage or
submodule that produced the line: ``ALIGN``, ``RPF``, ``ADAPTER``,
``PIPELINE``, ``PROFILE``, ``COVERAGE``, ``CODON``, ``CONTAM``,
``DEDUP``, ``VIS``, ``QC``, ``SUBSAMPLE``, ``COMPOSITION``, etc.

Rules (enforced by convention; no runtime check):

- **Never use bare** ``print(...)``, ``logging.info(...)``, ``sys.stderr.write``
  in package code. Every console line must flow through
  :func:`log_info`, :func:`log_warning`, :func:`log_error`, or
  :func:`log_progress` so output is uniformly decorated and also mirrored
  to the per-run ``mitoribopy.log`` file once
  :func:`configure_file_logging` has been called.
- **One component per log line.** Do not mix two tags like ``[ALIGN/RPF]``.
  If crossing a boundary, emit two lines.
- **Stages emit a ``Step K/N`` banner at each major step** via
  :func:`log_progress` or the pipeline's ``_emit_step_ok`` helper in
  :mod:`mitoribopy.pipeline.steps`. The banner format is:

      [PIPELINE] Step 3/7 OK: wrote unfiltered read-length QC outputs ...

- **Per-sample stages add a second-line detail** naming the sample:

      [ALIGN] sampleA: trim - adapter=TGGAATTCTCGGGTGCCAAGG, umi=0
      [ALIGN] sampleA: trim - kept 941,382/1,000,000 reads after cutadapt

- **WARNINGs start with "WARNING:"** (provided automatically by
  :func:`log_warning`), and ERRORs start with "ERROR:".

- **Progress bars** use :func:`format_progress_bar` via
  :func:`log_progress` / :func:`iter_with_progress`; do not invent new
  bar formats.

The one principled exception is CLI-time argparse validation errors,
which are allowed to go through ``print(..., file=sys.stderr)`` because
they happen before the file-logger is configured.
"""

from __future__ import annotations

import logging
import sys
from collections.abc import Callable, Iterator, Sequence
from pathlib import Path
from typing import TypeVar

import pandas as pd

_T = TypeVar("_T")
_LOGGER_NAME = "mitoribopy"
_FILE_HANDLER_NAME = "mitoribopy-file"


def _build_logger() -> logging.Logger:
    logger = logging.getLogger(_LOGGER_NAME)
    if logger.handlers:
        return logger

    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    logger.propagate = False
    return logger


LOGGER = _build_logger()


def configure_file_logging(log_path: str | Path) -> Path:
    """Attach a per-run file handler and return the resolved log path."""
    resolved_path = Path(log_path).expanduser().resolve()
    resolved_path.parent.mkdir(parents=True, exist_ok=True)

    for handler in list(LOGGER.handlers):
        if getattr(handler, "name", None) == _FILE_HANDLER_NAME:
            LOGGER.removeHandler(handler)
            handler.close()

    file_handler = logging.FileHandler(resolved_path, mode="w", encoding="utf-8")
    file_handler.name = _FILE_HANDLER_NAME
    file_handler.setFormatter(
        logging.Formatter(
            "%(asctime)s %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    )
    LOGGER.addHandler(file_handler)
    return resolved_path


def log_message(message: str) -> None:
    """Emit a raw log line without adding a component prefix."""
    LOGGER.info("%s", message)


def log_info(component: str, message: str) -> None:
    """Emit an informational console message."""
    LOGGER.info("[%s] %s", component, message)


def log_warning(component: str, message: str) -> None:
    """Emit a warning console message."""
    LOGGER.warning("[%s] WARNING: %s", component, message)


def log_error(component: str, message: str) -> None:
    """Emit an error console message."""
    LOGGER.error("[%s] ERROR: %s", component, message)


def format_progress_bar(current: int, total: int, width: int = 24) -> str:
    """Format a simple ASCII progress bar."""
    total = max(int(total), 1)
    current = min(max(int(current), 0), total)
    ratio = current / total
    filled = int(round(ratio * width))
    bar = "#" * filled + "-" * (width - filled)
    return f"[{bar}] {int(round(ratio * 100)):3d}%"


def log_progress(component: str, current: int, total: int, message: str) -> None:
    """Emit a progress-bar console message."""
    LOGGER.info("[%s] %s %s", component, format_progress_bar(current, total), message)


def iter_with_progress(
    items: Sequence[_T],
    *,
    component: str,
    noun: str,
    labeler: Callable[[_T], str] | None = None,
) -> Iterator[_T]:
    """Yield items while logging a simple per-item progress bar."""
    total = len(items)
    if total == 0:
        log_info(component, f"No {noun}s to process.")
        return

    for index, item in enumerate(items, start=1):
        detail = ""
        if labeler is not None:
            detail = f" ({labeler(item)})"
        log_progress(component, index, total, f"Processing {noun} {index}/{total}{detail}.")
        yield item


def log_dataframe_preview(
    component: str,
    title: str,
    dataframe: pd.DataFrame,
    *,
    max_rows: int = 12,
) -> None:
    """Emit a compact textual preview of a dataframe."""
    if dataframe.empty:
        log_info(component, f"{title}: no rows.")
        return

    preview = dataframe.head(max_rows)
    LOGGER.info("[%s] %s\n%s", component, title, preview.to_string(index=False))
    if len(dataframe) > max_rows:
        LOGGER.info(
            "[%s] Showing first %d of %d row(s).",
            component,
            max_rows,
            len(dataframe),
        )
