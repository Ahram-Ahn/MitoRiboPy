"""Console logging and progress helpers for package CLI runs."""

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
