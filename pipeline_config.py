"""Backward-compatible wrapper for legacy `pipeline_config` imports.

Phase II note:
The implementation moved to:
- ``mitoribopy.config.runtime``
"""

from __future__ import annotations

from pathlib import Path
import sys
import warnings


# Keep `python main.py ...` working from repository root without requiring install.
_src_dir = Path(__file__).resolve().parent / "src"
if _src_dir.exists() and str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))


warnings.warn(
    "pipeline_config is deprecated; import from mitoribopy.config instead.",
    DeprecationWarning,
    stacklevel=2,
)

from mitoribopy.config.runtime import (  # noqa: E402,F401
    DEFAULT_CONFIG,
    load_user_config,
    resolve_rpf_range,
)

__all__ = [
    "DEFAULT_CONFIG",
    "load_user_config",
    "resolve_rpf_range",
]
