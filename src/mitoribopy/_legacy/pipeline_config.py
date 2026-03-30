"""Legacy import wrapper for Phase II config module migration."""

from __future__ import annotations

import warnings


warnings.warn(
    "mitoribopy._legacy.pipeline_config is deprecated; use mitoribopy.config runtime helpers.",
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
