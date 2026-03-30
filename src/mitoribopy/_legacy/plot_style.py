"""Legacy import wrapper for Phase II plotting-style migration."""

from __future__ import annotations

import warnings


warnings.warn(
    "mitoribopy._legacy.plot_style is deprecated; use mitoribopy.plotting.style.",
    DeprecationWarning,
    stacklevel=2,
)

from mitoribopy.plotting.style import (  # noqa: E402,F401
    PUBLICATION_STYLE,
    apply_publication_style,
)

__all__ = ["PUBLICATION_STYLE", "apply_publication_style"]
