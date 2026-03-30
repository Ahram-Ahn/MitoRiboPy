"""Legacy import wrapper for Phase II BED subsampling tool migration."""

from __future__ import annotations

import warnings


warnings.warn(
    "mitoribopy._legacy.subsample_bed is deprecated; use mitoribopy.tools.subsample.",
    DeprecationWarning,
    stacklevel=2,
)

from mitoribopy.tools.subsample import (  # noqa: E402,F401
    main,
    parse_args,
    reservoir_sample_lines,
)

__all__ = ["parse_args", "reservoir_sample_lines", "main"]
