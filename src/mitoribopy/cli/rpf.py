"""``mitoribopy rpf`` subcommand.

Runs the Ribo-seq analysis pipeline from BED or BAM inputs. Internally
delegates to :func:`mitoribopy.pipeline.runner.run_pipeline_cli`; this
module peels the common flags (``--dry-run``, ``--threads``,
``--log-level``) off the argv before handing the remainder to the
pipeline parser so that ``--config``, ``-f``, ``-s``, ``-rpf``, and
the rest of the rpf flags continue to be owned by ``pipeline.runner``.
"""

from __future__ import annotations

import sys
from typing import Iterable

from . import common


RPF_SUBCOMMAND_HELP = (
    "Run the Ribo-seq analysis pipeline from BED or BAM inputs. "
    "See 'mitoribopy rpf --help' for the full flag list."
)


def run(argv: Iterable[str]) -> int:
    """Entry point for ``mitoribopy rpf <args>``.

    *argv* is the list of tokens that follow the ``rpf`` subcommand word.
    """
    # Avoid a top-level import of mitoribopy.cli to prevent an import cycle;
    # importing inside run() also lets tests monkeypatch
    # ``mitoribopy.cli.run_pipeline_cli`` between calls.
    from mitoribopy import cli as _cli_pkg

    argv_list = list(argv)

    # Peel the common --dry-run / --threads / --log-level flags; leave
    # --config in place because the rpf parser in pipeline.runner already
    # handles it (and now supports YAML/TOML via
    # config.runtime.load_user_config).
    common_ns, remaining = common.peel_common_arguments(argv_list)
    common.apply_common_arguments(common_ns)

    if common_ns.dry_run:
        return common.emit_dry_run(
            "rpf",
            [
                "load --config (if given) via mitoribopy.config.load_user_config",
                "initialize outputs, annotations, and RPF range",
                "run unfiltered read-length QC and filter BED/BAM inputs",
                "compute offset enrichment and select P-site/A-site offsets",
                "run translation-profile analysis, coverage plots, and "
                "optional structure-density / codon-correlation modules",
            ],
        )

    try:
        return _cli_pkg.run_pipeline_cli(remaining)
    except RuntimeError as exc:
        print(f"[mitoribopy rpf] {exc}", file=sys.stderr)
        return 2
