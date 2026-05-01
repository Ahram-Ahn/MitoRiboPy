"""Per-run summary aggregation (P1.6, P1.7, P1.8).

* :func:`build_summary_qc` aggregates per-sample QC metrics from the
  align / rpf / rnaseq stage outputs into a single roll-up.
* :func:`write_summary_qc` writes the aggregated rows as
  ``summary_qc.tsv`` with the standard schema-version header.
* :func:`render_summary_md` produces the human-readable ``SUMMARY.md``
  for a finished run.
"""

from .qc import (
    SUMMARY_QC_COLUMNS,
    build_summary_qc,
    write_summary_qc,
)
from .render import render_summary_md


__all__ = [
    "SUMMARY_QC_COLUMNS",
    "build_summary_qc",
    "render_summary_md",
    "write_summary_qc",
]
