"""Loader for per-sample per-gene RPF read counts from a prior rpf run.

The rpf pipeline emits ``<output>/rpf_counts.tsv`` with one row per
(sample, gene) and a ``count`` column. This module reads that file into
the structure the TE math wants: a nested dict
``{gene: {sample: count}}``.
"""

from __future__ import annotations

import csv
from pathlib import Path


def load_ribo_counts(path: Path) -> dict[str, dict[str, int]]:
    """Load ``rpf_counts.tsv`` into ``{gene: {sample: count}}``.

    The file is expected to be tab-delimited with a header row
    containing at least the columns ``sample``, ``gene``, ``count``.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"Ribo-seq counts file not found: {path}. "
            "This file is emitted by a MitoRiboPy >=0.3.0 rpf run as "
            "<output>/rpf_counts.tsv."
        )

    counts: dict[str, dict[str, int]] = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [col for col in ("sample", "gene", "count") if col not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(
                f"Ribo-seq counts file {path} is missing required column(s): "
                + ", ".join(missing)
            )
        for row in reader:
            sample = (row["sample"] or "").strip()
            gene = (row["gene"] or "").strip()
            try:
                count = int(row["count"])
            except (TypeError, ValueError):
                continue
            if not sample or not gene:
                continue
            counts.setdefault(gene, {})[sample] = count
    return counts
