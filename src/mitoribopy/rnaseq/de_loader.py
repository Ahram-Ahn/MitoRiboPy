"""DE-table loader with format auto-detection.

Supports DESeq2, Xtail, and Anota2Seq output shapes directly. For
unknown schemas the caller supplies an explicit column mapping via
``--de-gene-col / --de-log2fc-col / --de-padj-col / --de-basemean-col``.

File formats: CSV, TSV, or TXT. The delimiter is auto-detected by
sniffing the first data line.
"""

from __future__ import annotations

import csv
from pathlib import Path

from ._types import DE_COLUMN_ALIASES, DeColumnMap, DeFormat, DeTable


def _detect_delimiter(sample_line: str) -> str:
    """Return ``\\t`` if the header line has more tabs than commas, else ``,``."""
    tab_count = sample_line.count("\t")
    comma_count = sample_line.count(",")
    return "\t" if tab_count >= comma_count and tab_count > 0 else ","


def detect_de_format(header: list[str]) -> tuple[DeFormat, DeColumnMap]:
    """Identify the DE-table schema from its column headers.

    Detection is keyed on the log2FC column name because that is the
    most distinctive column across DESeq2 / Xtail / Anota2Seq.

    Returns
    -------
    (format, column_map):
        ``format`` is the literal for the matched schema, or ``"custom"``
        when no recognized pattern fires.
        ``column_map`` has its ``gene_id`` field overwritten to match
        the ACTUAL first column of the table, so downstream lookup
        works regardless of how the user named the gene-id column.
    """
    header_set = set(header)
    gene_id_col = header[0] if header else "gene_id"

    for fmt, alias in DE_COLUMN_ALIASES.items():
        if alias.log2fc and alias.log2fc in header_set:
            return fmt, DeColumnMap(
                gene_id=gene_id_col,
                log2fc=alias.log2fc if alias.log2fc in header_set else None,
                padj=alias.padj if (alias.padj and alias.padj in header_set) else None,
                basemean=(
                    alias.basemean
                    if (alias.basemean and alias.basemean in header_set)
                    else None
                ),
            )
    return "custom", DeColumnMap(
        gene_id=gene_id_col, log2fc=None, padj=None, basemean=None
    )


def _safe_float(value: str) -> float | None:
    if value in ("", "NA", "NaN", "na", "nan", None):
        return None
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    # DESeq2 / Xtail emit 'Inf' sometimes; keep as None for numerical sanity.
    if f != f or f == float("inf") or f == float("-inf"):
        return None
    return f


def load_de_table(
    path: Path,
    *,
    column_map: DeColumnMap | None = None,
) -> DeTable:
    """Load a DE table from *path* into a canonical :class:`DeTable`.

    Parameters
    ----------
    path:
        CSV / TSV / TXT. The delimiter is auto-detected from the header
        line.
    column_map:
        Optional explicit column mapping. When omitted,
        :func:`detect_de_format` is called on the header.

    Returns
    -------
    :class:`DeTable`
        with ``rows`` as a list of dicts keyed by canonical names.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"DE table not found: {path}")

    text = path.read_text(encoding="utf-8")
    if not text.strip():
        raise ValueError(f"DE table is empty: {path}")

    first_line = text.splitlines()[0]
    delimiter = _detect_delimiter(first_line)

    reader = csv.DictReader(text.splitlines(), delimiter=delimiter)
    header = list(reader.fieldnames or [])
    if not header:
        raise ValueError(f"DE table has no header row: {path}")

    if column_map is None:
        fmt, column_map = detect_de_format(header)
    else:
        # Honor user's explicit column map; still need a format tag for
        # the manifest. Re-run detect to tag it, but keep the user's map.
        fmt, _ = detect_de_format(header)

    canonical_rows: list[dict] = []
    for raw in reader:
        gene_id = (raw.get(column_map.gene_id) or "").strip()
        if not gene_id:
            continue
        canonical_rows.append(
            {
                "gene_id": gene_id,
                "log2fc": _safe_float(raw.get(column_map.log2fc))
                if column_map.log2fc
                else None,
                "padj": _safe_float(raw.get(column_map.padj))
                if column_map.padj
                else None,
                "basemean": _safe_float(raw.get(column_map.basemean))
                if column_map.basemean
                else None,
            }
        )

    return DeTable(format=fmt, column_map=column_map, rows=canonical_rows)
