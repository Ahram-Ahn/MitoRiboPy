"""Read-count table parsing and normalization helpers."""

from __future__ import annotations

import re
from typing import Iterable

import pandas as pd

from ..console import log_dataframe_preview, log_info, log_warning


def _normalize_column_name(name: object) -> str:
    """Normalize a column name for flexible matching."""
    return re.sub(r"[^a-z0-9]+", "", str(name).strip().lower())


def _resolve_column_name(
    columns: Iterable[str],
    requested_name: str | None,
    fallback_names: list[str],
    label: str,
    fallback_index: int | None = None,
) -> str:
    """Resolve a column by explicit name, flexible fallback names, or column order."""
    cols = list(columns)
    normalized_map = {_normalize_column_name(col): col for col in cols}

    if requested_name:
        match = normalized_map.get(_normalize_column_name(requested_name))
        if match is None:
            raise ValueError(
                f"{label} column '{requested_name}' not found. Available columns: {cols}"
            )
        return match

    for candidate in fallback_names:
        match = normalized_map.get(_normalize_column_name(candidate))
        if match is not None:
            return match

    if fallback_index is not None and 0 <= fallback_index < len(cols):
        chosen = cols[fallback_index]
        log_info(
            "COUNTS",
            f"Falling back to column index {fallback_index + 1} for {label}: '{chosen}'",
        )
        return chosen

    raise ValueError(
        f"Could not auto-detect {label} column. Available columns: {cols}. "
        "Please pass it explicitly."
    )


def _match_any_substring(value: object, patterns: Iterable[str]) -> bool:
    """Return True if any pattern appears in the provided value string."""
    text = str(value).lower()
    return any(str(pattern).lower() in text for pattern in patterns)


def _read_count_table(count_summary_path: str) -> pd.DataFrame:
    """Read a count-summary table from CSV, TSV, or generic delimited text."""
    counts_df = pd.read_csv(count_summary_path, sep=None, engine="python", comment="#")
    if len(counts_df.columns) == 1:
        counts_df = pd.read_csv(count_summary_path, sep=r"\s+", engine="python", comment="#")
    return counts_df


def compute_total_counts(
    count_summary_path: str,
    sample_col: str | None = None,
    reads_col: str | None = None,
    normalization_mode: str = "total",
    reference_col: str | None = None,
    mt_mrna_substring_patterns: list[str] | None = None,
    *,
    mrna_ref_patterns: list[str] | None = None,
) -> tuple[dict[str, int], pd.DataFrame]:
    """Aggregate read totals per sample from a count-summary table.

    The ``mrna_ref_patterns`` keyword is a deprecated alias for the
    canonical ``mt_mrna_substring_patterns`` and is accepted only so
    pinned third-party callers keep working; new code should use the
    new name.
    """
    if normalization_mode not in {"total", "mt_mrna"}:
        raise ValueError(
            f"Unsupported normalization_mode='{normalization_mode}'. "
            "Use 'total' or 'mt_mrna'."
        )

    if mrna_ref_patterns is not None and mt_mrna_substring_patterns is None:
        mt_mrna_substring_patterns = mrna_ref_patterns
    if not mt_mrna_substring_patterns:
        mt_mrna_substring_patterns = ["mt_genome", "mt-mrna", "mt_mrna"]

    counts_df = _read_count_table(count_summary_path)
    if counts_df.empty:
        log_warning("COUNTS", f"{count_summary_path} is empty or invalid.")
        return {}, pd.DataFrame()

    sample_col = _resolve_column_name(
        counts_df.columns,
        requested_name=sample_col,
        fallback_names=["Sample", "sample", "sample_name", "sampleid"],
        label="sample",
        fallback_index=0,
    )

    # Funnel-table fast path: the align stage's read_counts.tsv is one
    # row per sample with the mt-aligned post-dedup count already in a
    # named column. No reference filter needed — pull that column
    # directly. This is the layout `mitoribopy all` produces, so this
    # branch handles the orchestrated case.
    normalized_cols = {_normalize_column_name(c): c for c in counts_df.columns}
    funnel_columns = (
        "mtalignedafterdedup",  # mt_aligned_after_dedup
        "mtalignedaftermapq",   # mt_aligned_after_mapq (older layout)
        "mtaligned",            # mt_aligned (raw mt-tx alignments)
    )
    funnel_col = next(
        (normalized_cols[k] for k in funnel_columns if k in normalized_cols), None,
    )
    is_funnel_layout = funnel_col is not None and "reference" not in normalized_cols

    if is_funnel_layout:
        if normalization_mode == "mt_mrna":
            # Use the post-dedup mt-mRNA-aligned count as the RPM
            # denominator. Same biological meaning as
            # "sum(reads where reference matches mt-mRNA)" but read off
            # the funnel layout directly.
            selected_reads_col = funnel_col
        else:
            selected_reads_col = _resolve_column_name(
                counts_df.columns,
                requested_name=reads_col,
                fallback_names=["total_reads", "Reads", "reads", "read_count"],
                label="read-count",
                fallback_index=1,
            )
        work_df = counts_df[[sample_col, selected_reads_col]].copy()
        work_df.columns = ["Sample", "Total_reads"]
        work_df["Sample"] = (
            work_df["Sample"].astype(str).str.strip().str.replace(".bed", "", regex=False)
        )
        work_df["Total_reads"] = pd.to_numeric(work_df["Total_reads"], errors="coerce")
        work_df.dropna(inplace=True)
        if work_df.empty:
            log_warning(
                "COUNTS",
                f"Funnel layout: '{selected_reads_col}' column had no usable rows.",
            )
            return {}, pd.DataFrame(columns=["Sample", "Total_reads"])
        work_df["Total_reads"] = work_df["Total_reads"].astype(int)
        total_counts_map: dict[str, int] = {}
        for sample_name, total_reads in zip(work_df["Sample"], work_df["Total_reads"]):
            total_counts_map[str(sample_name)] = int(total_reads)
            total_counts_map[str(sample_name).lower()] = int(total_reads)
        log_info(
            "COUNTS",
            f"Funnel layout: using '{selected_reads_col}' column as the "
            f"RPM denominator (mode={normalization_mode}).",
        )
        log_dataframe_preview("COUNTS", "Read-total summary", work_df)
        return total_counts_map, work_df

    reads_col = _resolve_column_name(
        counts_df.columns,
        requested_name=reads_col,
        fallback_names=["Reads", "reads", "read_count", "read_counts", "count", "counts"],
        label="read-count",
        fallback_index=(2 if len(counts_df.columns) >= 3 else 1),
    )

    selected_cols = [sample_col, reads_col]
    if normalization_mode == "mt_mrna":
        reference_col = _resolve_column_name(
            counts_df.columns,
            requested_name=reference_col,
            fallback_names=["reference", "ref", "target", "gene", "feature"],
            label="reference",
            fallback_index=1,
        )
        selected_cols.append(reference_col)

    work_df = counts_df[selected_cols].copy()
    work_df[sample_col] = (
        work_df[sample_col].astype(str).str.strip().str.replace(".bed", "", regex=False)
    )
    work_df[reads_col] = pd.to_numeric(work_df[reads_col], errors="coerce")

    if normalization_mode == "mt_mrna":
        work_df = work_df[
            work_df[reference_col].astype(str).apply(
                lambda value: _match_any_substring(
                    value, mt_mrna_substring_patterns
                )
            )
        ].copy()

    work_df.dropna(subset=[sample_col, reads_col], inplace=True)
    if work_df.empty:
        log_warning(
            "COUNTS",
            "No rows remain after applying "
            f"normalization_mode='{normalization_mode}'.",
        )
        return {}, pd.DataFrame(columns=["Sample", "Total_reads"])

    total_counts_df = (
        work_df.groupby(sample_col, dropna=False)[reads_col]
        .sum()
        .reset_index(name="Total_reads")
    )
    total_counts_df.rename(columns={sample_col: "Sample"}, inplace=True)
    total_counts_df["Total_reads"] = total_counts_df["Total_reads"].astype(int)

    total_counts_map: dict[str, int] = {}
    for sample_name, total_reads in zip(
        total_counts_df["Sample"], total_counts_df["Total_reads"]
    ):
        total_counts_map[str(sample_name)] = int(total_reads)
        total_counts_map[str(sample_name).lower()] = int(total_reads)

    log_info("COUNTS", f"Computed read totals successfully (mode={normalization_mode}).")
    if normalization_mode == "mt_mrna":
        log_info(
            "COUNTS",
            "Reference filter substrings: "
            + ", ".join(str(pattern) for pattern in mt_mrna_substring_patterns),
        )
    log_dataframe_preview("COUNTS", "Read-total summary", total_counts_df)
    return total_counts_map, total_counts_df
