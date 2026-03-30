"""Read-count table parsing and normalization helpers."""

from __future__ import annotations

from typing import Iterable

import pandas as pd


def _resolve_column_name(
    columns: Iterable[str],
    requested_name: str | None,
    fallback_names: list[str],
    label: str,
) -> str:
    """Resolve a column by explicit name or from case-insensitive fallbacks."""
    cols = list(columns)
    lower_map = {str(col).lower(): col for col in cols}

    if requested_name:
        match = lower_map.get(str(requested_name).lower())
        if match is None:
            raise ValueError(
                f"{label} column '{requested_name}' not found. Available columns: {cols}"
            )
        return match

    for candidate in fallback_names:
        match = lower_map.get(str(candidate).lower())
        if match is not None:
            return match

    raise ValueError(
        f"Could not auto-detect {label} column. Available columns: {cols}. "
        "Please pass it explicitly."
    )


def _match_any_substring(value: object, patterns: Iterable[str]) -> bool:
    """Return True if any pattern appears in the provided value string."""
    text = str(value).lower()
    return any(str(pattern).lower() in text for pattern in patterns)


def compute_total_counts(
    count_summary_path: str,
    sample_col: str | None = None,
    reads_col: str | None = None,
    normalization_mode: str = "total",
    reference_col: str | None = None,
    mrna_ref_patterns: list[str] | None = None,
) -> tuple[dict[str, int], pd.DataFrame]:
    """Aggregate read totals per sample from a count-summary table."""
    if normalization_mode not in {"total", "mt_mrna"}:
        raise ValueError(
            f"Unsupported normalization_mode='{normalization_mode}'. "
            "Use 'total' or 'mt_mrna'."
        )

    if not mrna_ref_patterns:
        mrna_ref_patterns = ["mt_genome", "mt-mrna", "mt_mrna"]

    counts_df = pd.read_csv(count_summary_path, sep=None, engine="python")
    if counts_df.empty:
        print(f"[compute_total_counts] {count_summary_path} is empty or invalid.")
        return {}, pd.DataFrame()

    sample_col = _resolve_column_name(
        counts_df.columns,
        requested_name=sample_col,
        fallback_names=["Sample", "sample", "sample_name", "Sample_name"],
        label="sample",
    )
    reads_col = _resolve_column_name(
        counts_df.columns,
        requested_name=reads_col,
        fallback_names=["Reads", "reads", "read_count", "read_counts", "count"],
        label="read-count",
    )

    selected_cols = [sample_col, reads_col]
    if normalization_mode == "mt_mrna":
        reference_col = _resolve_column_name(
            counts_df.columns,
            requested_name=reference_col,
            fallback_names=["reference", "Reference", "ref", "target", "gene", "feature"],
            label="reference",
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
                lambda value: _match_any_substring(value, mrna_ref_patterns)
            )
        ].copy()

    work_df.dropna(subset=[sample_col, reads_col], inplace=True)
    if work_df.empty:
        print(
            "[compute_total_counts] No rows remain after applying "
            f"normalization_mode='{normalization_mode}'."
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

    print(
        "[compute_total_counts] Successfully computed read counts "
        f"(mode={normalization_mode}):"
    )
    if normalization_mode == "mt_mrna":
        print(
            "[compute_total_counts] Reference filter substrings: "
            + ", ".join(str(pattern) for pattern in mrna_ref_patterns)
        )
    print(total_counts_df)
    return total_counts_map, total_counts_df

