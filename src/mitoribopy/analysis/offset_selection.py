"""Offset selection utilities for read-length specific P/A-site offsets."""

from __future__ import annotations

import pandas as pd

from ..console import log_dataframe_preview, log_info, log_warning


def _validate_offset_mask_nt(offset_mask_nt: int) -> int:
    """Validate the near-anchor masking window."""
    offset_mask_nt = int(offset_mask_nt)
    if offset_mask_nt < 0:
        raise ValueError("offset_mask_nt must be zero or a positive integer.")
    return offset_mask_nt


def _convert_offsets_between_sites(
    offsets_df: pd.DataFrame, from_site: str, to_site: str
) -> pd.DataFrame:
    """Convert 5' and 3' offsets between P-site and A-site conventions."""
    from_site = str(from_site).lower()
    to_site = str(to_site).lower()
    if from_site == to_site:
        return offsets_df.copy()
    if from_site not in {"p", "a"} or to_site not in {"p", "a"}:
        raise ValueError(f"Unsupported site conversion: {from_site} -> {to_site}")

    converted = offsets_df.copy()
    if from_site == "p" and to_site == "a":
        converted["5' Offset"] = converted["5' Offset"] + 3
        converted["3' Offset"] = converted["3' Offset"] - 3
    elif from_site == "a" and to_site == "p":
        converted["5' Offset"] = converted["5' Offset"] - 3
        converted["3' Offset"] = converted["3' Offset"] + 3
    return converted


def _pick_most_enriched_with_tiebreak(
    counts_df: pd.DataFrame, offset_col: str
) -> pd.DataFrame:
    """Pick one enriched offset per read length and return diagnostic columns.

    Returns a DataFrame with one row per read length and these columns:

    * ``Read Length``
    * ``<offset_col>``           — the selected offset value (the picked one)
    * ``n_reads``                — total reads in the input window for this length
    * ``top_count``              — read count at the selected offset
    * ``second_best_count``      — read count at the runner-up offset (0 when none)
    * ``second_best_offset``     — the runner-up offset value (NaN when none)
    * ``delta_score``            — top_count - second_best_count (raw separation)
    * ``enrichment_score``       — top_count / n_reads (fraction of mass on the pick)

    Selection still uses the deterministic neighbour-support tie-break:
    ties on ``Count`` are broken by the sum of counts at offsets
    ``[x-1, x, x+1]`` (highest wins), then by ``abs(offset)``, then by
    raw offset value.
    """
    import math

    picked_rows = []
    for read_length, sub in counts_df.groupby("Read Length", sort=True):
        sub = sub.sort_values([offset_col]).copy()
        n_reads = int(sub["Count"].sum())
        top_count = int(sub["Count"].max()) if not sub.empty else 0
        tied = sub[sub["Count"] == top_count].copy()

        count_map = dict(zip(sub[offset_col], sub["Count"]))
        tied["NeighborSupport"] = tied[offset_col].apply(
            lambda x: count_map.get(x - 1, 0) + count_map.get(x, 0) + count_map.get(x + 1, 0)
        )
        tied["AbsOffset"] = tied[offset_col].abs()
        tied = tied.sort_values(
            ["NeighborSupport", "AbsOffset", offset_col],
            ascending=[False, True, True],
        )
        best_value = tied.iloc[0][offset_col]

        # Second-best by raw count, ignoring the chosen offset.
        runner = sub[sub[offset_col] != best_value]
        if runner.empty:
            second_best_count = 0
            second_best_offset = math.nan
        else:
            runner_top = runner["Count"].max()
            second_best_count = int(runner_top)
            # When the runner-up is itself a multi-offset tie, prefer the
            # one closest to the picked offset for a stable diagnostic.
            runner_top_rows = runner[runner["Count"] == runner_top]
            runner_top_rows = runner_top_rows.assign(
                _absdist=(runner_top_rows[offset_col] - best_value).abs()
            ).sort_values(["_absdist", offset_col])
            second_best_offset = runner_top_rows.iloc[0][offset_col]

        delta_score = top_count - second_best_count
        enrichment_score = (top_count / n_reads) if n_reads > 0 else 0.0

        picked_rows.append({
            "Read Length": read_length,
            offset_col: best_value,
            "n_reads": n_reads,
            "top_count": top_count,
            "second_best_count": second_best_count,
            "second_best_offset": second_best_offset,
            "delta_score": delta_score,
            "enrichment_score": enrichment_score,
        })

    return pd.DataFrame(picked_rows)


# Confidence-label thresholds. These are explicit numbers so the
# defensibility argument lives in code, not in tribal knowledge.
#
# enrichment_score = top_count / n_reads (fraction of mass at the pick).
# delta_ratio      = (top_count - second_best_count) / top_count
#                    (relative separation from the runner-up).
# n_reads          = total reads in the per-length window.
#
# A "high" call requires enough reads for a stable estimate, a clear
# enrichment peak, and meaningful separation from the runner-up. A
# "medium" call relaxes one of those. Anything weaker is "low".
_CONF_HIGH_ENRICH = 0.40
_CONF_HIGH_DELTA_RATIO = 0.30
_CONF_HIGH_MIN_READS = 200

_CONF_MEDIUM_ENRICH = 0.25
_CONF_MEDIUM_MIN_READS = 50


def _confidence_label(
    *, n_reads: int, top_count: int, second_best_count: int, enrichment_score: float
) -> str:
    """Classify a per-length offset pick as high / medium / low / insufficient.

    Special-case labels (``manual`` and ``fallback_combined``) are
    written by callers when those modes apply and bypass this function.
    """
    if n_reads <= 0 or top_count <= 0:
        return "insufficient"
    delta_ratio = (top_count - second_best_count) / top_count if top_count else 0.0
    if (
        enrichment_score >= _CONF_HIGH_ENRICH
        and delta_ratio >= _CONF_HIGH_DELTA_RATIO
        and n_reads >= _CONF_HIGH_MIN_READS
    ):
        return "high"
    if (
        enrichment_score >= _CONF_MEDIUM_ENRICH
        and n_reads >= _CONF_MEDIUM_MIN_READS
    ):
        return "medium"
    return "low"


def determine_p_site_offsets(
    offsets_df: pd.DataFrame,
    align_to: str,
    out_file: str,
    offset_min: int = 11,
    offset_max: int = 20,
    five_offset_min: int | None = None,
    five_offset_max: int | None = None,
    three_offset_min: int | None = None,
    three_offset_max: int | None = None,
    offset_mask_nt: int = 5,
    offset_site: str = "p",
    selection_reference: str = "p_site",
) -> pd.DataFrame | None:
    """Select most-enriched offsets per read length and write them to CSV."""
    if offsets_df.empty:
        return None

    offset_site = str(offset_site).lower()
    if offset_site not in {"p", "a"}:
        raise ValueError(f"Unsupported offset_site='{offset_site}'. Use 'p' or 'a'.")

    selection_reference = str(selection_reference).lower()
    # 'selected_site' is a deprecated alias for 'reported_site' (kept so
    # old YAML configs still parse). Canonicalise here so the downstream
    # checks only see the new name.
    if selection_reference == "selected_site":
        selection_reference = "reported_site"
    if selection_reference not in {"reported_site", "p_site"}:
        raise ValueError(
            "Unsupported selection_reference='"
            + selection_reference
            + "'. Use 'reported_site' or 'p_site' "
            "(legacy alias 'selected_site' also accepted)."
        )
    offset_mask_nt = _validate_offset_mask_nt(offset_mask_nt)
    five_offset_min = offset_min if five_offset_min is None else int(five_offset_min)
    five_offset_max = offset_max if five_offset_max is None else int(five_offset_max)
    three_offset_min = offset_min if three_offset_min is None else int(three_offset_min)
    three_offset_max = offset_max if three_offset_max is None else int(three_offset_max)

    reported_df = offsets_df.copy()
    working_df = reported_df.copy()
    if selection_reference == "p_site" and offset_site == "a":
        # Pick in canonical P-site space, then convert chosen offsets back to A-site.
        working_df = _convert_offsets_between_sites(reported_df, from_site="a", to_site="p")

    five_mask = (
        (reported_df["5' Offset"].abs() >= five_offset_min)
        & (reported_df["5' Offset"].abs() <= five_offset_max)
        & (reported_df["5' Offset"].abs() > offset_mask_nt)
    )
    three_mask = (
        (reported_df["3' Offset"].abs() >= three_offset_min)
        & (reported_df["3' Offset"].abs() <= three_offset_max)
        & (reported_df["3' Offset"].abs() > offset_mask_nt)
    )

    five_df = working_df.loc[five_mask].copy()
    three_df = working_df.loc[three_mask].copy()

    if five_df.empty and three_df.empty:
        log_warning("OFFSET", "No offsets were found inside the requested selection ranges.")
        return None

    five_diag_cols = [
        "n_reads_5", "top_count_5", "second_best_count_5",
        "second_best_offset_5", "delta_score_5", "enrichment_score_5",
        "confidence_5",
    ]
    three_diag_cols = [
        "n_reads_3", "top_count_3", "second_best_count_3",
        "second_best_offset_3", "delta_score_3", "enrichment_score_3",
        "confidence_3",
    ]

    if five_df.empty:
        five_pick = pd.DataFrame(
            columns=["Read Length", "Most Enriched 5' Offset"] + five_diag_cols
        )
    else:
        five_counts = (
            five_df.groupby(["Read Length", "5' Offset"]).size().reset_index(name="Count")
        )
        five_pick = _pick_most_enriched_with_tiebreak(five_counts, "5' Offset")
        five_pick["confidence_5"] = five_pick.apply(
            lambda r: _confidence_label(
                n_reads=int(r["n_reads"]),
                top_count=int(r["top_count"]),
                second_best_count=int(r["second_best_count"]),
                enrichment_score=float(r["enrichment_score"]),
            ),
            axis=1,
        )
        five_pick = five_pick.rename(columns={
            "5' Offset": "Most Enriched 5' Offset",
            "n_reads": "n_reads_5",
            "top_count": "top_count_5",
            "second_best_count": "second_best_count_5",
            "second_best_offset": "second_best_offset_5",
            "delta_score": "delta_score_5",
            "enrichment_score": "enrichment_score_5",
        })

    if three_df.empty:
        three_pick = pd.DataFrame(
            columns=["Read Length", "Most Enriched 3' Offset"] + three_diag_cols
        )
    else:
        three_counts = (
            three_df.groupby(["Read Length", "3' Offset"]).size().reset_index(name="Count")
        )
        three_pick = _pick_most_enriched_with_tiebreak(three_counts, "3' Offset")
        three_pick["confidence_3"] = three_pick.apply(
            lambda r: _confidence_label(
                n_reads=int(r["n_reads"]),
                top_count=int(r["top_count"]),
                second_best_count=int(r["second_best_count"]),
                enrichment_score=float(r["enrichment_score"]),
            ),
            axis=1,
        )
        three_pick = three_pick.rename(columns={
            "3' Offset": "Most Enriched 3' Offset",
            "n_reads": "n_reads_3",
            "top_count": "top_count_3",
            "second_best_count": "second_best_count_3",
            "second_best_offset": "second_best_offset_3",
            "delta_score": "delta_score_3",
            "enrichment_score": "enrichment_score_3",
        })

    selected_offsets = pd.merge(
        five_pick[["Read Length", "Most Enriched 5' Offset"] + five_diag_cols],
        three_pick[["Read Length", "Most Enriched 3' Offset"] + three_diag_cols],
        on="Read Length",
        how="outer",
    )

    if selection_reference == "p_site" and offset_site == "a":
        selected_offsets["Most Enriched 5' Offset"] = (
            selected_offsets["Most Enriched 5' Offset"] + 3
        )
        selected_offsets["Most Enriched 3' Offset"] = (
            selected_offsets["Most Enriched 3' Offset"] - 3
        )

    selected_offsets.to_csv(out_file, index=False)
    log_info(
        "OFFSET",
        "Selected offsets saved => "
        f"{out_file} (align_to={align_to}, offset_site={offset_site}, "
        f"selection_reference={selection_reference})",
    )
    log_dataframe_preview("OFFSET", "Selected offsets by read length", selected_offsets)
    return selected_offsets
