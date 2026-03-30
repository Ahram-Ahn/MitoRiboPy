"""Offset selection utilities for read-length specific P/A-site offsets."""

from __future__ import annotations

import pandas as pd


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
    """Pick one enriched offset per read length using deterministic tie-break rules."""
    picked_rows = []
    for read_length, sub in counts_df.groupby("Read Length", sort=True):
        sub = sub.sort_values([offset_col]).copy()
        top_count = sub["Count"].max()
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
        picked_rows.append({"Read Length": read_length, offset_col: best_value})

    return pd.DataFrame(picked_rows)


def determine_p_site_offsets(
    offsets_df: pd.DataFrame,
    align_to: str,
    out_file: str,
    offset_min: int = 11,
    offset_max: int = 20,
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
    if selection_reference not in {"selected_site", "p_site"}:
        raise ValueError(
            "Unsupported selection_reference='"
            + selection_reference
            + "'. Use 'selected_site' or 'p_site'."
        )

    working_df = offsets_df.copy()
    if selection_reference == "p_site" and offset_site == "a":
        # Pick in canonical P-site space, then convert chosen offsets back to A-site.
        working_df = _convert_offsets_between_sites(working_df, from_site="a", to_site="p")

    five_df = working_df[
        (working_df["5' Offset"].abs() >= offset_min)
        & (working_df["5' Offset"].abs() <= offset_max)
    ]
    three_df = working_df[
        (working_df["3' Offset"].abs() >= offset_min)
        & (working_df["3' Offset"].abs() <= offset_max)
    ]

    if five_df.empty and three_df.empty:
        print("[determine_p_site_offsets] No offsets found in range => skipping.")
        return None

    five_counts = (
        five_df.groupby(["Read Length", "5' Offset"]).size().reset_index(name="Count")
    )
    five_pick = _pick_most_enriched_with_tiebreak(five_counts, "5' Offset")
    five_pick.rename(columns={"5' Offset": "Most Enriched 5' Offset"}, inplace=True)

    three_counts = (
        three_df.groupby(["Read Length", "3' Offset"]).size().reset_index(name="Count")
    )
    three_pick = _pick_most_enriched_with_tiebreak(three_counts, "3' Offset")
    three_pick.rename(columns={"3' Offset": "Most Enriched 3' Offset"}, inplace=True)

    selected_offsets = pd.merge(
        five_pick[["Read Length", "Most Enriched 5' Offset"]],
        three_pick[["Read Length", "Most Enriched 3' Offset"]],
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
    print(
        "[determine_p_site_offsets] Selected offsets => "
        f"{out_file} (align_to={align_to}, offset_site={offset_site}, "
        f"selection_reference={selection_reference})"
    )
    return selected_offsets

