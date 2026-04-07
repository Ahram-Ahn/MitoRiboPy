from __future__ import annotations

from mitoribopy.data import (
    available_codon_table_names,
    human_annotation_df,
    load_annotation_table,
    load_codon_table,
    resolve_sequence_name,
)


def test_packaged_codon_tables_are_loaded_from_json() -> None:
    assert "standard" in available_codon_table_names()
    assert "vertebrate_mitochondrial" in available_codon_table_names()
    assert load_codon_table(table_name="vertebrate_mitochondrial")["AGA"] == "*"


def test_builtin_human_annotation_splits_bicistronic_pairs_with_default_baselines() -> None:
    assert {"ATP8", "ATP6", "ND4L", "ND4"}.issubset(set(human_annotation_df["transcript"]))
    assert "ATP86" not in set(human_annotation_df["transcript"])
    assert "ND4L4" not in set(human_annotation_df["transcript"])

    atp_rows = human_annotation_df[human_annotation_df["transcript"].isin(["ATP8", "ATP6"])]
    nd4_rows = human_annotation_df[human_annotation_df["transcript"].isin(["ND4L", "ND4"])]

    assert set(atp_rows["display_name"]) == {"ATP8/ATP6"}
    assert set(atp_rows["sequence_name"]) == {"ATP6"}
    assert set(nd4_rows["display_name"]) == {"ND4L/ND4"}
    assert set(nd4_rows["sequence_name"]) == {"ND4"}
    assert all("ATP86" in alias_text for alias_text in atp_rows["sequence_aliases"])
    assert all("ND4L4" in alias_text for alias_text in nd4_rows["sequence_aliases"])


def test_builtin_annotation_baselines_can_be_overridden() -> None:
    annotation_df = load_annotation_table(
        preset="h",
        atp8_atp6_baseline="ATP8",
        nd4l_nd4_baseline="ND4L",
    )

    atp_rows = annotation_df[annotation_df["transcript"].isin(["ATP8", "ATP6"])]
    nd4_rows = annotation_df[annotation_df["transcript"].isin(["ND4L", "ND4"])]

    assert set(atp_rows["sequence_name"]) == {"ATP8"}
    assert set(nd4_rows["sequence_name"]) == {"ND4L"}


def test_resolve_sequence_name_uses_legacy_bicistronic_aliases() -> None:
    atp8_row = human_annotation_df[human_annotation_df["transcript"] == "ATP8"].iloc[0]
    nd4_row = human_annotation_df[human_annotation_df["transcript"] == "ND4L"].iloc[0]

    assert resolve_sequence_name(atp8_row, {"ATP86"}) == "ATP86"
    assert resolve_sequence_name(nd4_row, {"ND4L4"}) == "ND4L4"
