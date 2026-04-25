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


def test_canonical_strain_resolves_short_aliases() -> None:
    from mitoribopy.data.reference_data import (
        BUILTIN_STRAIN_PRESETS,
        STRAIN_ALIASES,
        canonical_strain,
    )

    assert STRAIN_ALIASES == {"h": "h.sapiens", "y": "s.cerevisiae"}
    assert canonical_strain("h") == "h.sapiens"
    assert canonical_strain("y") == "s.cerevisiae"
    assert canonical_strain("h.sapiens") == "h.sapiens"
    assert canonical_strain("custom") == "custom"
    assert BUILTIN_STRAIN_PRESETS == frozenset({"h.sapiens", "s.cerevisiae"})


def test_load_annotation_table_canonicalises_short_alias() -> None:
    """A pre-existing pipeline_config.yaml that still uses `strain: h`
    must still load the built-in human annotation."""
    df_short = load_annotation_table(preset="h")
    df_long = load_annotation_table(preset="h.sapiens")
    assert list(df_short["transcript"]) == list(df_long["transcript"])


def test_load_annotation_table_rejects_dropped_strain() -> None:
    """`vm` and `ym` no longer exist; they are not in STRAIN_ALIASES so
    they fall through to the 'no built-in annotation' error path."""
    import pytest

    with pytest.raises(ValueError, match="built-in annotation table"):
        load_annotation_table(preset="vm")
    with pytest.raises(ValueError, match="built-in annotation table"):
        load_annotation_table(preset="ym")


def test_load_codon_table_picks_default_for_short_strain_alias() -> None:
    """The default codon-table picker must canonicalise the short alias
    so `preset='h'` still picks vertebrate_mitochondrial."""
    table_short = load_codon_table(preset="h")
    table_long = load_codon_table(preset="h.sapiens")
    assert table_short == table_long
    assert table_short["AGA"] == "*"  # vertebrate_mitochondrial signature
