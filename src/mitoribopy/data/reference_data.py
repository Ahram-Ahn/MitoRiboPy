"""Load built-in and user-supplied codon tables and transcript annotations."""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from importlib.resources import files
from pathlib import Path

import pandas as pd


DEFAULT_CODON_TABLE_BY_STRAIN = {
    # Built-in presets ship annotation + codon table out of the box.
    "y": "yeast_mitochondrial",
    "h": "vertebrate_mitochondrial",
    # "Any ___ mitochondria" presets share the codon table with their
    # anchor strain (vm = any vertebrate mito, ym = any yeast-mito-code
    # fungus). They require a user-supplied --annotation_file.
    "vm": "vertebrate_mitochondrial",
    "ym": "yeast_mitochondrial",
}

DEFAULT_START_CODONS_BY_STRAIN = {
    "y": ["ATG"],
    "h": ["ATG", "ATA"],
    "vm": ["ATG", "ATA", "ATT"],  # vertebrate mito: MT-ND2 uses ATT; be permissive
    "ym": ["ATG"],
    "custom": ["ATG"],
}

BUILTIN_ANNOTATION_FILES = {
    "y": "yeast_annotation.csv",
    "h": "human_annotation.csv",
    # vm / ym deliberately absent: user must supply --annotation_file.
}


# Presets that ship a ready-to-use annotation + RPF range with no user
# input required. Anything not in this set must have --annotation_file
# and an explicit -rpf MIN MAX.
BUILTIN_ANNOTATION_PRESETS = frozenset(BUILTIN_ANNOTATION_FILES)

BICISTRONIC_CONFIGS = (
    {
        "members": ("ATP8", "ATP6"),
        "label": "ATP8/ATP6",
        "default_baseline": "ATP6",
        "legacy_aliases": ["ATP86"],
    },
    {
        "members": ("ND4L", "ND4"),
        "label": "ND4L/ND4",
        "default_baseline": "ND4",
        "legacy_aliases": ["ND4L4"],
    },
)


def _package_data_path(filename: str) -> Path:
    return Path(files("mitoribopy.data").joinpath(filename))


def _normalize_codon_table(table: Mapping[str, object]) -> dict[str, str]:
    normalized: dict[str, str] = {}
    for codon, amino_acid in table.items():
        normalized[str(codon).upper()] = str(amino_acid).upper()
    return normalized


def load_codon_tables(table_file: str | None = None) -> dict[str, dict[str, str]]:
    """Load codon tables from package data or a user-supplied JSON file."""
    table_path = _package_data_path("codon_tables.json") if table_file is None else Path(table_file)
    with table_path.open("r", encoding="utf-8") as handle:
        raw = json.load(handle)

    if not isinstance(raw, dict):
        raise ValueError("Codon-table JSON must be an object.")

    if raw and all(isinstance(value, str) for value in raw.values()):
        return {"default": _normalize_codon_table(raw)}

    tables: dict[str, dict[str, str]] = {}
    for table_name, table_value in raw.items():
        if not isinstance(table_value, Mapping):
            raise ValueError(
                "Each named codon table must map codons to amino-acid codes."
            )
        tables[str(table_name)] = _normalize_codon_table(table_value)
    return tables


def load_codon_table(
    *,
    preset: str | None = None,
    table_name: str | None = None,
    table_file: str | None = None,
) -> dict[str, str]:
    """Load one codon table from built-ins or a user-supplied JSON file."""
    tables = load_codon_tables(table_file)

    resolved_table_name = table_name
    if resolved_table_name is None:
        if table_file and len(tables) == 1:
            return next(iter(tables.values()))
        if preset in DEFAULT_CODON_TABLE_BY_STRAIN:
            resolved_table_name = DEFAULT_CODON_TABLE_BY_STRAIN[preset]
        elif "standard" in tables:
            resolved_table_name = "standard"
        elif len(tables) == 1:
            return next(iter(tables.values()))
        else:
            raise ValueError(
                "Unable to choose a codon table automatically. "
                "Provide --codon_table_name."
            )

    if resolved_table_name not in tables:
        available = ", ".join(sorted(tables))
        raise ValueError(
            f"Codon table '{resolved_table_name}' was not found. Available tables: {available}"
        )
    return tables[resolved_table_name]


def available_codon_table_names(table_file: str | None = None) -> list[str]:
    """Return the available codon-table names."""
    return sorted(load_codon_tables(table_file))


def resolve_start_codons(strain: str, start_codons: Iterable[str] | None = None) -> list[str]:
    """Resolve the start-codon set from CLI input or preset defaults."""
    if start_codons:
        seen: set[str] = set()
        normalized: list[str] = []
        for codon in start_codons:
            codon_upper = str(codon).upper()
            if codon_upper not in seen:
                normalized.append(codon_upper)
                seen.add(codon_upper)
        return normalized
    return list(DEFAULT_START_CODONS_BY_STRAIN.get(strain, ["ATG"]))


def _split_aliases(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    text = str(value).strip()
    if not text:
        return []

    aliases: list[str] = []
    for chunk in text.replace(",", ";").split(";"):
        alias = chunk.strip()
        if alias:
            aliases.append(alias)
    return aliases


def _unique_preserving_order(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        if not value:
            continue
        if value in seen:
            continue
        seen.add(value)
        result.append(value)
    return result


def _pair_label_variants(label: str, legacy_aliases: Iterable[str]) -> list[str]:
    return _unique_preserving_order(
        [
            label,
            label.replace("/", ""),
            label.replace("/", "_"),
            label.replace("/", "-"),
            *legacy_aliases,
        ]
    )


def _apply_bicistronic_defaults(
    annotation_df: pd.DataFrame,
    *,
    atp8_atp6_baseline: str,
    nd4l_nd4_baseline: str,
) -> pd.DataFrame:
    baseline_by_pair = {
        ("ATP8", "ATP6"): atp8_atp6_baseline,
        ("ND4L", "ND4"): nd4l_nd4_baseline,
    }

    for config in BICISTRONIC_CONFIGS:
        members = config["members"]
        label = str(config["label"])
        member_mask = annotation_df["transcript"].isin(members)
        if member_mask.sum() != len(members):
            continue

        baseline = baseline_by_pair.get(members, config["default_baseline"])
        pair_aliases = _pair_label_variants(label, config["legacy_aliases"])

        for member in members:
            row_mask = annotation_df["transcript"] == member
            sequence_name = annotation_df.loc[row_mask, "sequence_name"].astype(str).str.strip()
            missing_sequence = sequence_name.eq("") | sequence_name.eq("nan")
            if missing_sequence.all():
                annotation_df.loc[row_mask, "sequence_name"] = baseline

            display_name = annotation_df.loc[row_mask, "display_name"].astype(str).str.strip()
            missing_display = display_name.eq("") | display_name.eq("nan")
            if missing_display.all():
                annotation_df.loc[row_mask, "display_name"] = label

            merged_aliases = _unique_preserving_order(
                [
                    *_split_aliases(annotation_df.loc[row_mask, "sequence_aliases"].iloc[0]),
                    baseline,
                    *members,
                    *pair_aliases,
                ]
            )
            annotation_df.loc[row_mask, "sequence_aliases"] = ";".join(merged_aliases)

    return annotation_df


def normalize_annotation_table(
    annotation_df: pd.DataFrame,
    *,
    atp8_atp6_baseline: str = "ATP6",
    nd4l_nd4_baseline: str = "ND4",
) -> pd.DataFrame:
    """Validate and normalize an annotation table for package use."""
    required_columns = {"transcript", "l_tr", "l_utr5", "l_utr3"}
    missing_columns = sorted(required_columns - set(annotation_df.columns))
    if missing_columns:
        raise ValueError(
            "Annotation CSV is missing required columns: " + ", ".join(missing_columns)
        )

    normalized = annotation_df.copy()
    normalized["transcript"] = normalized["transcript"].astype(str).str.strip()

    for column in ("l_tr", "l_utr5", "l_utr3"):
        normalized[column] = pd.to_numeric(normalized[column], errors="raise").astype(int)

    if "l_cds" not in normalized.columns:
        normalized["l_cds"] = normalized["l_tr"] - normalized["l_utr5"] - normalized["l_utr3"]
    else:
        normalized["l_cds"] = pd.to_numeric(normalized["l_cds"], errors="coerce")
        missing_cds = normalized["l_cds"].isna()
        normalized.loc[missing_cds, "l_cds"] = (
            normalized.loc[missing_cds, "l_tr"]
            - normalized.loc[missing_cds, "l_utr5"]
            - normalized.loc[missing_cds, "l_utr3"]
        )
        normalized["l_cds"] = normalized["l_cds"].astype(int)

    for column in ("sequence_name", "sequence_aliases", "display_name"):
        if column not in normalized.columns:
            normalized[column] = ""
        normalized[column] = normalized[column].fillna("").astype(str).str.strip()

    normalized = _apply_bicistronic_defaults(
        normalized,
        atp8_atp6_baseline=atp8_atp6_baseline,
        nd4l_nd4_baseline=nd4l_nd4_baseline,
    )

    missing_sequence = normalized["sequence_name"].eq("")
    normalized.loc[missing_sequence, "sequence_name"] = normalized.loc[missing_sequence, "transcript"]
    missing_display = normalized["display_name"].eq("")
    normalized.loc[missing_display, "display_name"] = normalized.loc[missing_display, "transcript"]

    normalized["start_codon"] = normalized["l_utr5"].astype(int)
    normalized["stop_codon"] = (
        normalized["l_tr"].astype(int) - normalized["l_utr3"].astype(int) - 3
    )

    ordered_columns = [
        "transcript",
        "sequence_name",
        "sequence_aliases",
        "display_name",
        "l_tr",
        "l_utr5",
        "l_cds",
        "l_utr3",
        "start_codon",
        "stop_codon",
    ]
    remaining_columns = [
        column for column in normalized.columns if column not in ordered_columns
    ]
    return normalized[ordered_columns + remaining_columns]


def load_annotation_table(
    *,
    preset: str | None = None,
    annotation_file: str | None = None,
    atp8_atp6_baseline: str = "ATP6",
    nd4l_nd4_baseline: str = "ND4",
) -> pd.DataFrame:
    """Load an annotation CSV from built-ins or from a user-provided path."""
    if annotation_file is None:
        if preset not in BUILTIN_ANNOTATION_FILES:
            raise ValueError(
                "A built-in annotation table is only available for strain presets 'y' and 'h'."
            )
        annotation_path = _package_data_path(BUILTIN_ANNOTATION_FILES[preset])
    else:
        annotation_path = Path(annotation_file)

    annotation_df = pd.read_csv(annotation_path)
    return normalize_annotation_table(
        annotation_df,
        atp8_atp6_baseline=atp8_atp6_baseline,
        nd4l_nd4_baseline=nd4l_nd4_baseline,
    )


def annotation_sequence_candidates(annotation_row: pd.Series | Mapping[str, object]) -> list[str]:
    """Return candidate FASTA/BED sequence IDs for one annotation row."""
    if isinstance(annotation_row, pd.Series):
        row = annotation_row.to_dict()
    else:
        row = dict(annotation_row)

    display_name = str(row.get("display_name", "")).strip()
    display_variants = []
    if display_name:
        display_variants = _pair_label_variants(display_name, [])

    candidates = _unique_preserving_order(
        [
            str(row.get("sequence_name", "")).strip(),
            str(row.get("transcript", "")).strip(),
            *display_variants,
            *_split_aliases(row.get("sequence_aliases", "")),
        ]
    )
    return candidates


def resolve_sequence_name(
    annotation_row: pd.Series | Mapping[str, object],
    available_sequence_names: Iterable[str],
) -> str | None:
    """Resolve one annotation row onto an available FASTA/BED sequence name."""
    available = list(available_sequence_names)
    exact_map = {name: name for name in available}
    lowered_map = {name.lower(): name for name in available}

    for candidate in annotation_sequence_candidates(annotation_row):
        if candidate in exact_map:
            return exact_map[candidate]
        lowered = candidate.lower()
        if lowered in lowered_map:
            return lowered_map[lowered]
    return None


def build_sequence_display_map(
    annotation_df: pd.DataFrame,
    available_sequence_names: Iterable[str],
) -> dict[str, str]:
    """Build a display-label map keyed by resolved sequence ID."""
    display_map: dict[str, list[str]] = {}
    for _, row in annotation_df.iterrows():
        resolved_name = resolve_sequence_name(row, available_sequence_names)
        if resolved_name is None:
            continue
        label = str(row.get("display_name", "")).strip() or str(row["transcript"]).strip()
        display_map.setdefault(resolved_name, [])
        if label not in display_map[resolved_name]:
            display_map[resolved_name].append(label)

    return {
        sequence_name: " / ".join(labels) if len(labels) > 1 else labels[0]
        for sequence_name, labels in display_map.items()
    }


def transcript_display_title(
    annotation_row: pd.Series | Mapping[str, object],
    *,
    include_transcript: bool = False,
) -> str:
    """Return a stable display title for one logical transcript row."""
    if isinstance(annotation_row, pd.Series):
        row = annotation_row.to_dict()
    else:
        row = dict(annotation_row)

    transcript = str(row.get("transcript", "")).strip()
    display_name = str(row.get("display_name", "")).strip() or transcript
    if include_transcript and display_name != transcript:
        return f"{display_name} [{transcript}]"
    return display_name
