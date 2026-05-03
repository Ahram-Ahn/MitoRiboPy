"""Reference data tables used by MitoRiboPy analyses."""

from .reference_data import (
    available_codon_table_names,
    annotation_sequence_candidates,
    build_sequence_display_map,
    load_annotation_table,
    load_codon_table,
    load_codon_tables,
    resolve_sequence_name,
    resolve_start_codons,
    transcript_display_title,
)
from .codon_tables import (
    human_mitochondrial_codon_table,
    standard_codon_table,
    yeast_mitochondrial_codon_table,
)
from .transcript_annotations import human_annotation_df, yeast_annotation_df

import json
from pathlib import Path


_RUN_MANIFEST_SCHEMA_PATH = Path(__file__).parent / "run_manifest.schema.json"


def load_run_manifest_schema() -> dict:
    """Return the JSON Schema (draft 2020-12) for ``run_manifest.json``.

    The schema lives next to the package's other data files so it ships
    with the wheel and can be looked up without a network round-trip.
    Use :func:`mitoribopy.data.validate_run_manifest` for one-shot
    validation against this schema.
    """
    return json.loads(_RUN_MANIFEST_SCHEMA_PATH.read_text(encoding="utf-8"))


def validate_run_manifest(manifest: dict) -> list[str]:
    """Validate a parsed ``run_manifest.json`` against the package schema.

    Returns a list of human-readable error strings; empty list means
    the manifest is valid. Uses ``jsonschema`` when available; falls
    back to a minimal hand-written required-key + type check otherwise
    so the package has no new hard dependency.
    """
    schema = load_run_manifest_schema()
    try:
        import jsonschema  # type: ignore
    except ImportError:
        return _fallback_manifest_validation(manifest, schema)
    validator = jsonschema.Draft202012Validator(schema)
    return [
        f"{'/'.join(str(p) for p in err.absolute_path)}: {err.message}"
        for err in sorted(validator.iter_errors(manifest), key=lambda e: e.path)
    ]


def _fallback_manifest_validation(manifest: dict, schema: dict) -> list[str]:
    """Minimal stand-in for jsonschema: required-key + top-level type checks.

    The fallback covers the publication-grade reproducibility minimum
    (top-level required keys + their primitive types). It does NOT
    validate nested structures or pattern constraints; install
    ``jsonschema`` for full coverage.
    """
    errors: list[str] = []
    if not isinstance(manifest, dict):
        return ["root: manifest is not a JSON object"]
    type_map = {
        "string": (str,),
        "object": (dict,),
        "array": (list,),
        "number": (int, float),
        "integer": (int,),
        "boolean": (bool,),
        "null": (type(None),),
    }
    for required in schema.get("required", []):
        if required not in manifest:
            errors.append(f"{required}: required field is missing")
            continue
        spec = schema.get("properties", {}).get(required, {})
        types = spec.get("type")
        if types is None:
            continue
        type_list = types if isinstance(types, list) else [types]
        py_types = tuple(t for s in type_list for t in type_map.get(s, ()))
        if py_types and not isinstance(manifest[required], py_types):
            errors.append(
                f"{required}: expected {types}, got {type(manifest[required]).__name__}"
            )
    return errors


__all__ = [
    "annotation_sequence_candidates",
    "available_codon_table_names",
    "build_sequence_display_map",
    "load_annotation_table",
    "load_codon_table",
    "load_codon_tables",
    "load_run_manifest_schema",
    "resolve_sequence_name",
    "resolve_start_codons",
    "standard_codon_table",
    "transcript_display_title",
    "validate_run_manifest",
    "yeast_mitochondrial_codon_table",
    "human_mitochondrial_codon_table",
    "yeast_annotation_df",
    "human_annotation_df",
]
