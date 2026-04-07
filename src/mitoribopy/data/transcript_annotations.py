"""Built-in transcript annotations loaded from packaged CSV data."""

from __future__ import annotations

from .reference_data import load_annotation_table, normalize_annotation_table


yeast_annotation_df = load_annotation_table(preset="y")
human_annotation_df = load_annotation_table(preset="h")

__all__ = [
    "human_annotation_df",
    "yeast_annotation_df",
    "load_annotation_table",
    "normalize_annotation_table",
]
