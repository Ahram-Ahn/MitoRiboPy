"""Built-in codon tables loaded from packaged JSON data."""

from __future__ import annotations

from .reference_data import available_codon_table_names, load_codon_table, load_codon_tables


standard_codon_table = load_codon_table(table_name="standard")
yeast_mitochondrial_codon_table = load_codon_table(table_name="yeast_mitochondrial")
human_mitochondrial_codon_table = load_codon_table(table_name="vertebrate_mitochondrial")

__all__ = [
    "available_codon_table_names",
    "load_codon_table",
    "load_codon_tables",
    "standard_codon_table",
    "yeast_mitochondrial_codon_table",
    "human_mitochondrial_codon_table",
]
