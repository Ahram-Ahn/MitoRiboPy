from __future__ import annotations

import itertools

from mitoribopy.data.codon_tables import (
    human_mitochondrial_codon_table,
    standard_codon_table,
    yeast_mitochondrial_codon_table,
)


def _all_codons() -> set[str]:
    return {"".join(parts) for parts in itertools.product("ATGC", repeat=3)}


def test_all_codon_tables_have_full_64_codon_space() -> None:
    expected = _all_codons()
    for table in (
        standard_codon_table,
        yeast_mitochondrial_codon_table,
        human_mitochondrial_codon_table,
    ):
        assert set(table.keys()) == expected
        assert len(table) == 64


def test_mito_specific_reassignments_are_present() -> None:
    assert standard_codon_table["TGA"] == "*"
    assert yeast_mitochondrial_codon_table["TGA"] == "W"
    assert human_mitochondrial_codon_table["TGA"] == "W"
    assert standard_codon_table["AGA"] == "R"
    assert human_mitochondrial_codon_table["AGA"] == "*"
