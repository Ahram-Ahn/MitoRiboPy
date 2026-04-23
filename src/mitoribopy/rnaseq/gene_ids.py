"""Gene-ID convention handling for mt-mRNAs.

Maps the user's chosen ``--gene-id-convention`` to the set of
identifiers we expect to see in the DE table, and matches them against
the mt-mRNA gene set defined in the annotation CSVs.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class MtGene:
    """One mt-mRNA's IDs across every supported convention.

    The human registry below is checked against NCBI / Ensembl / HGNC.
    Add a new vendor convention by extending this dataclass and
    populating the corresponding field on each gene.
    """

    hgnc: str  # canonical 'MT-...' symbol; identical to mt_prefixed
    bare: str  # without the 'MT-' prefix
    ensembl: str | None = None
    refseq: str | None = None


# Human mt-mRNAs, checked against Ensembl GRCh38.p14 / NCBI 2024.
# Order matches the canonical genomic layout so reviewers can spot
# omissions.
HUMAN_MT_MRNAS: tuple[MtGene, ...] = (
    MtGene(hgnc="MT-ND1", bare="ND1", ensembl="ENSG00000198888", refseq="YP_003024026.1"),
    MtGene(hgnc="MT-ND2", bare="ND2", ensembl="ENSG00000198763", refseq="YP_003024027.1"),
    MtGene(hgnc="MT-CO1", bare="CO1", ensembl="ENSG00000198804", refseq="YP_003024028.1"),
    MtGene(hgnc="MT-CO2", bare="CO2", ensembl="ENSG00000198712", refseq="YP_003024029.1"),
    MtGene(hgnc="MT-ATP8", bare="ATP8", ensembl="ENSG00000228253", refseq="YP_003024030.1"),
    MtGene(hgnc="MT-ATP6", bare="ATP6", ensembl="ENSG00000198899", refseq="YP_003024031.1"),
    MtGene(hgnc="MT-CO3", bare="CO3", ensembl="ENSG00000198938", refseq="YP_003024032.1"),
    MtGene(hgnc="MT-ND3", bare="ND3", ensembl="ENSG00000198840", refseq="YP_003024033.1"),
    MtGene(hgnc="MT-ND4L", bare="ND4L", ensembl="ENSG00000212907", refseq="YP_003024034.1"),
    MtGene(hgnc="MT-ND4", bare="ND4", ensembl="ENSG00000198886", refseq="YP_003024035.1"),
    MtGene(hgnc="MT-ND5", bare="ND5", ensembl="ENSG00000198786", refseq="YP_003024036.1"),
    MtGene(hgnc="MT-ND6", bare="ND6", ensembl="ENSG00000198695", refseq="YP_003024037.1"),
    MtGene(hgnc="MT-CYB", bare="CYB", ensembl="ENSG00000198727", refseq="YP_003024038.1"),
)


# Yeast mtDNA protein-coding genes. S. cerevisiae has 8 mt-mRNAs (COX1,
# COX2, COX3, COB, ATP6, ATP8, ATP9, VAR1). Yeast has no 'MT-' naming
# convention nor Ensembl Gene IDs for mt-DNA; we record the SGD-style
# bare symbols.
YEAST_MT_MRNAS: tuple[MtGene, ...] = (
    MtGene(hgnc="COX1", bare="COX1"),
    MtGene(hgnc="COX2", bare="COX2"),
    MtGene(hgnc="COX3", bare="COX3"),
    MtGene(hgnc="COB", bare="COB"),
    MtGene(hgnc="ATP6", bare="ATP6"),
    MtGene(hgnc="ATP8", bare="ATP8"),
    MtGene(hgnc="ATP9", bare="ATP9"),
    MtGene(hgnc="VAR1", bare="VAR1"),
)


def _registry_for_organism(organism: str) -> tuple[MtGene, ...]:
    organism = (organism or "").lower()
    if organism in {"y", "yeast", "s_cerevisiae"}:
        return YEAST_MT_MRNAS
    if organism in {"h", "human", "homo_sapiens"}:
        return HUMAN_MT_MRNAS
    raise ValueError(
        f"Unknown organism for mt-mRNA registry: {organism!r}. "
        "Supported: 'h', 'human', 'y', 'yeast'."
    )


def _id_for_convention(gene: MtGene, convention: str) -> str | None:
    if convention == "hgnc":
        return gene.hgnc
    if convention == "mt_prefixed":
        return gene.hgnc
    if convention == "bare":
        return gene.bare
    if convention == "ensembl":
        return gene.ensembl
    if convention == "refseq":
        return gene.refseq
    raise ValueError(
        f"Unknown --gene-id-convention: {convention!r}. "
        "Supported: ensembl, refseq, hgnc, mt_prefixed, bare."
    )


def normalize_gene_id(raw_id: str) -> str:
    """Case-normalize and strip whitespace; does NOT change convention.

    Used before matching so 'mt-nd1' equals 'MT-ND1'.
    """
    return (raw_id or "").strip().upper()


def expected_ids(convention: str, organism: str = "h") -> list[str]:
    """Return the expected gene IDs in *convention* for *organism*.

    Genes whose registry entry is missing that ID are skipped (e.g.,
    yeast has no Ensembl IDs for mtDNA so ``expected_ids('ensembl',
    'y')`` returns an empty list).
    """
    registry = _registry_for_organism(organism)
    ids: list[str] = []
    for gene in registry:
        candidate = _id_for_convention(gene, convention)
        if candidate:
            ids.append(candidate)
    return ids


def match_mt_mrnas(
    de_gene_ids: list[str],
    convention: str,
    organism: str = "h",
) -> dict[str, list[str]]:
    """Return the intersection and the missing gene list.

    Returns a dict with two keys:

    * ``matched``  sorted list of IDs from *de_gene_ids* that correspond
                   to known mt-mRNAs in the chosen convention.
    * ``missing``  sorted list of expected mt-mRNA IDs for *organism*
                   that are NOT present in *de_gene_ids*. Used to emit
                   the 'fewer than 13 mt-mRNAs matched' warning.
    """
    normalized = {normalize_gene_id(g) for g in de_gene_ids}
    expected = expected_ids(convention, organism)
    expected_norm = {normalize_gene_id(g): g for g in expected}

    matched = sorted(g for g in expected if normalize_gene_id(g) in normalized)
    missing = sorted(g for g in expected if normalize_gene_id(g) not in normalized)

    return {"matched": matched, "missing": missing}
