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
    aliases: tuple[str, ...] = ()


# Human mt-mRNAs, checked against Ensembl GRCh38.p14 / NCBI 2024.
# Order matches the canonical genomic layout so reviewers can spot
# omissions.
HUMAN_MT_MRNAS: tuple[MtGene, ...] = (
    MtGene(hgnc="MT-ND1", bare="ND1", ensembl="ENSG00000198888", refseq="YP_003024026.1"),
    MtGene(hgnc="MT-ND2", bare="ND2", ensembl="ENSG00000198763", refseq="YP_003024027.1"),
    MtGene(
        hgnc="MT-CO1",
        bare="CO1",
        ensembl="ENSG00000198804",
        refseq="YP_003024028.1",
        aliases=("COX1", "MT-COX1"),
    ),
    MtGene(
        hgnc="MT-CO2",
        bare="CO2",
        ensembl="ENSG00000198712",
        refseq="YP_003024029.1",
        aliases=("COX2", "MT-COX2"),
    ),
    MtGene(hgnc="MT-ATP8", bare="ATP8", ensembl="ENSG00000228253", refseq="YP_003024030.1"),
    MtGene(hgnc="MT-ATP6", bare="ATP6", ensembl="ENSG00000198899", refseq="YP_003024031.1"),
    MtGene(
        hgnc="MT-CO3",
        bare="CO3",
        ensembl="ENSG00000198938",
        refseq="YP_003024032.1",
        aliases=("COX3", "MT-COX3"),
    ),
    MtGene(hgnc="MT-ND3", bare="ND3", ensembl="ENSG00000198840", refseq="YP_003024033.1"),
    MtGene(hgnc="MT-ND4L", bare="ND4L", ensembl="ENSG00000212907", refseq="YP_003024034.1"),
    MtGene(hgnc="MT-ND4", bare="ND4", ensembl="ENSG00000198886", refseq="YP_003024035.1"),
    MtGene(hgnc="MT-ND5", bare="ND5", ensembl="ENSG00000198786", refseq="YP_003024036.1"),
    MtGene(hgnc="MT-ND6", bare="ND6", ensembl="ENSG00000198695", refseq="YP_003024037.1"),
    MtGene(
        hgnc="MT-CYB",
        bare="CYB",
        ensembl="ENSG00000198727",
        refseq="YP_003024038.1",
        aliases=("CYTB", "MT-CYTB"),
    ),
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
    if organism in {"y", "yeast", "s_cerevisiae", "s.cerevisiae"}:
        return YEAST_MT_MRNAS
    if organism in {"h", "human", "homo_sapiens", "h.sapiens"}:
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


def _all_gene_ids(gene: MtGene) -> tuple[str, ...]:
    """Return every recognized spelling for one registry gene."""
    ids = [gene.hgnc, gene.bare, gene.ensembl, gene.refseq, *gene.aliases]
    return tuple(str(value) for value in ids if value)


def _alias_lookup(organism: str) -> dict[str, str]:
    """Map every known alias to the registry canonical key.

    The canonical key is ``gene.hgnc`` for human and the SGD-style gene
    symbol for yeast. This lets external-DE identifiers and RPF count
    gene names meet in the middle without requiring both files to use
    identical spelling.
    """
    lookup: dict[str, str] = {}
    for gene in _registry_for_organism(organism):
        for alias in _all_gene_ids(gene):
            lookup[normalize_gene_id(alias)] = gene.hgnc
    return lookup


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
    alias_lookup = _alias_lookup(organism)
    present_canonical = {
        alias_lookup[norm]
        for g in de_gene_ids
        if (norm := normalize_gene_id(g)) in alias_lookup
    }
    expected = expected_ids(convention, organism)
    expected_canonical = [
        gene.hgnc for gene in _registry_for_organism(organism)
        if _id_for_convention(gene, convention)
    ]

    matched = sorted(
        expected_id
        for expected_id, canonical in zip(expected, expected_canonical, strict=False)
        if canonical in present_canonical
    )
    missing = sorted(
        expected_id
        for expected_id, canonical in zip(expected, expected_canonical, strict=False)
        if canonical not in present_canonical
    )

    return {"matched": matched, "missing": missing}


def map_de_gene_ids_to_ribo_ids(
    de_gene_ids: list[str],
    ribo_gene_ids: list[str],
    convention: str,
    organism: str = "h",
) -> dict[str, str]:
    """Return a DE-table gene-id -> RPF-counts gene-id mapping.

    The DE table is allowed to use the selected convention while
    ``rpf_counts.tsv`` may carry transcript/reference names from the
    Ribo-seq FASTA. Exact matches are preserved, and known mt-mRNA
    aliases are bridged. For example, a human DE row ``MT-CO1`` maps to
    an RPF counts gene named ``COX1``.
    """
    alias_lookup = _alias_lookup(organism)
    ribo_by_norm = {normalize_gene_id(g): g for g in ribo_gene_ids}
    ribo_by_canonical: dict[str, str] = {}
    for gene_id in ribo_gene_ids:
        canonical = alias_lookup.get(normalize_gene_id(gene_id))
        if canonical is not None:
            ribo_by_canonical.setdefault(canonical, gene_id)

    mapping: dict[str, str] = {}
    for gene_id in de_gene_ids:
        norm = normalize_gene_id(gene_id)
        if norm in ribo_by_norm:
            mapping[gene_id] = ribo_by_norm[norm]
            continue
        canonical = alias_lookup.get(norm)
        if canonical is not None and canonical in ribo_by_canonical:
            mapping[gene_id] = ribo_by_canonical[canonical]
    return mapping


def match_de_to_ribo_genes(
    de_gene_ids: list[str],
    ribo_gene_ids: list[str],
    convention: str,
    organism: str = "h",
) -> dict[str, object]:
    """Return DE/RPF mt-mRNA match diagnostics for run provenance."""
    alias_lookup = _alias_lookup(organism)
    de_canonical = {
        alias_lookup[norm]
        for gene_id in de_gene_ids
        if (norm := normalize_gene_id(gene_id)) in alias_lookup
    }
    ribo_canonical = {
        alias_lookup[norm]
        for gene_id in ribo_gene_ids
        if (norm := normalize_gene_id(gene_id)) in alias_lookup
    }
    expected_by_canonical = {
        gene.hgnc: _id_for_convention(gene, convention)
        for gene in _registry_for_organism(organism)
        if _id_for_convention(gene, convention)
    }
    de_to_ribo = map_de_gene_ids_to_ribo_ids(
        de_gene_ids, ribo_gene_ids, convention, organism=organism
    )
    return {
        "matched": sorted(de_to_ribo.values()),
        "de_to_ribo": de_to_ribo,
        "missing_de": sorted(
            value
            for canonical, value in expected_by_canonical.items()
            if canonical not in de_canonical and value is not None
        ),
        "missing_ribo": sorted(
            value
            for canonical, value in expected_by_canonical.items()
            if canonical not in ribo_canonical and value is not None
        ),
        "unmatched_ribo_genes": sorted(
            gene_id
            for gene_id in ribo_gene_ids
            if normalize_gene_id(gene_id) not in {
                normalize_gene_id(v) for v in de_to_ribo.values()
            }
        ),
    }
