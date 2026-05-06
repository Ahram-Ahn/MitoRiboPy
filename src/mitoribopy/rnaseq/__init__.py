"""RNA-seq integration for MitoRiboPy.

Consumes a pre-computed differential-expression table (DESeq2 / Xtail /
Anota2Seq) plus a prior ``mitoribopy rpf`` run and produces TE and
delta-TE tables + plots. A SHA256 reference-consistency gate ensures
Ribo-seq and RNA-seq were aligned to the SAME transcript set - a hard
requirement for TE comparisons.

The from-FASTQ flow (``--rna-fastq``) additionally runs cutadapt +
bowtie2 + pyDESeq2 inside the subcommand and is exploratory only;
publication-grade DE belongs in an external full-transcriptome run
fed back via ``--de-table``.
"""

from ._types import (
    DE_COLUMN_ALIASES,
    GENE_ID_CONVENTIONS,
    DeColumnMap,
    DeFormat,
    DeTable,
    DTeRow,
    GeneIdConvention,
    TeRow,
)
from .de_loader import detect_de_format, load_de_table
from .gene_ids import (
    HUMAN_MT_MRNAS,
    YEAST_MT_MRNAS,
    MtGene,
    map_de_gene_ids_to_ribo_ids,
    match_de_to_ribo_genes,
    match_mt_mrnas,
    normalize_gene_id,
)
from .reference_gate import (
    ReferenceMismatchError,
    compute_reference_checksum,
    verify_reference_consistency,
)
from .te import compute_delta_te, compute_te
from .counts import load_ribo_counts
from .fastq_pairing import FastqSample, detect_samples, enumerate_fastqs
from .umi_detect import UmiDetectionResult, detect_umi

__all__ = [
    "DE_COLUMN_ALIASES",
    "DTeRow",
    "DeColumnMap",
    "DeFormat",
    "DeTable",
    "FastqSample",
    "GENE_ID_CONVENTIONS",
    "GeneIdConvention",
    "HUMAN_MT_MRNAS",
    "MtGene",
    "ReferenceMismatchError",
    "TeRow",
    "UmiDetectionResult",
    "YEAST_MT_MRNAS",
    "compute_delta_te",
    "compute_reference_checksum",
    "compute_te",
    "detect_de_format",
    "detect_samples",
    "detect_umi",
    "enumerate_fastqs",
    "load_de_table",
    "load_ribo_counts",
    "map_de_gene_ids_to_ribo_ids",
    "match_de_to_ribo_genes",
    "match_mt_mrnas",
    "normalize_gene_id",
    "verify_reference_consistency",
]
