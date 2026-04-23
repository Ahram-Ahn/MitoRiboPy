"""Shared types and constants for the ``mitoribopy align`` pipeline.

Every inter-module interface in :mod:`mitoribopy.align` uses a frozen
dataclass from this module. This keeps step boundaries typed and keeps
all configuration constants (kit presets, strandedness literals,
deduplication strategies) in one place so a reviewer can audit them
without hunting across modules.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal


# ---------------------------------------------------------------------------
# Literals
# ---------------------------------------------------------------------------

Strandedness = Literal["forward", "reverse", "unstranded"]
"""Library strandedness.

``forward``  first read is same strand as the mRNA (e.g. NEBNext Small RNA,
             TruSeq Small RNA). Passed to bowtie2 as ``--norc``.
``reverse``  first read is opposite strand to the mRNA (dUTP-stranded).
             Passed to bowtie2 as ``--nofw``.
``unstranded`` no strand constraint; bowtie2 is left permissive.
               Triggers a WARNING on genome references; safe on
               Path A (transcriptome) because each transcript record
               is pre-oriented.
"""

DedupStrategy = Literal["auto", "umi-tools", "skip", "mark-duplicates"]
"""Deduplication strategy.

``auto``            resolves to ``umi-tools`` when ``--umi-length > 0`` else
                    ``skip``. Safe default for mt-Ribo-seq.
``umi-tools``       UMI-aware collapse (only correct choice when UMIs exist).
``skip``            no deduplication; required default for UMI-less,
                    low-complexity mt-Ribo-seq libraries.
``mark-duplicates`` picard MarkDuplicates. Destroys codon-occupancy signal
                    on mt-Ribo-seq; only enabled when the user also passes
                    the long confirmation flag.
"""

UmiPosition = Literal["5p", "3p"]


# ---------------------------------------------------------------------------
# Kit presets
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class KitPreset:
    """Library-prep kit defaults for adapter and UMI handling.

    The ``adapter`` may be ``None`` for ``custom``, in which case the user
    MUST pass ``--adapter`` explicitly. This is intentional: a silently
    wrong adapter is one of the worst failure modes for mt-Ribo-seq
    because it shifts the RPF length distribution outside the 15-45 nt
    filter and drops reads the user does not realize they lost.
    """

    name: str
    adapter: str | None
    umi_length: int
    umi_position: UmiPosition
    description: str


# The preset registry. Every addition here MUST be cross-checked against
# the kit vendor's documentation; a mismatch silently corrupts real data.
KIT_PRESETS: dict[str, KitPreset] = {
    "truseq_smallrna": KitPreset(
        name="truseq_smallrna",
        adapter="TGGAATTCTCGGGTGCCAAGG",
        umi_length=0,
        umi_position="5p",
        description="Illumina TruSeq Small RNA (3' adapter TGGAATTCTCGGGTGCCAAGG).",
    ),
    "nebnext_smallrna": KitPreset(
        name="nebnext_smallrna",
        adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        umi_length=0,
        umi_position="5p",
        description=(
            "NEB Next Multiplex Small RNA Library Prep "
            "(3' SRA adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA)."
        ),
    ),
    "nebnext_ultra_umi": KitPreset(
        name="nebnext_ultra_umi",
        adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        umi_length=8,
        umi_position="5p",
        description=(
            "NEBNext Ultra II with UMI Adapters (3' SRA adapter; 8 nt 5' UMI). "
            "Confirm against your kit's insert diagram before using."
        ),
    ),
    "qiaseq_mirna": KitPreset(
        name="qiaseq_mirna",
        adapter="AACTGTAGGCACCATCAAT",
        umi_length=12,
        umi_position="3p",
        description=(
            "QIAseq miRNA Library Kit (3' adapter AACTGTAGGCACCATCAAT; "
            "12 nt 3' UMI inside the insert)."
        ),
    ),
    "custom": KitPreset(
        name="custom",
        adapter=None,
        umi_length=0,
        umi_position="5p",
        description=(
            "No defaults applied. --adapter must be supplied explicitly. "
            "Use this preset whenever your kit is not one of the named ones."
        ),
    ),
}


@dataclass(frozen=True)
class ResolvedKit:
    """The effective (kit, adapter, UMI) tuple after CLI overrides.

    Produced by :func:`mitoribopy.align.trim.resolve_kit_settings` and
    written verbatim into the Phase 6 run manifest so a reader can
    reconstruct the exact trimming that happened.
    """

    kit: str
    adapter: str
    umi_length: int
    umi_position: UmiPosition


# ---------------------------------------------------------------------------
# Tool metadata + step results
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ToolInfo:
    """A single external tool's resolved location and version."""

    name: str
    path: str
    version: str


@dataclass(frozen=True)
class CutadaptResult:
    """Counts extracted from cutadapt's ``--json`` log."""

    input_reads: int
    reads_with_adapter: int
    reads_passing_filters: int
    log_json_path: Path


@dataclass(frozen=True)
class ContamResult:
    """Counts from the contaminant-subtraction bowtie2 step."""

    total_reads: int
    aligned_to_contam: int
    unaligned_reads: int  # post_rRNA_filter


@dataclass(frozen=True)
class AlignResult:
    """Result of the mt-transcriptome bowtie2 alignment step."""

    total_reads: int
    aligned: int
    bam_path: Path


@dataclass(frozen=True)
class DedupResult:
    """Result of the deduplication step."""

    strategy: str
    input_reads: int
    output_reads: int
    bam_path: Path


@dataclass(frozen=True)
class SampleCounts:
    """Per-sample stage counts written into ``read_counts.tsv``.

    This is the provenance spine: every downstream analysis that wants
    to know "how many reads made it to mt-aligned?" reads this table.
    Column order below is preserved verbatim in the TSV header so
    downstream consumers can index by position or by name.
    """

    sample: str
    total_reads: int
    post_trim: int
    rrna_aligned: int  # reads removed by the contam subtract step
    post_rrna_filter: int  # reads entering mt-transcriptome alignment
    mt_aligned: int  # reads aligned to the mt-transcriptome
    unaligned_to_mt: int  # post_rrna_filter minus mt_aligned
    mt_aligned_after_mapq: int
    mt_aligned_after_dedup: int
