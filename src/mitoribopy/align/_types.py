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

DedupStrategy = Literal["auto", "umi-tools", "skip"]
"""Deduplication strategy.

``auto``            resolves to ``umi-tools`` when ``--umi-length > 0`` else
                    ``skip``. Safe default for mt-Ribo-seq.
``umi-tools``       UMI-aware collapse (only correct choice when UMIs exist).
``skip``            no deduplication; required default for UMI-less,
                    low-complexity mt-Ribo-seq libraries.

The legacy ``mark-duplicates`` (picard) option is not supported because
coordinate-only deduplication flattens codon-occupancy signal on
mt-Ribo-seq libraries. See docs/validation/taco1_ko_regression.md for
the empirical evidence.
"""

UmiPosition = Literal["5p", "3p", "both"]
"""Where the UMI sits within the insert.

``5p``    UMI at the 5' end of the read; cutadapt extracts in pass 1.
``3p``    UMI at the 3' end of the read; cutadapt extracts in pass 2
          (after the adapter trim, since cutadapt's ``-u -N`` runs
          BEFORE adapter search in a single call).
``both``  Dual-end UMI (e.g. xGen Duplex, Twist) — pass 1 extracts
          the 5' UMI, pass 2 extracts the 3' UMI and APPENDS to the
          QNAME so umi_tools sees the concatenated token after the
          last ``_`` separator. The two UMIs can have different
          lengths (``umi_length_5p`` + ``umi_length_3p``).
"""


# ---------------------------------------------------------------------------
# Kit presets
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class KitPreset:
    """Library-prep kit signature used for adapter detection and reporting.

    Kit presets are an INTERNAL concept since v0.7.1: the user-facing
    ``--kit-preset`` flag was removed. Users supply ``--adapter`` or let
    auto-detection identify the kit; the names here are produced by the
    detector and surface only in ``kit_resolution.tsv`` for provenance.
    """

    name: str
    adapter: str | None
    umi_length: int
    umi_position: UmiPosition
    description: str
    # Per-end UMI lengths for ``umi_position == "both"`` libraries (e.g.
    # xGen Duplex, Twist). For single-end UMI presets these stay at 0
    # and ``umi_length`` continues to carry the canonical length.
    umi_length_5p: int = 0
    umi_length_3p: int = 0


# Canonical preset registry, organized by adapter family rather than
# kit brand. Most commercial library prep kits share one of three 3'
# adapter sequences; differentiating presets by adapter family (instead
# of by vendor name) keeps the registry small and makes detection
# unambiguous (no two presets share an adapter).
#
# Two non-detecting sentinels live alongside the real adapter families:
# ``pretrimmed`` (no adapter to trim — auto-inferred or declared via
# ``--pretrimmed``) and ``custom`` (the kit name reported when the user
# supplied an unknown ``--adapter`` sequence directly).
#
# Every adapter sequence here MUST be cross-checked against the kit
# vendor's documentation; a mismatch silently corrupts real data.
KIT_PRESETS: dict[str, KitPreset] = {
    # --- Sentinels / fallbacks (internal labels; not user-selectable) ---
    "pretrimmed": KitPreset(
        name="pretrimmed",
        adapter=None,
        umi_length=0,
        umi_position="5p",
        description=(
            "Already-trimmed FASTQs (e.g. SRA-deposited data, or output of "
            "a prior trim step) — the adapter has already been clipped, so "
            "cutadapt skips the -a flag and only enforces length and "
            "quality filtering. 'pretrimmed' describes the INPUT STATE, "
            "not an action this pipeline takes. Auto-inferred when adapter "
            "detection finds no known kit signature; users can also force "
            "it via ``--pretrimmed``."
        ),
    ),
    "custom": KitPreset(
        name="custom",
        adapter=None,
        umi_length=0,
        umi_position="5p",
        description=(
            "Internal label for runs where the user supplied ``--adapter "
            "<SEQ>`` and the sequence did not match any known adapter "
            "family. The reported adapter is the user's literal sequence."
        ),
    ),
    # --- Real adapter families (these are what the detector returns) ---
    "illumina_smallrna": KitPreset(
        name="illumina_smallrna",
        adapter="TGGAATTCTCGGGTGCCAAGG",
        umi_length=0,
        umi_position="5p",
        description=(
            "Illumina TruSeq Small RNA Library Prep adapter "
            "TGGAATTCTCGGGTGCCAAGG, no UMI. Covers Illumina TruSeq Small "
            "RNA and any compatible kit using the same 3' adapter."
        ),
    ),
    "illumina_truseq": KitPreset(
        name="illumina_truseq",
        adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        umi_length=0,
        umi_position="5p",
        description=(
            "Illumina TruSeq Read 1 adapter "
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA, no UMI. Covers most "
            "current Illumina-compatible total/stranded RNA preps: NEB "
            "Next Multiplex Small RNA, TruSeq Stranded Total RNA "
            "(Gold), Takara SMARTer Stranded Total RNA-Seq v3 Pico "
            "Input, Bio-Rad SEQuoia Express Standard, and any other "
            "prep using the standard Illumina R1 adapter without a UMI."
        ),
    ),
    "illumina_truseq_umi": KitPreset(
        name="illumina_truseq_umi",
        adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        umi_length=8,
        umi_position="5p",
        description=(
            "Illumina TruSeq Read 1 adapter with 8 nt 5' UMI. Covers "
            "NEBNext Ultra II with UMI Adapters, Bio-Rad SEQuoia "
            "Complete (UMI), and other UMI-bearing variants of the "
            "Illumina R1 adapter family. Confirm the UMI length against "
            "your kit's insert diagram before using."
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
}


@dataclass(frozen=True)
class ResolvedKit:
    """The effective (kit, adapter, UMI) tuple after CLI overrides.

    Produced by :func:`mitoribopy.align.trim.resolve_kit_settings` and
    written verbatim into the run manifest so a reader can reconstruct
    the exact trimming that happened.

    ``adapter`` is ``None`` only for the ``pretrimmed`` kit; in every
    other case it carries the resolved 3' adapter sequence. cutadapt
    skips the ``-a`` flag entirely when adapter is ``None``, falling
    back to length + quality filtering only.
    """

    kit: str
    adapter: str | None
    umi_length: int
    umi_position: UmiPosition
    # Per-end UMI lengths. For ``umi_position == "both"`` libraries
    # ``umi_length`` MUST equal ``umi_length_5p + umi_length_3p`` (the
    # final QNAME UMI token is the concatenated 5' + 3' string).
    umi_length_5p: int = 0
    umi_length_3p: int = 0


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
class SampleOverride:
    """Per-sample CLI/YAML override for adapter + UMI + dedup resolution.

    Populated from the ``align.samples:`` YAML list (via the ``mitoribopy
    all`` orchestrator) or from a TSV passed as ``--sample-overrides``
    to ``mitoribopy align``. Any field set to ``None`` falls through to
    the global CLI default for that field, so a user can override only
    ``umi_length`` for one sample and let everything else inherit.

    The ``sample`` field matches the FASTQ basename with the trailing
    ``.fq[.gz]`` / ``.fastq[.gz]`` suffix removed (the same naming
    convention :func:`mitoribopy.cli.align._sample_name` produces).

    ``pretrimmed=True`` declares already-trimmed input for that sample
    and is mutually exclusive with ``adapter`` (a sample is either
    pretrimmed or has an adapter to trim).
    """

    sample: str
    adapter: str | None = None
    pretrimmed: bool | None = None
    umi_length: int | None = None
    umi_position: UmiPosition | None = None
    dedup_strategy: DedupStrategy | None = None
    # Per-end UMI lengths for ``umi_position == "both"`` libraries.
    # Inherit the global default when ``None``.
    umi_length_5p: int | None = None
    umi_length_3p: int | None = None


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
    # ``post_trim`` = read count AFTER the cutadapt step has run. For a
    # ``pretrimmed`` kit (already-clipped FASTQ) cutadapt only enforces
    # length + quality, so this number can equal ``total_reads`` minus
    # the length-filter losses. The name does NOT mean "reads were
    # adapter-trimmed".
    post_trim: int
    rrna_aligned: int  # reads removed by the contam subtract step
    post_rrna_filter: int  # reads entering mt-transcriptome alignment
    mt_aligned: int  # reads aligned to the mt-transcriptome
    unaligned_to_mt: int  # post_rrna_filter minus mt_aligned
    mt_aligned_after_mapq: int
    mt_aligned_after_dedup: int
