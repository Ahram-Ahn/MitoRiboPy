"""FASTQ -> BAM + BED preprocessing pipeline for mt-Ribo-seq.

The ``mitoribopy align`` subcommand wraps external tools (cutadapt, bowtie2,
samtools, umi_tools) behind strict availability checks and produces outputs
that drop into the ``mitoribopy rpf`` subcommand unchanged (Path A -
transcriptome reference).

Modules:

* :mod:`mitoribopy.align._types`       shared dataclasses and literals
* :mod:`mitoribopy.align.tool_check`   PATH + version verification
* :mod:`mitoribopy.align.trim`         cutadapt wrapper with kit presets
* :mod:`mitoribopy.align.contam`       bowtie2 contaminant subtraction
* :mod:`mitoribopy.align.align`        bowtie2 mt-transcriptome alignment
* :mod:`mitoribopy.align.dedup`        umi_tools / skip / mark-duplicates
* :mod:`mitoribopy.align.bam_utils`    pysam-based MAPQ filter + BAM -> BED6
* :mod:`mitoribopy.align.read_counts`  per-sample provenance table
"""

from ._types import (
    KIT_PRESETS,
    AlignResult,
    ContamResult,
    CutadaptResult,
    DedupResult,
    DedupStrategy,
    KitPreset,
    ResolvedKit,
    SampleCounts,
    Strandedness,
    ToolInfo,
)
from .tool_check import (
    ToolNotFoundError,
    bowtie2_strand_flag,
    check_tool,
    ensure_tools_available,
    get_tool_version,
)

__all__ = [
    "KIT_PRESETS",
    "AlignResult",
    "ContamResult",
    "CutadaptResult",
    "DedupResult",
    "DedupStrategy",
    "KitPreset",
    "ResolvedKit",
    "SampleCounts",
    "Strandedness",
    "ToolInfo",
    "ToolNotFoundError",
    "bowtie2_strand_flag",
    "check_tool",
    "ensure_tools_available",
    "get_tool_version",
]
