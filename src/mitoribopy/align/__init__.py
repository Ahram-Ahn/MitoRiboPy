"""FASTQ -> BAM + BED preprocessing pipeline for mt-Ribo-seq.

The ``mitoribopy align`` subcommand wraps external tools (cutadapt, bowtie2,
samtools, umi_tools) behind strict availability checks and produces outputs
that drop into the ``mitoribopy rpf`` subcommand unchanged (Path A -
transcriptome reference).

Modules:

* :mod:`mitoribopy.align._types`         shared dataclasses and literals
* :mod:`mitoribopy.align.tool_check`     PATH + version verification
* :mod:`mitoribopy.align.adapter_detect` FASTQ-head adapter sanity check
* :mod:`mitoribopy.align.trim`           cutadapt wrapper with kit presets
* :mod:`mitoribopy.align.contam`         bowtie2 contaminant subtraction
* :mod:`mitoribopy.align.align`          bowtie2 mt-transcriptome alignment
* :mod:`mitoribopy.align.dedup`          umi_tools / skip / mark-duplicates
* :mod:`mitoribopy.align.bam_utils`      pysam-based MAPQ filter + BAM -> BED6
* :mod:`mitoribopy.align.read_counts`    per-sample provenance table
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
from .adapter_detect import (
    DetectionResult,
    detect_adapter,
    format_per_kit_rates,
)
from .align import align_mt
from .bam_utils import (
    bam_to_bed6,
    count_mapped_reads,
    filter_bam_mapq,
)
from .contam import (
    parse_bowtie2_stderr,
    subtract_contaminants,
)
from .dedup import (
    CONFIRM_MARK_DUPLICATES_FLAG,
    resolve_dedup_strategy,
    run_dedup,
    run_mark_duplicates,
    run_umi_tools_dedup,
    skip_dedup,
)
from .read_counts import (
    assemble_sample_counts,
    format_row,
    read_counts_columns,
    write_read_counts_table,
)
from .trim import (
    parse_cutadapt_json,
    resolve_kit_settings,
    run_cutadapt,
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
    "CONFIRM_MARK_DUPLICATES_FLAG",
    "DetectionResult",
    "align_mt",
    "assemble_sample_counts",
    "bam_to_bed6",
    "count_mapped_reads",
    "detect_adapter",
    "filter_bam_mapq",
    "format_per_kit_rates",
    "format_row",
    "parse_bowtie2_stderr",
    "parse_cutadapt_json",
    "read_counts_columns",
    "resolve_dedup_strategy",
    "resolve_kit_settings",
    "run_cutadapt",
    "run_dedup",
    "run_mark_duplicates",
    "run_umi_tools_dedup",
    "skip_dedup",
    "subtract_contaminants",
    "write_read_counts_table",
]
