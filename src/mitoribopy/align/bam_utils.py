"""BAM-level utilities used by ``mitoribopy align`` - MAPQ filtering,
BAM -> BED6 conversion, and primary-mapped read counting. Implemented
with pysam per the Phase 3.6 decision so we do not depend on samtools or
bedtools being on PATH.

Why BED6 and not BED3
---------------------
The rpf pipeline's BED ingest already tolerates BED6 (strand column is
column 6). Producing BED6 preserves the ``--library-strandedness``
decision made at alignment time. On Path A (transcriptome reference)
every transcript is a separate FASTA record so the ``+`` / ``-`` column
records the read orientation relative to the mRNA; that lets downstream
QC flag any reverse-complement contamination if the ``--norc``/``--nofw``
strand guard ever failed.

Why 0-based half-open start/end
-------------------------------
Standard BED semantics (UCSC) and pysam's ``reference_start`` is already
0-based; ``reference_end`` is already half-open. We emit those values
verbatim without adjustment.
"""

from __future__ import annotations

from pathlib import Path

import pysam


# ---------------------------------------------------------------------------
# Primary-mapped read counting
# ---------------------------------------------------------------------------


def count_mapped_reads(bam_path: Path) -> int:
    """Return the number of primary mapped records in *bam_path*.

    Excludes unmapped records (FLAG 0x4), secondary alignments (0x100),
    and supplementary alignments (0x800). This is the definition used
    across :mod:`mitoribopy.align` wherever "how many reads made it to
    this stage?" needs a single honest answer.
    """
    bam_path = Path(bam_path)
    with pysam.AlignmentFile(str(bam_path), "rb") as handle:
        return sum(
            1
            for read in handle
            if not read.is_unmapped
            and not read.is_secondary
            and not read.is_supplementary
        )


# ---------------------------------------------------------------------------
# MAPQ filtering
# ---------------------------------------------------------------------------


def filter_bam_mapq(
    *,
    bam_in: Path,
    bam_out: Path,
    mapq_threshold: int,
    indexer=pysam.index,
) -> int:
    """Write a new BAM containing only records with ``mapq >= threshold``.

    Unmapped, secondary, and supplementary records are dropped
    regardless of their MAPQ, matching :func:`count_mapped_reads` so the
    stage counts stay consistent.

    The primary use of MAPQ filtering on mt-Ribo-seq is NUMT cross-talk
    suppression: pseudogenized copies of mt-mRNAs in the nuclear genome
    accept a small fraction of reads at low MAPQ, and filtering at
    ``--mapq 10`` (the Phase 3 default) drops most of them without
    losing genuine mt reads.

    Returns
    -------
    int
        Number of records written to ``bam_out``.
    """
    bam_in = Path(bam_in)
    bam_out = Path(bam_out)
    bam_out.parent.mkdir(parents=True, exist_ok=True)

    kept = 0
    with pysam.AlignmentFile(str(bam_in), "rb") as source:
        with pysam.AlignmentFile(str(bam_out), "wb", template=source) as sink:
            for read in source:
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if (read.mapping_quality or 0) < int(mapq_threshold):
                    continue
                sink.write(read)
                kept += 1

    indexer(str(bam_out))
    return kept


# ---------------------------------------------------------------------------
# BAM -> BED6 conversion
# ---------------------------------------------------------------------------


def bam_to_bed6(
    *,
    bam_in: Path,
    bed_out: Path,
) -> int:
    """Convert a BAM to a BED6 file, preserving strand.

    Column layout:

    =========  =================================================
    chrom      reference name (for Path A: the mt-transcript id)
    start      0-based start (= ``read.reference_start``)
    end        half-open end (= ``read.reference_end``)
    name       read name (includes the UMI tail after '_' when
               cutadapt extracted one)
    score      MAPQ, clamped into [0, 1000] per BED spec
    strand     ``"-"`` if ``read.is_reverse`` else ``"+"``
    =========  =================================================

    Unmapped / secondary / supplementary records are skipped. Returns
    the number of BED rows written so the caller can reconcile it
    against the stage counts.
    """
    bam_in = Path(bam_in)
    bed_out = Path(bed_out)
    bed_out.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    with pysam.AlignmentFile(str(bam_in), "rb") as source, bed_out.open(
        "w", encoding="utf-8"
    ) as sink:
        for read in source:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.reference_name is None or read.reference_end is None:
                continue  # defensive; should not occur for primary mapped reads

            strand = "-" if read.is_reverse else "+"
            score = max(0, min(1000, int(read.mapping_quality or 0)))
            row = (
                f"{read.reference_name}\t"
                f"{int(read.reference_start)}\t"
                f"{int(read.reference_end)}\t"
                f"{read.query_name or '.'}\t"
                f"{score}\t"
                f"{strand}\n"
            )
            sink.write(row)
            written += 1

    return written
