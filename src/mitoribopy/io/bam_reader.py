"""BAM input preprocessor for ``mitoribopy rpf`` (Phase 4).

Converts BAM alignment files to BED6 so they can flow through the
existing BED-centric rpf pipeline unchanged. Implemented with pysam
(htslib) via :mod:`mitoribopy.align.bam_utils` - no samtools or bedtools
PATH dependency.

BAM -> BED6 field mapping
-------------------------

Every primary mapped record (excluding FLAG 0x4 unmapped, 0x100
secondary, 0x800 supplementary) produces one BED6 line:

========  ===========================================================
BED col   Source
========  ===========================================================
chrom     ``read.reference_name`` (on Path A, the mt-transcript id
          that matches the annotation CSV's ``sequence_name`` column)
start     ``read.reference_start`` - already 0-based in htslib, no
          adjustment
end       ``read.reference_end``   - already half-open in htslib
name      ``read.query_name`` (carries the UMI tail after ``_`` when
          cutadapt extracted one in Phase 3 Step B)
score     ``max(0, min(1000, read.mapping_quality))`` - clamped into
          the BED spec's [0, 1000] range; MAPQ is [0, 255] in practice
strand    ``"-"`` if ``read.is_reverse`` else ``"+"``
========  ===========================================================

Strand preservation matters even though Phase 3's Path A
(transcriptome reference) already enforces polarity via bowtie2
``--norc``/``--nofw``: any reverse-strand rows in the resulting BED6
flag a library-prep mismatch that QC can surface downstream. Stripping
the strand column (emitting BED3) silently erases that signal.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

from ..align.bam_utils import bam_to_bed6, filter_bam_mapq
from ..console import log_info, log_warning


def convert_bam_to_bed(
    bam_in: Path,
    bed_out: Path,
    *,
    mapq_threshold: int = 0,
) -> int:
    """Convert a BAM file to a BED6 file, optionally MAPQ-filtering first.

    Parameters
    ----------
    bam_in:
        Input BAM (coordinate-sorted or unsorted; pysam reads both).
    bed_out:
        Destination BED6 path. Parent directory is created if needed.
    mapq_threshold:
        MAPQ cutoff applied before conversion. ``0`` disables filtering.
        Default matches ``mitoribopy align --mapq`` (10 in the orchestrator),
        which is the NUMT-suppression level recommended in Phase 3.1.

    Returns
    -------
    int
        Number of BED6 rows written.
    """
    bam_in = Path(bam_in)
    bed_out = Path(bed_out)
    bed_out.parent.mkdir(parents=True, exist_ok=True)

    if mapq_threshold and mapq_threshold > 0:
        with tempfile.TemporaryDirectory(
            prefix="mitoribopy_bam_mapq_", dir=bed_out.parent
        ) as tmp:
            filtered_bam = Path(tmp) / "mapq_filtered.bam"
            filter_bam_mapq(
                bam_in=bam_in,
                bam_out=filtered_bam,
                mapq_threshold=mapq_threshold,
            )
            return bam_to_bed6(bam_in=filtered_bam, bed_out=bed_out)

    return bam_to_bed6(bam_in=bam_in, bed_out=bed_out)


def prepare_bam_inputs(
    input_dir: Path,
    converted_dir: Path,
    *,
    mapq_threshold: int = 0,
) -> list[Path]:
    """Convert every ``*.bam`` file in *input_dir* to BED6 under *converted_dir*.

    Returns the list of produced BED paths in sorted order. Each BAM
    ``<sample>.bam`` becomes ``<converted_dir>/<sample>.bed``; the sample
    name is the BAM filename stem.

    If both ``<sample>.bam`` and ``<sample>.bed`` exist in
    ``input_dir`` (name conflict), the BAM is skipped with a warning so
    the user's explicit BED wins. This matches the Phase 4 decision
    documented in the CHANGELOG.

    An empty list is returned when ``input_dir`` contains no BAM files;
    the caller can then fall back to the unchanged BED-only path.
    """
    input_dir = Path(input_dir)
    converted_dir = Path(converted_dir)

    if not input_dir.is_dir():
        return []

    bam_files = sorted(p for p in input_dir.iterdir() if p.suffix == ".bam")
    if not bam_files:
        return []

    bed_stems = {p.stem for p in input_dir.iterdir() if p.suffix == ".bed"}

    converted_dir.mkdir(parents=True, exist_ok=True)
    produced: list[Path] = []
    for bam_path in bam_files:
        sample = bam_path.stem
        if sample in bed_stems:
            log_warning(
                "BAM",
                f"Both {sample}.bed and {sample}.bam found in {input_dir}; "
                "using the explicit BED and skipping the BAM.",
            )
            continue

        bed_path = converted_dir / f"{sample}.bed"
        log_info(
            "BAM",
            f"Converting {bam_path.name} -> {bed_path} "
            f"(MAPQ filter: "
            f"{'off' if mapq_threshold <= 0 else f'>= {mapq_threshold}'}).",
        )
        convert_bam_to_bed(
            bam_in=bam_path,
            bed_out=bed_path,
            mapq_threshold=mapq_threshold,
        )
        produced.append(bed_path)

    return produced
