"""bowtie2 alignment to the mt-transcriptome reference.

Aligns the non-contaminant FASTQ to a bowtie2 index built over the
mt-mRNA transcriptome. The output is a coordinate-sorted, indexed BAM
that drops into downstream MAPQ filtering and BAM -> BED conversion
without further shell glue.

Why a transcriptome reference and not the whole mtDNA: aligning to
per-transcript FASTA records sidesteps the ND5 / ND6 antisense
overlap at the genome level (each mt-mRNA is a separate FASTA record);
strand polarity is enforced at alignment time via ``--norc`` / ``--nofw``.

sort/index use pysam.sort and pysam.index rather than shelling out to
samtools. pysam wheels bundle htslib, so mt-transcriptome sort + index
do not add a PATH dependency.

Historical design notes (Path A vs Path B, ToT decision tables) live in
``docs/developer/architecture_history.md``.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
from typing import Iterable

import pysam

from ._types import AlignResult, Strandedness
from .contam import parse_bowtie2_stderr
from .tool_check import bowtie2_strand_flag


_INDEX_SUFFIXES = (".1.bt2", ".1.bt2l")


def _validate_index_prefix(mt_index: Path) -> None:
    mt_index = Path(mt_index)
    for suffix in _INDEX_SUFFIXES:
        if Path(str(mt_index) + suffix).exists():
            return
    raise FileNotFoundError(
        f"bowtie2 mt-transcriptome index not found at prefix {mt_index!s}. "
        "Build one first, for example:\n"
        f"  bowtie2-build mt_transcripts.fa {mt_index}\n"
        "The FASTA should contain each mt-mRNA as a separate record whose "
        "header matches the 'sequence_name' column of the annotation CSV."
    )


def _build_bowtie2_align_command(
    *,
    mt_index: Path,
    fastq_in: Path,
    sam_out: Path,
    strandedness: Strandedness,
    threads: int,
    seed: int,
    extra_flags: Iterable[str] = (),
) -> list[str]:
    cmd: list[str] = [
        "bowtie2",
        "-x",
        str(mt_index),
        "-U",
        str(fastq_in),
        "--end-to-end",
        "--very-sensitive",
        "-L",
        "18",
        "--no-unal",
        "-p",
        str(max(int(threads), 1)),
        "--seed",
        str(seed),
    ]
    cmd += bowtie2_strand_flag(strandedness)
    cmd += list(extra_flags)
    cmd += ["-S", str(sam_out)]
    return cmd


def align_mt(
    *,
    fastq_in: Path,
    mt_index: Path,
    bam_out: Path,
    strandedness: Strandedness = "forward",
    threads: int = 1,
    seed: int = 42,
    extra_flags: Iterable[str] = (),
    runner=subprocess.run,
    sorter=pysam.sort,
    indexer=pysam.index,
) -> AlignResult:
    """Align reads to the mt-transcriptome and return a sorted, indexed BAM.

    Parameters
    ----------
    fastq_in:
        Trimmed, contaminant-subtracted FASTQ from Step D.
    mt_index:
        bowtie2 index prefix of the mt-transcriptome (per-transcript
        FASTA, headers matching annotation ``sequence_name``).
    bam_out:
        Destination for the sorted, indexed BAM.
    strandedness:
        Library strandedness. ``forward`` -> ``--norc``, ``reverse`` ->
        ``--nofw``, ``unstranded`` -> neither. The choice is recorded
        here because ND5 / ND6, ATP8 / ATP6, and ND4L / ND4 overlap in
        the parent genome, and even on a transcriptome reference a
        wrongly-stranded library would reverse-complement reads onto
        the wrong transcript record.
    threads, seed, extra_flags:
        Plumbing pass-throughs.
    runner, sorter, indexer:
        Injection points for tests; defaults run real subprocess /
        pysam.sort / pysam.index.
    """
    fastq_in = Path(fastq_in)
    mt_index = Path(mt_index)
    bam_out = Path(bam_out)

    _validate_index_prefix(mt_index)
    bam_out.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(
        prefix="mitoribopy_align_", dir=bam_out.parent
    ) as tmp:
        tmp_sam = Path(tmp) / "aligned.sam"

        cmd = _build_bowtie2_align_command(
            mt_index=mt_index,
            fastq_in=fastq_in,
            sam_out=tmp_sam,
            strandedness=strandedness,
            threads=threads,
            seed=seed,
            extra_flags=extra_flags,
        )
        completed = runner(cmd, check=False, capture_output=True, text=True)
        if completed.returncode != 0:
            stderr = (getattr(completed, "stderr", "") or "").strip()
            raise RuntimeError(
                "bowtie2 (mt-transcriptome alignment) failed with exit "
                f"code {completed.returncode}: "
                f"{stderr or '<no stderr captured>'}"
            )

        total, unaligned = parse_bowtie2_stderr(
            getattr(completed, "stderr", "") or ""
        )

        sorter(
            "-O",
            "bam",
            "-o",
            str(bam_out),
            str(tmp_sam),
            "-@",
            str(max(int(threads), 1)),
        )

    indexer(str(bam_out))

    return AlignResult(
        total_reads=total,
        aligned=total - unaligned,
        bam_path=bam_out,
    )
