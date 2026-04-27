"""Per-sample bowtie2 alignment for the from-FASTQ rnaseq mode.

Reuses the SE machinery in :mod:`mitoribopy.align` (cutadapt wrapper,
adapter auto-detect, kit-preset registry) and adds the bits the SE
pipeline does not need:

* a transcriptome-FASTA → bowtie2 index builder with content-addressed
  caching so repeated runs over the same FASTA do not rebuild,
* a paired-end cutadapt + bowtie2 wrapper (the SE align module is
  single-end only),
* a per-transcript count helper that reads a sorted, indexed BAM and
  returns ``{ref_name: count}`` of primary-mapped reads (and read1 only
  for paired records, so each fragment contributes 1).

PE + UMI is currently :class:`NotImplementedError`. Reason: cutadapt's
PE UMI handling needs ``--rename`` against both R1 and R2 with matching
linkers and we have not validated the rename interaction; better to
fail loudly than to silently corrupt UMI dedup.
"""

from __future__ import annotations

import hashlib
import re
import subprocess
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import pysam

from ..align._types import KIT_PRESETS, ResolvedKit
from ..align.adapter_detect import DetectionResult, detect_adapter
from ..align.trim import run_cutadapt
from .fastq_pairing import FastqSample
from .umi_detect import UmiDetectionResult, detect_umi


_BOWTIE2_TOTAL_RE = re.compile(r"^\s*(\d+)\s+reads;\s+of\s+these", re.MULTILINE)
_BOWTIE2_UNALIGNED_RE = re.compile(
    r"^\s*(\d+)\s+\([\d.]+%\)\s+aligned\s+(?:0\s+times|concordantly\s+0\s+times)",
    re.MULTILINE,
)


@dataclass(frozen=True)
class SampleAlignmentResult:
    """Per-sample alignment outcome."""

    sample: str
    bam_path: Path
    counts: dict[str, int]
    paired: bool
    total_reads: int
    aligned_reads: int
    resolved_kit: ResolvedKit


# ---------------------------------------------------------------------------
# Index build (cached)
# ---------------------------------------------------------------------------


def _hash_fasta_prefix(fasta: Path) -> str:
    digest = hashlib.sha256()
    with Path(fasta).open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()[:12]


def build_bowtie2_index(
    reference_fasta: Path,
    index_dir: Path,
    *,
    runner=subprocess.run,
) -> Path:
    """Return a bowtie2 index prefix for ``reference_fasta``.

    The index is cached at ``<index_dir>/transcriptome_<sha12>``; if a
    sidecar ``.1.bt2`` or ``.1.bt2l`` already exists for that prefix the
    function returns the prefix without rebuilding.
    """
    fasta = Path(reference_fasta)
    if not fasta.is_file():
        raise FileNotFoundError(f"Reference FASTA not found: {fasta}")

    index_dir = Path(index_dir)
    index_dir.mkdir(parents=True, exist_ok=True)
    prefix = index_dir / f"transcriptome_{_hash_fasta_prefix(fasta)}"

    for suffix in (".1.bt2", ".1.bt2l"):
        if Path(str(prefix) + suffix).exists():
            return prefix

    cmd = ["bowtie2-build", str(fasta), str(prefix)]
    completed = runner(cmd, check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        stderr = (getattr(completed, "stderr", "") or "").strip()
        raise RuntimeError(
            f"bowtie2-build failed with exit code {completed.returncode}: "
            f"{stderr or '<no stderr captured>'}"
        )
    return prefix


# ---------------------------------------------------------------------------
# Per-sample kit + UMI resolution
# ---------------------------------------------------------------------------


def resolve_sample_kit(
    sample: FastqSample,
    *,
    detect_umis: bool = True,
    detector=detect_adapter,
    umi_detector=detect_umi,
) -> tuple[ResolvedKit, DetectionResult, UmiDetectionResult]:
    """Resolve adapter + UMI for one sample by scanning its R1.

    Decisions:

    * ambiguous adapter call OR pretrimmed → ``"pretrimmed"`` kit
      (``adapter=None``, ``umi_length=0``).
    * unambiguous detection → use the named preset's adapter sequence.
    * preset's own UMI length wins when > 0 (kit knows best); otherwise
      :func:`detect_umi` may add a length when ``detect_umis`` is true.
    """
    detection = detector(sample.r1)

    if detection.preset_name is None or detection.ambiguous or detection.pretrimmed:
        preset = KIT_PRESETS["pretrimmed"]
        umi = (
            umi_detector(sample.r1)
            if detect_umis
            else UmiDetectionResult(
                length=0,
                position="5p",
                entropy_5p=tuple(),
                entropy_3p=tuple(),
                n_reads_scanned=0,
            )
        )
        # In pretrimmed mode the preset has no UMI metadata of its own,
        # so the entropy detector has the floor.
        resolved = ResolvedKit(
            kit=preset.name,
            adapter=None,
            umi_length=umi.length,
            umi_position=umi.position,
        )
        return resolved, detection, umi

    preset = KIT_PRESETS[detection.preset_name]
    umi = (
        umi_detector(sample.r1)
        if detect_umis
        else UmiDetectionResult(
            length=0,
            position="5p",
            entropy_5p=tuple(),
            entropy_3p=tuple(),
            n_reads_scanned=0,
        )
    )

    if preset.umi_length > 0:
        umi_length = preset.umi_length
        umi_position = preset.umi_position
    else:
        umi_length = umi.length
        umi_position = umi.position

    resolved = ResolvedKit(
        kit=preset.name,
        adapter=preset.adapter,
        umi_length=umi_length,
        umi_position=umi_position,
    )
    return resolved, detection, umi


# ---------------------------------------------------------------------------
# Counting helpers
# ---------------------------------------------------------------------------


def count_per_transcript(bam: Path) -> dict[str, int]:
    """Return ``{ref_name: count}`` of primary-mapped reads in ``bam``.

    Skips records where any of ``is_unmapped``, ``is_secondary``,
    ``is_supplementary`` is true. For paired records additionally
    requires ``is_read1`` so each fragment contributes 1, not 2.

    Every reference present in the BAM header is included with a
    starting value of zero so downstream matrices have stable columns
    even when a reference saw no reads.
    """
    counts: dict[str, int] = {}
    with pysam.AlignmentFile(str(bam), "rb") as handle:
        for ref in handle.references:
            counts[ref] = 0
        for record in handle.fetch(until_eof=True):
            if record.is_unmapped or record.is_secondary or record.is_supplementary:
                continue
            if record.is_paired and not record.is_read1:
                continue
            ref = record.reference_name
            if ref is None:
                continue
            counts[ref] = counts.get(ref, 0) + 1
    return counts


def write_counts_matrix(
    results: Sequence[SampleAlignmentResult],
    path: Path,
) -> None:
    """Write a wide gene-by-sample TSV with zero-fill for missing genes."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    samples = [r.sample for r in results]
    genes: set[str] = set()
    for r in results:
        genes.update(r.counts.keys())
    sorted_genes = sorted(genes)

    with path.open("w", encoding="utf-8") as handle:
        handle.write("gene\t" + "\t".join(samples) + "\n")
        for gene in sorted_genes:
            row = [gene]
            for r in results:
                row.append(str(r.counts.get(gene, 0)))
            handle.write("\t".join(row) + "\n")


def write_long_counts(
    results: Sequence[SampleAlignmentResult],
    path: Path,
) -> None:
    """Write a long-format TSV (``sample\\tgene\\tcount``).

    Compatible with :func:`mitoribopy.rnaseq.counts.load_ribo_counts` so
    the existing TE / delta-TE path can consume it unchanged.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("sample\tgene\tcount\n")
        for r in results:
            for gene in sorted(r.counts.keys()):
                handle.write(f"{r.sample}\t{gene}\t{r.counts[gene]}\n")


# ---------------------------------------------------------------------------
# Subprocess wrappers
# ---------------------------------------------------------------------------


def _invoke(runner, cmd: list[str], *, label: str) -> "subprocess.CompletedProcess[str]":
    completed = runner(cmd, check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        stderr = (getattr(completed, "stderr", "") or "").strip()
        raise RuntimeError(
            f"{label} failed with exit code {completed.returncode}: "
            f"{stderr or '<no stderr captured>'}"
        )
    return completed


def _parse_bowtie2_pair_counts(stderr: str) -> tuple[int, int]:
    total_match = _BOWTIE2_TOTAL_RE.search(stderr)
    unaligned_match = _BOWTIE2_UNALIGNED_RE.search(stderr)
    if total_match is None or unaligned_match is None:
        raise ValueError(
            "Could not parse bowtie2 stderr summary. Got:\n" + stderr
        )
    total = int(total_match.group(1))
    unaligned = int(unaligned_match.group(1))
    return total, total - unaligned


def _run_bowtie2_se(
    *,
    bt2_index: Path,
    fastq_in: Path,
    bam_out: Path,
    threads: int,
    seed: int,
    runner,
) -> tuple[int, int]:
    sam_out = bam_out.with_suffix(".sam")
    cmd = [
        "bowtie2",
        "-x", str(bt2_index),
        "-U", str(fastq_in),
        "--end-to-end", "--very-sensitive",
        "-L", "18",
        "--no-unal",
        "-p", str(max(int(threads), 1)),
        "--seed", str(seed),
        "-S", str(sam_out),
    ]
    completed = _invoke(runner, cmd, label="bowtie2 (SE rnaseq alignment)")
    total, aligned = _parse_bowtie2_pair_counts(getattr(completed, "stderr", "") or "")
    pysam.sort("-o", str(bam_out), str(sam_out))
    pysam.index(str(bam_out))
    sam_out.unlink(missing_ok=True)
    return total, aligned


def _run_cutadapt_pe(
    *,
    fastq_r1: Path,
    fastq_r2: Path,
    out_r1: Path,
    out_r2: Path,
    adapter: str | None,
    threads: int,
    min_length: int,
    quality: int,
    runner,
) -> None:
    cmd: list[str] = ["cutadapt"]
    if adapter is not None:
        cmd += ["-a", adapter, "-A", adapter]
    cmd += [
        "--minimum-length", str(min_length),
        "-q", str(quality),
    ]
    if threads and threads > 1:
        cmd += ["--cores", str(threads)]
    cmd += [
        "-o", str(out_r1),
        "-p", str(out_r2),
        str(fastq_r1), str(fastq_r2),
    ]
    _invoke(runner, cmd, label="cutadapt (PE rnaseq trim)")


def _run_bowtie2_pe(
    *,
    bt2_index: Path,
    fastq_r1: Path,
    fastq_r2: Path,
    bam_out: Path,
    threads: int,
    seed: int,
    runner,
) -> tuple[int, int]:
    sam_out = bam_out.with_suffix(".sam")
    cmd = [
        "bowtie2",
        "-x", str(bt2_index),
        "-1", str(fastq_r1),
        "-2", str(fastq_r2),
        "--end-to-end", "--very-sensitive",
        "-L", "18",
        "--no-unal", "--no-mixed", "--no-discordant",
        "-p", str(max(int(threads), 1)),
        "--seed", str(seed),
        "-S", str(sam_out),
    ]
    completed = _invoke(runner, cmd, label="bowtie2 (PE rnaseq alignment)")
    total, aligned = _parse_bowtie2_pair_counts(getattr(completed, "stderr", "") or "")
    pysam.sort("-o", str(bam_out), str(sam_out))
    pysam.index(str(bam_out))
    sam_out.unlink(missing_ok=True)
    return total, aligned


# ---------------------------------------------------------------------------
# Per-sample orchestrator
# ---------------------------------------------------------------------------


def align_sample(
    sample: FastqSample,
    *,
    bt2_index: Path,
    workdir: Path,
    threads: int = 4,
    seed: int = 42,
    detect_umis: bool = True,
    runner=subprocess.run,
) -> SampleAlignmentResult:
    """Trim + align one sample, return counts + provenance.

    SE path reuses :func:`mitoribopy.align.trim.run_cutadapt`. PE path
    runs cutadapt + bowtie2 directly with paired flags. PE + UMI raises
    :class:`NotImplementedError`.
    """
    workdir = Path(workdir)
    sample_dir = workdir / sample.sample
    sample_dir.mkdir(parents=True, exist_ok=True)

    resolved, _detection, _umi = resolve_sample_kit(
        sample, detect_umis=detect_umis
    )

    bam_out = sample_dir / f"{sample.sample}.bam"

    if not sample.paired:
        trimmed = sample_dir / f"{sample.sample}.trimmed.fq.gz"
        log_json = sample_dir / f"{sample.sample}.cutadapt.json"
        run_cutadapt(
            fastq_in=sample.r1,
            fastq_out=trimmed,
            resolved=resolved,
            min_length=15,
            max_length=10_000,  # transcript reads are not RPF-length-bounded
            quality=20,
            threads=threads,
            log_json=log_json,
            runner=runner,
        )
        total, aligned = _run_bowtie2_se(
            bt2_index=bt2_index,
            fastq_in=trimmed,
            bam_out=bam_out,
            threads=threads,
            seed=seed,
            runner=runner,
        )
    else:
        if resolved.umi_length > 0:
            raise NotImplementedError(
                f"Sample {sample.sample!r}: paired-end + UMI is not yet "
                "supported. Either preprocess UMIs into the read name "
                "before invoking from-FASTQ mode, or pass --de-table with "
                "your own pre-computed DE results."
            )
        trimmed_r1 = sample_dir / f"{sample.sample}.trimmed_R1.fq.gz"
        trimmed_r2 = sample_dir / f"{sample.sample}.trimmed_R2.fq.gz"
        _run_cutadapt_pe(
            fastq_r1=sample.r1,
            fastq_r2=sample.r2,
            out_r1=trimmed_r1,
            out_r2=trimmed_r2,
            adapter=resolved.adapter,
            threads=threads,
            min_length=15,
            quality=20,
            runner=runner,
        )
        total, aligned = _run_bowtie2_pe(
            bt2_index=bt2_index,
            fastq_r1=trimmed_r1,
            fastq_r2=trimmed_r2,
            bam_out=bam_out,
            threads=threads,
            seed=seed,
            runner=runner,
        )

    counts = count_per_transcript(bam_out)
    return SampleAlignmentResult(
        sample=sample.sample,
        bam_path=bam_out,
        counts=counts,
        paired=sample.paired,
        total_reads=total,
        aligned_reads=aligned,
        resolved_kit=resolved,
    )
