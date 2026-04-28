"""Stream-split a FASTQ into two pseudo-replicates by record parity.

Used by from-FASTQ mode when a condition has only one biological sample
on the RNA-seq or Ribo-seq side: pyDESeq2 needs at least two replicates
per condition to estimate dispersion, and the existing delta-TE math
falls back to a "single_replicate_no_statistics" note when there is no
internal Ribo log2FC available.

Splitting by record parity (even-indexed records → rep1, odd-indexed →
rep2) keeps the position / fragment-length distribution identical
across the two halves, which is what we want when the goal is "let the
math run" rather than "estimate biological variance" — the latter
genuinely needs more libraries, not more software.

Both pseudo-replicates are obviously NOT biological replicates: they
share the exact same read population minus the order. DESeq2's
dispersion estimates will be artificially low. The orchestrator emits
a stderr WARNING every time a split happens so the user is never
surprised.

Implementation note: large RNA-seq FASTQs (1–4 GB compressed,
50–110M reads) take ~5–10 min EACH through pure-Python text-mode
gzip. We therefore shell out to native ``gunzip`` / ``gzip`` for the
codec (typically 5–10× faster than Python's gzip module) and keep
only the line-routing loop in Python — operating on raw bytes,
which avoids UTF-8 decoding overhead. When ``gunzip`` / ``gzip`` are
not on ``$PATH`` we fall back to the Python ``gzip`` module so the
splitter still works on minimal environments (it is just slower
there).
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
from pathlib import Path
from typing import IO

from .fastq_pairing import FastqSample


def _open_input_bytes(
    path: Path,
) -> tuple[IO[bytes], "subprocess.Popen | None"]:
    """Open *path* for reading bytes; native ``gunzip -c`` if available."""
    p = str(path)
    if not p.endswith(".gz"):
        return open(p, "rb"), None
    gunzip_bin = shutil.which("gunzip") or shutil.which("gzip")
    if gunzip_bin is None:
        return gzip.open(p, "rb"), None
    proc = subprocess.Popen(
        [gunzip_bin, "-c", p], stdout=subprocess.PIPE
    )
    assert proc.stdout is not None
    return proc.stdout, proc


def _open_output_bytes(
    path: Path,
) -> tuple[IO[bytes], "subprocess.Popen | None", "IO[bytes] | None"]:
    """Open *path* for writing bytes; native ``gzip -c`` if available.

    Returns a 3-tuple ``(writable, gzip_proc, out_file)``. The latter
    two are non-``None`` only when we shelled out; the caller must
    close them in the same order.
    """
    p = str(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not p.endswith(".gz"):
        return open(p, "wb"), None, None
    gzip_bin = shutil.which("gzip")
    if gzip_bin is None:
        return gzip.open(p, "wb"), None, None
    out_file = open(p, "wb")
    proc = subprocess.Popen(
        [gzip_bin, "-c"], stdin=subprocess.PIPE, stdout=out_file
    )
    assert proc.stdin is not None
    return proc.stdin, proc, out_file


def stream_split_by_parity(src: Path, even_dst: Path, odd_dst: Path) -> int:
    """Split a FASTQ into two halves by record parity.

    Even-indexed records (0, 2, 4, ...) go to ``even_dst``; odd-indexed
    records go to ``odd_dst``. Returns the total number of records
    written across both halves.

    Performance: shells out to ``gunzip`` / ``gzip`` for the codec when
    they are on ``$PATH`` (typical conda / homebrew install). On a
    50M-read FASTQ this is ~5–10× faster than the pure-Python gzip
    module.
    """
    src = Path(src)
    even_dst = Path(even_dst)
    odd_dst = Path(odd_dst)

    src_fh, src_proc = _open_input_bytes(src)
    even_fh, even_proc, even_out = _open_output_bytes(even_dst)
    odd_fh, odd_proc, odd_out = _open_output_bytes(odd_dst)

    n = 0
    truncated = False
    try:
        record: list[bytes] = []
        for line in src_fh:
            record.append(line)
            if len(record) == 4:
                target = even_fh if (n % 2 == 0) else odd_fh
                target.writelines(record)
                record.clear()
                n += 1
        if record:
            truncated = True
    finally:
        # Close in pipeline order: source first (stops gunzip), then
        # the gzip stdin pipes (lets each gzip flush + finalize), then
        # wait + close output files.
        src_fh.close()
        if src_proc is not None:
            src_proc.wait()

        even_fh.close()
        if even_proc is not None:
            even_proc.wait()
        if even_out is not None:
            even_out.close()

        odd_fh.close()
        if odd_proc is not None:
            odd_proc.wait()
        if odd_out is not None:
            odd_out.close()

    if truncated:
        raise ValueError(
            f"Truncated FASTQ record at end of {src} "
            f"(got {len(record)} lines after the last complete record)"
        )
    return n


def split_sample_into_pseudo_replicates(
    sample: FastqSample, workdir: Path
) -> tuple[FastqSample, FastqSample]:
    """Stream-split ``sample`` into two pseudo-replicates by record parity.

    For paired-end samples both R1 and R2 are split; record N's R1
    always lands on the same parity as record N's R2 so the pairs stay
    matched.

    Output FASTQs are written under ``workdir/<rep_name>/`` preserving
    the original filename so users can recognise the source.
    """
    workdir = Path(workdir)
    rep1_name = f"{sample.sample}_rep1"
    rep2_name = f"{sample.sample}_rep2"
    rep1_dir = workdir / rep1_name
    rep2_dir = workdir / rep2_name

    rep1_r1 = rep1_dir / sample.r1.name
    rep2_r1 = rep2_dir / sample.r1.name
    stream_split_by_parity(sample.r1, rep1_r1, rep2_r1)

    rep1_r2: Path | None = None
    rep2_r2: Path | None = None
    if sample.r2 is not None:
        rep1_r2 = rep1_dir / sample.r2.name
        rep2_r2 = rep2_dir / sample.r2.name
        stream_split_by_parity(sample.r2, rep1_r2, rep2_r2)

    return (
        FastqSample(sample=rep1_name, r1=rep1_r1, r2=rep1_r2),
        FastqSample(sample=rep2_name, r1=rep2_r1, r2=rep2_r2),
    )
