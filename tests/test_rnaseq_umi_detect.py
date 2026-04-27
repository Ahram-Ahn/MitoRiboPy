"""Per-position-entropy UMI detection on synthetic FASTQs."""

from __future__ import annotations

import gzip
import random
from pathlib import Path

from mitoribopy.rnaseq.umi_detect import detect_umi


_BASES = ("A", "C", "G", "T")


def _biased_body(rng: random.Random, length: int) -> str:
    """Return a low-entropy (codon-biased) body sequence."""
    # Mostly A with the occasional C — yields entropy < 1.0 across all positions.
    return "".join(rng.choices("AAAAAACG", k=length))


def _uniform(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(_BASES) for _ in range(length))


def _write_fastq(path: Path, reads: list[str], gzipped: bool = False) -> None:
    lines: list[str] = []
    for i, seq in enumerate(reads):
        lines.append(f"@read{i}")
        lines.append(seq)
        lines.append("+")
        lines.append("I" * len(seq))
    blob = "\n".join(lines) + "\n"
    if gzipped:
        with gzip.open(str(path), "wt") as h:
            h.write(blob)
    else:
        path.write_text(blob)


def test_detect_8nt_5p_umi(tmp_path: Path) -> None:
    rng = random.Random(0)
    reads = [_uniform(rng, 8) + _biased_body(rng, 50) for _ in range(200)]
    fq = tmp_path / "umi5p.fastq"
    _write_fastq(fq, reads)
    result = detect_umi(fq)
    assert result.length == 8
    assert result.position == "5p"
    assert result.n_reads_scanned == 200


def test_detect_6nt_3p_umi(tmp_path: Path) -> None:
    rng = random.Random(0)
    reads = [_biased_body(rng, 50) + _uniform(rng, 6) for _ in range(200)]
    fq = tmp_path / "umi3p.fastq"
    _write_fastq(fq, reads)
    result = detect_umi(fq)
    assert result.length == 6
    assert result.position == "3p"


def test_no_umi_returns_zero_length(tmp_path: Path) -> None:
    rng = random.Random(0)
    reads = [_biased_body(rng, 50) for _ in range(200)]
    fq = tmp_path / "noumi.fastq"
    _write_fastq(fq, reads)
    result = detect_umi(fq)
    assert result.length == 0


def test_detect_works_on_gzipped_input(tmp_path: Path) -> None:
    rng = random.Random(0)
    reads = [_uniform(rng, 8) + _biased_body(rng, 50) for _ in range(200)]
    fq = tmp_path / "umi5p.fastq.gz"
    _write_fastq(fq, reads, gzipped=True)
    result = detect_umi(fq)
    assert result.length == 8
    assert result.position == "5p"
