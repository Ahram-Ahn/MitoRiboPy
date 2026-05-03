"""Generate the smoke-test FASTQs and bowtie2 index from human_mt_tiny.fa.

Deterministic (seed 42). Re-running the script overwrites the previous
outputs. Reads are simulated as perfectly periodic 32-nt mt-Ribo-seq
footprints with the canonical 5' P-site offset of 12 nt — i.e. each
read starts at a codon-aligned position 12 nt before its CDS coordinate
of interest. The KO sample additionally carries a 5-fold pile-up at
the 25th codon of MT-CO1 to mimic a stalling event so the smoke run
exercises the codon-correlation panel.

Usage:
    cd examples/smoke
    python generate_smoke_fastqs.py
"""

from __future__ import annotations

import gzip
import random
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

# Per-sample read budget. Small enough that bowtie2 finishes in seconds,
# large enough that periodicity / coverage panels look populated.
READS_PER_SAMPLE = 3000
READ_LENGTH = 32
P_SITE_OFFSET = 12       # 5' end → P-site, canonical mt offset
ADAPTER_SEQ = "CTGTAGGCACCATCAAT"  # Illumina sRNA-like; cutadapt --trim removes it.
SEED = 42


def _read_fasta(path: Path) -> dict[str, str]:
    """Tiny FASTA reader (Biopython is overkill for 3 contigs)."""
    out: dict[str, str] = {}
    name = None
    seq: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(seq)
            name = line[1:].split()[0].strip()
            seq = []
        else:
            seq.append(line.strip())
    if name is not None:
        out[name] = "".join(seq)
    return out


def _qual(read_length: int) -> str:
    """Phred 30 across the whole read — uniform high quality so the
    base-quality filter never kicks in on the smoke fixture."""
    return "?" * read_length  # ASCII '?' = phred 30


def _simulate_reads(
    refs: dict[str, str],
    *,
    sample: str,
    n_reads: int,
    rng: random.Random,
    add_stall: bool = False,
) -> list[tuple[str, str, str]]:
    """Return a list of (qname, seq, qual) tuples.

    Reads are placed so the 5' end is `P_SITE_OFFSET` upstream of a
    codon-aligned position. If `add_stall=True`, 5x extra reads land at
    codon 25 of MT-CO1 (mimics a polyproline-style stall).
    """
    rows: list[tuple[str, str, str]] = []
    qual = _qual(READ_LENGTH)
    contigs = list(refs.items())
    n_per_contig = n_reads // len(contigs)
    counter = 0
    for contig, seq in contigs:
        max_start = len(seq) - READ_LENGTH
        for _ in range(n_per_contig):
            # Pick a codon-aligned downstream position; place the 5'
            # end P_SITE_OFFSET upstream so the P-site lands on the
            # codon (matches the real mt-Ribo-seq geometry).
            codon_index = rng.randint(1, max(2, (max_start - P_SITE_OFFSET) // 3))
            start = codon_index * 3 - P_SITE_OFFSET
            start = max(0, min(start, max_start))
            read_seq = seq[start:start + READ_LENGTH]
            qname = f"@smoke:{sample}:{contig}:{counter}"
            rows.append((qname, read_seq, qual))
            counter += 1

    if add_stall and "MT-CO1" in refs:
        seq = refs["MT-CO1"]
        max_start = len(seq) - READ_LENGTH
        # Codon 25 → P-site at position 75; 5' end at position 75 - 12.
        stall_start = max(0, min(75 - P_SITE_OFFSET, max_start))
        read_seq = seq[stall_start:stall_start + READ_LENGTH]
        for _ in range(150):
            qname = f"@smoke:{sample}:MT-CO1:stall:{counter}"
            rows.append((qname, read_seq, qual))
            counter += 1
    return rows


def _add_adapter_and_write(
    rows: list[tuple[str, str, str]], out_path: Path,
) -> None:
    """Append a real adapter sequence + write gzipped FASTQ.

    Cutadapt's auto-detection should pick this adapter up; the smoke
    config relies on the default auto-detection path.
    """
    with gzip.open(out_path, "wt", encoding="utf-8") as h:
        for qname, seq, qual in rows:
            full_seq = seq + ADAPTER_SEQ
            full_qual = qual + ("?" * len(ADAPTER_SEQ))
            h.write(f"{qname}\n{full_seq}\n+\n{full_qual}\n")


def _build_bowtie2_index(reference: Path, prefix: Path) -> None:
    """Run bowtie2-build. Errors out if the tool is missing."""
    if shutil.which("bowtie2-build") is None:
        sys.exit(
            "ERROR: bowtie2-build not on PATH. Install bowtie2 first "
            "(e.g. `conda install -c bioconda bowtie2`)."
        )
    subprocess.run(
        ["bowtie2-build", str(reference), str(prefix)],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def main() -> int:
    reference = HERE / "human_mt_tiny.fa"
    if not reference.is_file():
        sys.exit(f"ERROR: reference {reference} missing.")
    refs = _read_fasta(reference)
    if not refs:
        sys.exit("ERROR: reference FASTA is empty.")

    fastq_dir = HERE / "fastqs"
    fastq_dir.mkdir(exist_ok=True)

    rng = random.Random(SEED)
    samples = (
        ("WT_smoke_1", False),
        ("KO_smoke_1", True),  # KO carries the synthetic stall.
    )
    for sample, add_stall in samples:
        rows = _simulate_reads(
            refs,
            sample=sample,
            n_reads=READS_PER_SAMPLE,
            rng=rng,
            add_stall=add_stall,
        )
        out = fastq_dir / f"{sample}.fastq.gz"
        _add_adapter_and_write(rows, out)
        print(f"  wrote {out.relative_to(HERE)}  ({len(rows)} reads)")

    index_dir = HERE / "bowtie2_index"
    index_dir.mkdir(exist_ok=True)
    index_prefix = index_dir / "human_mt_tiny"
    _build_bowtie2_index(reference, index_prefix)
    print(f"  wrote bowtie2 index under {index_dir.relative_to(HERE)}")

    print("\nNext step:")
    print("  mitoribopy all --config pipeline_config.smoke.yaml --output results/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
