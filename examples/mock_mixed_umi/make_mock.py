"""Generate a tiny mixed-UMI mock dataset for end-to-end smoke tests.

Produces (under the script's directory):

* ``human-mt-mRNA.mini.fasta`` - 3 mt-mRNA transcripts (ND1, COX1, ATP86)
  copied verbatim from the project's full reference, kept as a subset to
  keep ``bowtie2-build`` and the alignment loop fast.
* ``human_rrna.mini.fa``       - a small synthetic rRNA decoy
  (one ~600 nt random sequence with a fixed seed).
* ``sample_UMI.fq.gz``         - 2,000 reads with an 8 nt 5' UMI prepended
  followed by a 30 nt mt-mRNA fragment + the Illumina TruSeq R1 adapter
  (``AGATCGGAAGAGCACACGTCTGAACTCCAGTCA``).
* ``sample_NOUMI.fq.gz``       - 2,000 already-trimmed (SRA-style) reads
  with no adapter and no UMI, just the 30 nt fragment.

Both FASTQs deliberately overlap the same set of fragments so we can
sanity-check that dedup vs no-dedup produces the expected post-dedup
counts (UMI sample collapses identical fragments by their 8 nt tag;
no-UMI sample keeps every duplicate because dedup is skipped).

Run: ``python make_mock.py`` (idempotent; overwrites existing outputs).
"""

from __future__ import annotations

import gzip
import os
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[2]
FULL_FASTA = PROJECT_ROOT / "input_data" / "human-mt-mRNA.fasta"

ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
UMI_LEN = 8
FRAGMENT_LEN = 30
N_READS = 2000
SEED = 42

KEEP_TRANSCRIPTS = {"ND1", "COX1", "ATP86"}


def parse_fasta(path: Path) -> dict[str, str]:
    """Tiny FASTA parser (no biopython dep)."""
    records: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            records[name] = "".join(chunks)
    return records


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("w") as handle:
        for name, seq in records.items():
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                handle.write(seq[i : i + 60] + "\n")


def random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choices("ACGT", k=length))


def build_fragments(
    records: dict[str, str], n_reads: int, rng: random.Random
) -> list[tuple[str, str]]:
    """Pick (transcript, fragment) pairs uniformly across the kept transcripts.

    Deliberately reuses a small pool of fragment offsets so the UMI dedup
    step has duplicates to collapse (different UMI on the same coordinate
    -> kept; same UMI on the same coordinate -> collapsed).
    """
    transcripts = sorted(records)
    fragments: list[tuple[str, str]] = []
    for _ in range(n_reads):
        tx = rng.choice(transcripts)
        seq = records[tx]
        max_start = max(0, len(seq) - FRAGMENT_LEN)
        # Bias toward a small set of starts so we get duplicates.
        start = rng.choice(range(0, max_start, 5))
        fragments.append((tx, seq[start : start + FRAGMENT_LEN]))
    return fragments


def fastq_record(name: str, seq: str) -> str:
    qual = "I" * len(seq)
    return f"@{name}\n{seq}\n+\n{qual}\n"


def write_umi_fastq(
    out_path: Path, fragments: list[tuple[str, str]], rng: random.Random
) -> None:
    """Each read = 8nt UMI + 30nt fragment + adapter, gzipped."""
    with gzip.open(out_path, "wt") as handle:
        for index, (tx, fragment) in enumerate(fragments, start=1):
            umi = random_dna(UMI_LEN, rng)
            seq = umi + fragment + ADAPTER
            handle.write(fastq_record(f"umi_{index}_{tx}", seq))


def write_no_umi_fastq(
    out_path: Path, fragments: list[tuple[str, str]]
) -> None:
    """Each read = bare 30nt fragment (already trimmed, no UMI), gzipped."""
    with gzip.open(out_path, "wt") as handle:
        for index, (tx, fragment) in enumerate(fragments, start=1):
            handle.write(fastq_record(f"noumi_{index}_{tx}", fragment))


def main() -> None:
    rng = random.Random(SEED)

    if not FULL_FASTA.exists():
        raise SystemExit(
            f"Reference FASTA not found at {FULL_FASTA}. Adjust the path "
            "or re-run from a checkout that includes input_data/."
        )
    full = parse_fasta(FULL_FASTA)
    mini = {name: seq for name, seq in full.items() if name in KEEP_TRANSCRIPTS}
    if set(mini) != KEEP_TRANSCRIPTS:
        missing = KEEP_TRANSCRIPTS - set(mini)
        raise SystemExit(
            f"Reference FASTA is missing expected transcripts: {missing}. "
            "Update KEEP_TRANSCRIPTS to match what your reference ships."
        )

    write_fasta(HERE / "human-mt-mRNA.mini.fasta", mini)

    rrna_seq = random_dna(600, rng)
    write_fasta(HERE / "human_rrna.mini.fa", {"mock_rrna_1": rrna_seq})

    fragments = build_fragments(mini, N_READS, rng)
    write_umi_fastq(HERE / "sample_UMI.fq.gz", fragments, rng)
    write_no_umi_fastq(HERE / "sample_NOUMI.fq.gz", fragments)

    print("Wrote mock dataset under", HERE)
    for path in sorted(HERE.iterdir()):
        if path.suffix in {".fasta", ".fa", ".gz"}:
            print(f"  {path.name}: {os.path.getsize(path)} bytes")


if __name__ == "__main__":
    main()
