"""End-to-end integration test for ``mitoribopy align``.

This test exercises the real cutadapt / bowtie2 / samtools / umi_tools
tools on tiny synthetic inputs. It is gated behind
``@pytest.mark.requires_tools`` and auto-skipped by ``conftest.py`` when
any of those tools is missing from ``PATH`` (which is the case in the
sandbox this repository was developed in).

The fixtures are generated deterministically at test time (no binary
files committed to the repo):

* ``mt_transcripts.fa``:  two synthetic 300-nt mt-mRNA records whose
                          headers match the naming convention used in
                          the annotation CSV (sequence_name).
* ``contam.fa``:          one synthetic 300-nt rRNA-like sequence.
* ``sample_A.fq.gz``:     120 reads:
                           - 60 substrings sampled from the mt FASTA
                             (should survive contam + align to mt),
                           - 30 substrings from the contam FASTA
                             (should be subtracted at contam step),
                           - 30 random sequences (should fail contam
                             and then fail mt align).
                          Every read has a 3' TruSeq small-RNA adapter
                          appended; cutadapt is expected to trim it.

The assertions verify the pipeline mechanics end-to-end: each stage of
the read-counts table has the expected ordering, and the BED6 output
exists with strand-aware rows.
"""

from __future__ import annotations

import gzip
import random
import subprocess
from pathlib import Path

import pytest


TRUSEQ_ADAPTER = "TGGAATTCTCGGGTGCCAAGG"


# ---------- fixture generators ----------------------------------------------


def _seeded_sequence(seed: int, length: int) -> str:
    """Return a deterministic ACGT string of *length* driven by *seed*."""
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, seq in records.items():
            handle.write(f">{name}\n{seq}\n")


def _write_fastq_gz(
    path: Path,
    *,
    mt_refs: dict[str, str],
    contam_refs: dict[str, str],
    rng_seed: int,
    n_mt: int,
    n_contam: int,
    n_random: int,
    read_length: int = 28,
) -> None:
    rng = random.Random(rng_seed)
    reads: list[tuple[str, str]] = []

    def _sample(seq: str, count: int, prefix: str) -> None:
        for i in range(count):
            start = rng.randrange(0, len(seq) - read_length)
            body = seq[start : start + read_length]
            reads.append((f"{prefix}_{i}", body + TRUSEQ_ADAPTER))

    for name, seq in mt_refs.items():
        _sample(seq, n_mt // len(mt_refs), f"mt_{name}")
    for name, seq in contam_refs.items():
        _sample(seq, n_contam // len(contam_refs), f"contam_{name}")
    for i in range(n_random):
        body = _seeded_sequence(rng_seed + 1000 + i, read_length)
        reads.append((f"random_{i}", body + TRUSEQ_ADAPTER))

    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for read_id, seq in reads:
            qual = "I" * len(seq)
            handle.write(f"@{read_id}\n{seq}\n+\n{qual}\n")


def _build_bowtie2_index(fasta: Path, prefix: Path) -> None:
    result = subprocess.run(
        ["bowtie2-build", "--quiet", str(fasta), str(prefix)],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"bowtie2-build failed on {fasta}: {result.stderr}"
        )


# ---------- the actual test --------------------------------------------------


@pytest.mark.requires_tools
def test_mitoribopy_align_end_to_end_on_tiny_synthetic_library(tmp_path):
    """Run 'mitoribopy align' against a tiny synthetic library end-to-end."""
    # 1) Build synthetic references.
    mt_fa = tmp_path / "mt_transcripts.fa"
    contam_fa = tmp_path / "contam.fa"
    mt_refs = {
        "ND1": _seeded_sequence(seed=1, length=300),
        "COX1": _seeded_sequence(seed=2, length=300),
    }
    contam_refs = {"RNR1_like": _seeded_sequence(seed=3, length=300)}
    _write_fasta(mt_fa, mt_refs)
    _write_fasta(contam_fa, contam_refs)

    # 2) Build bowtie2 indexes.
    idx_dir = tmp_path / "idx"
    idx_dir.mkdir()
    mt_index = idx_dir / "mt"
    contam_index = idx_dir / "contam"
    _build_bowtie2_index(mt_fa, mt_index)
    _build_bowtie2_index(contam_fa, contam_index)

    # 3) Generate a synthetic FASTQ sample.
    fq_dir = tmp_path / "fq"
    fq_dir.mkdir()
    sample_fq = fq_dir / "sample_A.fq.gz"
    _write_fastq_gz(
        sample_fq,
        mt_refs=mt_refs,
        contam_refs=contam_refs,
        rng_seed=42,
        n_mt=60,
        n_contam=30,
        n_random=30,
    )

    out_dir = tmp_path / "out"

    # 4) Invoke the CLI through its Python entry point (no subprocess so
    #    we get a clear traceback on failure).
    from mitoribopy import cli

    exit_code = cli.main(
        [
            "align",
            # v0.7.1 removed --kit-preset; pass the small-RNA adapter
            # explicitly (this is the same sequence the synthetic reads
            # were constructed with above).
            "--adapter",
            TRUSEQ_ADAPTER,
            "--library-strandedness",
            "unstranded",  # the synthetic reads are not strand-constrained
            "--fastq",
            str(sample_fq),
            "--contam-index",
            str(contam_index),
            "--mt-index",
            str(mt_index),
            "--output",
            str(out_dir),
            "--seed",
            "42",
        ]
    )
    assert exit_code == 0, "align orchestrator exited non-zero"

    # 5) Verify output directory layout.
    # `deduped/` is intentionally NOT in this list: when --dedup-strategy
    # auto resolves to 'skip' (the case for non-UMI libraries like this
    # synthetic one), the orchestrator does not write a duplicate
    # `deduped/<sample>.dedup.bam` — the upstream mapq.bam goes straight
    # into BED conversion. Documented in the README and validated
    # separately in test_align_cli.py.
    for subdir in ("trimmed", "contam_filtered", "aligned", "bed"):
        assert (out_dir / subdir).is_dir(), f"missing {subdir}/"
    # `deduped/` should NOT exist on a non-UMI synthetic library because
    # auto -> skip drops the redundant copy of the BAM.
    assert not (out_dir / "deduped").exists(), (
        "deduped/ should not be written when dedup auto-resolves to skip"
    )
    assert (out_dir / "read_counts.tsv").is_file()
    assert (out_dir / "run_settings.json").is_file()

    # 6) read_counts.tsv has a row for sample_A with the stage invariants.
    # The first line is a `# schema_version: ...` comment (introduced in
    # the v0.6+ schema-version banner); skip lines starting with '#' to
    # land on the column header.
    lines = [
        ln for ln in (out_dir / "read_counts.tsv").read_text().splitlines()
        if ln and not ln.startswith("#")
    ]
    header = lines[0].split("\t")
    data = dict(zip(header, lines[1].split("\t")))
    assert data["sample"] == "sample_A"
    total = int(data["total_reads"])
    post_trim = int(data["post_trim"])
    rrna = int(data["rrna_aligned"])
    post_rrna = int(data["post_rrna_filter"])
    mt_aligned = int(data["mt_aligned"])
    mapq = int(data["mt_aligned_after_mapq"])
    dedup = int(data["mt_aligned_after_dedup"])

    assert total == 120
    # Cutadapt should recover most reads (length 28 + adapter trimmed).
    assert post_trim > 0
    assert rrna + post_rrna == post_trim
    assert mt_aligned <= post_rrna
    assert dedup <= mapq <= mt_aligned
    # Expect most mt-sampled reads to align to the mt index.
    assert mt_aligned >= 30

    # 7) BED6 output is populated and strand-aware.
    bed_path = out_dir / "bed" / "sample_A.bed"
    assert bed_path.is_file()
    bed_rows = [row.split("\t") for row in bed_path.read_text().splitlines() if row]
    assert bed_rows, "BED file is empty"
    for row in bed_rows:
        assert len(row) == 6, f"BED row not BED6: {row}"
        assert row[0] in mt_refs, f"unexpected chrom: {row[0]}"
        assert row[5] in {"+", "-"}, f"invalid strand: {row[5]}"
