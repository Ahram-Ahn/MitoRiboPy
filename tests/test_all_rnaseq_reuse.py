"""P0.2 — `mitoribopy all` reuses rpf_counts.tsv in rnaseq from-FASTQ mode.

Today, the rnaseq from-FASTQ flow could re-align Ribo FASTQs even
though the rpf stage in the same `all` invocation just produced
rpf_counts.tsv from the same library. The reuse-by-default change
threads `--upstream-rpf-counts <path>` from the orchestrator into the
rnaseq subcommand; the in-tree `_run_from_fastq` synthesizes per-sample
counts from that file and skips bowtie2 entirely on the Ribo side.
Users opt out with `rnaseq.recount_ribo_fastq: true`.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pytest

from mitoribopy import cli
from mitoribopy.cli import all_ as all_cli
from mitoribopy.cli import rnaseq as rnaseq_cli


# ---------- _auto_wire_paths sets upstream_rpf_counts ----------------------


def test_auto_wire_sets_upstream_rpf_counts_in_from_fastq_mode(
    tmp_path: Path,
) -> None:
    cfg: dict = {
        "rpf": {"strain": "h", "fasta": "/tmp/tx.fa"},
        "rnaseq": {
            "rna_fastq": ["rna1.fq.gz"],
            "gene_id_convention": "hgnc",
            "condition_a": "WT",
            "condition_b": "KO",
        },
    }
    all_cli._auto_wire_paths(
        cfg,
        run_root=tmp_path / "results",
        has_align=False,
        has_rpf=True,
        has_rnaseq=True,
        has_de_table=False,
        has_fastq_mode=True,
    )
    expected = str(tmp_path / "results" / "rpf" / "rpf_counts.tsv")
    assert cfg["rnaseq"]["upstream_rpf_counts"] == expected


def test_auto_wire_skips_upstream_rpf_counts_when_recount_set(
    tmp_path: Path,
) -> None:
    """User opt-out: `rnaseq.recount_ribo_fastq: true` disables reuse wiring."""
    cfg: dict = {
        "rpf": {"strain": "h", "fasta": "/tmp/tx.fa"},
        "rnaseq": {
            "rna_fastq": ["rna1.fq.gz"],
            "recount_ribo_fastq": True,
            "gene_id_convention": "hgnc",
            "condition_a": "WT",
            "condition_b": "KO",
        },
    }
    all_cli._auto_wire_paths(
        cfg,
        run_root=tmp_path / "results",
        has_align=False,
        has_rpf=True,
        has_rnaseq=True,
        has_de_table=False,
        has_fastq_mode=True,
    )
    assert "upstream_rpf_counts" not in cfg["rnaseq"]


def test_auto_wire_does_not_set_upstream_for_de_table_flow(
    tmp_path: Path,
) -> None:
    """The de_table flow has its own --ribo-dir wiring; no upstream_rpf_counts."""
    cfg: dict = {
        "rpf": {"strain": "h", "fasta": "/tmp/tx.fa"},
        "rnaseq": {
            "de_table": "de.tsv",
            "gene_id_convention": "hgnc",
            "reference_gtf": "ref.fa",
        },
    }
    all_cli._auto_wire_paths(
        cfg,
        run_root=tmp_path / "results",
        has_align=False,
        has_rpf=True,
        has_rnaseq=True,
        has_de_table=True,
        has_fastq_mode=False,
    )
    assert "upstream_rpf_counts" not in cfg["rnaseq"]
    assert cfg["rnaseq"]["ribo-dir"] == str(tmp_path / "results" / "rpf")


def test_auto_wire_skips_upstream_when_no_rpf_stage(tmp_path: Path) -> None:
    cfg: dict = {
        "rnaseq": {
            "rna_fastq": ["rna1.fq.gz"],
            "gene_id_convention": "hgnc",
            "condition_a": "WT",
            "condition_b": "KO",
        },
    }
    all_cli._auto_wire_paths(
        cfg,
        run_root=tmp_path / "results",
        has_align=False,
        has_rpf=False,
        has_rnaseq=True,
        has_de_table=False,
        has_fastq_mode=True,
    )
    assert "upstream_rpf_counts" not in cfg["rnaseq"]


# ---------- argv passthrough ------------------------------------------------


def test_all_passes_upstream_rpf_counts_argv_to_rnaseq(
    tmp_path: Path, monkeypatch
) -> None:
    """End-to-end: the orchestrator must hand --upstream-rpf-counts to rnaseq."""
    captured: dict[str, list[str]] = {}

    def fake_align(argv):
        out = Path(tmp_path / "results" / "align")
        out.mkdir(parents=True, exist_ok=True)
        (out / "read_counts.tsv").write_text("sample\n")
        (out / "run_settings.json").write_text("{}")
        return 0

    def fake_rpf(argv):
        out = Path(tmp_path / "results" / "rpf")
        out.mkdir(parents=True, exist_ok=True)
        (out / "rpf_counts.tsv").write_text("sample\tgene\tcount\n")
        (out / "run_settings.json").write_text("{}")
        return 0

    def fake_rnaseq(argv):
        captured["rnaseq"] = list(argv)
        out = Path(tmp_path / "results" / "rnaseq")
        out.mkdir(parents=True, exist_ok=True)
        (out / "delta_te.tsv").write_text("gene\n")
        (out / "run_settings.json").write_text("{}")
        return 0

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rpf as rpf_cli
    monkeypatch.setattr(align_cli, "run", fake_align)
    monkeypatch.setattr(rpf_cli, "run", fake_rpf)
    monkeypatch.setattr(rnaseq_cli, "run", fake_rnaseq)

    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "align:\n  adapter: TGGAATTCTCGGGTGCCAAGG\n"
        "rpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
        "rnaseq:\n"
        "  rna_fastq:\n    - rna1.fq.gz\n"
        "  reference_fasta: /tmp/ref.fa\n"
        "  gene_id_convention: hgnc\n"
        "  condition_a: WT\n  condition_b: KO\n"
    )
    exit_code = cli.main([
        "all",
        "--config", str(cfg),
        "--output", str(tmp_path / "results"),
    ])
    assert exit_code == 0
    rnaseq_argv = captured["rnaseq"]
    assert "--upstream-rpf-counts" in rnaseq_argv
    idx = rnaseq_argv.index("--upstream-rpf-counts")
    expected = str(tmp_path / "results" / "rpf" / "rpf_counts.tsv")
    assert rnaseq_argv[idx + 1] == expected


def test_all_recount_opt_out_drops_upstream_argv(
    tmp_path: Path, monkeypatch
) -> None:
    captured: dict[str, list[str]] = {}

    def fake_align(argv):
        out = Path(tmp_path / "results" / "align")
        out.mkdir(parents=True, exist_ok=True)
        (out / "read_counts.tsv").write_text("sample\n")
        (out / "run_settings.json").write_text("{}")
        return 0

    def fake_rpf(argv):
        out = Path(tmp_path / "results" / "rpf")
        out.mkdir(parents=True, exist_ok=True)
        (out / "rpf_counts.tsv").write_text("sample\tgene\tcount\n")
        (out / "run_settings.json").write_text("{}")
        return 0

    def fake_rnaseq(argv):
        captured["rnaseq"] = list(argv)
        out = Path(tmp_path / "results" / "rnaseq")
        out.mkdir(parents=True, exist_ok=True)
        (out / "delta_te.tsv").write_text("gene\n")
        (out / "run_settings.json").write_text("{}")
        return 0

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rpf as rpf_cli
    monkeypatch.setattr(align_cli, "run", fake_align)
    monkeypatch.setattr(rpf_cli, "run", fake_rpf)
    monkeypatch.setattr(rnaseq_cli, "run", fake_rnaseq)

    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "align:\n  adapter: TGGAATTCTCGGGTGCCAAGG\n"
        "rpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
        "rnaseq:\n"
        "  rna_fastq:\n    - rna1.fq.gz\n"
        "  reference_fasta: /tmp/ref.fa\n"
        "  gene_id_convention: hgnc\n"
        "  recount_ribo_fastq: true\n"
        "  condition_a: WT\n  condition_b: KO\n"
    )
    exit_code = cli.main([
        "all",
        "--config", str(cfg),
        "--output", str(tmp_path / "results"),
    ])
    assert exit_code == 0
    rnaseq_argv = captured["rnaseq"]
    assert "--upstream-rpf-counts" not in rnaseq_argv
    assert "--recount-ribo-fastq" in rnaseq_argv


# ---------- _synth_ribo_results_from_upstream_counts ------------------------


def test_synth_ribo_results_pivots_counts_per_sample(tmp_path: Path) -> None:
    counts = tmp_path / "rpf_counts.tsv"
    counts.write_text(
        "sample\tgene\tcount\n"
        "A1\tMT-ND1\t100\n"
        "A1\tMT-CO1\t200\n"
        "B1\tMT-ND1\t50\n"
    )
    results = rnaseq_cli._synth_ribo_results_from_upstream_counts(counts)
    assert [r.sample for r in results] == ["A1", "B1"]
    a1 = next(r for r in results if r.sample == "A1")
    assert a1.counts == {"MT-ND1": 100, "MT-CO1": 200}
    assert a1.total_reads == 300
    assert a1.aligned_reads == 300
    assert a1.paired is False
    # Sentinel kit must be detectable in provenance.
    assert a1.resolved_kit.kit == "upstream-rpf-counts"


# ---------- end-to-end via _run_from_fastq with reuse path ------------------


def test_run_from_fastq_reuse_path_skips_alignment(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """When upstream_rpf_counts is set and recount is False, _run_from_fastq
    must NOT call align_sample on Ribo samples and must NOT write
    <rnaseq>/rpf_counts.tsv (the upstream file is reused in place)."""
    pytest.importorskip("pandas")
    pytest.importorskip("pydeseq2")

    # Tiny synthetic counts so build_sample_sheet + run_deseq2 have data.
    # Two replicates per condition on each assay so pyDESeq2 can fit the
    # dispersion (n=1 per condition trips the singleton gate before the
    # reuse path is exercised).
    upstream = tmp_path / "rpf" / "rpf_counts.tsv"
    upstream.parent.mkdir(parents=True, exist_ok=True)
    upstream.write_text(
        "sample\tgene\tcount\n"
        "WT_Ribo_1\tMT-ND1\t100\n"
        "WT_Ribo_1\tMT-CO1\t250\n"
        "WT_Ribo_2\tMT-ND1\t105\n"
        "WT_Ribo_2\tMT-CO1\t245\n"
        "KO_Ribo_1\tMT-ND1\t110\n"
        "KO_Ribo_1\tMT-CO1\t300\n"
        "KO_Ribo_2\tMT-ND1\t108\n"
        "KO_Ribo_2\tMT-CO1\t295\n"
    )

    # Stub out align/index/heavy I/O. The reuse path should skip Ribo
    # alignment entirely; we still need RNA alignment to be fake to
    # keep the test free of bowtie2/cutadapt.
    from mitoribopy.rnaseq import alignment as align_mod
    from mitoribopy.rnaseq.alignment import SampleAlignmentResult
    from mitoribopy.align._types import ResolvedKit

    align_calls: list[str] = []

    fake_kit = ResolvedKit(
        kit="pretrimmed", adapter=None, umi_length=0, umi_position="5p"
    )

    def fake_align_sample(sample, *, bt2_index, workdir, threads):
        align_calls.append(sample.sample)
        return SampleAlignmentResult(
            sample=sample.sample,
            bam_path=Path("/dev/null"),
            counts={"MT-ND1": 80, "MT-CO1": 200},
            paired=False,
            total_reads=280,
            aligned_reads=280,
            resolved_kit=fake_kit,
        )

    def fake_build_index(ref, cache):
        return Path("/tmp/fake_index")

    monkeypatch.setattr(rnaseq_cli, "_run_from_fastq", rnaseq_cli._run_from_fastq)
    monkeypatch.setattr(align_mod, "align_sample", fake_align_sample)
    monkeypatch.setattr(align_mod, "build_bowtie2_index", fake_build_index)

    # Prepare RNA FASTQs (mutually-distinct conditions; need at least 2 to
    # satisfy DESeq2's minimum two-condition requirement). The files
    # don't actually have to contain data because align_sample is faked.
    # Filenames must avoid the bcl2fastq/_R1/_1 mate suffixes that
    # detect_samples canonicalises away — otherwise the sample names
    # the orchestrator sees will not match the condition_map keys.
    rna_dir = tmp_path / "rna"
    rna_dir.mkdir()
    (rna_dir / "WT_RNA_a.fq.gz").write_bytes(b"")
    (rna_dir / "WT_RNA_b.fq.gz").write_bytes(b"")
    (rna_dir / "KO_RNA_a.fq.gz").write_bytes(b"")
    (rna_dir / "KO_RNA_b.fq.gz").write_bytes(b"")
    cmap = tmp_path / "cmap.tsv"
    cmap.write_text(
        "sample\tcondition\n"
        "WT_RNA_a\tWT\nWT_RNA_b\tWT\n"
        "KO_RNA_a\tKO\nKO_RNA_b\tKO\n"
        "WT_Ribo_1\tWT\nWT_Ribo_2\tWT\n"
        "KO_Ribo_1\tKO\nKO_Ribo_2\tKO\n"
    )

    args = argparse.Namespace(
        rna_fastq=[str(rna_dir)],
        ribo_fastq=None,
        reference_fasta=str(tmp_path / "ref.fa"),
        bowtie2_index=None,
        workdir=str(tmp_path / "work"),
        align_threads=1,
        condition_map=str(cmap),
        condition_a="WT",
        condition_b="KO",
        no_auto_pseudo_replicate=False,
        allow_pseudo_replicates=False,
        upstream_rpf_counts=str(upstream),
        recount_ribo_fastq=False,
    )
    # _run_from_fastq writes args.de_table, args.ribo_counts, etc.
    (tmp_path / "ref.fa").write_bytes(b">MT-ND1\nA\n>MT-CO1\nC\n")

    out = tmp_path / "rnaseq_out"
    rc = rnaseq_cli._run_from_fastq(args, out)
    assert rc == 0  # ribo_results synthesised => continues past short-circuit

    # Side 1: align_sample was NOT called for Ribo samples — the reuse
    # path skips bowtie2 entirely. RNA samples were aligned (fake).
    ribo_samples_aligned = [s for s in align_calls if "Ribo" in s]
    assert ribo_samples_aligned == []

    # Side 2: <out>/rpf_counts.tsv must NOT exist. The upstream file
    # is the source of truth and we point args.ribo_counts at it.
    assert not (out / "rpf_counts.tsv").exists()
    assert args.ribo_counts == str(upstream)

    # Side 3: provenance reflects the reuse mode.
    prov = args._fastq_provenance
    assert prov["mode"] == "from-fastq-rna-with-upstream-rpf-counts"
    assert prov["upstream_rpf_counts"] == str(upstream)

    # Side 4: stderr banner mentions reuse so the user is never surprised.
    err = capsys.readouterr().err
    assert "reusing upstream rpf counts" in err
