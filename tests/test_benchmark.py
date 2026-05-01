"""Tests for `mitoribopy benchmark` (P2.13)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from mitoribopy import cli
from mitoribopy.cli import benchmark as benchmark_cli


# ---------- _measure_max_rss_mb --------------------------------------------


def test_measure_max_rss_mb_returns_positive_float() -> None:
    rss = benchmark_cli._measure_max_rss_mb()
    assert isinstance(rss, float)
    assert rss > 0


# ---------- _directory_size_mb ---------------------------------------------


def test_directory_size_mb_sums_file_bytes(tmp_path: Path) -> None:
    (tmp_path / "a.bin").write_bytes(b"x" * (1024 * 1024))  # 1 MB
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "b.bin").write_bytes(b"y" * (512 * 1024))  # 0.5 MB
    size = benchmark_cli._directory_size_mb(tmp_path)
    assert 1.4 < size < 1.6  # roughly 1.5 MB


# ---------- _subsample_fastq_inputs ----------------------------------------


def test_subsample_fastq_inputs_writes_subsampled_copies(
    tmp_path: Path,
) -> None:
    """Real-FASTQ subsample: 100 records in, 10 out per file."""
    src_dir = tmp_path / "raw"
    src_dir.mkdir()
    for fname in ("a.fq", "b.fq"):
        with (src_dir / fname).open("w") as out:
            for i in range(100):
                out.write(f"@r{i}\nACGT\n+\nIIII\n")

    cfg = {"align": {"fastq_dir": str(src_dir)}}
    work = tmp_path / "subsamples"
    out_cfg = benchmark_cli._subsample_fastq_inputs(
        cfg, n_reads=10, seed=1, work_dir=work
    )
    fastqs = out_cfg["align"]["fastq"]
    assert len(fastqs) == 2
    for path in fastqs:
        p = Path(path)
        assert p.exists()
        records = sum(1 for line in p.read_text().splitlines() if line.startswith("@"))
        assert records == 10
    # fastq_dir was nulled (replaced by the explicit list).
    assert out_cfg["align"]["fastq_dir"] is None


def test_subsample_fastq_inputs_handles_explicit_list(tmp_path: Path) -> None:
    src = tmp_path / "raw.fq"
    with src.open("w") as out:
        for i in range(50):
            out.write(f"@r{i}\nACGT\n+\nIIII\n")
    cfg = {"align": {"fastq": [str(src)]}}
    out_cfg = benchmark_cli._subsample_fastq_inputs(
        cfg, n_reads=5, seed=1, work_dir=tmp_path / "out"
    )
    sampled = Path(out_cfg["align"]["fastq"][0])
    records = sum(1 for line in sampled.read_text().splitlines() if line.startswith("@"))
    assert records == 5


# ---------- benchmark CLI ---------------------------------------------------


def test_benchmark_cli_writes_tsv_and_summary_md(
    tmp_path: Path, monkeypatch
) -> None:
    cfg_path = tmp_path / "c.yaml"
    cfg_path.write_text("rpf:\n  strain: h\n  fasta: /tmp/tx.fa\n")

    def fake_rpf(argv):
        out = Path(tmp_path / "results" / "rpf")
        out.mkdir(parents=True, exist_ok=True)
        (out / "rpf_counts.tsv").write_text(
            "# schema_version: 1.0\nsample\tgene\tcount\n"
        )
        (out / "run_settings.json").write_text("{}")
        return 0

    from mitoribopy.cli import rpf as rpf_cli
    monkeypatch.setattr(rpf_cli, "run", fake_rpf)

    rc = cli.main([
        "benchmark",
        "--config", str(cfg_path),
        "--output", str(tmp_path / "results"),
        "--threads", "2",
    ])
    assert rc == 0

    tsv_path = tmp_path / "results" / "benchmark.tsv"
    md_path = tmp_path / "results" / "benchmark_summary.md"
    assert tsv_path.is_file()
    assert md_path.is_file()
    raw = tsv_path.read_text().splitlines()
    assert raw[0].startswith("# schema_version:")
    rows = [r for r in raw if not r.startswith("#")]
    header = rows[0].split("\t")
    assert header[:3] == ["stage", "status", "wall_time_sec"]
    # 4 rows: align, rpf, rnaseq, total.
    assert len(rows) - 1 == 4
    total_row = rows[-1].split("\t")
    assert total_row[0] == "total"
    # threads column populated.
    threads_idx = header.index("threads")
    assert total_row[threads_idx] == "2"

    md = md_path.read_text()
    assert "MitoRiboPy benchmark summary" in md
    assert "Per-stage timing" in md
    assert "Resource totals" in md


def test_benchmark_cli_propagates_inner_failure_exit_code(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    cfg_path = tmp_path / "c.yaml"
    cfg_path.write_text("rpf:\n  strain: h\n  fasta: /tmp/tx.fa\n")

    def fake_rpf(argv):
        return 7

    from mitoribopy.cli import rpf as rpf_cli
    monkeypatch.setattr(rpf_cli, "run", fake_rpf)

    rc = cli.main([
        "benchmark",
        "--config", str(cfg_path),
        "--output", str(tmp_path / "results"),
    ])
    assert rc == 7
    err = capsys.readouterr().err
    assert "exited 7" in err
    # benchmark.tsv is still written even when the inner run failed.
    assert (tmp_path / "results" / "benchmark.tsv").is_file()


def test_benchmark_cli_listed_in_top_level_help(capsys) -> None:
    rc = cli.main(["--help"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "benchmark" in out


def test_benchmark_cli_subsample_actually_subsamples(
    tmp_path: Path, monkeypatch
) -> None:
    """End-to-end: --subsample N rewrites the canonical config to point
    at smaller FASTQs under <output>/.benchmark_subsamples/."""
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    src = raw_dir / "s1.fq"
    with src.open("w") as out:
        for i in range(200):
            out.write(f"@r{i}\nACGT\n+\nIIII\n")

    cfg_path = tmp_path / "c.yaml"
    cfg_path.write_text(
        "align:\n  kit_preset: auto\n"
        f"  fastq_dir: {raw_dir}\n"
        "rpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
    )

    def fake_align(argv):
        out = Path(tmp_path / "results" / "align")
        out.mkdir(parents=True, exist_ok=True)
        (out / "read_counts.tsv").write_text("sample\n")
        (out / "run_settings.json").write_text("{}")
        return 0

    def fake_rpf(argv):
        out = Path(tmp_path / "results" / "rpf")
        out.mkdir(parents=True, exist_ok=True)
        (out / "rpf_counts.tsv").write_text(
            "# schema_version: 1.0\nsample\tgene\tcount\n"
        )
        (out / "run_settings.json").write_text("{}")
        return 0

    from mitoribopy.cli import align as align_cli
    from mitoribopy.cli import rpf as rpf_cli
    monkeypatch.setattr(align_cli, "run", fake_align)
    monkeypatch.setattr(rpf_cli, "run", fake_rpf)

    rc = cli.main([
        "benchmark",
        "--config", str(cfg_path),
        "--output", str(tmp_path / "results"),
        "--subsample", "20",
        "--seed", "7",
    ])
    assert rc == 0
    sub_dir = tmp_path / "results" / ".benchmark_subsamples"
    assert sub_dir.is_dir()
    sampled = sub_dir / "s1.fq"
    assert sampled.exists()
    records = sum(
        1 for line in sampled.read_text().splitlines() if line.startswith("@")
    )
    assert records == 20

    # benchmark.tsv records the subsample read count.
    tsv = (tmp_path / "results" / "benchmark.tsv").read_text()
    assert "\t20\t" in tsv or tsv.endswith("\t20\n") or "\t20\n" in tsv
