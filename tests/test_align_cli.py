"""Unit tests for ``mitoribopy.cli.align``.

The orchestrator's external-tool calls (cutadapt, bowtie2, umi_tools,
picard) and BAM operations (pysam) are replaced via monkeypatch so the
tests run without any bioinformatics tools installed and without real
BAM I/O.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from mitoribopy import cli
from mitoribopy.align._types import (
    AlignResult,
    ContamResult,
    CutadaptResult,
    DedupResult,
    ResolvedKit,
    SampleCounts,
    ToolInfo,
)
from mitoribopy.cli import align as align_cli
from mitoribopy.cli.align import _enumerate_fastqs, _sample_name, build_parser


# ---------- build_parser / help --------------------------------------------


def test_build_parser_help_mentions_new_flags() -> None:
    parser = build_parser()
    help_text = parser.format_help()
    assert "--kit-preset" in help_text
    assert "--adapter" in help_text
    assert "--umi-length" in help_text
    assert "--library-strandedness" in help_text
    assert "--contam-index" in help_text
    assert "--mt-index" in help_text
    assert "--dedup-strategy" in help_text
    assert "--i-understand-mark-duplicates-destroys-mt-ribo-seq-signal" in help_text


def test_align_subcommand_help_via_cli(capsys) -> None:
    with pytest.raises(SystemExit) as exc:
        cli.main(["align", "--help"])
    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert "mitoribopy align" in captured.out
    assert "--kit-preset" in captured.out


# ---------- _enumerate_fastqs + _sample_name --------------------------------


def test_enumerate_fastqs_globs_directory_and_deduplicates(tmp_path) -> None:
    d = tmp_path / "fq"
    d.mkdir()
    (d / "a.fq.gz").write_text("")
    (d / "b.fastq").write_text("")
    (d / "c.fq").write_text("")
    (d / "nope.txt").write_text("")

    # Pass d twice via --fastq to test dedup (absolute + relative combos
    # resolve to the same path).
    paths = _enumerate_fastqs(str(d), [str(d / "a.fq.gz")])
    names = sorted(p.name for p in paths)
    assert names == ["a.fq.gz", "b.fastq", "c.fq"]


def test_enumerate_fastqs_raises_for_missing_dir(tmp_path) -> None:
    with pytest.raises(FileNotFoundError):
        _enumerate_fastqs(str(tmp_path / "missing"), None)


@pytest.mark.parametrize(
    "filename, expected",
    [
        ("sample.fq.gz", "sample"),
        ("sample.fastq.gz", "sample"),
        ("sample.fq", "sample"),
        ("sample.fastq", "sample"),
        ("has.dots.in.name.fq.gz", "has.dots.in.name"),
    ],
)
def test_sample_name_strips_fastq_suffixes(filename, expected) -> None:
    assert _sample_name(Path(f"/tmp/{filename}")) == expected


# ---------- dry-run ---------------------------------------------------------


def test_dry_run_exits_zero_and_prints_plan_with_resolved_settings(
    capsys, tmp_path, monkeypatch
) -> None:
    # Sentinel to prove the orchestrator does NOT enter real execution.
    def explode(*args, **kwargs):
        raise AssertionError("tool_check must not run in --dry-run")

    monkeypatch.setattr(align_cli.tool_check, "ensure_tools_available", explode)

    exit_code = cli.main(
        [
            "align",
            "--dry-run",
            "--kit-preset",
            "illumina_smallrna",
            "--fastq-dir",
            str(tmp_path),
        ]
    )

    assert exit_code == 0
    out = capsys.readouterr().out
    assert "dry-run" in out
    assert "illumina_smallrna" in out
    assert "TGGAATTCTCGGGTGCCAAGG" in out  # resolved adapter
    assert "strand=forward" in out


def test_dry_run_fails_fast_on_custom_kit_without_adapter(capsys, tmp_path) -> None:
    exit_code = cli.main(
        [
            "align",
            "--dry-run",
            "--kit-preset",
            "custom",  # no --adapter
            "--fastq-dir",
            str(tmp_path),
        ]
    )
    assert exit_code == 2
    assert "--adapter" in capsys.readouterr().err


def test_dry_run_fails_fast_on_mark_duplicates_without_confirm_flag(
    capsys, tmp_path
) -> None:
    exit_code = cli.main(
        [
            "align",
            "--dry-run",
            "--kit-preset",
            "illumina_smallrna",
            "--dedup-strategy",
            "mark-duplicates",
            "--fastq-dir",
            str(tmp_path),
        ]
    )
    assert exit_code == 2
    err = capsys.readouterr().err
    assert "i-understand" in err


# ---------- non-dry-run validation ------------------------------------------


def test_non_dry_run_requires_output_contam_index_mt_index_inputs(
    capsys,
) -> None:
    exit_code = cli.main(["align", "--kit-preset", "illumina_smallrna"])
    err = capsys.readouterr().err
    assert exit_code == 2
    assert "--output" in err
    assert "--contam-index" in err
    assert "--mt-index" in err


# ---------- end-to-end orchestration (mocked) -------------------------------


@pytest.fixture
def mocked_align_step_functions(monkeypatch, tmp_path):
    """Replace every external-tool and pysam-backed call with a fake.

    The fakes write plausible on-disk artifacts so the orchestrator's
    path logic still exercises (directories created, files referenced).
    """
    # Tool check returns populated ToolInfo for every required tool.
    def fake_ensure(required, optional=()):
        return {
            name: ToolInfo(name=name, path=f"/usr/local/bin/{name}", version="1.0")
            for name in list(required) + list(optional)
        }

    monkeypatch.setattr(align_cli.tool_check, "ensure_tools_available", fake_ensure)

    def fake_run_cutadapt(*, fastq_in, fastq_out, resolved, log_json, **kwargs):
        Path(fastq_out).parent.mkdir(parents=True, exist_ok=True)
        Path(fastq_out).write_text("fake trimmed")
        Path(log_json).parent.mkdir(parents=True, exist_ok=True)
        Path(log_json).write_text("{}")
        return CutadaptResult(
            input_reads=1000,
            reads_with_adapter=900,
            reads_passing_filters=850,
            log_json_path=Path(log_json),
        )

    monkeypatch.setattr(align_cli.trim_step, "run_cutadapt", fake_run_cutadapt)

    def fake_subtract(*, fastq_in, contam_index, fastq_out_unaligned, **kwargs):
        Path(fastq_out_unaligned).parent.mkdir(parents=True, exist_ok=True)
        Path(fastq_out_unaligned).write_text("fake uncontam")
        return ContamResult(
            total_reads=850, aligned_to_contam=350, unaligned_reads=500
        )

    monkeypatch.setattr(align_cli.contam_step, "subtract_contaminants", fake_subtract)

    def fake_align_mt(*, fastq_in, mt_index, bam_out, **kwargs):
        Path(bam_out).parent.mkdir(parents=True, exist_ok=True)
        Path(bam_out).write_bytes(b"")
        return AlignResult(total_reads=500, aligned=480, bam_path=Path(bam_out))

    monkeypatch.setattr(align_cli.align_step, "align_mt", fake_align_mt)

    def fake_filter_mapq(*, bam_in, bam_out, mapq_threshold, **kwargs):
        Path(bam_out).parent.mkdir(parents=True, exist_ok=True)
        Path(bam_out).write_bytes(b"")
        return 470

    monkeypatch.setattr(align_cli.bam_utils, "filter_bam_mapq", fake_filter_mapq)

    def fake_run_dedup(*, bam_in, bam_out, **kwargs):
        Path(bam_out).parent.mkdir(parents=True, exist_ok=True)
        Path(bam_out).write_bytes(b"")
        return DedupResult(
            strategy="skip",
            input_reads=470,
            output_reads=470,
            bam_path=Path(bam_out),
        )

    monkeypatch.setattr(align_cli.dedup_step, "run_dedup", fake_run_dedup)

    def fake_bam_to_bed6(*, bam_in, bed_out, **kwargs):
        Path(bed_out).parent.mkdir(parents=True, exist_ok=True)
        Path(bed_out).write_text(
            "ND1\t0\t30\tr1\t60\t+\nCOX1\t10\t42\tr2\t55\t-\n"
        )
        return 2

    monkeypatch.setattr(align_cli.bam_utils, "bam_to_bed6", fake_bam_to_bed6)


def _make_fake_index(prefix: Path) -> Path:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"):
        (Path(str(prefix) + suffix)).write_text("idx")
    return prefix


def test_end_to_end_orchestration_produces_expected_outputs(
    tmp_path, mocked_align_step_functions
) -> None:
    fq_dir = tmp_path / "fq"
    fq_dir.mkdir()
    (fq_dir / "sampleA.fq.gz").write_text("@r1\nAAA\n+\nIII\n")
    (fq_dir / "sampleB.fq.gz").write_text("@r1\nAAA\n+\nIII\n")

    contam_idx = _make_fake_index(tmp_path / "idx" / "contam")
    mt_idx = _make_fake_index(tmp_path / "idx" / "mt")
    out_dir = tmp_path / "out"

    exit_code = cli.main(
        [
            "align",
            "--kit-preset",
            "illumina_smallrna",
            "--fastq-dir",
            str(fq_dir),
            "--contam-index",
            str(contam_idx),
            "--mt-index",
            str(mt_idx),
            "--output",
            str(out_dir),
        ]
    )

    assert exit_code == 0

    # Directory layout:
    for subdir in ("trimmed", "contam_filtered", "aligned", "deduped", "bed"):
        assert (out_dir / subdir).is_dir(), f"missing {subdir}/"

    # Per-sample artifacts:
    for sample in ("sampleA", "sampleB"):
        assert (out_dir / "trimmed" / f"{sample}.trimmed.fq.gz").exists()
        assert (out_dir / "trimmed" / f"{sample}.cutadapt.json").exists()
        assert (out_dir / "contam_filtered" / f"{sample}.nocontam.fq.gz").exists()
        assert (out_dir / "aligned" / f"{sample}.bam").exists()
        assert (out_dir / "aligned" / f"{sample}.mapq.bam").exists()
        assert (out_dir / "deduped" / f"{sample}.dedup.bam").exists()
        assert (out_dir / "bed" / f"{sample}.bed").exists()

    # Provenance files:
    counts_path = out_dir / "read_counts.tsv"
    assert counts_path.exists()
    counts_lines = counts_path.read_text().splitlines()
    assert counts_lines[0].split("\t")[0] == "sample"
    # Samples sorted by name.
    assert counts_lines[1].split("\t")[0] == "sampleA"
    assert counts_lines[2].split("\t")[0] == "sampleB"

    settings_path = out_dir / "run_settings.json"
    assert settings_path.exists()
    settings = json.loads(settings_path.read_text())
    assert settings["subcommand"] == "align"
    assert settings["kit_preset_default"] == "illumina_smallrna"
    assert settings["library_strandedness"] == "forward"
    assert settings["dedup_strategy"] == "skip"  # auto -> skip because umi_length==0
    assert settings["mapq_threshold"] == 10
    assert "tool_versions" in settings
    assert settings["tool_versions"]["cutadapt"] == "1.0"
    # Per-sample resolution table is the new spine for kit/dedup decisions.
    by_sample = {row["sample"]: row for row in settings["per_sample"]}
    for sample in ("sampleA", "sampleB"):
        row = by_sample[sample]
        assert row["applied_kit"] == "illumina_smallrna"
        assert row["adapter"] == "TGGAATTCTCGGGTGCCAAGG"
        assert row["umi_length"] == 0
        assert row["dedup_strategy"] == "skip"

    # The per-sample TSV is the human-friendly view of the same data.
    kit_tsv = out_dir / "kit_resolution.tsv"
    assert kit_tsv.exists()
    header = kit_tsv.read_text().splitlines()[0].split("\t")
    assert header[0] == "sample"
    assert "applied_kit" in header
    assert "dedup_strategy" in header


def test_end_to_end_uses_explicit_fastq_list(
    tmp_path, mocked_align_step_functions
) -> None:
    f1 = tmp_path / "s1.fq.gz"
    f1.write_text("")
    contam_idx = _make_fake_index(tmp_path / "c")
    mt_idx = _make_fake_index(tmp_path / "m")
    out_dir = tmp_path / "out2"

    exit_code = cli.main(
        [
            "align",
            "--kit-preset",
            "illumina_smallrna",
            "--fastq",
            str(f1),
            "--contam-index",
            str(contam_idx),
            "--mt-index",
            str(mt_idx),
            "--output",
            str(out_dir),
        ]
    )

    assert exit_code == 0
    assert (out_dir / "bed" / "s1.bed").exists()


def test_end_to_end_picks_up_umi_dedup_when_umi_present(
    tmp_path, mocked_align_step_functions
) -> None:
    fq = tmp_path / "s.fq.gz"
    fq.write_text("")
    out_dir = tmp_path / "out3"

    exit_code = cli.main(
        [
            "align",
            "--kit-preset",
            "illumina_truseq_umi",  # 8 nt 5' UMI
            "--fastq",
            str(fq),
            "--contam-index",
            str(_make_fake_index(tmp_path / "c")),
            "--mt-index",
            str(_make_fake_index(tmp_path / "m")),
            "--output",
            str(out_dir),
        ]
    )

    assert exit_code == 0
    settings = json.loads((out_dir / "run_settings.json").read_text())
    assert settings["dedup_strategy"] == "umi-tools"
    assert settings["per_sample"][0]["umi_length"] == 8
    assert settings["per_sample"][0]["applied_kit"] == "illumina_truseq_umi"
