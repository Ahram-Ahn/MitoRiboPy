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
    # The legacy --i-understand-mark-duplicates-destroys-mt-ribo-seq-signal
    # flag was removed in v0.4.5 along with the mark-duplicates strategy.
    assert "i-understand-mark-duplicates" not in help_text
    assert "mark-duplicates" not in help_text or "removed in v0.4.5" in help_text


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


def test_mark_duplicates_strategy_is_rejected_by_argparse(capsys, tmp_path) -> None:
    """The picard mark-duplicates dedup strategy was removed in v0.4.5;
    argparse should reject it as an invalid choice."""
    with pytest.raises(SystemExit) as exc:
        cli.main(
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
    assert exc.value.code == 2
    err = capsys.readouterr().err
    assert "invalid choice" in err
    assert "mark-duplicates" in err


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


# ---------- --config flag (flat YAML / JSON / TOML) -------------------------


def test_align_config_flag_loads_yaml_values_into_args(tmp_path, capsys) -> None:
    """A flat YAML passed via --config should populate args for every
    flag whose dest matches a top-level key. Catches the silent no-op
    where --config was registered on the parser but never consumed."""
    fastq_dir = tmp_path / "fq"
    fastq_dir.mkdir()
    cfg = tmp_path / "align.yaml"
    cfg.write_text(
        "kit_preset: illumina_truseq_umi\n"
        "min_length: 99\n"
        "max_length: 50\n"
        "quality: 42\n"
        "mapq: 31\n"
        "seed: 7\n"
    )

    exit_code = cli.main([
        "align",
        "--config", str(cfg),
        "--fastq-dir", str(fastq_dir),
        "--contam-index", str(tmp_path / "ci"),
        "--mt-index", str(tmp_path / "mi"),
        "--dry-run",
    ])
    out = capsys.readouterr().out
    assert exit_code == 0
    # The dry-run plan prints the resolved kit + cutadapt window;
    # both must reflect the YAML values, not the argparse defaults.
    assert "kit=illumina_truseq_umi" in out
    assert "umi_length=8" in out                 # comes from the kit preset
    assert "cutadapt trim per sample (min=99 max=50 nt)" in out
    assert "MAPQ filter at q >= 31" in out


def test_align_cli_flag_wins_over_config_value(tmp_path, capsys) -> None:
    """An explicit CLI flag must override the same key set in --config."""
    fastq_dir = tmp_path / "fq"
    fastq_dir.mkdir()
    cfg = tmp_path / "align.yaml"
    cfg.write_text("min_length: 99\nmapq: 31\n")

    exit_code = cli.main([
        "align",
        "--config", str(cfg),
        "--min-length", "25",                    # overrides 99
        "--fastq-dir", str(fastq_dir),
        "--contam-index", str(tmp_path / "ci"),
        "--mt-index", str(tmp_path / "mi"),
        "--dry-run",
    ])
    out = capsys.readouterr().out
    assert exit_code == 0
    assert "(min=25 max=45 nt)" in out           # CLI win, max stays at default
    assert "MAPQ filter at q >= 31" in out       # mapq still from config


def test_align_config_warns_on_unknown_key(tmp_path, caplog) -> None:
    """Unknown YAML keys are logged but do not abort the run."""
    import logging

    fastq_dir = tmp_path / "fq"
    fastq_dir.mkdir()
    cfg = tmp_path / "align.yaml"
    cfg.write_text("kit_preset: illumina_smallrna\nzz_typo_key: 1\n")

    # The package logger sets propagate=False so caplog (which hooks the
    # root logger) would not see the record. Re-enable propagation for
    # the duration of the test, restore it afterwards.
    pkg_logger = logging.getLogger("mitoribopy")
    saved_propagate = pkg_logger.propagate
    pkg_logger.propagate = True
    try:
        with caplog.at_level(logging.WARNING, logger="mitoribopy"):
            exit_code = cli.main([
                "align",
                "--config", str(cfg),
                "--fastq-dir", str(fastq_dir),
                "--contam-index", str(tmp_path / "ci"),
                "--mt-index", str(tmp_path / "mi"),
                "--dry-run",
            ])
    finally:
        pkg_logger.propagate = saved_propagate

    assert exit_code == 0
    assert any(
        "Ignoring unknown --config keys" in rec.message
        and "zz_typo_key" in rec.message
        for rec in caplog.records
    )


def test_align_config_missing_file_errors_with_clear_message(
    tmp_path, capsys
) -> None:
    exit_code = cli.main([
        "align",
        "--config", str(tmp_path / "does_not_exist.yaml"),
        "--fastq-dir", str(tmp_path),
        "--contam-index", "x",
        "--mt-index", "y",
        "--dry-run",
    ])
    err = capsys.readouterr().err
    assert exit_code == 2
    assert "ERROR" in err
    assert "Config file not found" in err


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

    # The CLI bypasses dedup_step.skip_dedup for the no-UMI ('skip')
    # path and instead calls dedup_step.skip_dedup_in_place. Stub the
    # whole helper so the placeholder mapq.bam (empty bytes) does not
    # blow up pysam.
    def fake_skip_in_place(bam_path, **kwargs):
        return DedupResult(
            strategy="skip",
            input_reads=470,
            output_reads=470,
            bam_path=Path(bam_path),
        )

    monkeypatch.setattr(
        align_cli.dedup_step, "skip_dedup_in_place", fake_skip_in_place
    )

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

    # Directory layout — only the dirs that hold a kept artefact remain
    # populated. trimmed/ and contam_filtered/ exist because cutadapt /
    # contam-subtract created them, but the FASTQ files inside them
    # have been deleted as soon as the next step consumed them.
    for subdir in ("trimmed", "contam_filtered", "aligned", "bed"):
        assert (out_dir / subdir).is_dir(), f"missing {subdir}/"

    # Per-sample artifacts. Intermediate FASTQs and the pre-MAPQ BAM
    # are removed by the orchestrator unless --keep-intermediates is
    # passed; only the cutadapt JSON log, the post-MAPQ BAM, and the
    # final BED survive. There is no separate dedup.bam for the no-UMI
    # path: the orchestrator wires mapq.bam straight into bam_to_bed6.
    for sample in ("sampleA", "sampleB"):
        assert not (out_dir / "trimmed" / f"{sample}.trimmed.fq.gz").exists()
        assert (out_dir / "trimmed" / f"{sample}.cutadapt.json").exists()
        assert not (out_dir / "contam_filtered" / f"{sample}.nocontam.fq.gz").exists()
        assert not (out_dir / "aligned" / f"{sample}.bam").exists()
        assert (out_dir / "aligned" / f"{sample}.mapq.bam").exists()
        assert not (out_dir / "deduped" / f"{sample}.dedup.bam").exists()
        assert (out_dir / "bed" / f"{sample}.bed").exists()

    # Provenance files:
    counts_path = out_dir / "read_counts.tsv"
    assert counts_path.exists()
    raw_counts = counts_path.read_text().splitlines()
    assert raw_counts[0].startswith("# schema_version:")  # P1.12
    counts_lines = [line for line in raw_counts if not line.startswith("#")]
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
    raw = kit_tsv.read_text().splitlines()
    assert raw[0].startswith("# schema_version:")  # P1.12
    header = [line for line in raw if not line.startswith("#")][0].split("\t")
    assert header[0] == "sample"
    assert "applied_kit" in header
    assert "dedup_strategy" in header
    assert "umi_source" in header  # P1.11


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


# ---------- --max-parallel-samples concurrency ----------------------------


def test_max_parallel_samples_default_is_auto() -> None:
    """The publication-readiness refactor (ref-5) made 'auto' the
    default: the resolved plan is computed at runtime, so the parser
    sees ``None`` and the resource_plan layer treats that as 'auto'."""
    parser = build_parser()
    ns = parser.parse_args([])
    assert ns.max_parallel_samples is None
    assert ns.single_sample_mode is False


def test_max_parallel_samples_parses_value() -> None:
    parser = build_parser()
    ns = parser.parse_args(["--max-parallel-samples", "4"])
    assert ns.max_parallel_samples == "4"


def test_single_sample_mode_flag_parses() -> None:
    parser = build_parser()
    ns = parser.parse_args(["--single-sample-mode"])
    assert ns.single_sample_mode is True


def test_max_parallel_samples_accepts_auto_string() -> None:
    parser = build_parser()
    ns = parser.parse_args(["--max-parallel-samples", "auto"])
    assert ns.max_parallel_samples == "auto"


def test_max_parallel_samples_help_mentions_threads_division() -> None:
    parser = build_parser()
    help_text = parser.format_help()
    assert "--max-parallel-samples" in help_text
    # The help should explain that --threads is divided across workers,
    # so a user setting both flags doesn't double-allocate CPU.
    assert "T // N" in help_text or "threads // N" in help_text


@pytest.mark.parametrize(
    "threads, max_parallel, expected",
    [
        (8, 4, 2),
        (8, 1, 8),
        (1, 4, 1),  # threads < max_parallel -> floor at 1
        (7, 4, 1),  # integer-division floor
        (None, 4, 1),  # threads=None falls back to 1
        (None, 1, 1),
        (16, 4, 4),
        (3, 2, 1),
    ],
)
def test_per_worker_threads_divides_budget_across_workers(
    threads, max_parallel, expected
) -> None:
    from mitoribopy.cli.align import _per_worker_threads

    assert _per_worker_threads(threads, max_parallel) == expected


def test_align_parallel_completes_all_samples(
    tmp_path, mocked_align_step_functions
) -> None:
    """End-to-end run with 4 mocked samples and --max-parallel-samples 2.

    All per-sample BED files and the aggregated read_counts.tsv must
    match the serial reference run; only execution order changes.
    """
    fq_dir = tmp_path / "fq"
    fq_dir.mkdir()
    sample_names = ["sA", "sB", "sC", "sD"]
    for name in sample_names:
        (fq_dir / f"{name}.fq.gz").write_text("@r1\nAAA\n+\nIII\n")

    contam_idx = _make_fake_index(tmp_path / "idx" / "c")
    mt_idx = _make_fake_index(tmp_path / "idx" / "m")
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
            "--threads",
            "8",
            "--max-parallel-samples",
            "2",
        ]
    )

    assert exit_code == 0
    for name in sample_names:
        assert (out_dir / "bed" / f"{name}.bed").exists(), f"missing bed/{name}.bed"
        assert (out_dir / "aligned" / f"{name}.mapq.bam").exists()

    raw_counts = (out_dir / "read_counts.tsv").read_text().splitlines()
    counts_lines = [line for line in raw_counts if not line.startswith("#")]
    # Header + one row per sample, sorted by sample name.
    assert counts_lines[0].split("\t")[0] == "sample"
    sample_col = [line.split("\t")[0] for line in counts_lines[1:]]
    assert sample_col == sorted(sample_names)


def test_align_parallel_passes_divided_threads_to_tools(
    tmp_path, monkeypatch, mocked_align_step_functions
) -> None:
    """With --threads 8 --max-parallel-samples 4, each tool sees threads=2."""
    seen_threads: list[int] = []

    real_cutadapt = align_cli.trim_step.run_cutadapt

    def capturing_cutadapt(*, threads, **kwargs):
        seen_threads.append(threads)
        return real_cutadapt(threads=threads, **kwargs)

    monkeypatch.setattr(
        align_cli.trim_step, "run_cutadapt", capturing_cutadapt
    )

    fq_dir = tmp_path / "fq"
    fq_dir.mkdir()
    for name in ("a", "b", "c", "d"):
        (fq_dir / f"{name}.fq.gz").write_text("")

    out_dir = tmp_path / "out"
    exit_code = cli.main(
        [
            "align",
            "--kit-preset",
            "illumina_smallrna",
            "--fastq-dir",
            str(fq_dir),
            "--contam-index",
            str(_make_fake_index(tmp_path / "c")),
            "--mt-index",
            str(_make_fake_index(tmp_path / "m")),
            "--output",
            str(out_dir),
            "--threads",
            "8",
            "--max-parallel-samples",
            "4",
        ]
    )

    assert exit_code == 0
    assert len(seen_threads) == 4
    # 8 // 4 = 2 per worker.
    assert all(t == 2 for t in seen_threads), seen_threads


def test_align_parallel_failure_propagates(
    tmp_path, monkeypatch, mocked_align_step_functions
) -> None:
    """If one sample raises, the run exits non-zero (fail-fast)."""
    real_align_mt = align_cli.align_step.align_mt

    def selective_align(*, fastq_in, **kwargs):
        # Sample 'b' blows up in alignment; others succeed.
        if "b.nocontam" in str(fastq_in) or str(fastq_in).endswith("b.fq.gz"):
            raise RuntimeError("synthetic alignment failure for sample b")
        return real_align_mt(fastq_in=fastq_in, **kwargs)

    monkeypatch.setattr(align_cli.align_step, "align_mt", selective_align)

    fq_dir = tmp_path / "fq"
    fq_dir.mkdir()
    for name in ("a", "b", "c", "d"):
        (fq_dir / f"{name}.fq.gz").write_text("")

    out_dir = tmp_path / "out_fail"
    with pytest.raises(RuntimeError, match="synthetic alignment failure"):
        cli.main(
            [
                "align",
                "--kit-preset",
                "illumina_smallrna",
                "--fastq-dir",
                str(fq_dir),
                "--contam-index",
                str(_make_fake_index(tmp_path / "c")),
                "--mt-index",
                str(_make_fake_index(tmp_path / "m")),
                "--output",
                str(out_dir),
                "--max-parallel-samples",
                "2",
            ]
        )

    # The aggregated read_counts.tsv must NOT be written when the run
    # fails -- callers rely on its presence to gate downstream stages.
    assert not (out_dir / "read_counts.tsv").exists()


# ---------- per-stage timing + summary --------------------------------------


def test_align_logs_per_stage_timing_and_emits_summary(
    tmp_path, mocked_align_step_functions
) -> None:
    """Each stage line carries a duration; a summary table follows the run."""
    fq_dir = tmp_path / "fq"
    fq_dir.mkdir()
    for name in ("sA", "sB"):
        (fq_dir / f"{name}.fq.gz").write_text("@r1\nAAA\n+\nIII\n")

    out_dir = tmp_path / "out"
    exit_code = cli.main(
        [
            "align",
            "--kit-preset",
            "illumina_smallrna",
            "--fastq-dir",
            str(fq_dir),
            "--contam-index",
            str(_make_fake_index(tmp_path / "c")),
            "--mt-index",
            str(_make_fake_index(tmp_path / "m")),
            "--output",
            str(out_dir),
        ]
    )
    assert exit_code == 0

    log_text = (out_dir / "mitoribopy.log").read_text()

    # Compact per-stage lines: "[ALIGN] sample: stage <pad> <duration> — detail"
    # The em-dash separator is the marker for the new format.
    for sample in ("sA", "sB"):
        for stage in ("trim", "contam-filter", "mt-align", "mapq-filter",
                       "dedup", "bam-to-bed"):
            assert f"{sample}: {stage}" in log_text, f"missing stage line {sample}: {stage}"
        assert f"{sample}: ✓ done in" in log_text, f"missing per-sample done line for {sample}"

    # End-of-run summary table lines.
    assert "Timing summary (2 sample(s)" in log_text
    assert "stage" in log_text and "total" in log_text and "max" in log_text
    assert "wall:" in log_text


def test_align_summary_skipped_when_all_samples_resumed(
    tmp_path, mocked_align_step_functions
) -> None:
    """If --resume reloads every sample from cache, no stages run -> no summary."""
    from mitoribopy.align._types import SampleCounts

    fq_dir = tmp_path / "fq"
    fq_dir.mkdir()
    (fq_dir / "sA.fq.gz").write_text("")

    out_dir = tmp_path / "out_resumed"
    sample_done_dir = out_dir / ".sample_done"
    sample_done_dir.mkdir(parents=True)
    # Pre-seed a completion marker that matches the SampleCounts schema
    # so resume can reload it.
    cached = SampleCounts(
        sample="sA",
        total_reads=1000,
        post_trim=850,
        rrna_aligned=350,
        post_rrna_filter=500,
        mt_aligned=480,
        unaligned_to_mt=20,
        mt_aligned_after_mapq=470,
        mt_aligned_after_dedup=470,
    )
    (sample_done_dir / "sA.json").write_text(
        json.dumps({k: getattr(cached, k) for k in SampleCounts.__dataclass_fields__})
    )

    exit_code = cli.main(
        [
            "align",
            "--kit-preset",
            "illumina_smallrna",
            "--fastq-dir",
            str(fq_dir),
            "--contam-index",
            str(_make_fake_index(tmp_path / "c")),
            "--mt-index",
            str(_make_fake_index(tmp_path / "m")),
            "--output",
            str(out_dir),
            "--resume",
        ]
    )
    assert exit_code == 0

    log_text = (out_dir / "mitoribopy.log").read_text()
    # No stage timings recorded -> no summary block emitted.
    assert "Timing summary" not in log_text
    # But the resume notice is still there.
    assert "resumed from" in log_text
