"""Tests for the summary subsystem (P1.6, P1.7, P1.8).

* `build_summary_qc` aggregates per-sample metrics from align +
  rpf outputs.
* `write_summary_qc` writes the schema-versioned TSV.
* `render_summary_md` produces a Markdown report.
* `mitoribopy summarize <run-dir>` regenerates both.
* `mitoribopy all` auto-invokes summarize after the manifest write.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from mitoribopy import cli
from mitoribopy.summary.qc import (
    QC_THRESHOLDS,
    SUMMARY_QC_COLUMNS,
    build_summary_qc,
    write_summary_qc,
)
from mitoribopy.summary.render import render_summary_md


# ---------- build_summary_qc -----------------------------------------------


def _seed_align_outputs(run_root: Path, *, samples: list[tuple[str, str]]) -> None:
    """Write a tiny align/read_counts.tsv and align/kit_resolution.tsv."""
    align = run_root / "align"
    align.mkdir(parents=True, exist_ok=True)
    rc_lines = [
        "# schema_version: 1.0",
        "sample\ttotal_reads\tpost_trim\trrna_aligned\tpost_rrna_filter\t"
        "mt_aligned\tunaligned_to_mt\tmt_aligned_after_mapq\tmt_aligned_after_dedup",
    ]
    # Numbers chosen to exceed every QC threshold so qc_status='pass'.
    # post_trim=200_000, mt_mrna_fraction=400/800=0.5 (>0.05),
    # dedup_removed_fraction=20/380≈0.05 (<0.95).
    for sample_id, _ in samples:
        rc_lines.append(
            f"{sample_id}\t250000\t200000\t100\t800\t400\t400\t380\t360"
        )
    (align / "read_counts.tsv").write_text("\n".join(rc_lines) + "\n")

    kit_lines = [
        "# schema_version: 1.1",
        "sample\tfastq\tapplied_kit\tadapter\tumi_length\tumi_position\t"
        "dedup_strategy\tdetected_kit\tdetection_match_rate\t"
        "detection_ambiguous\tsource\tumi_source",
    ]
    for sample_id, _ in samples:
        kit_lines.append(
            f"{sample_id}\t{sample_id}.fq.gz\tilllumina_truseq_umi\t"
            f"AGATCGGAAGAGC\t8\t5p\tumi-tools\tilllumina_truseq_umi\t"
            f"0.92\tfalse\tdetected\tdeclared"
        )
    (align / "kit_resolution.tsv").write_text("\n".join(kit_lines) + "\n")


def _seed_sample_sheet(run_root: Path, samples: list[tuple[str, str]]) -> Path:
    sheet = run_root / "samples.tsv"
    lines = ["sample_id\tassay\tcondition\tfastq_1"]
    for sample_id, condition in samples:
        lines.append(f"{sample_id}\tribo\t{condition}\t{sample_id}.fq.gz")
    sheet.write_text("\n".join(lines) + "\n")
    return sheet


def test_build_summary_qc_aggregates_align_metrics(tmp_path: Path) -> None:
    samples = [("WT_Ribo_1", "WT"), ("KO_Ribo_1", "KO")]
    sheet = _seed_sample_sheet(tmp_path, samples)
    _seed_align_outputs(tmp_path, samples=samples)
    manifest = {"sample_sheet": str(sheet)}

    rows = build_summary_qc(tmp_path, manifest)
    assert len(rows) == 2
    by_id = {r["sample_id"]: r for r in rows}

    wt = by_id["WT_Ribo_1"]
    assert wt["assay"] == "ribo"
    assert wt["condition"] == "WT"
    assert wt["kit_applied"] == "illlumina_truseq_umi"
    assert wt["adapter_match_rate"] == "0.92"
    assert wt["umi_source"] == "declared"
    assert wt["post_trim_reads"] == "200000"
    # contam = (post_trim - post_rrna) / post_trim = (200000-800)/200000
    assert abs(float(wt["contam_fraction"]) - (200000 - 800) / 200000) < 1e-3
    # mt_mrna = mt_aligned / post_rrna = 400 / 800 = 0.5
    assert float(wt["mt_mrna_fraction"]) == 0.5
    # mapq_kept = 380/400 = 0.95
    assert abs(float(wt["mapq_kept_fraction"]) - 0.95) < 1e-3
    # dedup_removed = 1 - 360/380 ≈ 0.0526
    assert abs(float(wt["dedup_removed_fraction"]) - (1 - 360 / 380)) < 1e-3
    assert wt["qc_status"] == "pass"
    assert wt["qc_notes"] == ""


def test_build_summary_qc_flags_low_post_trim_reads(tmp_path: Path) -> None:
    samples = [("S1", "WT")]
    sheet = _seed_sample_sheet(tmp_path, samples)
    align = tmp_path / "align"
    align.mkdir()
    (align / "read_counts.tsv").write_text(
        "# schema_version: 1.0\n"
        "sample\ttotal_reads\tpost_trim\trrna_aligned\tpost_rrna_filter\t"
        "mt_aligned\tunaligned_to_mt\tmt_aligned_after_mapq\tmt_aligned_after_dedup\n"
        # post_trim < 100,000 threshold
        "S1\t1000\t500\t100\t400\t100\t300\t90\t85\n"
    )
    rows = build_summary_qc(tmp_path, {"sample_sheet": str(sheet)})
    assert rows[0]["qc_status"] == "warn"
    assert "post_trim_reads" in rows[0]["qc_notes"]


def test_build_summary_qc_falls_back_when_no_sample_sheet(tmp_path: Path) -> None:
    """Without a sheet, use align/read_counts.tsv as the identity source."""
    samples = [("S1", "")]
    _seed_align_outputs(tmp_path, samples=samples)
    rows = build_summary_qc(tmp_path, {})
    assert rows[0]["sample_id"] == "S1"
    assert rows[0]["assay"] == "ribo"


def test_build_summary_qc_handles_missing_align_dir(tmp_path: Path) -> None:
    """No align outputs at all — every align metric column is empty."""
    sheet = _seed_sample_sheet(tmp_path, [("S1", "WT")])
    rows = build_summary_qc(tmp_path, {"sample_sheet": str(sheet)})
    assert rows[0]["sample_id"] == "S1"
    assert rows[0]["post_trim_reads"] == ""
    assert rows[0]["kit_applied"] == ""


# ---------- write_summary_qc -----------------------------------------------


def test_write_summary_qc_emits_schema_header_and_columns(tmp_path: Path) -> None:
    rows = [
        {col: ("S1" if col == "sample_id" else "ribo" if col == "assay" else "")
         for col in SUMMARY_QC_COLUMNS}
    ]
    out = write_summary_qc(rows, tmp_path / "summary_qc.tsv")
    raw = out.read_text().splitlines()
    assert raw[0].startswith("# schema_version:")
    header = [line for line in raw if not line.startswith("#")][0].split("\t")
    assert header == list(SUMMARY_QC_COLUMNS)


# ---------- render_summary_md -----------------------------------------------


def test_render_summary_md_includes_key_sections(tmp_path: Path) -> None:
    manifest = {
        "mitoribopy_version": "0.5.1",
        "schema_version": "1.1.0",
        "git_commit": "abc1234",
        "command": "mitoribopy all --config c.yaml --output results/",
        "config_source": "c.yaml",
        "stages": {
            "align": {"status": "completed", "runtime_seconds": 12.3},
            "rpf": {"status": "completed", "runtime_seconds": 4.5},
            "rnaseq": {"status": "skipped", "reason": "no rnaseq section"},
        },
        "output_schemas": {"read_counts.tsv": "1.0", "te.tsv": "1.0"},
    }
    rows = [
        {
            "sample_id": "S1", "assay": "ribo", "condition": "WT",
            "kit_applied": "illumina_smallrna",
            "post_trim_reads": "900",
            "mt_mrna_fraction": "0.5",
            "qc_status": "pass", "qc_notes": "",
        }
    ]
    md = render_summary_md(
        run_root=tmp_path / "results",
        manifest=manifest,
        summary_qc_rows=rows,
        warnings=[],
    )
    assert md.startswith("# MitoRiboPy run summary")
    assert "## Stages" in md
    assert "**align**" in md and "completed" in md
    assert "**rnaseq**" in md and "skipped" in md
    assert "## Samples" in md
    assert "## Per-sample QC (subset)" in md
    assert "S1" in md and "illumina_smallrna" in md
    assert "## Outputs" in md
    assert "## Output schema versions" in md
    assert "no structured warnings" in md


def test_render_summary_md_includes_periodicity_confidence_section(tmp_path: Path) -> None:
    """v0.9.0+ SUMMARY.md surfaces the headline Fourier CI + permutation p."""
    run_root = tmp_path / "results"
    qc_dir = run_root / "rpf" / "qc"
    qc_dir.mkdir(parents=True)
    score_path = qc_dir / "fourier_period3_score_combined.tsv"
    score_path.write_text(
        "sample\tread_length\tgene_set\tregion\tn_genes\tn_sites_total\tn_nt"
        "\tspectral_ratio_3nt\tspectral_ratio_3nt_ci_low"
        "\tspectral_ratio_3nt_ci_high\tpermutation_p\tsnr_call\n"
        "WT_R1\t32\tcombined\torf_start\t9\t12345\t99"
        "\t8.42\t6.10\t10.55\t0.0005\thealthy\n"
        "KO_R1\t32\tcombined\torf_start\t9\t11500\t99"
        "\t1.42\t1.05\t1.85\t0.412\tbroken\n"
    )
    manifest = {
        "mitoribopy_version": "0.9.0",
        "schema_version": "1.3.0",
        "command": "mitoribopy all",
        "config_source": "c.yaml",
        "stages": {},
    }
    md = render_summary_md(
        run_root=run_root,
        manifest=manifest,
        summary_qc_rows=[],
        warnings=[],
    )
    assert "## Periodicity statistical confidence" in md
    assert "WT_R1" in md and "8.42x" in md
    assert "[6.10, 10.55]" in md
    assert "<0.001" in md  # The Laplace-smoothed perm p formatter.
    assert "KO_R1" in md and "broken" in md


def test_render_summary_md_skips_periodicity_section_when_table_missing(tmp_path: Path) -> None:
    """No score TSV → no spurious empty section."""
    run_root = tmp_path / "results"
    md = render_summary_md(
        run_root=run_root,
        manifest={
            "mitoribopy_version": "0.9.0",
            "schema_version": "1.3.0",
            "command": "mitoribopy all",
            "stages": {},
        },
        summary_qc_rows=[],
        warnings=[],
    )
    assert "## Periodicity statistical confidence" not in md


def test_render_summary_md_lists_warnings(tmp_path: Path) -> None:
    manifest = {
        "mitoribopy_version": "0.5.1",
        "schema_version": "1.1.0",
        "command": "mitoribopy all",
        "config_source": "c.yaml",
        "stages": {},
    }
    warnings = [
        {
            "component": "RNASEQ", "code": "UMI_INFERRED_NO_DECLARATION",
            "sample": "S1", "message": "umi inferred",
        }
    ]
    md = render_summary_md(
        run_root=tmp_path / "results",
        manifest=manifest,
        summary_qc_rows=[],
        warnings=warnings,
    )
    assert "UMI_INFERRED_NO_DECLARATION" in md
    assert "umi inferred" in md


# ---------- mitoribopy summarize CLI ----------------------------------------


def _make_finished_run(tmp_path: Path) -> Path:
    """Build a minimal "finished run" directory with a manifest."""
    run_root = tmp_path / "results"
    run_root.mkdir(parents=True, exist_ok=True)
    samples = [("WT_Ribo_1", "WT"), ("KO_Ribo_1", "KO")]
    sheet = _seed_sample_sheet(run_root, samples)
    _seed_align_outputs(run_root, samples=samples)
    (run_root / "rpf").mkdir()
    (run_root / "rpf" / "rpf_counts.tsv").write_text(
        "# schema_version: 1.0\nsample\tgene\tcount\n"
    )
    manifest = {
        "schema_version": "1.1.0",
        "mitoribopy_version": "0.5.1",
        "command": "mitoribopy all",
        "config_source": "c.yaml",
        "sample_sheet": str(sheet),
        "stages": {
            "align": {"status": "completed", "runtime_seconds": 1.0},
            "rpf": {"status": "completed", "runtime_seconds": 1.0},
        },
        "output_schemas": {"read_counts.tsv": "1.0"},
        "warnings": [],
    }
    (run_root / "run_manifest.json").write_text(json.dumps(manifest))
    return run_root


def test_summarize_cli_emits_summary_md_and_summary_qc_tsv(
    tmp_path: Path, capsys
) -> None:
    run_root = _make_finished_run(tmp_path)
    rc = cli.main(["summarize", str(run_root)])
    assert rc == 0
    summary_md = run_root / "SUMMARY.md"
    summary_qc = run_root / "summary_qc.tsv"
    assert summary_md.is_file()
    assert summary_qc.is_file()
    md = summary_md.read_text()
    assert "MitoRiboPy run summary" in md
    assert "WT_Ribo_1" in md or "Per-sample QC" in md
    qc = summary_qc.read_text()
    assert "# schema_version:" in qc
    assert "WT_Ribo_1" in qc


def test_summarize_cli_missing_directory_exits_2(tmp_path: Path, capsys) -> None:
    rc = cli.main(["summarize", str(tmp_path / "nope")])
    assert rc == 2
    assert "ERROR" in capsys.readouterr().err


def test_summarize_cli_missing_manifest_exits_2(tmp_path: Path, capsys) -> None:
    rc = cli.main(["summarize", str(tmp_path)])
    assert rc == 2
    assert "manifest not found" in capsys.readouterr().err


def test_summarize_cli_listed_in_top_level_help(capsys) -> None:
    rc = cli.main(["--help"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "summarize" in out


# ---------- mitoribopy all auto-emits SUMMARY ------------------------------


def test_mitoribopy_all_auto_emits_summary_md(tmp_path: Path, monkeypatch) -> None:
    """A real `mitoribopy all` invocation should produce SUMMARY.md and
    summary_qc.tsv at the run root without an extra `summarize` call."""
    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "rpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
    )

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
        "all", "--config", str(cfg), "--output", str(tmp_path / "results"),
    ])
    assert rc == 0
    assert (tmp_path / "results" / "SUMMARY.md").is_file()
    assert (tmp_path / "results" / "summary_qc.tsv").is_file()
