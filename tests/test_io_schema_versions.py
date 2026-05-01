"""Tests for the output-schema-version registry (P1.12).

Every advertised TSV output prepends a ``# schema_version: X.Y`` line
so consumers (pandas with ``comment='#'``, jq pipelines, this package's
own loaders) can detect drift without parsing the manifest.

The manifest itself records the same values under ``output_schemas`` so
a downstream script can fail fast on a major bump.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from mitoribopy.io.schema_versions import (
    OUTPUT_SCHEMA_VERSIONS,
    schema_header_line,
)


def test_registered_files_have_semver_versions() -> None:
    for name, version in OUTPUT_SCHEMA_VERSIONS.items():
        assert isinstance(name, str) and name.endswith(".tsv")
        # Loose semver: major.minor (no pre-release, no metadata).
        major, minor = version.split(".")
        assert major.isdigit() and minor.isdigit(), name


def test_schema_header_line_format() -> None:
    line = schema_header_line("read_counts.tsv")
    assert line.startswith("# schema_version: ")
    assert line.endswith("\n")


def test_schema_header_line_unknown_name_raises() -> None:
    with pytest.raises(KeyError):
        schema_header_line("not_a_real_output.tsv")


def test_read_counts_tsv_carries_schema_header(tmp_path: Path) -> None:
    """Real writer round-trip: write_read_counts_table prepends the header."""
    from mitoribopy.align._types import SampleCounts
    from mitoribopy.align.read_counts import write_read_counts_table

    out = tmp_path / "read_counts.tsv"
    write_read_counts_table(
        [
            SampleCounts(
                sample="A",
                total_reads=100,
                post_trim=90,
                rrna_aligned=10,
                post_rrna_filter=80,
                mt_aligned=70,
                unaligned_to_mt=10,
                mt_aligned_after_mapq=65,
                mt_aligned_after_dedup=60,
            )
        ],
        out,
    )
    text = out.read_text()
    assert text.startswith(f"# schema_version: {OUTPUT_SCHEMA_VERSIONS['read_counts.tsv']}\n")


def test_read_counts_tsv_parses_with_pandas_comment_arg(tmp_path: Path) -> None:
    """A consumer that opens the file with ``comment='#'`` (pandas /
    csv.DictReader with manual filter) must see the column header on
    the first non-comment line and the data rows after."""
    pd = pytest.importorskip("pandas")
    from mitoribopy.align._types import SampleCounts
    from mitoribopy.align.read_counts import write_read_counts_table

    out = tmp_path / "read_counts.tsv"
    write_read_counts_table(
        [
            SampleCounts(
                sample="A", total_reads=100, post_trim=90, rrna_aligned=10,
                post_rrna_filter=80, mt_aligned=70, unaligned_to_mt=10,
                mt_aligned_after_mapq=65, mt_aligned_after_dedup=60,
            )
        ],
        out,
    )
    df = pd.read_csv(out, sep="\t", comment="#")
    assert list(df.columns)[:3] == ["sample", "total_reads", "post_trim"]
    assert df["sample"].iloc[0] == "A"


def test_kit_resolution_tsv_includes_umi_source_column(tmp_path: Path) -> None:
    """P1.11 added umi_source as the trailing column at v1.1; the
    schema-version header must reflect the bump."""
    assert OUTPUT_SCHEMA_VERSIONS["kit_resolution.tsv"] == "1.1"


def test_kit_resolution_tsv_writer_emits_schema_header(tmp_path: Path) -> None:
    """End-to-end check: write_kit_resolution_tsv prepends the header
    line and includes umi_source as the trailing column."""
    from mitoribopy.align._types import KIT_PRESETS, ResolvedKit
    from mitoribopy.align.sample_resolve import (
        SampleResolution,
        write_kit_resolution_tsv,
    )

    fastq = tmp_path / "x.fq.gz"
    fastq.write_bytes(b"")
    res = SampleResolution(
        sample="A",
        fastq=fastq,
        kit=ResolvedKit(
            kit="illumina_truseq_umi",
            adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            umi_length=8,
            umi_position="5p",
        ),
        dedup_strategy="umi-tools",
        detected_kit="illumina_truseq_umi",
        detection_match_rate=0.92,
        detection_ambiguous=False,
        source="detected",
        umi_source="declared",
    )
    out = tmp_path / "kit_resolution.tsv"
    write_kit_resolution_tsv([res], out)
    raw = out.read_text().splitlines()
    assert raw[0].startswith(
        f"# schema_version: {OUTPUT_SCHEMA_VERSIONS['kit_resolution.tsv']}"
    )
    header = raw[1].split("\t")
    assert header[-1] == "umi_source"
    body = raw[2].split("\t")
    assert body[-1] == "declared"


def test_load_ribo_counts_skips_schema_header_line(tmp_path: Path) -> None:
    """`load_ribo_counts` must transparently skip the comment line."""
    from mitoribopy.rnaseq.counts import load_ribo_counts

    p = tmp_path / "rpf_counts.tsv"
    p.write_text(
        "# schema_version: 1.0\n"
        "sample\tgene\tcount\n"
        "A1\tMT-ND1\t100\n"
    )
    counts = load_ribo_counts(p)
    assert counts == {"MT-ND1": {"A1": 100}}


def test_te_and_delta_te_writers_emit_schema_header(tmp_path: Path) -> None:
    from mitoribopy.cli.rnaseq import _write_delta_te_table, _write_te_table
    from mitoribopy.rnaseq._types import DTeRow, TeRow

    te_path = tmp_path / "te.tsv"
    _write_te_table(
        [TeRow(sample="A", gene="MT-ND1", rpf_count=10, mrna_abundance=5.0, te=2.0)],
        te_path,
    )
    text = te_path.read_text()
    assert text.startswith(
        f"# schema_version: {OUTPUT_SCHEMA_VERSIONS['te.tsv']}\n"
    )

    dte_path = tmp_path / "delta_te.tsv"
    _write_delta_te_table(
        [
            DTeRow(
                gene="MT-ND1",
                mrna_log2fc=0.1,
                rpf_log2fc=0.5,
                delta_te_log2=0.4,
                padj=0.01,
                note="ok",
            )
        ],
        dte_path,
    )
    text = dte_path.read_text()
    assert text.startswith(
        f"# schema_version: {OUTPUT_SCHEMA_VERSIONS['delta_te.tsv']}\n"
    )


def test_manifest_records_output_schemas(tmp_path: Path, monkeypatch) -> None:
    """A real `mitoribopy all` run must include the output_schemas map."""
    from mitoribopy import cli
    from mitoribopy.cli import all_ as all_cli

    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "rpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
    )

    def fake_rpf(argv):
        out = Path(tmp_path / "results" / "rpf")
        out.mkdir(parents=True, exist_ok=True)
        (out / "rpf_counts.tsv").write_text("# schema_version: 1.0\nsample\tgene\tcount\n")
        (out / "run_settings.json").write_text("{}")
        return 0

    from mitoribopy.cli import rpf as rpf_cli
    monkeypatch.setattr(rpf_cli, "run", fake_rpf)

    rc = cli.main([
        "all", "--config", str(cfg), "--output", str(tmp_path / "results")
    ])
    assert rc == 0
    manifest = json.loads(
        (tmp_path / "results" / "run_manifest.json").read_text()
    )
    schemas = manifest.get("output_schemas")
    assert isinstance(schemas, dict)
    # Spot-check a few entries.
    assert schemas["read_counts.tsv"] == OUTPUT_SCHEMA_VERSIONS["read_counts.tsv"]
    assert schemas["te.tsv"] == OUTPUT_SCHEMA_VERSIONS["te.tsv"]
    assert schemas["kit_resolution.tsv"] == "1.1"
    # Manifest layout itself bumped to 1.1 to record the new field.
    assert manifest["schema_version"].startswith("1.")
