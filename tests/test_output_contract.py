"""§9 P0 — toy ``mitoribopy all`` run produces every documented output.

We mock the three per-stage runners (align, rpf, rnaseq) so the test
stays hermetic, then assert that the orchestrator wrote every file
listed in the §8 output contract:

* ``run_manifest.json``
* ``warnings.tsv``
* ``outputs_index.tsv``
* ``align/read_counts.tsv``
* ``align/kit_resolution.tsv``
* ``align/run_settings.json``
* ``rpf/rpf_counts.tsv``
* ``rpf/run_settings.json``
* ``rnaseq/te.tsv``        (when rnaseq enabled)
* ``rnaseq/delta_te.tsv``  (when rnaseq enabled)
* ``rnaseq/run_settings.json`` (when rnaseq enabled)

The test only inspects the file system and the manifest's
``outputs`` block; it does NOT validate column schemas (that lives
in tests/test_io_schema_versions.py and tests/test_outputs_index.py).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest


def _seed_align_outputs(run_root: Path) -> None:
    align_dir = run_root / "align"
    align_dir.mkdir(parents=True, exist_ok=True)
    (align_dir / "read_counts.tsv").write_text(
        "# schema_version: 1.0\nsample\ttotal_reads\n"
        "S1\t1000\n"
    )
    (align_dir / "kit_resolution.tsv").write_text(
        "# schema_version: 1.2\nsample\tapplied_kit\n"
        "S1\tillumina_truseq\n"
    )
    (align_dir / "run_settings.json").write_text("{}")
    bed_dir = align_dir / "bed"
    bed_dir.mkdir(exist_ok=True)
    (bed_dir / "S1.mt.bed").write_text("")


def _seed_rpf_outputs(run_root: Path) -> None:
    rpf_dir = run_root / "rpf"
    rpf_dir.mkdir(parents=True, exist_ok=True)
    (rpf_dir / "rpf_counts.tsv").write_text(
        "# schema_version: 1.0\nsample\tgene\tcount\nS1\tMT-ND1\t100\n"
    )
    (rpf_dir / "run_settings.json").write_text(
        json.dumps({"reference_checksum": "abc123"})
    )


def _seed_rnaseq_outputs(run_root: Path) -> None:
    rnaseq_dir = run_root / "rnaseq"
    rnaseq_dir.mkdir(parents=True, exist_ok=True)
    (rnaseq_dir / "te.tsv").write_text(
        "# schema_version: 2.0\n"
        "sample_id\tcondition\tassay\tgene\trpf_count\trna_abundance\tte\tlog2_te\tnote\n"
        "S1\tWT\tribo\tMT-ND1\t100\t1000\t0.1\t-3.32\tpublication_grade\n"
    )
    (rnaseq_dir / "delta_te.tsv").write_text(
        "# schema_version: 2.0\n"
        "gene\tbase_condition\tcompare_condition\tmrna_log2fc\trpf_log2fc\t"
        "delta_te_log2\tpadj_mrna\tpadj_rpf\tpadj_delta_te\tmethod\tnote\n"
        "MT-ND1\tWT\tKO\t0.1\t1.0\t0.9\t0.05\t\t\texternal_de_table\tpublication_grade\n"
    )
    (rnaseq_dir / "run_settings.json").write_text("{}")


@pytest.fixture
def cfg_align_rpf(tmp_path: Path) -> Path:
    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "align:\n  kit_preset: auto\n"
        "rpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
    )
    return cfg


@pytest.fixture
def cfg_align_rpf_rnaseq(tmp_path: Path) -> Path:
    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "align:\n  kit_preset: auto\n"
        "rpf:\n  strain: h\n  fasta: /tmp/tx.fa\n"
        "rnaseq:\n  de_table: /tmp/de.tsv\n"
    )
    return cfg


def _make_fake(name: str, run_root: Path):
    """Return a fake stage runner that seeds the matching outputs."""
    seeders = {
        "align": _seed_align_outputs,
        "rpf": _seed_rpf_outputs,
        "rnaseq": _seed_rnaseq_outputs,
    }

    def fake_run(argv):
        seeders[name](run_root)
        return 0

    return fake_run


def test_align_rpf_run_writes_every_promised_output(
    tmp_path: Path, cfg_align_rpf: Path, monkeypatch
) -> None:
    from mitoribopy import cli
    from mitoribopy.cli import align as align_cli, rpf as rpf_cli

    run_root = tmp_path / "results"
    monkeypatch.setattr(align_cli, "run", _make_fake("align", run_root))
    monkeypatch.setattr(rpf_cli, "run", _make_fake("rpf", run_root))

    rc = cli.main(
        ["all", "--config", str(cfg_align_rpf), "--output", str(run_root)]
    )
    assert rc == 0

    expected = [
        "run_manifest.json",
        "warnings.tsv",
        "outputs_index.tsv",
        "align/read_counts.tsv",
        "align/kit_resolution.tsv",
        "align/run_settings.json",
        "rpf/rpf_counts.tsv",
        "rpf/run_settings.json",
    ]
    missing = [p for p in expected if not (run_root / p).exists()]
    assert not missing, f"missing output(s): {missing}"


def test_rnaseq_outputs_present_when_enabled(
    tmp_path: Path, cfg_align_rpf_rnaseq: Path, monkeypatch
) -> None:
    from mitoribopy import cli
    from mitoribopy.cli import (
        align as align_cli,
        rnaseq as rnaseq_cli,
        rpf as rpf_cli,
    )

    run_root = tmp_path / "results"
    monkeypatch.setattr(align_cli, "run", _make_fake("align", run_root))
    monkeypatch.setattr(rpf_cli, "run", _make_fake("rpf", run_root))
    monkeypatch.setattr(rnaseq_cli, "run", _make_fake("rnaseq", run_root))

    rc = cli.main(
        [
            "all",
            "--config",
            str(cfg_align_rpf_rnaseq),
            "--output",
            str(run_root),
        ]
    )
    assert rc == 0

    for path in (
        "rnaseq/te.tsv",
        "rnaseq/delta_te.tsv",
        "rnaseq/run_settings.json",
    ):
        assert (run_root / path).exists(), f"rnaseq output missing: {path}"


def test_manifest_outputs_block_lists_existing_files(
    tmp_path: Path, cfg_align_rpf: Path, monkeypatch
) -> None:
    from mitoribopy import cli
    from mitoribopy.cli import align as align_cli, rpf as rpf_cli

    run_root = tmp_path / "results"
    monkeypatch.setattr(align_cli, "run", _make_fake("align", run_root))
    monkeypatch.setattr(rpf_cli, "run", _make_fake("rpf", run_root))

    rc = cli.main(
        ["all", "--config", str(cfg_align_rpf), "--output", str(run_root)]
    )
    assert rc == 0

    manifest = json.loads((run_root / "run_manifest.json").read_text())
    outputs = manifest.get("outputs")
    assert isinstance(outputs, list)
    output_types = {row["output_type"] for row in outputs}
    # The toy run produced these stage outputs; the manifest's
    # outputs block should reflect them.
    assert "read_counts" in output_types
    assert "rpf_counts" in output_types
    # And it must NOT advertise files we did not seed.
    assert "te_table" not in output_types

    # Top-level §8 fields are populated.
    for key in ("runtime_seconds", "platform", "python_version"):
        assert key in manifest, f"manifest missing top-level key: {key}"


def test_warnings_tsv_present_even_when_no_warnings(
    tmp_path: Path, cfg_align_rpf: Path, monkeypatch
) -> None:
    from mitoribopy import cli
    from mitoribopy.cli import align as align_cli, rpf as rpf_cli

    run_root = tmp_path / "results"
    monkeypatch.setattr(align_cli, "run", _make_fake("align", run_root))
    monkeypatch.setattr(rpf_cli, "run", _make_fake("rpf", run_root))

    rc = cli.main(
        ["all", "--config", str(cfg_align_rpf), "--output", str(run_root)]
    )
    assert rc == 0

    warnings_tsv = run_root / "warnings.tsv"
    assert warnings_tsv.exists()
    lines = warnings_tsv.read_text().splitlines()
    # At minimum: schema header comment + column header.
    assert any(l.startswith("# schema_version:") for l in lines)
    header = [l for l in lines if not l.startswith("#")][0].split("\t")
    assert header == [
        "stage",
        "sample_id",
        "severity",
        "code",
        "message",
        "suggested_action",
    ]
