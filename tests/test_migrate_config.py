"""Tests for the legacy-key migrator (P1.10).

Coverage:

* `mitoribopy.config.migrate.migrate` — pure-data rewrites for every
  legacy key the codebase historically accepted.
* `mitoribopy migrate-config <path>` — CLI: stdout (canonical YAML),
  stderr (per-rewrite log lines), exit codes.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from mitoribopy import cli
from mitoribopy.config.migrate import (
    LEGACY_ALIGN_TOPLEVEL,
    LEGACY_RPF_TOPLEVEL,
    STRAIN_SHORTCUTS,
    migrate,
)


# ---------- pure-data migrate() ---------------------------------------------


def test_migrate_renames_legacy_align_keys() -> None:
    cfg = {"align": {"fastq_dir": "input/"}}
    canonical, log = migrate(cfg)
    assert "fastq_dir" not in canonical["align"]
    assert canonical["align"]["fastq"] == "input/"
    assert any("fastq_dir" in line and "fastq" in line for line in log)


def test_migrate_renames_legacy_rpf_keys() -> None:
    cfg = {
        "rpf": {
            "merge_density": True,
            "mrna_ref_patterns": ["mt_genome"],
            "strain": "h",
        }
    }
    canonical, log = migrate(cfg)
    assert "merge_density" not in canonical["rpf"]
    assert canonical["rpf"]["codon_density_window"] is True
    assert "mrna_ref_patterns" not in canonical["rpf"]
    assert canonical["rpf"]["mt_mrna_substring_patterns"] == ["mt_genome"]
    assert canonical["rpf"]["strain"] == "h.sapiens"


def test_migrate_strain_shortcuts() -> None:
    canonical, _ = migrate({"rpf": {"strain": "y"}})
    assert canonical["rpf"]["strain"] == "s.cerevisiae"
    canonical, _ = migrate({"rpf": {"strain": "h"}})
    assert canonical["rpf"]["strain"] == "h.sapiens"
    # Already-canonical values are left alone.
    canonical, _ = migrate({"rpf": {"strain": "h.sapiens"}})
    assert canonical["rpf"]["strain"] == "h.sapiens"


def test_migrate_kit_preset_dropped_with_log() -> None:
    """v0.7.1: kit_preset removed from user-facing input. The migrator
    drops the key with a logged warning so the user knows to add an
    explicit ``adapter:`` / ``pretrimmed:`` value if needed."""
    canonical, log = migrate({"align": {"kit_preset": "truseq_smallrna"}})
    assert "kit_preset" not in canonical["align"]
    assert any("kit_preset" in line and "REMOVED" in line for line in log)


def test_migrate_offset_pick_reference_value_rewrite() -> None:
    canonical, log = migrate(
        {"rpf": {"offset_pick_reference": "selected_site"}}
    )
    assert canonical["rpf"]["offset_pick_reference"] == "reported_site"
    assert any("selected_site" in line for line in log)


def test_migrate_rnaseq_mode_value_rewrite() -> None:
    canonical, log = migrate(
        {"rnaseq": {"mode": "from-fastq"}}
    )
    assert canonical["rnaseq"]["mode"] == "from_fastq"
    assert any("from-fastq" in line and "from_fastq" in line for line in log)


def test_migrate_per_sample_kit_preset_dropped() -> None:
    """`align.samples[*].kit_preset` must also be stripped (v0.7.1)."""
    cfg = {
        "align": {
            "samples": [
                {"name": "A", "kit_preset": "truseq_smallrna"},
                {"name": "B", "kit_preset": "illumina_smallrna"},
            ]
        }
    }
    canonical, log = migrate(cfg)
    samples = canonical["align"]["samples"]
    assert "kit_preset" not in samples[0]
    assert "kit_preset" not in samples[1]
    assert any(
        "samples[0].kit_preset" in line and "REMOVED" in line for line in log
    )


def test_migrate_drops_legacy_when_canonical_already_set() -> None:
    """If the user wrote BOTH the legacy and canonical keys, drop the
    legacy one and log it (no silent shadow)."""
    cfg = {"align": {"fastq_dir": "old/", "fastq": "new/"}}
    canonical, log = migrate(cfg)
    # The canonical 'fastq' value wins.
    assert canonical["align"]["fastq"] == "new/"
    assert "fastq_dir" not in canonical["align"]
    assert any("dropped" in line for line in log)


def test_migrate_flat_config_without_section_wrapper() -> None:
    """Flat (per-stage) configs (no align:/rpf:/rnaseq: wrapper) must
    have rewrites applied at the top level."""
    canonical, log = migrate(
        {"strain": "h", "merge_density": True, "kit_preset": "truseq_smallrna"}
    )
    assert canonical["strain"] == "h.sapiens"
    assert "merge_density" not in canonical
    assert canonical["codon_density_window"] is True
    assert "kit_preset" not in canonical


def test_migrate_no_op_for_canonical_config() -> None:
    cfg = {
        "align": {"fastq": "input/", "adapter": "AGATCGGAAGAGC"},
        "rpf": {"strain": "h.sapiens", "codon_density_window": True},
    }
    canonical, log = migrate(cfg)
    assert canonical == cfg
    assert log == []


def test_migrate_does_not_mutate_input() -> None:
    cfg = {"align": {"fastq_dir": "input/"}}
    snapshot = json.dumps(cfg, sort_keys=True)
    migrate(cfg)
    assert json.dumps(cfg, sort_keys=True) == snapshot


# ---------- CLI: mitoribopy migrate-config ----------------------------------


def _write_yaml(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    return path


def test_migrate_config_cli_stdout_emits_canonical_yaml(
    tmp_path: Path, capsys
) -> None:
    cfg = _write_yaml(
        tmp_path / "old.yaml",
        "align:\n  fastq_dir: input/\n  kit_preset: truseq_smallrna\n"
        "rpf:\n  strain: h\n  merge_density: true\n",
    )
    rc = cli.main(["migrate-config", str(cfg)])
    assert rc == 0
    captured = capsys.readouterr()
    # stdout: canonical YAML (parseable).
    import yaml
    canonical = yaml.safe_load(captured.out)
    assert canonical["align"]["fastq"] == "input/"
    assert "fastq_dir" not in canonical["align"]
    # v0.7.1: kit_preset is dropped, not rewritten.
    assert "kit_preset" not in canonical["align"]
    assert canonical["rpf"]["strain"] == "h.sapiens"
    assert canonical["rpf"]["codon_density_window"] is True
    # stderr: human-readable change log.
    assert "fastq_dir" in captured.err
    assert "kit_preset" in captured.err
    assert "merge_density" in captured.err


def test_migrate_config_cli_no_op_message_when_already_canonical(
    tmp_path: Path, capsys
) -> None:
    cfg = _write_yaml(
        tmp_path / "good.yaml",
        "align:\n  fastq: input/\n",
    )
    rc = cli.main(["migrate-config", str(cfg)])
    assert rc == 0
    captured = capsys.readouterr()
    assert "no legacy keys" in captured.err.lower()


def test_migrate_config_cli_missing_file_exits_2(tmp_path: Path, capsys) -> None:
    rc = cli.main(["migrate-config", str(tmp_path / "nope.yaml")])
    assert rc == 2
    err = capsys.readouterr().err
    assert "ERROR" in err


def test_migrate_config_cli_in_place_preserves_backup(
    tmp_path: Path, capsys
) -> None:
    cfg = _write_yaml(
        tmp_path / "old.yaml",
        "align:\n  fastq_dir: input/\n",
    )
    original = cfg.read_text()
    rc = cli.main(["migrate-config", str(cfg), "--in-place"])
    assert rc == 0
    assert (tmp_path / "old.yaml.bak").exists()
    assert (tmp_path / "old.yaml.bak").read_text() == original
    new = cfg.read_text()
    assert "fastq_dir" not in new
    assert "fastq:" in new


def test_migrate_config_listed_in_top_level_help(capsys) -> None:
    """The new subcommand must appear in `mitoribopy --help`."""
    rc = cli.main(["--help"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "migrate-config" in out
