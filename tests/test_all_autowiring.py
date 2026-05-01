"""§9 P0 — verify the orchestrator's path-handoff between stages.

Three contracts:

1. ``align/bed/`` is wired into ``rpf.directory`` after align runs.
2. ``align/read_counts.tsv`` is wired into ``rpf.read_counts_file``.
3. The de_table flow's rnaseq stage receives ``ribo_dir`` pointing
   at the rpf run directory, so it consumes the rpf stage's
   ``rpf_counts.tsv`` rather than re-counting the Ribo FASTQs.

The check is performed against the post-canonicalisation config dict
the orchestrator builds via ``_auto_wire_paths`` rather than against
a real run, keeping the test hermetic and fast.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import pytest


def _autowire(config: dict, *, run_root: Path, **flags) -> dict:
    """Helper: invoke the private :func:`_auto_wire_paths` entry point.

    ``has_align`` / ``has_rpf`` / ``has_rnaseq`` default to True when
    the key is *present in the config* (even with an empty-dict
    value), matching the orchestrator's intent. The empty-dict case
    is the "stage configured but with all defaults" path the
    orchestrator hits whenever a stage section is declared but not
    populated.
    """
    from mitoribopy.cli import all_ as all_cli

    defaults = dict(
        has_align="align" in config,
        has_rpf="rpf" in config,
        has_rnaseq="rnaseq" in config,
        has_de_table=False,
        has_fastq_mode=False,
    )
    defaults.update(flags)
    all_cli._auto_wire_paths(config, run_root=run_root, **defaults)
    return config


class TestAlignToRpfWiring:
    def test_bed_dir_default_points_to_align_bed(self, tmp_path: Path) -> None:
        cfg = {"align": {}, "rpf": {"strain": "h.sapiens"}}
        out = _autowire(cfg, run_root=tmp_path)
        assert out["rpf"]["directory"] == str(tmp_path / "align" / "bed")

    def test_read_counts_file_default_points_to_align(
        self, tmp_path: Path
    ) -> None:
        cfg = {"align": {}, "rpf": {"strain": "h.sapiens"}}
        out = _autowire(cfg, run_root=tmp_path)
        assert (
            out["rpf"]["read_counts_file"]
            == str(tmp_path / "align" / "read_counts.tsv")
        )

    def test_user_set_rpf_directory_wins_over_default(
        self, tmp_path: Path
    ) -> None:
        explicit = tmp_path / "explicit_bed"
        cfg = {
            "align": {},
            "rpf": {"strain": "h.sapiens", "directory": str(explicit)},
        }
        out = _autowire(cfg, run_root=tmp_path)
        # User's explicit directory must NOT be overridden.
        assert out["rpf"]["directory"] == str(explicit)

    def test_skip_align_does_not_inject_default_rpf_directory(
        self, tmp_path: Path
    ) -> None:
        cfg = {"rpf": {"strain": "h.sapiens"}}
        out = _autowire(cfg, run_root=tmp_path, has_align=False)
        assert "directory" not in out["rpf"]


class TestRpfToRnaseqWiring:
    def test_de_table_flow_points_at_rpf_run_dir(self, tmp_path: Path) -> None:
        cfg = {
            "align": {},
            "rpf": {"strain": "h.sapiens"},
            "rnaseq": {"de_table": "de.tsv"},
        }
        out = _autowire(cfg, run_root=tmp_path, has_de_table=True)
        # Note: hyphen-form key is the orchestrator's canonical write
        # site (matches the CLI flag style); both are accepted.
        assert (
            out["rnaseq"].get("ribo-dir")
            or out["rnaseq"].get("ribo_dir")
        ) == str(tmp_path / "rpf")

    def test_rnaseq_output_defaults_to_run_root(self, tmp_path: Path) -> None:
        cfg = {
            "align": {},
            "rpf": {"strain": "h.sapiens"},
            "rnaseq": {"de_table": "de.tsv"},
        }
        out = _autowire(cfg, run_root=tmp_path, has_de_table=True)
        assert out["rnaseq"]["output"] == str(tmp_path / "rnaseq")

    def test_from_fastq_inherits_reference_fasta_from_rpf(
        self, tmp_path: Path
    ) -> None:
        rpf_fa = tmp_path / "tx.fa"
        rpf_fa.write_text(">x\nACGT\n")
        cfg = {
            "align": {},
            "rpf": {"strain": "h.sapiens", "fasta": str(rpf_fa)},
            "rnaseq": {"rna_fastq": ["a.fq.gz"]},
        }
        out = _autowire(cfg, run_root=tmp_path, has_fastq_mode=True)
        # rnaseq.reference_fasta should fall back to rpf.fasta when not set.
        assert out["rnaseq"]["reference_fasta"] == str(rpf_fa)


class TestRnaseqReusesUpstreamRpfCounts:
    """Section 9 spec: the de_table flow must NOT re-process Ribo
    FASTQs when an upstream rpf_counts.tsv already exists. A user
    explicitly opts in via ``rnaseq.recount_ribo_fastq: true`` /
    ``--recount-ribo-fastq``.

    The autowiring step is the only place the orchestrator imposes
    that default; the from-FASTQ orchestrator inside
    ``mitoribopy.cli.rnaseq`` honours it. We assert here that the
    autowiring sets up the upstream-counts pointer when both stages
    are configured.
    """

    def test_from_fastq_with_rpf_routes_to_upstream_counts(
        self, tmp_path: Path
    ) -> None:
        rpf_fa = tmp_path / "tx.fa"
        rpf_fa.write_text(">x\nACGT\n")
        cfg = {
            "align": {},
            "rpf": {"strain": "h.sapiens", "fasta": str(rpf_fa)},
            "rnaseq": {
                "rna_fastq": ["a.fq.gz"],
                "reference_fasta": str(rpf_fa),
            },
        }
        out = _autowire(cfg, run_root=tmp_path, has_fastq_mode=True)
        # The upstream rpf_counts pointer is the explicit handoff that
        # makes the from-FASTQ flow skip Ribo re-alignment.
        upstream_pointer = (
            out["rnaseq"].get("upstream_rpf_counts")
            or out["rnaseq"].get("upstream-rpf-counts")
        )
        assert upstream_pointer == str(tmp_path / "rpf" / "rpf_counts.tsv"), (
            "rnaseq stage was not auto-wired to consume the rpf "
            "stage's rpf_counts.tsv; the from-FASTQ flow would "
            "re-align Ribo FASTQs unnecessarily."
        )
