"""Tests for :mod:`mitoribopy.config.canonical` (§10 foundations)."""

from __future__ import annotations

import pytest

from mitoribopy.config import (
    CanonicalConfig,
    ConfigChange,
    canonicalize_config,
)


class TestCanonicalizeConfig:
    def test_already_canonical_no_changes(self) -> None:
        cfg = {
            "align": {"adapter": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "fastq": "fq/"},
            "rpf": {"strain": "h.sapiens"},
            "rnaseq": {"mode": "de_table"},
        }
        result = canonicalize_config(cfg)
        assert isinstance(result, CanonicalConfig)
        assert result.changes == []
        assert result.config["rpf"]["strain"] == "h.sapiens"

    def test_legacy_strain_shortcut_rewritten(self) -> None:
        result = canonicalize_config({"rpf": {"strain": "h"}})
        assert result.config["rpf"]["strain"] == "h.sapiens"
        assert any("strain" in c.message for c in result.changes)
        assert all(c.severity == "info" for c in result.changes)

    def test_legacy_align_fastq_dir_rewritten(self) -> None:
        result = canonicalize_config(
            {"align": {"fastq_dir": "/data/fq"}, "rpf": {}}
        )
        assert "fastq_dir" not in result.config["align"]
        assert result.config["align"]["fastq"] == "/data/fq"

    def test_legacy_kit_preset_dropped_with_log(self) -> None:
        """v0.7.1: kit_preset removed; the canonicaliser drops the key
        and logs the removal so the user sees what changed."""
        result = canonicalize_config(
            {"align": {"kit_preset": "truseq_smallrna"}}
        )
        assert "kit_preset" not in result.config["align"]
        assert any(
            "kit_preset" in c.message for c in result.changes
        )

    def test_unknown_top_level_section_emits_warning(self) -> None:
        result = canonicalize_config(
            {
                "align": {},
                "rpf": {},
                "totally_made_up_section": {"foo": 1},
            }
        )
        assert result.has_warnings
        warns = result.warning_messages()
        assert any("totally_made_up_section" in w for w in warns)
        # The unknown section is preserved in config, not silently dropped.
        assert "totally_made_up_section" in result.config

    def test_flat_per_stage_config_no_section_warning(self) -> None:
        # A flat rpf-style config has no align/rpf/rnaseq wrappers; the
        # unknown-section check should NOT fire.
        result = canonicalize_config(
            {"strain": "h", "directory": ".", "range": 20}
        )
        assert result.has_warnings is False
        # The strain shortcut DID rewrite, so info-level changes can still exist.
        assert any(c.severity == "info" for c in result.changes)

    def test_rejects_non_dict(self) -> None:
        with pytest.raises(TypeError):
            canonicalize_config("not a dict")  # type: ignore[arg-type]

    def test_rnaseq_mode_value_rewrite(self) -> None:
        result = canonicalize_config(
            {"align": {}, "rpf": {}, "rnaseq": {"mode": "from-fastq"}}
        )
        assert result.config["rnaseq"]["mode"] == "from_fastq"


class TestErrorCodes:
    def test_error_codes_are_string_constants(self) -> None:
        from mitoribopy import errors as err

        for name in (
            "E_CONFIG_UNKNOWN_KEY",
            "E_CONFIG_MUTUALLY_EXCLUSIVE",
            "E_SAMPLE_SHEET_MISSING_COLUMN",
            "E_REFERENCE_HASH_MISMATCH",
            "E_OFFSET_LOW_CONFIDENCE",
            "E_RNASEQ_INSUFFICIENT_REPLICATES",
            "E_TOOL_NOT_FOUND",
        ):
            value = getattr(err, name)
            assert isinstance(value, str)
            assert value == name, (
                f"Error code {name} value drifted from its name: {value!r}"
            )


class TestStageResults:
    def test_align_result_round_trip(self, tmp_path) -> None:
        from mitoribopy.pipeline import AlignResult

        bed = tmp_path / "bed"
        bed.mkdir()
        rc = tmp_path / "read_counts.tsv"
        rc.touch()
        kit = tmp_path / "kit_resolution.tsv"
        kit.touch()
        rs = tmp_path / "run_settings.json"
        rs.touch()

        r = AlignResult(
            bed_dir=bed, read_counts=rc, kit_resolution=kit, run_settings=rs
        )
        d = r.as_dict()
        assert d["bed_dir"] == str(bed)
        assert d["bam_dir"] is None

    def test_rpf_result_carries_reference_checksum(self, tmp_path) -> None:
        from mitoribopy.pipeline import RpfResult

        rpf = tmp_path / "rpf_counts.tsv"
        rpf.touch()
        diag = tmp_path / "diag"
        diag.mkdir()
        rs = tmp_path / "run_settings.json"
        rs.touch()
        r = RpfResult(
            rpf_counts=rpf,
            offset_diagnostics_dir=diag,
            run_settings=rs,
            reference_checksum="abc123",
        )
        assert r.reference_checksum == "abc123"
        assert r.as_dict()["reference_checksum"] == "abc123"

    def test_rnaseq_result_optional_fields(self, tmp_path) -> None:
        from mitoribopy.pipeline import RnaseqResult

        rs = tmp_path / "run_settings.json"
        rs.touch()
        r = RnaseqResult(mode="rna_only", run_settings=rs)
        assert r.te_table is None
        assert r.delta_te_table is None
        assert r.as_dict()["mode"] == "rna_only"


class TestCanonicalStages:
    def test_canonical_stage_constants(self) -> None:
        from mitoribopy.console import (
            CANONICAL_STAGES,
            STAGE_ALIGN,
            STAGE_RNASEQ,
            STAGE_RPF,
        )

        assert STAGE_ALIGN == "align"
        assert STAGE_RPF == "rpf"
        assert STAGE_RNASEQ == "rnaseq"
        assert "config" in CANONICAL_STAGES
        assert "sample-sheet" in CANONICAL_STAGES
        assert "resume" in CANONICAL_STAGES
