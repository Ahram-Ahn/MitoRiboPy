"""Tests for §7 rpf refinements (P5.7).

Covers:

* The new ``entropy`` and ``fallback_used`` columns produced by the
  offset selection diagnostics.
* The Shannon-entropy helper returning expected values for delta and
  uniform distributions.
* The ``mitoribopy validate-reference`` subcommand happy-path and
  failure modes (length mismatch, CDS not divisible by 3, bad start
  codon, missing FASTA record).
"""

from __future__ import annotations

import math
import textwrap
from pathlib import Path

import pandas as pd
import pytest


class TestShannonEntropy:
    def test_delta_distribution_zero(self) -> None:
        from mitoribopy.analysis.offset_selection import _shannon_entropy

        # All mass on a single offset -> entropy 0.
        assert _shannon_entropy([0, 0, 100, 0, 0]) == pytest.approx(0.0)

    def test_uniform_distribution_log2_n(self) -> None:
        from mitoribopy.analysis.offset_selection import _shannon_entropy

        # Uniform across 4 offsets -> entropy log2(4) = 2 bits.
        assert _shannon_entropy([10, 10, 10, 10]) == pytest.approx(2.0)

    def test_empty_distribution_zero(self) -> None:
        from mitoribopy.analysis.offset_selection import _shannon_entropy

        assert _shannon_entropy([]) == pytest.approx(0.0)
        assert _shannon_entropy([0, 0, 0]) == pytest.approx(0.0)


class TestOffsetSelectionEnrichment:
    def test_per_length_emits_entropy_and_fallback_columns(
        self, tmp_path: Path
    ) -> None:
        from mitoribopy.analysis.offset_selection import determine_p_site_offsets

        # Build a tiny offsets dataframe with a clear peak at 5'=12 for
        # length 30 (1000 reads) and a flatter distribution at length 31.
        rows = []
        # Length 30: 950 reads at offset 12, 30 at 11, 20 at 13.
        rows.extend([{"Read Length": 30, "5' Offset": 12, "3' Offset": -18}] * 950)
        rows.extend([{"Read Length": 30, "5' Offset": 11, "3' Offset": -19}] * 30)
        rows.extend([{"Read Length": 30, "5' Offset": 13, "3' Offset": -17}] * 20)
        # Length 31: bimodal — 250 at 12 and 240 at 14.
        rows.extend([{"Read Length": 31, "5' Offset": 12, "3' Offset": -19}] * 250)
        rows.extend([{"Read Length": 31, "5' Offset": 14, "3' Offset": -17}] * 240)
        df = pd.DataFrame(rows)

        out = tmp_path / "offsets.csv"
        result = determine_p_site_offsets(
            df,
            align_to="start",
            out_file=str(out),
            offset_min=10,
            offset_max=18,
            offset_mask_nt=5,
            offset_site="p",
            selection_reference="reported_site",
        )
        assert result is not None
        cols = list(result.columns)
        assert "entropy_5" in cols
        assert "fallback_used_5" in cols
        # Length 30 should have low entropy (peaky), length 31 high
        # entropy (bimodal).
        by_length = result.set_index("Read Length")
        assert by_length.loc[30, "entropy_5"] < by_length.loc[31, "entropy_5"]
        # fallback_used defaults to False for every per-length row.
        assert all(v is False for v in result["fallback_used_5"])


def _seed_reference(tmp_path: Path) -> tuple[Path, Path]:
    """Write a minimal valid (FASTA, annotation) pair under *tmp_path*.

    Layout per transcript (l_tr=18, l_utr5=7, l_utr3=5, l_cds=6):
        positions 0-6   UTR5  (7 nt of A or T)
        positions 7-9   START (ATG)
        positions 10-12 STOP  (TAG)
        positions 13-17 UTR3  (5 nt)
    """
    fasta = tmp_path / "ref.fa"
    fasta.write_text(
        textwrap.dedent(
            """\
            >GENE1
            AAAAAAAATGTAGCCCCC
            >GENE2
            TTTTTTTATGTAGCCCCC
            """
        )
    )
    annotation = tmp_path / "ann.csv"
    annotation.write_text(
        "transcript,sequence_name,l_tr,l_utr5,l_utr3\n"
        "GENE1,GENE1,18,7,5\n"
        "GENE2,GENE2,18,7,5\n"
    )
    return fasta, annotation


class TestValidateReference:
    def test_clean_pair_passes(self, tmp_path: Path) -> None:
        from mitoribopy.cli.validate_reference import validate_reference

        fasta, ann = _seed_reference(tmp_path)
        errors, warnings = validate_reference(
            fasta_path=fasta,
            annotation_path=ann,
            codon_table="vertebrate_mitochondrial",
        )
        assert errors == []

    def test_missing_fasta_record_flagged(self, tmp_path: Path) -> None:
        from mitoribopy.cli.validate_reference import validate_reference

        fasta, ann = _seed_reference(tmp_path)
        # Drop GENE2 from the FASTA.
        fasta.write_text(">GENE1\nAAAAAAAATGTAGCCCCC\n")
        errors, _ = validate_reference(
            fasta_path=fasta,
            annotation_path=ann,
            codon_table="vertebrate_mitochondrial",
        )
        assert any("GENE2" in e for e in errors)

    def test_length_mismatch_flagged(self, tmp_path: Path) -> None:
        from mitoribopy.cli.validate_reference import validate_reference

        fasta, ann = _seed_reference(tmp_path)
        # Annotation says 18 but FASTA gives 17 for GENE1.
        fasta.write_text(
            ">GENE1\nAAAAAAATGTAGCCCCC\n>GENE2\nTTTTTTTATGTAGCCCCC\n"
        )
        errors, _ = validate_reference(
            fasta_path=fasta,
            annotation_path=ann,
            codon_table="vertebrate_mitochondrial",
        )
        assert any("l_tr=18" in e and "GENE1" in e for e in errors)

    def test_cds_not_divisible_by_three_flagged(self, tmp_path: Path) -> None:
        from mitoribopy.cli.validate_reference import validate_reference

        fasta, ann = _seed_reference(tmp_path)
        # Make CDS length 7 by editing l_utr3=4 (so 18 - 7 - 4 = 7).
        ann.write_text(
            "transcript,sequence_name,l_tr,l_utr5,l_utr3\n"
            "GENE1,GENE1,18,7,4\n"
            "GENE2,GENE2,18,7,5\n"
        )
        errors, _ = validate_reference(
            fasta_path=fasta,
            annotation_path=ann,
            codon_table="vertebrate_mitochondrial",
        )
        assert any("not divisible by 3" in e for e in errors)

    def test_bad_start_codon_under_standard_table(self, tmp_path: Path) -> None:
        from mitoribopy.cli.validate_reference import validate_reference

        fasta, ann = _seed_reference(tmp_path)
        # GTG is allowed as a start under vertebrate_mitochondrial but
        # NOT under the standard table.
        fasta.write_text(
            ">GENE1\nAAAAAAAGTGTAGCCCCC\n>GENE2\nTTTTTTTATGTAGCCCCC\n"
        )
        errors, _ = validate_reference(
            fasta_path=fasta,
            annotation_path=ann,
            codon_table="standard",
        )
        assert any("start codon" in e for e in errors)

    def test_orphan_fasta_record_warns_only(self, tmp_path: Path) -> None:
        from mitoribopy.cli.validate_reference import validate_reference

        fasta, ann = _seed_reference(tmp_path)
        fasta.write_text(
            ">GENE1\nAAAAAAAATGTAGCCCCC\n"
            ">GENE2\nTTTTTTTATGTAGCCCCC\n"
            ">GENE_EXTRA\nAAAAAAAATGTAGCCCCC\n"
        )
        errors, warnings = validate_reference(
            fasta_path=fasta,
            annotation_path=ann,
            codon_table="vertebrate_mitochondrial",
        )
        assert errors == []
        assert any("GENE_EXTRA" in w for w in warnings)


class TestValidateReferenceCli:
    def test_command_dispatchable_via_main(self, tmp_path: Path, capsys) -> None:
        from mitoribopy.cli import main

        fasta, ann = _seed_reference(tmp_path)
        rc = main(
            [
                "validate-reference",
                "--fasta",
                str(fasta),
                "--annotation",
                str(ann),
            ]
        )
        assert rc == 0
        captured = capsys.readouterr()
        assert "looks consistent" in captured.err
