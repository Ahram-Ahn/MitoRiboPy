"""Tests for §6 align refinements (P5.6).

Covers:

* :func:`mitoribopy.align.read_counts.check_count_invariants` /
  :func:`enforce_count_invariants` rule set.
* :func:`mitoribopy.align.sample_resolve._second_best_from_rates`
  helper and the new kit_resolution.tsv columns.
* The :func:`mitoribopy.cli.align._strict_publication_mode_errors`
  gate.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mitoribopy.align._types import KIT_PRESETS, ResolvedKit, SampleCounts
from mitoribopy.align.read_counts import (
    CountInvariantError,
    CountInvariantViolation,
    check_count_invariants,
    enforce_count_invariants,
)
from mitoribopy.align.sample_resolve import (
    SampleResolution,
    _second_best_from_rates,
    write_kit_resolution_tsv,
)


def _good_counts(sample: str = "A") -> SampleCounts:
    return SampleCounts(
        sample=sample,
        total_reads=1000,
        post_trim=900,
        rrna_aligned=200,
        post_rrna_filter=700,
        mt_aligned=600,
        unaligned_to_mt=100,
        mt_aligned_after_mapq=500,
        mt_aligned_after_dedup=400,
    )


class TestCheckCountInvariants:
    def test_clean_rows_no_violations(self) -> None:
        assert check_count_invariants([_good_counts()]) == []

    def test_rrna_split_violation(self) -> None:
        bad = _good_counts()
        bad = SampleCounts(
            **{**bad.__dict__, "post_rrna_filter": 600},  # 200 + 600 != 900
        )
        violations = check_count_invariants([bad])
        rules = [v.rule for v in violations]
        assert "rrna_split" in rules

    def test_mt_split_violation(self) -> None:
        bad = SampleCounts(
            **{**_good_counts().__dict__, "unaligned_to_mt": 50},
        )  # 600 + 50 != 700
        violations = check_count_invariants([bad])
        assert any(v.rule == "mt_split" for v in violations)

    def test_mapq_le_mt_violation(self) -> None:
        bad = SampleCounts(
            **{**_good_counts().__dict__, "mt_aligned_after_mapq": 700},
        )
        violations = check_count_invariants([bad])
        assert any(v.rule == "mapq_le_mt" for v in violations)

    def test_dedup_le_mapq_violation(self) -> None:
        bad = SampleCounts(
            **{**_good_counts().__dict__, "mt_aligned_after_dedup": 600},
        )
        violations = check_count_invariants([bad])
        assert any(v.rule == "dedup_le_mapq" for v in violations)


class TestEnforceCountInvariants:
    def test_clean_rows_returns_empty_list(self) -> None:
        out = enforce_count_invariants([_good_counts()])
        assert out == []

    def test_default_raises_on_violation(self) -> None:
        bad = SampleCounts(
            **{**_good_counts().__dict__, "post_rrna_filter": 600},
        )
        with pytest.raises(CountInvariantError) as exc_info:
            enforce_count_invariants([bad])
        assert exc_info.value.violations
        assert exc_info.value.violations[0].rule == "rrna_split"

    def test_allow_warning_records_warnings_log_and_returns(self) -> None:
        from mitoribopy.io import warnings_log

        warnings_log.clear()
        bad = SampleCounts(
            **{**_good_counts().__dict__, "post_rrna_filter": 600},
        )
        out = enforce_count_invariants([bad], allow_warning=True)
        assert out  # the violation list, not raised
        records = warnings_log.collected()
        assert any(
            r.code == "E_OUTPUT_COUNT_INVARIANT" for r in records
        )
        warnings_log.clear()


class TestSecondBestFromRates:
    def test_returns_runner_up_and_margin(self) -> None:
        rates = {"alpha": 0.95, "beta": 0.40, "gamma": 0.10}
        sec, sec_rate, margin = _second_best_from_rates(
            rates, best_name="alpha"
        )
        assert sec == "beta"
        assert sec_rate == pytest.approx(0.40)
        assert margin == pytest.approx(0.55)

    def test_empty_rates_returns_zeros(self) -> None:
        sec, sec_rate, margin = _second_best_from_rates({}, best_name=None)
        assert sec is None
        assert sec_rate == 0.0
        assert margin == 0.0

    def test_single_kit_no_second_best(self) -> None:
        rates = {"only": 0.90}
        sec, sec_rate, margin = _second_best_from_rates(
            rates, best_name="only"
        )
        assert sec is None
        assert sec_rate == 0.0
        # margin falls back to the best rate when there is no runner-up
        assert margin == pytest.approx(0.90)


class TestKitResolutionTsvColumns:
    def test_new_columns_present_in_writer_output(self, tmp_path: Path) -> None:
        fastq = tmp_path / "x.fq.gz"
        fastq.write_bytes(b"")
        res = SampleResolution(
            sample="A",
            fastq=fastq,
            kit=ResolvedKit(
                kit="illumina_truseq",
                adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                umi_length=0,
                umi_position="5p",
            ),
            dedup_strategy="skip",
            detected_kit="illumina_truseq",
            detection_match_rate=0.92,
            detection_ambiguous=False,
            source="detected",
            umi_source="none",
            second_best_kit="illumina_smallrna",
            second_best_match_rate=0.04,
            confidence_margin=0.88,
        )
        out = tmp_path / "kit_resolution.tsv"
        write_kit_resolution_tsv([res], out)
        lines = [
            l for l in out.read_text().splitlines() if not l.startswith("#")
        ]
        header = lines[0].split("\t")
        for col in (
            "best_adapter",
            "best_match_rate",
            "second_best_kit",
            "second_best_match_rate",
            "confidence_margin",
        ):
            assert col in header, f"missing column: {col}"
        body = dict(zip(header, lines[1].split("\t")))
        assert body["second_best_kit"] == "illumina_smallrna"
        assert body["second_best_match_rate"].startswith("0.04")
        assert body["confidence_margin"].startswith("0.88")


class TestStrictPublicationMode:
    def _resolution(self, **overrides) -> SampleResolution:
        from pathlib import Path

        defaults = dict(
            sample="A",
            fastq=Path("/tmp/A.fq.gz"),
            kit=ResolvedKit(
                kit="illumina_truseq",
                adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                umi_length=0,
                umi_position="5p",
            ),
            dedup_strategy="skip",
            detected_kit="illumina_truseq",
            detection_match_rate=0.95,
            detection_ambiguous=False,
            source="detected",
            umi_source="none",
            second_best_kit="illumina_smallrna",
            second_best_match_rate=0.05,
            confidence_margin=0.90,
        )
        defaults.update(overrides)
        return SampleResolution(**defaults)

    def test_clean_resolution_passes(self) -> None:
        from mitoribopy.cli.align import _strict_publication_mode_errors

        assert _strict_publication_mode_errors([self._resolution()]) == []

    def test_inferred_pretrimmed_rejected(self) -> None:
        from mitoribopy.cli.align import _strict_publication_mode_errors

        r = self._resolution(source="inferred_pretrimmed")
        errs = _strict_publication_mode_errors([r])
        assert errs and "pretrimmed" in errs[0]

    def test_ambiguous_detection_rejected(self) -> None:
        from mitoribopy.cli.align import _strict_publication_mode_errors

        r = self._resolution(detection_ambiguous=True, confidence_margin=0.02)
        errs = _strict_publication_mode_errors([r])
        assert errs and "ambiguous" in errs[0]

    def test_low_confidence_margin_rejected(self) -> None:
        from mitoribopy.cli.align import _strict_publication_mode_errors

        r = self._resolution(
            detection_ambiguous=False,
            detection_match_rate=0.55,
            confidence_margin=0.05,  # below 0.10 threshold
        )
        errs = _strict_publication_mode_errors([r])
        assert errs and "ambiguous" in errs[0]

    def test_inferred_umi_rejected(self) -> None:
        from mitoribopy.cli.align import _strict_publication_mode_errors

        r = self._resolution(
            kit=ResolvedKit(
                kit="illumina_truseq_umi",
                adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                umi_length=8,
                umi_position="5p",
            ),
            umi_source="preset_default",  # not declared
        )
        errs = _strict_publication_mode_errors([r])
        assert errs and "UMI" in errs[0]
