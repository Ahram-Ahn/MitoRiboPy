"""Unit tests for ``mitoribopy.align.sample_resolve``."""

from __future__ import annotations

from pathlib import Path

import pytest

from mitoribopy.align.adapter_detect import DetectionResult
from mitoribopy.align.sample_resolve import (
    SampleResolutionError,
    required_dedup_tools,
    resolve_sample_resolutions,
    resolution_summary_lines,
    write_kit_resolution_tsv,
)


def _detection(name: str | None, rate: float = 0.95, ambiguous: bool = False) -> DetectionResult:
    return DetectionResult(
        preset_name=name,
        match_rate=rate,
        per_kit_rates={"truseq_smallrna": rate, "nebnext_smallrna": 0.0},
        n_reads_scanned=5000,
        ambiguous=ambiguous,
    )


def _detector(map_: dict[str, DetectionResult]):
    """Return a detector callable that looks up DetectionResult by filename."""

    def detector(path: Path) -> DetectionResult:
        return map_[Path(path).name]

    return detector


# ---------- happy path -------------------------------------------------------


def test_auto_detection_picks_per_sample_kit_independently(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    b = tmp_path / "b.fq.gz"
    a.write_text("")
    b.write_text("")

    detector = _detector(
        {
            "a.fq.gz": _detection("truseq_smallrna"),
            "b.fq.gz": _detection("nebnext_ultra_umi"),
        }
    )

    resolutions = resolve_sample_resolutions(
        [("a", a), ("b", b)],
        kit_preset="auto",
        adapter=None,
        umi_length=None,
        umi_position=None,
        dedup_strategy="auto",
        confirm_mark_duplicates=False,
        adapter_detection_mode="auto",
        detector=detector,
    )

    assert resolutions[0].kit.kit == "truseq_smallrna"
    assert resolutions[0].dedup_strategy == "skip"  # no UMI
    assert resolutions[1].kit.kit == "nebnext_ultra_umi"
    assert resolutions[1].dedup_strategy == "umi-tools"  # 8 nt UMI
    # Both samples have a recorded detection source.
    assert all(r.source == "detected" for r in resolutions)


def test_auto_falls_back_to_user_kit_when_detection_fails(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")

    detector = _detector({"a.fq.gz": _detection(None, rate=0.0)})

    resolutions = resolve_sample_resolutions(
        [("a", a)],
        kit_preset="truseq_smallrna",  # explicit fallback
        adapter=None,
        umi_length=None,
        umi_position=None,
        dedup_strategy="auto",
        confirm_mark_duplicates=False,
        adapter_detection_mode="auto",
        detector=detector,
    )
    assert resolutions[0].kit.kit == "truseq_smallrna"
    assert resolutions[0].source == "user_fallback"


def test_auto_hard_fails_when_no_detection_and_no_fallback(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")

    detector = _detector({"a.fq.gz": _detection(None, rate=0.0)})

    with pytest.raises(SampleResolutionError) as exc:
        resolve_sample_resolutions(
            [("a", a)],
            kit_preset="auto",  # no fallback
            adapter=None,
            umi_length=None,
            umi_position=None,
            dedup_strategy="auto",
            confirm_mark_duplicates=False,
            adapter_detection_mode="auto",
            detector=detector,
        )
    assert "auto-detection found no known kit" in str(exc.value)


# ---------- strict mode ------------------------------------------------------


def test_strict_hard_fails_on_disagreement(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")

    detector = _detector({"a.fq.gz": _detection("nebnext_smallrna", rate=0.9)})

    with pytest.raises(SampleResolutionError) as exc:
        resolve_sample_resolutions(
            [("a", a)],
            kit_preset="truseq_smallrna",
            adapter=None,
            umi_length=None,
            umi_position=None,
            dedup_strategy="auto",
            confirm_mark_duplicates=False,
            adapter_detection_mode="strict",
            detector=detector,
        )
    assert "strict" in str(exc.value).lower()
    assert "nebnext_smallrna" in str(exc.value)


def test_strict_uses_detected_when_user_left_auto(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")

    detector = _detector({"a.fq.gz": _detection("truseq_smallrna")})

    resolutions = resolve_sample_resolutions(
        [("a", a)],
        kit_preset="auto",
        adapter=None,
        umi_length=None,
        umi_position=None,
        dedup_strategy="auto",
        confirm_mark_duplicates=False,
        adapter_detection_mode="strict",
        detector=detector,
    )
    assert resolutions[0].kit.kit == "truseq_smallrna"
    assert resolutions[0].source == "detected"


# ---------- off mode --------------------------------------------------------


def test_off_mode_requires_explicit_preset(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")

    with pytest.raises(SampleResolutionError) as exc:
        resolve_sample_resolutions(
            [("a", a)],
            kit_preset="auto",
            adapter=None,
            umi_length=None,
            umi_position=None,
            dedup_strategy="auto",
            confirm_mark_duplicates=False,
            adapter_detection_mode="off",
        )
    assert "--adapter-detection off requires" in str(exc.value)


def test_off_mode_uses_explicit_preset_without_scanning(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")

    def explode(_path):
        raise AssertionError("detector must NOT run in off mode")

    resolutions = resolve_sample_resolutions(
        [("a", a)],
        kit_preset="truseq_smallrna",
        adapter=None,
        umi_length=None,
        umi_position=None,
        dedup_strategy="auto",
        confirm_mark_duplicates=False,
        adapter_detection_mode="off",
        detector=explode,
    )
    assert resolutions[0].kit.kit == "truseq_smallrna"
    assert resolutions[0].source == "explicit_off"


# ---------- mixed-kit error aggregation -------------------------------------


def test_resolver_aggregates_per_sample_errors(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    b = tmp_path / "b.fq.gz"
    a.write_text("")
    b.write_text("")

    detector = _detector(
        {
            "a.fq.gz": _detection(None, rate=0.0),
            "b.fq.gz": _detection(None, rate=0.0),
        }
    )

    with pytest.raises(SampleResolutionError) as exc:
        resolve_sample_resolutions(
            [("a", a), ("b", b)],
            kit_preset="auto",
            adapter=None,
            umi_length=None,
            umi_position=None,
            dedup_strategy="auto",
            confirm_mark_duplicates=False,
            adapter_detection_mode="auto",
            detector=detector,
        )
    msg = str(exc.value)
    # Both failing samples are named in the aggregated message so the
    # user can fix the config in one pass.
    assert "a:" in msg
    assert "b:" in msg


# ---------- helpers ---------------------------------------------------------


def test_required_dedup_tools_unions_per_sample_strategies(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    b = tmp_path / "b.fq.gz"
    a.write_text("")
    b.write_text("")

    detector = _detector(
        {
            "a.fq.gz": _detection("truseq_smallrna"),  # no UMI -> skip
            "b.fq.gz": _detection("nebnext_ultra_umi"),  # UMI -> umi-tools
        }
    )
    resolutions = resolve_sample_resolutions(
        [("a", a), ("b", b)],
        kit_preset="auto",
        adapter=None,
        umi_length=None,
        umi_position=None,
        dedup_strategy="auto",
        confirm_mark_duplicates=False,
        adapter_detection_mode="auto",
        detector=detector,
    )
    assert required_dedup_tools(resolutions) == {"umi_tools"}


def test_write_kit_resolution_tsv_writes_header_and_one_row_per_sample(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")
    detector = _detector({"a.fq.gz": _detection("truseq_smallrna")})
    resolutions = resolve_sample_resolutions(
        [("a", a)],
        kit_preset="auto",
        adapter=None,
        umi_length=None,
        umi_position=None,
        dedup_strategy="auto",
        confirm_mark_duplicates=False,
        adapter_detection_mode="auto",
        detector=detector,
    )
    out = tmp_path / "kit_resolution.tsv"
    write_kit_resolution_tsv(resolutions, out)
    lines = out.read_text().splitlines()
    assert lines[0].startswith("sample\tfastq\tapplied_kit")
    assert "truseq_smallrna" in lines[1]
    assert "skip" in lines[1]


def test_resolution_summary_lines_includes_detected_when_overridden(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")
    # Detector returns NEB but user explicitly chose TruSeq → user_fallback.
    detector = _detector({"a.fq.gz": _detection("nebnext_smallrna", rate=0.9)})
    resolutions = resolve_sample_resolutions(
        [("a", a)],
        kit_preset="truseq_smallrna",
        adapter=None,
        umi_length=None,
        umi_position=None,
        dedup_strategy="auto",
        confirm_mark_duplicates=False,
        adapter_detection_mode="auto",
        detector=detector,
    )
    line = resolution_summary_lines(resolutions)[0]
    assert "kit=truseq_smallrna" in line
    assert "detected=nebnext_smallrna" in line
    assert "user_fallback" in line
