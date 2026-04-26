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


def _detection(
    name: str | None,
    rate: float = 0.95,
    ambiguous: bool = False,
    pretrimmed: bool = False,
) -> DetectionResult:
    return DetectionResult(
        preset_name=name,
        match_rate=rate,
        per_kit_rates={"illumina_smallrna": rate, "illumina_truseq": 0.0},
        n_reads_scanned=5000,
        ambiguous=ambiguous,
        pretrimmed=pretrimmed,
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
            "a.fq.gz": _detection("illumina_smallrna"),
            "b.fq.gz": _detection("illumina_truseq_umi"),
        }
    )

    resolutions = resolve_sample_resolutions(
        [("a", a), ("b", b)],
        kit_preset="auto",
        adapter=None,
        umi_length=None,
        umi_position=None,
        dedup_strategy="auto",
        adapter_detection_mode="auto",
        detector=detector,
    )

    assert resolutions[0].kit.kit == "illumina_smallrna"
    assert resolutions[0].dedup_strategy == "skip"  # no UMI
    assert resolutions[1].kit.kit == "illumina_truseq_umi"
    assert resolutions[1].dedup_strategy == "umi-tools"  # 8 nt UMI
    # Both samples have a recorded detection source.
    assert all(r.source == "detected" for r in resolutions)


def test_auto_falls_back_to_user_kit_when_detection_fails(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")

    detector = _detector({"a.fq.gz": _detection(None, rate=0.0)})

    resolutions = resolve_sample_resolutions(
        [("a", a)],
        kit_preset="illumina_smallrna",  # explicit fallback
        adapter=None,
        umi_length=None,
        umi_position=None,
        dedup_strategy="auto",
        adapter_detection_mode="auto",
        detector=detector,
    )
    assert resolutions[0].kit.kit == "illumina_smallrna"
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
            adapter_detection_mode="auto",
            detector=detector,
        )
    assert "auto-detection found no known kit" in str(exc.value)


# ---------- strict mode ------------------------------------------------------


def test_strict_hard_fails_on_disagreement(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")

    detector = _detector({"a.fq.gz": _detection("illumina_truseq", rate=0.9)})

    with pytest.raises(SampleResolutionError) as exc:
        resolve_sample_resolutions(
            [("a", a)],
            kit_preset="illumina_smallrna",
            adapter=None,
            umi_length=None,
            umi_position=None,
            dedup_strategy="auto",
            adapter_detection_mode="strict",
            detector=detector,
        )
    assert "strict" in str(exc.value).lower()
    assert "illumina_truseq" in str(exc.value)


def test_strict_uses_detected_when_user_left_auto(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")

    detector = _detector({"a.fq.gz": _detection("illumina_smallrna")})

    resolutions = resolve_sample_resolutions(
        [("a", a)],
        kit_preset="auto",
        adapter=None,
        umi_length=None,
        umi_position=None,
        dedup_strategy="auto",
        adapter_detection_mode="strict",
        detector=detector,
    )
    assert resolutions[0].kit.kit == "illumina_smallrna"
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
        kit_preset="illumina_smallrna",
        adapter=None,
        umi_length=None,
        umi_position=None,
        dedup_strategy="auto",
        adapter_detection_mode="off",
        detector=explode,
    )
    assert resolutions[0].kit.kit == "illumina_smallrna"
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
            "a.fq.gz": _detection("illumina_smallrna"),  # no UMI -> skip
            "b.fq.gz": _detection("illumina_truseq_umi"),  # UMI -> umi-tools
        }
    )
    resolutions = resolve_sample_resolutions(
        [("a", a), ("b", b)],
        kit_preset="auto",
        adapter=None,
        umi_length=None,
        umi_position=None,
        dedup_strategy="auto",
        adapter_detection_mode="auto",
        detector=detector,
    )
    assert required_dedup_tools(resolutions) == {"umi_tools"}


def test_write_kit_resolution_tsv_writes_header_and_one_row_per_sample(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")
    detector = _detector({"a.fq.gz": _detection("illumina_smallrna")})
    resolutions = resolve_sample_resolutions(
        [("a", a)],
        kit_preset="auto",
        adapter=None,
        umi_length=None,
        umi_position=None,
        dedup_strategy="auto",
        adapter_detection_mode="auto",
        detector=detector,
    )
    out = tmp_path / "kit_resolution.tsv"
    write_kit_resolution_tsv(resolutions, out)
    lines = out.read_text().splitlines()
    assert lines[0].startswith("sample\tfastq\tapplied_kit")
    assert "illumina_smallrna" in lines[1]
    assert "skip" in lines[1]


def test_auto_infers_pretrimmed_when_no_kit_matches_and_no_fallback(tmp_path) -> None:
    """The fix for SRA-deposited / already-trimmed FASTQs: detection
    finds 0% across every kit, the scanner reports ``pretrimmed=True``,
    and the resolver falls through to the ``pretrimmed`` kit instead of
    erroring."""
    a = tmp_path / "a.fq.gz"
    a.write_text("")
    detector = _detector(
        {"a.fq.gz": _detection(None, rate=0.0, pretrimmed=True)}
    )

    resolutions = resolve_sample_resolutions(
        [("a", a)],
        kit_preset="auto",
        adapter=None,
        umi_length=None,
        umi_position=None,
        dedup_strategy="auto",
        adapter_detection_mode="auto",
        detector=detector,
    )
    assert resolutions[0].kit.kit == "pretrimmed"
    assert resolutions[0].kit.adapter is None
    assert resolutions[0].source == "inferred_pretrimmed"


def test_auto_pretrimmed_inference_can_be_disabled(tmp_path) -> None:
    """``allow_pretrimmed_inference=False`` reverts to the previous
    behaviour: detection failure with no fallback raises."""
    a = tmp_path / "a.fq.gz"
    a.write_text("")
    detector = _detector(
        {"a.fq.gz": _detection(None, rate=0.0, pretrimmed=True)}
    )

    with pytest.raises(SampleResolutionError):
        resolve_sample_resolutions(
            [("a", a)],
            kit_preset="auto",
            adapter=None,
            umi_length=None,
            umi_position=None,
            dedup_strategy="auto",
            adapter_detection_mode="auto",
            allow_pretrimmed_inference=False,
            detector=detector,
        )


def test_explicit_pretrimmed_kit_preset_resolves_with_no_adapter(tmp_path) -> None:
    """``--kit-preset pretrimmed`` resolves to adapter=None so cutadapt
    skips the -a flag entirely."""
    a = tmp_path / "a.fq.gz"
    a.write_text("")
    # Detector returns a real kit, but the user explicitly chose
    # pretrimmed → user's choice wins via the user_fallback path.
    detector = _detector({"a.fq.gz": _detection("illumina_smallrna")})

    resolutions = resolve_sample_resolutions(
        [("a", a)],
        kit_preset="pretrimmed",
        adapter=None,
        umi_length=None,
        umi_position=None,
        dedup_strategy="auto",
        adapter_detection_mode="auto",
        detector=detector,
    )
    assert resolutions[0].kit.kit == "pretrimmed"
    assert resolutions[0].kit.adapter is None
    # Dedup falls through to skip (no UMI on the pretrimmed preset).
    assert resolutions[0].dedup_strategy == "skip"


def test_mixed_batch_pretrimmed_and_raw_samples(tmp_path) -> None:
    """A batch with one pre-trimmed FASTQ and one raw FASTQ resolves
    each sample independently — no global decision contaminates the
    other sample."""
    a = tmp_path / "a.fq.gz"
    b = tmp_path / "b.fq.gz"
    a.write_text("")
    b.write_text("")
    detector = _detector(
        {
            "a.fq.gz": _detection(None, rate=0.0, pretrimmed=True),  # SRA-trimmed
            "b.fq.gz": _detection("illumina_smallrna"),                # raw
        }
    )

    resolutions = resolve_sample_resolutions(
        [("a", a), ("b", b)],
        kit_preset="auto",
        adapter=None,
        umi_length=None,
        umi_position=None,
        dedup_strategy="auto",
        adapter_detection_mode="auto",
        detector=detector,
    )
    by_sample = {res.sample: res for res in resolutions}
    assert by_sample["a"].kit.kit == "pretrimmed"
    assert by_sample["a"].source == "inferred_pretrimmed"
    assert by_sample["b"].kit.kit == "illumina_smallrna"
    assert by_sample["b"].source == "detected"


def test_legacy_kit_aliases_resolve_to_canonical_presets(tmp_path) -> None:
    """v0.4.0 vendor names (truseq_smallrna, nebnext_smallrna,
    nebnext_ultra_umi, …) must keep working as --kit-preset choices in
    YAML configs and CLI invocations even after the v0.4.1
    consolidation; they translate to the canonical adapter-family
    preset name through ``resolve_kit_alias``.
    """
    a = tmp_path / "a.fq.gz"
    a.write_text("")
    detector = _detector({"a.fq.gz": _detection("illumina_smallrna")})

    # User passes the legacy alias; resolver must translate it to the
    # canonical illumina_smallrna preset.
    resolutions = resolve_sample_resolutions(
        [("a", a)],
        kit_preset="truseq_smallrna",  # legacy alias
        adapter=None,
        umi_length=None,
        umi_position=None,
        dedup_strategy="auto",
        adapter_detection_mode="off",  # no scan; trust the alias
    )
    assert resolutions[0].kit.kit == "illumina_smallrna"
    assert resolutions[0].kit.adapter == "TGGAATTCTCGGGTGCCAAGG"


def test_resolution_summary_lines_includes_detected_when_overridden(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")
    # Detector returns NEB but user explicitly chose TruSeq → user_fallback.
    detector = _detector({"a.fq.gz": _detection("illumina_truseq", rate=0.9)})
    resolutions = resolve_sample_resolutions(
        [("a", a)],
        kit_preset="illumina_smallrna",
        adapter=None,
        umi_length=None,
        umi_position=None,
        dedup_strategy="auto",
        adapter_detection_mode="auto",
        detector=detector,
    )
    line = resolution_summary_lines(resolutions)[0]
    assert "kit=illumina_smallrna" in line
    assert "detected=illumina_truseq" in line
    assert "user_fallback" in line
