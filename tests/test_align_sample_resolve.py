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


# v0.7.1: kit names are an internal label only — users supply
# ``--adapter <SEQ>`` (the literal 3' adapter sequence) or ``--pretrimmed``.
ILLUMINA_SMALLRNA_ADAPTER = "TGGAATTCTCGGGTGCCAAGG"
ILLUMINA_TRUSEQ_ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"


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


def _resolve(samples, **kwargs):
    """Tiny helper to fill defaults that every test needs."""
    base: dict = {
        "adapter": None,
        "pretrimmed": False,
        "umi_length": None,
        "umi_position": None,
        "dedup_strategy": "auto",
        "adapter_detection_mode": "auto",
    }
    base.update(kwargs)
    return resolve_sample_resolutions(samples, **base)


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

    resolutions = _resolve(
        [("a", a), ("b", b)],
        detector=detector,
    )

    assert resolutions[0].kit.kit == "illumina_smallrna"
    assert resolutions[0].dedup_strategy == "skip"  # no UMI
    assert resolutions[1].kit.kit == "illumina_truseq_umi"
    assert resolutions[1].dedup_strategy == "umi-tools"  # 8 nt UMI
    # Both samples have a recorded detection source.
    assert all(r.source == "detected" for r in resolutions)


def test_auto_falls_back_to_user_adapter_when_detection_fails(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")

    detector = _detector({"a.fq.gz": _detection(None, rate=0.0)})

    resolutions = _resolve(
        [("a", a)],
        adapter=ILLUMINA_SMALLRNA_ADAPTER,
        detector=detector,
    )
    # Even though detection failed, the supplied adapter happens to match
    # a known family so the kit name is reported.
    assert resolutions[0].kit.kit == "illumina_smallrna"
    assert resolutions[0].kit.adapter == ILLUMINA_SMALLRNA_ADAPTER
    assert resolutions[0].source == "user_adapter"


def test_auto_hard_fails_when_no_detection_and_no_fallback(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")

    detector = _detector({"a.fq.gz": _detection(None, rate=0.0)})

    with pytest.raises(SampleResolutionError) as exc:
        _resolve(
            [("a", a)],
            detector=detector,
        )
    assert "auto-detection found no known kit" in str(exc.value)


# ---------- strict mode ------------------------------------------------------


def test_strict_hard_fails_on_disagreement(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")

    detector = _detector({"a.fq.gz": _detection("illumina_truseq", rate=0.9)})

    with pytest.raises(SampleResolutionError) as exc:
        _resolve(
            [("a", a)],
            adapter=ILLUMINA_SMALLRNA_ADAPTER,
            adapter_detection_mode="strict",
            detector=detector,
        )
    assert "strict" in str(exc.value).lower()
    assert "illumina_truseq" in str(exc.value)


def test_strict_uses_detected_when_no_user_adapter(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")

    detector = _detector({"a.fq.gz": _detection("illumina_smallrna")})

    resolutions = _resolve(
        [("a", a)],
        adapter_detection_mode="strict",
        detector=detector,
    )
    assert resolutions[0].kit.kit == "illumina_smallrna"
    assert resolutions[0].source == "detected"


# ---------- off mode --------------------------------------------------------


def test_off_mode_requires_explicit_adapter(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")

    with pytest.raises(SampleResolutionError) as exc:
        _resolve(
            [("a", a)],
            adapter_detection_mode="off",
        )
    assert "--adapter-detection off requires" in str(exc.value)


def test_off_mode_uses_explicit_adapter_without_scanning(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")

    def explode(_path):
        raise AssertionError("detector must NOT run in off mode")

    resolutions = _resolve(
        [("a", a)],
        adapter=ILLUMINA_SMALLRNA_ADAPTER,
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
        _resolve(
            [("a", a), ("b", b)],
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
    resolutions = _resolve(
        [("a", a), ("b", b)],
        detector=detector,
    )
    assert required_dedup_tools(resolutions) == {"umi_tools"}


def test_write_kit_resolution_tsv_writes_header_and_one_row_per_sample(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")
    detector = _detector({"a.fq.gz": _detection("illumina_smallrna")})
    resolutions = _resolve(
        [("a", a)],
        detector=detector,
    )
    out = tmp_path / "kit_resolution.tsv"
    write_kit_resolution_tsv(resolutions, out)
    raw = out.read_text().splitlines()
    # P1.12: schema-version comment line precedes the column header.
    assert raw[0].startswith("# schema_version:")
    lines = [line for line in raw if not line.startswith("#")]
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

    resolutions = _resolve(
        [("a", a)],
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
        _resolve(
            [("a", a)],
            allow_pretrimmed_inference=False,
            detector=detector,
        )


def test_explicit_pretrimmed_resolves_with_no_adapter(tmp_path) -> None:
    """``--pretrimmed`` resolves to adapter=None so cutadapt skips
    the -a flag entirely; detection is short-circuited."""
    a = tmp_path / "a.fq.gz"
    a.write_text("")

    def explode(_path):
        raise AssertionError("detector must NOT run when --pretrimmed is set")

    resolutions = _resolve(
        [("a", a)],
        pretrimmed=True,
        detector=explode,
    )
    assert resolutions[0].kit.kit == "pretrimmed"
    assert resolutions[0].kit.adapter is None
    assert resolutions[0].source == "explicit_pretrimmed"
    # Dedup falls through to skip (no UMI on the pretrimmed preset).
    assert resolutions[0].dedup_strategy == "skip"


def test_pretrimmed_and_adapter_are_mutually_exclusive(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")
    with pytest.raises(SampleResolutionError) as exc:
        _resolve(
            [("a", a)],
            pretrimmed=True,
            adapter=ILLUMINA_SMALLRNA_ADAPTER,
        )
    assert "mutually exclusive" in str(exc.value)


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

    resolutions = _resolve(
        [("a", a), ("b", b)],
        detector=detector,
    )
    by_sample = {res.sample: res for res in resolutions}
    assert by_sample["a"].kit.kit == "pretrimmed"
    assert by_sample["a"].source == "inferred_pretrimmed"
    assert by_sample["b"].kit.kit == "illumina_smallrna"
    assert by_sample["b"].source == "detected"


def test_user_adapter_overrides_detected_kit(tmp_path) -> None:
    """When detection succeeds AND --adapter is supplied, the user's
    literal sequence wins; the detected kit name is preserved for
    provenance reporting."""
    a = tmp_path / "a.fq.gz"
    a.write_text("")
    # Detector says the adapter is illumina_truseq, but user passes
    # the small-RNA adapter — we honour the user's literal sequence.
    detector = _detector({"a.fq.gz": _detection("illumina_truseq", rate=0.9)})
    resolutions = _resolve(
        [("a", a)],
        adapter=ILLUMINA_SMALLRNA_ADAPTER,
        detector=detector,
    )
    res = resolutions[0]
    # Provenance: kit name follows the detector (so the report is
    # honest about what was seen) but the adapter sequence is the
    # user's.
    assert res.detected_kit == "illumina_truseq"
    assert res.kit.adapter == ILLUMINA_SMALLRNA_ADAPTER
    assert res.source == "user_adapter"


def test_resolution_summary_lines_includes_detected_when_overridden(tmp_path) -> None:
    a = tmp_path / "a.fq.gz"
    a.write_text("")
    # Detector returns truseq but user explicitly supplied small-RNA
    # adapter → user_adapter source, with detected name preserved.
    detector = _detector({"a.fq.gz": _detection("illumina_truseq", rate=0.9)})
    resolutions = _resolve(
        [("a", a)],
        adapter=ILLUMINA_SMALLRNA_ADAPTER,
        detector=detector,
    )
    line = resolution_summary_lines(resolutions)[0]
    assert "kit=illumina_truseq" in line
    assert "user_adapter" in line
