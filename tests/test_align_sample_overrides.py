"""Tests for per-sample adapter / UMI overrides in ``mitoribopy align``.

Covers the full path from YAML ``align.samples:`` -> sidecar TSV ->
``--sample-overrides`` -> resolver. The resolver itself is also exercised
directly through :func:`resolve_sample_resolutions` so the override
precedence rules (per-sample beats global) are tested in isolation.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mitoribopy.align._types import SampleOverride
from mitoribopy.align.adapter_detect import DetectionResult
from mitoribopy.align.sample_resolve import (
    SampleResolutionError,
    read_sample_overrides_tsv,
    resolve_sample_resolutions,
    write_sample_overrides_tsv,
)
from mitoribopy.cli import all_ as all_cli


ILLUMINA_SMALLRNA_ADAPTER = "TGGAATTCTCGGGTGCCAAGG"
ILLUMINA_TRUSEQ_ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
QIASEQ_ADAPTER = "AACTGTAGGCACCATCAAT"


def _detection(
    name: str | None,
    rate: float = 0.95,
    pretrimmed: bool = False,
) -> DetectionResult:
    return DetectionResult(
        preset_name=name,
        match_rate=rate,
        per_kit_rates={"illumina_smallrna": rate, "illumina_truseq": 0.0},
        n_reads_scanned=5000,
        ambiguous=False,
        pretrimmed=pretrimmed,
    )


def _detector(map_):
    def detector(path):
        return map_[Path(path).name]
    return detector


def _resolve(samples, **kwargs):
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


# ---------- TSV I/O round-trip ----------------------------------------------


def test_sample_overrides_tsv_round_trips(tmp_path) -> None:
    overrides = [
        SampleOverride(
            sample="sampleA",
            adapter=ILLUMINA_TRUSEQ_ADAPTER,
            umi_length=8,
            umi_position="5p",
        ),
        SampleOverride(
            sample="sampleB",
            adapter=QIASEQ_ADAPTER,
            umi_length=12,
            umi_position="3p",
        ),
        SampleOverride(sample="sampleC", pretrimmed=True, umi_length=0),
    ]
    tsv_path = tmp_path / "sample_overrides.tsv"
    write_sample_overrides_tsv(overrides, tsv_path)

    loaded = read_sample_overrides_tsv(tsv_path)
    assert set(loaded) == {"sampleA", "sampleB", "sampleC"}
    assert loaded["sampleA"].umi_length == 8
    assert loaded["sampleA"].umi_position == "5p"
    assert loaded["sampleA"].adapter == ILLUMINA_TRUSEQ_ADAPTER
    assert loaded["sampleB"].umi_position == "3p"
    assert loaded["sampleC"].umi_length == 0
    assert loaded["sampleC"].umi_position is None
    assert loaded["sampleC"].pretrimmed is True


def test_sample_overrides_tsv_handles_blank_cells(tmp_path) -> None:
    """Empty / NA / None cells fall through to global defaults (None)."""
    tsv_path = tmp_path / "partial.tsv"
    tsv_path.write_text(
        "sample\tadapter\tpretrimmed\tumi_length\tumi_position\tdedup_strategy\n"
        f"sampleA\t{ILLUMINA_TRUSEQ_ADAPTER}\t\t8\t5p\t\n"
        "sampleB\t\t\tNA\t\tnone\n"
    )
    loaded = read_sample_overrides_tsv(tsv_path)
    a = loaded["sampleA"]
    assert a.adapter == ILLUMINA_TRUSEQ_ADAPTER
    assert a.umi_length == 8
    assert a.dedup_strategy is None

    b = loaded["sampleB"]
    assert b.adapter is None
    assert b.umi_length is None
    assert b.dedup_strategy is None
    assert b.pretrimmed is None


def test_sample_overrides_tsv_rejects_legacy_kit_preset_column(tmp_path) -> None:
    """v0.7.1: the legacy ``kit_preset`` column raises a clear error."""
    tsv_path = tmp_path / "legacy.tsv"
    tsv_path.write_text(
        "sample\tkit_preset\tumi_length\nsampleA\tillumina_truseq_umi\t8\n"
    )
    with pytest.raises(SampleResolutionError, match="kit_preset"):
        read_sample_overrides_tsv(tsv_path)


def test_sample_overrides_tsv_rejects_missing_sample_column(tmp_path) -> None:
    tsv_path = tmp_path / "bad.tsv"
    tsv_path.write_text("adapter\tumi_length\nfoo\t8\n")
    with pytest.raises(SampleResolutionError, match="sample"):
        read_sample_overrides_tsv(tsv_path)


def test_sample_overrides_tsv_rejects_duplicate_sample(tmp_path) -> None:
    tsv_path = tmp_path / "dup.tsv"
    tsv_path.write_text(
        "sample\tumi_length\nsampleA\t8\nsampleA\t12\n"
    )
    with pytest.raises(SampleResolutionError, match="duplicate"):
        read_sample_overrides_tsv(tsv_path)


def test_sample_overrides_tsv_rejects_bad_umi_length(tmp_path) -> None:
    tsv_path = tmp_path / "badlen.tsv"
    tsv_path.write_text("sample\tumi_length\nsampleA\teight\n")
    with pytest.raises(SampleResolutionError, match="non-integer"):
        read_sample_overrides_tsv(tsv_path)


def test_sample_overrides_tsv_rejects_bad_umi_position(tmp_path) -> None:
    tsv_path = tmp_path / "badpos.tsv"
    tsv_path.write_text("sample\tumi_position\nsampleA\tmiddle\n")
    with pytest.raises(SampleResolutionError, match="umi_position"):
        read_sample_overrides_tsv(tsv_path)


# ---------- resolver: per-sample override beats global ----------------------


def test_per_sample_override_beats_global_kit(tmp_path) -> None:
    """Sample with override picks its kit; sample without falls back to global."""
    a = tmp_path / "sampleA.fq.gz"
    b = tmp_path / "sampleB.fq.gz"
    a.write_text("")
    b.write_text("")

    detector = _detector(
        {
            "sampleA.fq.gz": _detection("illumina_smallrna"),
            "sampleB.fq.gz": _detection("illumina_smallrna"),
        }
    )
    overrides = {
        "sampleA": SampleOverride(
            sample="sampleA",
            adapter=ILLUMINA_TRUSEQ_ADAPTER,
            umi_length=8,
            umi_position="5p",
        )
    }

    resolutions = _resolve(
        [("sampleA", a), ("sampleB", b)],
        detector=detector,
        sample_overrides=overrides,
    )

    by_name = {r.sample: r for r in resolutions}
    # sampleA's override pins the truseq adapter + UMI; the kit
    # name follows the family lookup.
    assert by_name["sampleA"].kit.adapter == ILLUMINA_TRUSEQ_ADAPTER
    assert by_name["sampleA"].kit.umi_length == 8
    assert by_name["sampleA"].kit.umi_position == "5p"
    assert by_name["sampleA"].dedup_strategy == "umi-tools"
    assert by_name["sampleA"].source.startswith("per_sample_override:")

    # sampleB has no override -> falls back to detected (auto) kit.
    assert by_name["sampleB"].kit.kit == "illumina_smallrna"
    assert by_name["sampleB"].kit.umi_length == 0
    assert by_name["sampleB"].dedup_strategy == "skip"
    assert "per_sample_override" not in by_name["sampleB"].source


def test_mixed_umi_and_no_umi_batch_resolves_correctly(tmp_path) -> None:
    """A real mixed batch: one UMI sample, one already-trimmed SRA-style sample."""
    umi_fq = tmp_path / "sampleU.fq.gz"
    no_fq = tmp_path / "sampleN.fq.gz"
    umi_fq.write_text("")
    no_fq.write_text("")

    detector = _detector(
        {
            "sampleU.fq.gz": _detection("illumina_truseq"),  # adapter detected
            "sampleN.fq.gz": _detection(None, rate=0.0, pretrimmed=True),
        }
    )
    overrides = {
        "sampleU": SampleOverride(
            sample="sampleU",
            adapter=ILLUMINA_TRUSEQ_ADAPTER,
            umi_length=8,
            umi_position="5p",
        ),
        "sampleN": SampleOverride(
            sample="sampleN", pretrimmed=True, umi_length=0
        ),
    }

    resolutions = _resolve(
        [("sampleU", umi_fq), ("sampleN", no_fq)],
        detector=detector,
        sample_overrides=overrides,
    )

    by_name = {r.sample: r for r in resolutions}
    # UMI sample -> umi_tools dedup.
    assert by_name["sampleU"].kit.umi_length == 8
    assert by_name["sampleU"].dedup_strategy == "umi-tools"
    # No-UMI / pretrimmed sample -> skip dedup, no -a flag at trim time.
    assert by_name["sampleN"].kit.kit == "pretrimmed"
    assert by_name["sampleN"].kit.adapter is None
    assert by_name["sampleN"].kit.umi_length == 0
    assert by_name["sampleN"].dedup_strategy == "skip"


def test_unknown_sample_in_overrides_raises(tmp_path) -> None:
    a = tmp_path / "sampleA.fq.gz"
    a.write_text("")

    overrides = {
        "sampleZZZ": SampleOverride(sample="sampleZZZ", umi_length=8)
    }
    with pytest.raises(SampleResolutionError, match="unknown sample"):
        _resolve(
            [("sampleA", a)],
            adapter=ILLUMINA_SMALLRNA_ADAPTER,
            adapter_detection_mode="off",
            sample_overrides=overrides,
        )


# ---------- mitoribopy all -> sidecar TSV materialisation -------------------


def test_normalize_align_inputs_materialises_samples_block(tmp_path) -> None:
    cfg = {
        "fastq": "input_data/",
        "samples": [
            {
                "name": "sampleA",
                "adapter": ILLUMINA_TRUSEQ_ADAPTER,
                "umi_length": 8,
                "umi_position": "5p",
            },
            {
                "name": "sampleB",
                "pretrimmed": True,
                "umi_length": 0,
            },
        ],
    }
    out = all_cli._normalize_align_inputs(cfg, run_root=tmp_path)
    # samples key consumed; sample_overrides points at the materialised TSV.
    assert "samples" not in out
    overrides_path = Path(out["sample_overrides"])
    assert overrides_path.is_file()

    loaded = read_sample_overrides_tsv(overrides_path)
    assert loaded["sampleA"].umi_length == 8
    assert loaded["sampleA"].umi_position == "5p"
    assert loaded["sampleA"].adapter == ILLUMINA_TRUSEQ_ADAPTER
    assert loaded["sampleB"].pretrimmed is True


def test_normalize_align_inputs_dry_run_drops_samples(tmp_path) -> None:
    """Dry-run path: no run_root -> no TSV is written, key is dropped."""
    cfg = {
        "fastq": "input_data/",
        "samples": [{"name": "sampleA", "umi_length": 8}],
    }
    out = all_cli._normalize_align_inputs(cfg, run_root=None)
    assert out.get("sample_overrides") is None
    assert "samples" not in out


def test_samples_block_requires_name(tmp_path) -> None:
    cfg = {"samples": [{"umi_length": 8}]}  # missing name
    with pytest.raises(ValueError, match="name"):
        all_cli._normalize_align_inputs(cfg, run_root=tmp_path)


def test_samples_block_rejects_legacy_kit_preset(tmp_path) -> None:
    cfg = {
        "samples": [
            {"name": "sampleA", "kit_preset": "illumina_truseq_umi"},
        ],
    }
    with pytest.raises(ValueError, match="kit_preset"):
        all_cli._normalize_align_inputs(cfg, run_root=tmp_path)
