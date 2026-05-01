"""Unit tests for the unified sample sheet parser."""

from __future__ import annotations

from pathlib import Path

import pytest

from mitoribopy.sample_sheet import (
    SampleSheet,
    SampleSheetError,
    load_sample_sheet,
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_minimum_required_columns_round_trip(tmp_path: Path) -> None:
    sheet = _write(
        tmp_path / "samples.tsv",
        "sample_id\tassay\tcondition\tfastq_1\n"
        "WT_Ribo_1\tribo\tWT\tfq/WT_Ribo_1.fq.gz\n"
        "KO_Ribo_1\tribo\tKO\tfq/KO_Ribo_1.fq.gz\n",
    )
    s = load_sample_sheet(sheet)
    assert isinstance(s, SampleSheet)
    assert len(s.rows) == 2
    row = s.rows[0]
    assert row.sample_id == "WT_Ribo_1"
    assert row.assay == "ribo"
    assert row.condition == "WT"
    assert row.fastq_1 == Path("fq/WT_Ribo_1.fq.gz")
    # Optional fields default to None / False:
    assert row.fastq_2 is None
    assert row.umi_length is None
    assert row.exclude is False
    assert s.condition_map() == {"WT_Ribo_1": "WT", "KO_Ribo_1": "KO"}


def test_paired_end_and_per_sample_overrides(tmp_path: Path) -> None:
    sheet = _write(
        tmp_path / "samples.tsv",
        "sample_id\tassay\tcondition\treplicate\tfastq_1\tfastq_2\t"
        "kit_preset\tumi_length\tumi_position\tstrandedness\n"
        "WT_RNA_1\trna\tWT\t1\trna/WT_R1.fq.gz\trna/WT_R2.fq.gz\t"
        "pretrimmed\t0\t5p\tforward\n"
        "WT_Ribo_1\tribo\tWT\t1\tribo/WT_Ribo_1.fq.gz\t\t"
        "illumina_truseq_umi\t8\t5p\tforward\n",
    )
    s = load_sample_sheet(sheet)
    rna = s.by_assay("rna")
    ribo = s.by_assay("ribo")
    assert len(rna) == 1 and len(ribo) == 1
    assert rna[0].fastq_2 == Path("rna/WT_R2.fq.gz")
    assert rna[0].umi_length == 0
    assert ribo[0].umi_length == 8
    assert ribo[0].umi_position == "5p"
    assert ribo[0].strandedness == "forward"


def test_exclude_flag_drops_row_from_active_view(tmp_path: Path) -> None:
    sheet = _write(
        tmp_path / "samples.tsv",
        "sample_id\tassay\tcondition\tfastq_1\texclude\n"
        "WT_Ribo_1\tribo\tWT\tfq/A.fq.gz\tfalse\n"
        "WT_Ribo_2\tribo\tWT\tfq/B.fq.gz\ttrue\n"
        "KO_Ribo_1\tribo\tKO\tfq/C.fq.gz\t\n",
    )
    s = load_sample_sheet(sheet)
    assert len(s.rows) == 3  # raw row count unchanged
    active_ids = [row.sample_id for row in s.active()]
    assert active_ids == ["WT_Ribo_1", "KO_Ribo_1"]  # B excluded
    # Excluded rows do not appear in the condition map either:
    assert "WT_Ribo_2" not in s.condition_map()


def test_blank_lines_and_comments_are_ignored(tmp_path: Path) -> None:
    sheet = _write(
        tmp_path / "samples.tsv",
        "# project XYZ\n"
        "\n"
        "sample_id\tassay\tcondition\tfastq_1\n"
        "# the following are pilot samples\n"
        "WT_Ribo_1\tribo\tWT\tA.fq.gz\n"
        "\n"
        "KO_Ribo_1\tribo\tKO\tB.fq.gz\n",
    )
    s = load_sample_sheet(sheet)
    assert len(s.rows) == 2


def test_null_tokens_become_none(tmp_path: Path) -> None:
    """Empty / NA / None / - / null all read as None."""
    sheet = _write(
        tmp_path / "samples.tsv",
        "sample_id\tassay\tcondition\tfastq_1\tadapter\tumi_length\tnotes\n"
        "A\tribo\tWT\tA.fq.gz\tNA\t\tnull\n"
        "B\tribo\tWT\tB.fq.gz\t-\tNone\t-\n",
    )
    s = load_sample_sheet(sheet)
    for row in s.rows:
        assert row.adapter is None
        assert row.umi_length is None
        assert row.notes is None


def test_missing_required_column_lists_all(tmp_path: Path) -> None:
    sheet = _write(
        tmp_path / "samples.tsv",
        "sample_id\tassay\tfastq_1\n"
        "WT_Ribo_1\tribo\tA.fq.gz\n",
    )
    with pytest.raises(SampleSheetError, match="missing required column"):
        load_sample_sheet(sheet)


def test_unknown_column_is_rejected(tmp_path: Path) -> None:
    sheet = _write(
        tmp_path / "samples.tsv",
        "sample_id\tassay\tcondition\tfastq_1\tmystery\n"
        "A\tribo\tWT\tA.fq.gz\tboo\n",
    )
    with pytest.raises(SampleSheetError, match="unknown column"):
        load_sample_sheet(sheet)


def test_duplicate_sample_id_is_an_error(tmp_path: Path) -> None:
    sheet = _write(
        tmp_path / "samples.tsv",
        "sample_id\tassay\tcondition\tfastq_1\n"
        "A\tribo\tWT\tA.fq.gz\n"
        "B\tribo\tWT\tB.fq.gz\n"
        "A\trna\tKO\tA_rna.fq.gz\n",
    )
    with pytest.raises(SampleSheetError, match="duplicate sample_id 'A'"):
        load_sample_sheet(sheet)


def test_unknown_assay_is_an_error(tmp_path: Path) -> None:
    sheet = _write(
        tmp_path / "samples.tsv",
        "sample_id\tassay\tcondition\tfastq_1\n"
        "A\tprotein\tWT\tA.fq.gz\n",
    )
    with pytest.raises(SampleSheetError, match="assay 'protein' must be one of"):
        load_sample_sheet(sheet)


def test_malformed_umi_length_is_an_error(tmp_path: Path) -> None:
    sheet = _write(
        tmp_path / "samples.tsv",
        "sample_id\tassay\tcondition\tfastq_1\tumi_length\n"
        "A\tribo\tWT\tA.fq.gz\teight\n",
    )
    with pytest.raises(SampleSheetError, match="umi_length 'eight' is not an integer"):
        load_sample_sheet(sheet)


def test_negative_umi_length_is_an_error(tmp_path: Path) -> None:
    sheet = _write(
        tmp_path / "samples.tsv",
        "sample_id\tassay\tcondition\tfastq_1\tumi_length\n"
        "A\tribo\tWT\tA.fq.gz\t-5\n",
    )
    with pytest.raises(SampleSheetError, match="umi_length must be >= 0"):
        load_sample_sheet(sheet)


def test_invalid_strandedness_is_an_error(tmp_path: Path) -> None:
    sheet = _write(
        tmp_path / "samples.tsv",
        "sample_id\tassay\tcondition\tfastq_1\tstrandedness\n"
        "A\tribo\tWT\tA.fq.gz\tweird\n",
    )
    with pytest.raises(SampleSheetError, match="strandedness 'weird'"):
        load_sample_sheet(sheet)


def test_invalid_dedup_strategy_is_an_error(tmp_path: Path) -> None:
    sheet = _write(
        tmp_path / "samples.tsv",
        "sample_id\tassay\tcondition\tfastq_1\tdedup_strategy\n"
        "A\tribo\tWT\tA.fq.gz\tmark-duplicates\n",
    )
    with pytest.raises(SampleSheetError, match="dedup_strategy 'mark-duplicates'"):
        load_sample_sheet(sheet)


def test_all_errors_reported_in_one_pass(tmp_path: Path) -> None:
    """Multiple errors across rows should all be reported at once so the
    user does not have to fix-and-retry one row at a time."""
    sheet = _write(
        tmp_path / "samples.tsv",
        "sample_id\tassay\tcondition\tfastq_1\tumi_length\n"
        "A\tribo\tWT\tA.fq.gz\teight\n"        # bad umi_length
        "A\trna\tWT\tA_rna.fq.gz\t8\n"          # duplicate sample_id
        "B\tprotein\tKO\tB.fq.gz\t8\n",         # bad assay
    )
    with pytest.raises(SampleSheetError) as exc:
        load_sample_sheet(sheet)
    msg = str(exc.value)
    assert "umi_length 'eight'" in msg
    assert "duplicate sample_id 'A'" in msg
    assert "assay 'protein'" in msg


def test_empty_file_is_an_error(tmp_path: Path) -> None:
    sheet = _write(tmp_path / "samples.tsv", "")
    with pytest.raises(SampleSheetError, match="file is empty"):
        load_sample_sheet(sheet)


def test_header_only_file_is_an_error(tmp_path: Path) -> None:
    sheet = _write(
        tmp_path / "samples.tsv",
        "sample_id\tassay\tcondition\tfastq_1\n",
    )
    with pytest.raises(SampleSheetError, match="no data rows"):
        load_sample_sheet(sheet)


def test_fastq_paths_helper_yields_r1_then_r2(tmp_path: Path) -> None:
    sheet = _write(
        tmp_path / "samples.tsv",
        "sample_id\tassay\tcondition\tfastq_1\tfastq_2\n"
        "A\trna\tWT\tA_R1.fq.gz\tA_R2.fq.gz\n"
        "B\trna\tKO\tB.fq.gz\t\n"
        "C\tribo\tWT\tC.fq.gz\t\n",
    )
    s = load_sample_sheet(sheet)
    assert s.fastq_paths("rna") == [
        Path("A_R1.fq.gz"),
        Path("A_R2.fq.gz"),
        Path("B.fq.gz"),
    ]
    assert s.fastq_paths("ribo") == [Path("C.fq.gz")]
