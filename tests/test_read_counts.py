from __future__ import annotations

from pathlib import Path

from mitoribopy.io.read_counts import compute_total_counts


def test_compute_total_counts_accepts_flexible_column_names(tmp_path: Path) -> None:
    count_file = tmp_path / "counts.csv"
    count_file.write_text(
        "Sample,Reference,Counts\nS1,mt_genome.aligned,10\nS1,other,5\n",
        encoding="utf-8",
    )

    total_counts_map, total_counts_df = compute_total_counts(
        str(count_file),
        normalization_mode="mt_mrna",
    )

    assert total_counts_map["S1"] == 10
    assert total_counts_df["Total_reads"].tolist() == [10]


def test_compute_total_counts_falls_back_to_column_order_for_txt(tmp_path: Path) -> None:
    count_file = tmp_path / "counts.txt"
    count_file.write_text(
        "col_a\tcol_b\tcol_c\nalpha1\tmt_genome.aligned\t7\nalpha1\tother\t3\n",
        encoding="utf-8",
    )

    total_counts_map, total_counts_df = compute_total_counts(
        str(count_file),
        normalization_mode="mt_mrna",
    )

    assert total_counts_map["alpha1"] == 7
    assert total_counts_df["Total_reads"].tolist() == [7]
