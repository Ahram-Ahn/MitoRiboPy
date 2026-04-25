from __future__ import annotations

from types import SimpleNamespace

import pandas as pd

from mitoribopy.plotting.coverage_profile_plots import run_coverage_profile_plots
from mitoribopy.plotting.visualization import _selected_offset_guides


def test_selected_offset_guides_use_plot_axis_positions_but_keep_true_offset_labels() -> None:
    row = pd.Series(
        {
            "Read Length": 30,
            "Most Enriched 5' Offset": 13,
            "Most Enriched 3' Offset": 17,
        }
    )

    guides = _selected_offset_guides(
        row,
        {-13, 17},
    )

    assert [guide["axis_offset"] for guide in guides] == [-13, 17]
    assert [guide["label"] for guide in guides] == [
        "Selected 5' offset (13 nt)",
        "Selected 3' offset (17 nt)",
    ]


def test_selected_offset_guides_skip_offsets_outside_plotted_axis_range() -> None:
    row = pd.Series(
        {
            "Read Length": 31,
            "Most Enriched 5' Offset": 13,
            "Most Enriched 3' Offset": 17,
        }
    )

    guides = _selected_offset_guides(
        row,
        {17},
    )

    assert [guide["axis_offset"] for guide in guides] == [17]


def test_run_coverage_profile_plots_writes_codon_binned_outputs(tmp_path) -> None:
    fasta_path = tmp_path / "tiny.fa"
    fasta_path.write_text(">ND1\n" + ("A" * 120) + "\n", encoding="utf-8")

    annotation_df = pd.DataFrame(
        [
            {
                "transcript": "ND1",
                "sequence_name": "ND1",
                "l_tr": 120,
                "l_utr5": 9,
                "l_utr3": 9,
            }
        ]
    )
    filtered_bed_df = pd.DataFrame(
        [
            {
                "chrom": "ND1",
                "start": 9,
                "end": 39,
                "read_length": 30,
                "sample_name": "sample1",
            },
            {
                "chrom": "ND1",
                "start": 12,
                "end": 42,
                "read_length": 30,
                "sample_name": "sample1",
            },
        ]
    )
    selected_offsets_df = pd.DataFrame(
        [{"Read Length": 30, "Most Enriched 5' Offset": 12, "Most Enriched 3' Offset": 18}]
    )
    args = SimpleNamespace(
        plot_format="svg",
        offset_site="p",
        total_mrna_map={"sample1": 100},
        cap_percentile=0.999,
    )

    run_coverage_profile_plots(
        sample_dirs=[str(tmp_path / "sample1")],
        selected_offsets_df=selected_offsets_df,
        offset_type="5",
        fasta_file=fasta_path,
        output_dir=tmp_path / "coverage",
        args=args,
        annotation_df=annotation_df,
        filtered_bed_df=filtered_bed_df,
    )

    # Read-coverage (site-independent) sits at the top of output_dir
    # since this single-call test does not pass `read_coverage_dir`.
    assert (tmp_path / "coverage" / "read_coverage_rpm" / "ND1_read_coverage_(rpm).svg").is_file()
    assert (
        tmp_path
        / "coverage"
        / "read_coverage_rpm_codon"
        / "ND1_read_coverage_(rpm,_codon-binned_cds).svg"
    ).is_file()
    # Frame-coloured P-site plots, renamed from p_site_coverage_*_frame
    # to site_density_*_frame in v0.4.x; the parent dir names the site.
    assert (
        tmp_path
        / "coverage"
        / "site_density_rpm_frame"
        / "ND1_p-site_density_(rpm,_cds_frame_coloring).svg"
    ).is_file()
    assert (
        tmp_path
        / "coverage"
        / "site_density_raw_frame"
        / "ND1_p-site_density_(raw_counts,_cds_frame_coloring).svg"
    ).is_file()


def test_cds_nt_slice_map_excludes_utr_and_aligns_on_start_codon() -> None:
    """Sanity: the frame-colored slice must start at CDS start, exclude UTRs,
    and be trimmed to a multiple of 3."""
    import numpy as np
    from mitoribopy.plotting.coverage_profile_plots import _build_cds_nt_slice_map

    # Transcript length 30: 5' UTR 3 nt, CDS 24 nt (8 codons), 3' UTR 3 nt.
    arr = np.arange(30, dtype=float)
    cov_map = {"sample1": {"TXA": arr}}
    cds_lookup = {"TXA": (3, 27)}  # CDS covers indices 3..26 (end-exclusive 27)

    out = _build_cds_nt_slice_map(cov_map, cds_lookup)
    slice_arr, cds_start = out["sample1"]["TXA"]

    assert cds_start == 3
    # Exactly 24 nt, all from the CDS, no UTR leakage.
    assert len(slice_arr) == 24
    assert slice_arr[0] == 3.0
    assert slice_arr[-1] == 26.0
    # Length is a multiple of 3 (usable codons).
    assert len(slice_arr) % 3 == 0
