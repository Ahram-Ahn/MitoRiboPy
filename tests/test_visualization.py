from __future__ import annotations

import pandas as pd

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
