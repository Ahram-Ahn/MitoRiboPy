"""Reference transcript annotation tables for yeast and human mitochondrial coding regions."""

from __future__ import annotations

import pandas as pd


yeast_annotation_df = pd.DataFrame(
    {
        "transcript": ["COX1", "ATP8", "ATP6", "COB", "ATP9", "VAR1", "COX2", "COX3"],
        "l_tr": [2147, 898, 1167, 2223, 960, 1952, 886, 1528],
        "l_utr5": [460, 364, 287, 954, 630, 162, 54, 604],
        "l_utr3": [82, 387, 100, 111, 99, 593, 76, 114],
    }
)
yeast_annotation_df["l_cds"] = (
    yeast_annotation_df["l_tr"] - yeast_annotation_df["l_utr5"] - yeast_annotation_df["l_utr3"]
)
yeast_annotation_df = yeast_annotation_df[["transcript", "l_tr", "l_utr5", "l_cds", "l_utr3"]]


human_annotation_df = pd.DataFrame(
    {
        "transcript": [
            "ND1",
            "ND2",
            "COX1",
            "COX2",
            "ATP86",
            "COX3",
            "ND3",
            "ND4L4",
            "ND5",
            "ND6",
            "CYTB",
        ],
        "l_tr": [958, 1042, 1617, 708, 843, 784, 346, 1668, 2379, 538, 1141],
        "l_utr5": [2, 0, 3, 0, 162, 0, 0, 290, 0, 0, 0],
        "l_utr3": [0, 0, 72, 24, 0, 0, 0, 0, 567, 13, 0],
    }
)
human_annotation_df["l_cds"] = (
    human_annotation_df["l_tr"] - human_annotation_df["l_utr5"] - human_annotation_df["l_utr3"]
)
human_annotation_df = human_annotation_df[["transcript", "l_tr", "l_utr5", "l_cds", "l_utr3"]]

