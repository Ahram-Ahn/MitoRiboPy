"""Plot writers for translation-profile analysis outputs."""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from .style import apply_publication_style

apply_publication_style()


def plot_codon_usage_dataframe(
    df_rows: pd.DataFrame,
    output_path: str,
    title: str,
    codon_label_order: list[str],
) -> None:
    """Render a codon-usage bar plot."""
    plot_df = df_rows.copy()
    plot_df["Label"] = plot_df["Codon"] + "-" + plot_df["AA"]

    plt.figure(figsize=(20, 6))
    sns.barplot(
        x="Label",
        y="CoverageDivFreq",
        hue="Category",
        data=plot_df,
        dodge=False,
        order=codon_label_order,
    )
    plt.xticks(rotation=90)
    y_max = plot_df["CoverageDivFreq"].max() * 1.2 if not plot_df.empty else 1
    plt.ylim(0, y_max)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_site_depth_profile(
    footprint_df: pd.DataFrame,
    output_path: str,
    *,
    site_column: str,
    site_label: str,
    sample_name: str,
    transcript_name: str,
    color: str,
) -> None:
    """Render a selected-site depth profile for a transcript."""
    plt.figure(figsize=(12, 3))
    plt.bar(
        footprint_df["Position"],
        footprint_df[site_column],
        color=color,
        width=3.0,
    )
    plt.title(f"{transcript_name} - {site_label} ({sample_name})")
    plt.xlabel("Position")
    plt.ylabel(f"{site_label} Depth")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_frame_usage_total(
    frame_df: pd.DataFrame,
    output_path: str,
    *,
    primary_site_label: str,
    sample_name: str,
) -> None:
    """Render the per-sample total frame-usage plot."""
    plt.figure(figsize=(6, 6))
    sns.barplot(
        x="Frame",
        y="Proportion",
        data=frame_df,
        hue="Frame",
        dodge=False,
        palette="Set2",
        legend=False,
    )
    plt.title(f"Total Frame Usage Proportion ({primary_site_label}) - {sample_name}")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_frame_usage_by_transcript(
    frame_df: pd.DataFrame,
    output_path: str,
    *,
    primary_site_label: str,
    sample_name: str,
) -> None:
    """Render the per-transcript frame-usage plot."""
    if frame_df.empty:
        return

    plt.figure(figsize=(12, 6))
    sns.barplot(
        x="Transcript",
        y="Proportion",
        hue="Frame",
        data=frame_df,
        palette="Set3",
    )
    plt.title(f"Frame Usage Proportion by Transcript ({primary_site_label}) - {sample_name}")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
