# This script parses the user's CSV, computes enrichment of "mRNA (mt_genome) vs contaminants",
# and plots either a horizontal bar chart or a pie chart (or both) using matplotlib.
#
# You can customize behavior via the PARAMETERS section below or by passing arguments to the
# `summarize_and_plot()` function directly.
#
# Notes (adjustable assumptions):
# - "mRNA" is defined by references that contain the substring in MRNA_REFS (default: "mt_genome.aligned").
# - "Contaminants" are everything else unless CONTAMINANT_REFS is explicitly provided.
# - You can either show each contaminant separately or collapse them into one bucket.
#
# Chart rules (per your environment/tooling requirements):
# - Uses matplotlib (not seaborn)
# - One chart per figure (no subplots)
# - No explicit color settings (matplotlib defaults)
#
# Files are saved to /mnt/data/.

import io
from dataclasses import dataclass, asdict
from typing import Iterable, Optional, Dict, Tuple, List, Union
import pandas as pd
import matplotlib.pyplot as plt
import os

# -------------------- PARAMETERS (edit these) --------------------
CSV_TEXT = """
Control,mt-mRNA,8353378
Control,cyto-tRNA,3255
Control,cyto-rRNA,2782245
Control,mt-rRNA,18971518
Control,mt-tRNA,271935
Oligo,mt-mRNA,3829325
Oligo,cyto-tRNA,6091
Oligo,cyto-rRNA,4792064
Oligo,mt-rRNA,24505010
Oligo,mt-tRNA,238302
"""

SAMPLE: Optional[str] = "Oligo"          # If None, aggregate across all samples (here there's only E1)
MRNA_REFS: Iterable[str] = ("mt-mRNA",)  # Substrings that define "mRNA"
CONTAMINANT_REFS: Optional[Iterable[str]] = None   # If None, all non-mRNA rows are contaminants
COMBINE_CONTAMINANTS: bool = False    # If True, collapse all contaminants into a single "Contaminants" bucket
LABEL_MIN_FRAC: float = 0.001          # Minimum fraction to display data labels (e.g., 0.01 = 1%)
SAVE_DIR: str = "./pie_chart"           # Where to save charts
MAKE_BARH: bool = False                # Produce horizontal bar chart
MAKE_PIE: bool = True                 # Produce pie chart as well
BARH_FILENAME: str = "Oligo_enrichment_barh.svg"
PIE_FILENAME: str = "Oligo_enrichment_pie.svg" 
# ---------------------------------------------------------------


@dataclass
class EnrichmentConfig:
    sample: Optional[str] = SAMPLE
    mrna_refs: Iterable[str] = MRNA_REFS
    contaminant_refs: Optional[Iterable[str]] = CONTAMINANT_REFS
    combine_contaminants: bool = COMBINE_CONTAMINANTS
    label_min_frac: float = LABEL_MIN_FRAC

def _load_csv(csv_source: Union[str, io.StringIO, os.PathLike]) -> pd.DataFrame:
    """
    Load CSV from a string containing newline-separated values, a file path, or a file-like object.
    Expected columns: sample, reference, reads
    """
    if isinstance(csv_source, (str, os.PathLike)) and os.path.exists(str(csv_source)):
        df = pd.read_csv(csv_source, header=None, names=["sample", "reference", "reads"])
    else:
        # Assume it's a raw CSV string content or file-like
        df = pd.read_csv(io.StringIO(str(csv_source)), header=None, names=["sample", "reference", "reads"])
    # Coerce 'reads' to integers
    df["reads"] = pd.to_numeric(df["reads"], errors="coerce").fillna(0).astype(int)
    return df

def _match_any(text: str, patterns: Iterable[str]) -> bool:
    return any(pat in text for pat in patterns)

def compute_enrichment_table(
    df: pd.DataFrame,
    cfg: EnrichmentConfig
) -> Tuple[pd.DataFrame, Dict[str, int], Dict[str, float]]:
    """
    Returns:
      - summary_df: rows = categories (e.g., mRNA, contaminants...), columns: reads, frac
      - counts: dict of absolute counts
      - fracs: dict of fractions (0-1)
    """
    if cfg.sample is not None:
        df = df.loc[df["sample"] == cfg.sample].copy()

    # Determine which references are mRNA
    is_mrna = df["reference"].apply(lambda r: _match_any(str(r), cfg.mrna_refs))
    if cfg.contaminant_refs is not None:
        is_contaminant = df["reference"].apply(lambda r: _match_any(str(r), cfg.contaminant_refs))
    else:
        is_contaminant = ~is_mrna

    mrna_df = df.loc[is_mrna].copy()
    cont_df = df.loc[is_contaminant].copy()

    # Aggregate
    mrna_total = int(mrna_df["reads"].sum())
    cont_by_ref = cont_df.groupby("reference", as_index=False)["reads"].sum().sort_values("reads", ascending=False)
    cont_total = int(cont_by_ref["reads"].sum())
    grand_total = mrna_total + cont_total if (mrna_total + cont_total) > 0 else 1

    # Build summary rows
    rows: List[Tuple[str, int, float]] = []
    rows.append(("mt-mRNA", mrna_total, mrna_total / grand_total))

    if cfg.combine_contaminants:
        rows.append(("Contaminants (all)", cont_total, cont_total / grand_total))
    else:
        for _, rec in cont_by_ref.iterrows():
            rows.append((str(rec["reference"]), int(rec["reads"]), int(rec["reads"]) / grand_total))

    summary_df = pd.DataFrame(rows, columns=["category", "reads", "frac"])
    counts = {row[0]: row[1] for row in rows}
    fracs = {row[0]: row[2] for row in rows}
    return summary_df, counts, fracs

def plot_barh(summary_df: pd.DataFrame, label_min_frac: float, title: str, save_path: Optional[str] = None):
    fig, ax = plt.subplots(figsize=(8, 4 + 0.3 * len(summary_df)))
    # Horizontal bar: x = percent, y = category
    percents = (summary_df["frac"] * 100).values
    ax.barh(summary_df["category"].values, percents)
    ax.set_xlabel("Percent of aligned reads (%)")
    ax.set_ylabel("Category")
    ax.set_title(title)
    # Add labels
    for i, v in enumerate(percents):
        if summary_df.loc[i, "frac"] >= label_min_frac:
            ax.text(v, i, f"{v:.2f}%", va="center", ha="left", fontsize=9)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.show()

def plot_pie(summary_df: pd.DataFrame, label_min_frac: float, title: str, save_path: Optional[str] = None):
    fig, ax = plt.subplots(figsize=(6, 6))
    sizes = summary_df["frac"].values
    # Pie chart with labels only for slices above threshold
    labels = [
        f"{cat} ({frac*100:.1f}%)" if frac >= label_min_frac else ""
        for cat, frac in zip(summary_df["category"].values, sizes)
    ]
    ax.pie(sizes, labels=labels, autopct=lambda p: f"{p:.1f}%" if p/100 >= label_min_frac else "")
    ax.set_title(title)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.show()

def summarize_and_plot(
    csv_source: Union[str, io.StringIO, os.PathLike, pd.DataFrame],
    cfg: Optional[EnrichmentConfig] = None,
    make_barh: bool = True,
    make_pie: bool = False,
    save_dir: Optional[str] = None,
    barh_filename: str = "enrichment_barh.svg",
    pie_filename: str = "enrichment_pie.svg",
    return_summary: bool = True,
):
    # Load
    if isinstance(csv_source, pd.DataFrame):
        df = csv_source.copy()
    else:
        df = _load_csv(csv_source)

    # Default config
    cfg = cfg or EnrichmentConfig()

    # Compute
    summary_df, counts, fracs = compute_enrichment_table(df, cfg)

    # Prepare titles and save paths
    sample_label = cfg.sample if cfg.sample is not None else "ALL"
    title = f"{sample_label}"
    barh_path = os.path.join(save_dir, barh_filename) if save_dir else None
    pie_path = os.path.join(save_dir, pie_filename) if save_dir else None

    # Plot
    if make_barh:
        plot_barh(summary_df, cfg.label_min_frac, title, save_path=barh_path)
    if make_pie:
        plot_pie(summary_df, cfg.label_min_frac, title, save_path=pie_path)

    # Compute enrichment summary
    mrna_key = next((k for k in counts.keys() if k.startswith("mRNA")), None)
    mrna_reads = counts.get(mrna_key, 0)
    contam_reads = sum(v for k, v in counts.items() if k != mrna_key)
    total = mrna_reads + contam_reads if (mrna_reads + contam_reads) > 0 else 1
    enrichment_pct = 100.0 * mrna_reads / total
    ratio = (mrna_reads / contam_reads) if contam_reads > 0 else float('inf')

    # Show a compact table to the user
    from caas_jupyter_tools import display_dataframe_to_user
    display_dataframe_to_user(
        name=f"Read composition (Sample {sample_label})",
        dataframe=summary_df.assign(percent=(summary_df["frac"]*100).round(2)).drop(columns=["frac"])
    )

    print("CONFIG:", asdict(cfg))
    print(f"Sample: {sample_label}")
    print(f"mRNA reads: {mrna_reads:,}")
    print(f"Contaminant reads: {contam_reads:,}")
    print(f"mRNA % of aligned reads: {enrichment_pct:.2f}%")
    print(f"mRNA : contaminants ratio = {ratio:.3f}:1" if ratio != float('inf') else "All reads are mRNA (no contaminants).")

    # Return for programmatic use
    if return_summary:
        return {
            "summary_df": summary_df,
            "counts": counts,
            "fractions": fracs,
            "mrna_percent": enrichment_pct,
            "ratio_mrna_to_contaminants": ratio,
            "barh_path": barh_path,
            "pie_path": pie_path,
        }

def main() -> dict:
    """Run read-composition summary with current module-level defaults."""
    os.makedirs(SAVE_DIR, exist_ok=True)

    cfg = EnrichmentConfig(
        sample=SAMPLE,
        mrna_refs=MRNA_REFS,
        contaminant_refs=CONTAMINANT_REFS,
        combine_contaminants=COMBINE_CONTAMINANTS,
        label_min_frac=LABEL_MIN_FRAC,
    )

    return summarize_and_plot(
        CSV_TEXT,
        cfg=cfg,
        make_barh=MAKE_BARH,
        make_pie=MAKE_PIE,
        save_dir=SAVE_DIR,
        barh_filename=BARH_FILENAME,
        pie_filename=PIE_FILENAME,
    )


if __name__ == "__main__":
    main()
