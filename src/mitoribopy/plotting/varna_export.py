# varna_plot.py
"""
Generates VARNA color files from footprint density CSVs, storing them in a separate directory.
Additionally includes both 'VARNA_log2' and 'VARNA_Scaled' columns in the final CSV.
"""

import os
import numpy as np
import pandas as pd
from math import log2

def run_varna_plot(
    inframe_out_dir,
    varna_norm_perc=0.98,
    varna_out_dir="analysis_results/varna",
    varna_site_col="A_site",
):
    """
    Searches `inframe_out_dir` for '*_footprint_density.csv'.
    For each CSV:
      - read varna_site_col,
      - coverage=0 => 0, coverage>0 => log2(coverage+1),
      - store that in a 'VARNA_log2' column,
      - cap coverage at 'varna_norm_perc' among nonzero => scale => 'VARNA_Scaled',
      - produce 2 new columns => 'VARNA_log2' + 'VARNA_Scaled'.
    Write output to `varna_out_dir`, creating:
      {transcript}_varna_color.txt (one scaled value per line)
      {transcript}_varna_scaled.csv (original + 'VARNA_log2','VARNA_Scaled' columns)
    """

    print(
        f"[VARNA] Searching in '{inframe_out_dir}' for *_footprint_density.csv files "
        f"using column '{varna_site_col}'."
    )
    os.makedirs(varna_out_dir, exist_ok=True)

    # Gather all footprint_density CSV files
    footprint_csvs = []
    for root, dirs, files in os.walk(inframe_out_dir):
        for f in files:
            if f.endswith("_footprint_density.csv"):
                footprint_csvs.append(os.path.join(root, f))

    if not footprint_csvs:
        print("[VARNA] No footprint_density CSV found. Skipping varna plot.")
        return

    print(f"[VARNA] Found {len(footprint_csvs)} files. Using varna_norm_perc={varna_norm_perc}")

    for csv_file in footprint_csvs:
        base_name = os.path.basename(csv_file)  # e.g. "COX1_footprint_density.csv"
        transcript = base_name.replace("_footprint_density.csv","")
        # Expected path: <inframe_out_dir>/<sample>/footprint_density/<file>
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(csv_file)))
        sample_varna_dir = os.path.join(varna_out_dir, sample_name)
        os.makedirs(sample_varna_dir, exist_ok=True)

        try:
            df = pd.read_csv(csv_file)
        except Exception as e:
            print(f"[VARNA] Error reading {csv_file}: {e}")
            continue

        if varna_site_col not in df.columns:
            print(f"[VARNA] '{varna_site_col}' not found in {csv_file}. Skipping.")
            continue

        coverage = df[varna_site_col].values
        log2_vals = []
        for c_ in coverage:
            if c_ <= 0:
                log2_vals.append(0.0)
            else:
                log2_vals.append(log2(c_+1))

        df["VARNA_log2"] = log2_vals

        arr = np.array(log2_vals)
        nonzero = arr[arr>0]
        if len(nonzero)>0:
            cap_val = np.percentile(nonzero, varna_norm_perc*100)
        else:
            cap_val = 1.0
        if cap_val < 1e-9:
            cap_val = 1e-9

        scaled_vals = []
        for val in arr:
            if val> cap_val:
                val = cap_val
            scaled = val / cap_val
            scaled_vals.append(scaled)

        df["VARNA_Scaled"] = scaled_vals

        # Output #1 => varna_color.txt
        varna_txt = os.path.join(sample_varna_dir, f"{transcript}_varna_color.txt")
        with open(varna_txt, "w") as fout:
            for sc_ in scaled_vals:
                fout.write(f"{sc_:.4f}\n")

        # Output #2 => CSV
        varna_csv = os.path.join(sample_varna_dir, f"{transcript}_varna_scaled.csv")
        df.to_csv(varna_csv, index=False)

        print(f"[VARNA] {sample_name}/{transcript} => wrote {varna_txt}, {varna_csv}")

    print("[VARNA] All transcripts done. Varna outputs in", varna_out_dir)
