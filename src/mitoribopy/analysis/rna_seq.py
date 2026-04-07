# rna_seq_analysis.py
"""
RNA-seq analysis from a single directory with .bed files named after the RNA sample names,
plus a one-to-one mapping to Ribo-seq sample names for ratio.

Example usage:
  run_rna_seq_analysis(
      rna_seq_dir="total_bed_files",
      rna_order=["WT_total","PET_total","MSS116_total","MSS51_total"],
      ribo_order=["WT","PET","MSS116","MSS51"],    # to match indices
      annotation_df=...,
      fasta_file=...,
      output_dir="rna_seq_results",
      translation_profile_dir="analysis_results",
      do_merge_with_ribo=True
  )
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from ..console import iter_with_progress, log_info, log_warning
from ..plotting.style import apply_publication_style

apply_publication_style()

def run_rna_seq_analysis(
    rna_seq_dir,
    rna_order,
    annotation_df,
    fasta_file,
    output_dir,
    translation_profile_dir=None,
    do_merge_with_ribo=False,
    ribo_order=None
):
    """
    :param rna_seq_dir: single folder containing e.g. WT_total.bed, PET_total.bed, ...
    :param rna_order: list of RNA sample names (the .bed base names).
    :param annotation_df: transcripts info
    :param fasta_file: path to .fasta
    :param output_dir: where to write coverage, summary, ratio
    :param translation_profile_dir: where Ribo footprint outputs are stored, e.g. 'analysis_results'
    :param do_merge_with_ribo: produce ratio => p-site / (CoveragePerLength)
    :param ribo_order: list of Ribo sample names, e.g. ["WT","PET","MSS116","MSS51"]
                       must have the same length as rna_order for index matching.
    """

    log_info("RNA-SEQ", "Starting RNA-seq analysis with paired RNA and Ribo sample ordering.")
    os.makedirs(output_dir, exist_ok=True)

    # parse FASTA
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    coverage_map = {}       # coverage_map[rna_name][transcript] = coverage array
    transcript_counts = {}  # transcript_counts[rna_name][transcript] = total coverage

    # 1) Read each RNA sample => e.g. rna_name = "WT_total" => "WT_total.bed"
    for rna_name in iter_with_progress(
        list(rna_order),
        component="RNA-SEQ",
        noun="RNA sample",
        labeler=str,
    ):
        bedfile = f"{rna_name}.bed"
        bedpath = os.path.join(rna_seq_dir, bedfile)
        if not os.path.isfile(bedpath):
            log_warning("RNA-SEQ", f"BED file missing => {bedpath}; skipping sample {rna_name}.")
            continue

        log_info("RNA-SEQ", f"Processing RNA sample => {rna_name} from {bedfile}")
        coverage_map[rna_name] = {}
        transcript_counts[rna_name] = {}

        # each sample => subfolder
        sample_out = os.path.join(output_dir, rna_name)
        os.makedirs(sample_out, exist_ok=True)
        coverage_subdir = os.path.join(sample_out, "rna_coverage")
        os.makedirs(coverage_subdir, exist_ok=True)

        # read bed
        try:
            bed_df = pd.read_csv(bedpath, sep="\t", header=None)
        except Exception as e:
            log_warning("RNA-SEQ", f"Error reading {bedfile}: {e}")
            continue

        # minimal columns
        if bed_df.shape[1]>=3:
            std_cols = ["chrom","start","end","name","score","strand"][:bed_df.shape[1]]
        else:
            log_warning("RNA-SEQ", f"Unexpected BED format => {bedfile}; skipping.")
            continue
        bed_df.columns = std_cols
        bed_df["start"] = pd.to_numeric(bed_df["start"], errors="coerce")
        bed_df["end"]   = pd.to_numeric(bed_df["end"], errors="coerce")
        bed_df.dropna(subset=["start","end"], inplace=True)
        bed_df["start"] = bed_df["start"].astype(int)
        bed_df["end"]   = bed_df["end"].astype(int)

        for chrom_ in bed_df["chrom"].unique():
            if chrom_ not in fasta_dict:
                continue
            seq_len = len(fasta_dict[chrom_].seq)
            if chrom_ not in coverage_map[rna_name]:
                coverage_map[rna_name][chrom_] = np.zeros(seq_len, dtype=int)

            sub_ = bed_df[bed_df["chrom"]== chrom_]
            for _, row_ in sub_.iterrows():
                # BED coordinates are 0-based, end-exclusive.
                st_ = max(int(row_["start"]), 0)
                ed_ = min(int(row_["end"]), seq_len)
                if st_ >= ed_:
                    continue
                coverage_map[rna_name][chrom_][st_:ed_] += 1

        # write coverage + plots
        for chrom_ in coverage_map[rna_name]:
            arr_ = coverage_map[rna_name][chrom_]
            seq_len = len(arr_)
            tot_cov = arr_.sum()
            transcript_counts[rna_name][chrom_] = tot_cov

            df_cov = pd.DataFrame({
                "Position": range(1, seq_len+1),
                "Coverage": arr_
            })
            cov_csv = os.path.join(coverage_subdir, f"{chrom_}_rna_coverage.csv")
            df_cov.to_csv(cov_csv, index=False)

            plt.figure(figsize=(14,4))
            plt.bar(df_cov["Position"], df_cov["Coverage"], color='blue', width=1.0)
            plt.title(f"{chrom_} RNA-seq Coverage ({rna_name})")
            plt.xlabel("Position (1-based)")
            plt.ylabel("Coverage")
            plt.tight_layout()
            out_png = os.path.join(coverage_subdir, f"{chrom_}_rna_coverage.png")
            plt.savefig(out_png)
            plt.close()

    # 2) transcript-level summary
    summary_data = []
    for rna_name in coverage_map:
        for tr_ in coverage_map[rna_name]:
            length_ = len(coverage_map[rna_name][tr_])
            tc_     = transcript_counts[rna_name][tr_]
            cpl_    = tc_/length_ if length_>0 else 0
            summary_data.append({
                "Sample": rna_name,
                "Transcript": tr_,
                "Length": length_,
                "TotalCoverage": tc_,
                "CoveragePerLength": cpl_
            })
    summary_df = pd.DataFrame(summary_data)
    out_summary_csv = os.path.join(output_dir, "rna_seq_transcript_summary.csv")
    summary_df.to_csv(out_summary_csv, index=False)
    log_info("RNA-SEQ", f"Wrote transcript-level summary => {out_summary_csv}")

    # 3) Merge with Ribo => ratio
    if not do_merge_with_ribo:
        log_info("RNA-SEQ", "RNA-only analysis complete; skipping Ribo merge.")
        return

    if not translation_profile_dir or not os.path.isdir(translation_profile_dir):
        log_warning(
            "RNA-SEQ",
            f"translation_profile_dir not found => {translation_profile_dir}; skipping ratio.",
        )
        return

    # We need ribo_order to match rna_order by index
    if (ribo_order is None) or (len(ribo_order)!=len(rna_order)):
        log_warning("RNA-SEQ", "Mismatch in length between rna_order and ribo_order; skipping ratio.")
        return

    log_info("RNA-SEQ", "Merging Ribo footprints into RNA/Ribo ratio outputs.")
    ratio_dir = os.path.join(output_dir, "ribo_ratio")
    os.makedirs(ratio_dir, exist_ok=True)

    # Build coverage_per_len_map => coverage_per_len_map[rna_name][transcript]
    coverage_per_len_map = {}
    for rna_name in coverage_map:
        coverage_per_len_map[rna_name] = {}
    for _, row_ in summary_df.iterrows():
        s_  = row_["Sample"]   # e.g. WT_total
        tr_ = row_["Transcript"]
        cpl_= row_["CoveragePerLength"]
        coverage_per_len_map[s_][tr_] = cpl_

    # We will collect total p-site coverage for each transcript in each sample
    # then compute ratio = total_p_site / coveragePerLen.
    all_ratio_records = []

    # For each index i:
    #   rna_name = rna_order[i]
    #   ribo_name= ribo_order[i]
    # Then footprints are at translation_profile_dir/ribo_name/footprint_density
    # ratio => (sum of p-site coverage) / coveragePerLen(rna_name).
    for i in range(len(rna_order)):
        rna_name  = rna_order[i]   # e.g. "WT_total"
        ribo_name = ribo_order[i]  # e.g. "WT"
        log_info("RNA-SEQ", f"Pairing Ribo '{ribo_name}' with RNA '{rna_name}'")

        foot_dir = os.path.join(translation_profile_dir, ribo_name, "footprint_density")
        if not os.path.isdir(foot_dir):
            log_warning(
                "RNA-SEQ",
                f"Ribo footprint directory not found => {foot_dir}; skipping {ribo_name}.",
            )
            continue

        sample_ratio_dir = os.path.join(ratio_dir, rna_name)
        os.makedirs(sample_ratio_dir, exist_ok=True)

        foot_csvs = [f for f in os.listdir(foot_dir) if f.endswith("_footprint_density.csv")]
        for fc in foot_csvs:
            transcript = fc.replace("_footprint_density.csv","")
            foot_path = os.path.join(foot_dir, fc)

            try:
                df_foot = pd.read_csv(foot_path)
            except Exception as exc:
                log_warning("RNA-SEQ", f"Failed reading footprint CSV {foot_path}: {exc}")
                continue
            if "P_site" not in df_foot.columns:
                continue

            # total p-site coverage for this transcript
            total_p_site = df_foot["P_site"].sum()

            # coveragePerLength from the RNA summary
            cpl_ = coverage_per_len_map[rna_name].get(transcript, 0)

            # ratio for this transcript
            ratio_val = 0
            if cpl_ != 0:
                ratio_val = total_p_site / cpl_

            # store per-position ratio in df_foot if desired
            df_foot["NormalizedPsite"] = 0
            if cpl_ != 0:
                df_foot["NormalizedPsite"] = df_foot["P_site"]/ cpl_

            out_ratio_csv = os.path.join(sample_ratio_dir, f"{transcript}_rna_ribo_ratio.csv")
            df_foot.to_csv(out_ratio_csv, index=False)

            # Save the position-wise bar plot
            plt.figure(figsize=(14,4))
            plt.bar(df_foot["Position"], df_foot["NormalizedPsite"], color='green', width=1.0)
            plt.title(f"{transcript} - ratio => Ribo:{ribo_name} / RNA:{rna_name}")
            plt.xlabel("Position (1-based)")
            plt.ylabel("P-site / RNA_CovPerLen")
            plt.tight_layout()
            out_png = os.path.join(sample_ratio_dir, f"{transcript}_rna_ribo_ratio.png")
            plt.savefig(out_png)
            plt.close()

            # Collect the total ratio for the final summary plot
            all_ratio_records.append({
                "Transcript": transcript,
                "Sample": ribo_name,  # or rna_name if you prefer matching
                "Ratio": ratio_val
            })

    # -- Now produce a single CSV and a grouped barplot showing ratio vs. transcript vs. sample
    ratio_summary_df = pd.DataFrame(all_ratio_records)
    combined_ratio_csv = os.path.join(ratio_dir, "rna_ribo_ratio_summary.csv")
    ratio_summary_df.to_csv(combined_ratio_csv, index=False)
    log_info("RNA-SEQ", f"Wrote combined ratio summary => {combined_ratio_csv}")

    if not ratio_summary_df.empty:
        plt.figure(figsize=(12, 6))
        sns.barplot(data=ratio_summary_df, x="Transcript", y="Ratio", hue="Sample")
        plt.title("Per-transcript Ribo/RNA Ratio")
        plt.xticks(rotation=90)
        plt.tight_layout()
        combined_ratio_png = os.path.join(ratio_dir, "rna_ribo_ratio_per_transcript.png")
        plt.savefig(combined_ratio_png)
        plt.close()
        log_info("RNA-SEQ", f"Wrote ratio barplot => {combined_ratio_png}")
    else:
        log_warning("RNA-SEQ", "No ratio data available to plot.")

    log_info("RNA-SEQ", f"Done merging; ratio outputs stored in {ratio_dir}")
