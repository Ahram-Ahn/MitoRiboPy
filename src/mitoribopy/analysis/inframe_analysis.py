"""In-frame footprint, frame-usage, and codon-usage analysis."""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from math import log2
from ..plotting.style import apply_publication_style

from Bio import SeqIO
from ..data import (
    yeast_mitochondrial_codon_table,
    human_mitochondrial_codon_table
)

apply_publication_style()

def run_inframe_codon_analysis(
    sample_dirs,
    a_site_offsets_df,
    offset_type,
    fasta_file,
    output_dir,
    args,
    annotation_df
):
    """
    Analyzes coverage in frames 0,1,2, producing codon usage CSV/plots (both P-site & A-site)
    and frame usage plots. We respect each transcript’s start_codon when iterating over codons
    to ensure proper in-frame alignment.

    If args.merge_density == True:
        We add coverage from the adjacent nucleotides (i - 1, i + 1) to the central position i.
    """

    print("Starting in-frame codon analysis (frames 0/1/2). merge_density =", args.merge_density)
    os.makedirs(output_dir, exist_ok=True)

    # 1) Load FASTA
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    # 2) Pick codon table
    if args.strain == "y":
        codon_table = yeast_mitochondrial_codon_table
        start_codons = {"ATG"}
        hydrophobic = {"A", "F", "L", "I", "M", "V", "W", "Y"}
        hydrophilic = {"S", "T", "N", "Q"}
        ionic = {"K", "R", "H", "D", "E"}
    else:
        codon_table = human_mitochondrial_codon_table
        start_codons = {"ATG", "ATA"}
        hydrophobic = {"A", "F", "L", "I", "M", "V", "W", "Y"}
        hydrophilic = {"S", "T", "N", "Q"}
        ionic = {"K", "R", "H", "D", "E"}


    # 3) Define codon classification


    def classify_codon(codon_str):
        c_ = codon_str.upper()
        aa_ = codon_table.get(c_, "X")
        if c_ in start_codons:
            category = "Start"
        elif aa_ == "*":
            category = "Stop"
        elif aa_ in hydrophobic:
            category = "Hydrophobic"
        elif aa_ in hydrophilic:
            category = "Hydrophilic"
        elif aa_ in ionic:
            category = "Ionic"
        else:
            category = "Other"
        return aa_, category

    # 4) Sort codons for consistent output
    aa_order = [
        "M", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
        "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*"
    ]
    all_codon_keys = []
    for cod in codon_table.keys():
        aa, cat = classify_codon(cod)
        all_codon_keys.append((cod, aa, cat))

    def sort_key(item):
        c_, aa_, cat_ = item
        try:
            i_aa = aa_order.index(aa_)
        except ValueError:
            i_aa = len(aa_order)
        return (i_aa, c_)

    all_codon_keys = sorted(all_codon_keys, key=sort_key)
    codon_label_order = [f"{cod}-{aa}" for cod, aa, _ in all_codon_keys]

    # 5) Precompute reference codon frequency in frame0
    transcript_codon_counts = {}
    for t_name in annotation_df["transcript"].unique():
        if t_name not in fasta_dict:
            continue
        seq_obj = fasta_dict[t_name]
        seq_str = str(seq_obj.seq).upper()

        ann_row = annotation_df[annotation_df["transcript"] == t_name].iloc[0]
        start_pos = ann_row["start_codon"]
        stop_pos = ann_row["stop_codon"] + 3  # typically +3 to include last codon
        start_pos = max(start_pos, 0)
        stop_pos = min(stop_pos, len(seq_str))
        print(f"{t_name}", start_pos, stop_pos)
        freq_map = {}
        for i in range(start_pos, stop_pos - 2, 3):
            cod3 = seq_str[i:i + 3]
            aa_, cat_ = classify_codon(cod3)
            freq_map.setdefault((cod3, aa_, cat_), 0)
            freq_map[(cod3, aa_, cat_)] += 1
        transcript_codon_counts[t_name] = freq_map

    # Overall structures
    coverage_usage_per_sample = {}
    frame_usage_map = {}

    offset_site = str(getattr(args, "offset_site", "p")).lower()
    if offset_site not in {"p", "a"}:
        offset_site = "p"
    primary_site_key = "P_site" if offset_site == "p" else "A_site"
    primary_site_label = "P-site" if offset_site == "p" else "A-site"
    primary_site_plot_color = "forestgreen" if offset_site == "p" else "darkorange"

    # -----------
    # Process Each Sample
    # -----------
    for sample_dir in sample_dirs:
        sample_name = os.path.basename(sample_dir)
        print(f"\nProcessing sample => {sample_name}")

        sample_out = os.path.join(output_dir, sample_name)
        os.makedirs(sample_out, exist_ok=True)

        # subdirs
        foot_dir = os.path.join(sample_out, "footprint_density")
        os.makedirs(foot_dir, exist_ok=True)

        codon_dir = os.path.join(sample_out, "codon_usage")
        os.makedirs(codon_dir, exist_ok=True)

        frame_dir = os.path.join(sample_out, "translating_frame")
        os.makedirs(frame_dir, exist_ok=True)

        debug_dir = os.path.join(sample_out, "debug_csv")
        os.makedirs(debug_dir, exist_ok=True)

        bed_files = [f for f in os.listdir(sample_dir) if f.endswith(".bed")]
        coverage_usage_per_sample[sample_name] = {}
        frame_usage_map[sample_name] = {}

        # ----------
        # read_cov_map:
        #   { transcript => { "A_site": [...], "P_site": [...], "E_site": [...], "seq": <string> } }
        # ----------
        read_cov_map = {}
        for bf in bed_files:
            rl_match = re.search(r"_(\d+)nt\.bed$", bf)
            if not rl_match:
                continue
            rl = int(rl_match.group(1))

            offset_row = a_site_offsets_df[a_site_offsets_df["Read Length"] == rl]
            if offset_row.empty:
                continue

            if offset_type == "5":
                if "Most Enriched 5' Offset" in offset_row:
                    offset_ = offset_row["Most Enriched 5' Offset"].values[0]
                else:
                    continue
            else:  # offset_type == "3"
                if "Most Enriched 3' Offset" in offset_row:
                    offset_ = offset_row["Most Enriched 3' Offset"].values[0]
                else:
                    continue

            if pd.isna(offset_):
                continue
            offset_ = int(offset_)

            bed_path = os.path.join(sample_dir, bf)
            bed_df = pd.read_csv(bed_path, sep="\t", header=None)
            if bed_df.shape[1] >= 2:
                std_cols = ["chrom", "start", "end", "name", "score", "strand"]
                if bed_df.shape[1] <= len(std_cols):
                    bed_df.columns = std_cols[:bed_df.shape[1]]
                else:
                    extras = [f"extra_{i+1}" for i in range(bed_df.shape[1] - len(std_cols))]
                    bed_df.columns = std_cols + extras
            else:
                continue
            # I should give a filtered_bed not the bed_df
            bed_df["start"] = pd.to_numeric(bed_df["start"], errors="coerce")
            bed_df["end"] = pd.to_numeric(bed_df.get("end"), errors="coerce")
            bed_df.dropna(subset=["start", "end"], inplace=True)
            bed_df["start"] = bed_df["start"].astype(int)
            bed_df["end"] = bed_df["end"].astype(int)
            bed_df = bed_df[bed_df["end"] > bed_df["start"]].copy()
            if bed_df.empty:
                continue

            # Compute A/P/E sites
            if offset_type == "5":
                if offset_site == "p":
                    bed_df["P_site"] = bed_df["start"] + offset_
                else:
                    # Offset points to A-site in this mode; convert to P-site.
                    bed_df["P_site"] = bed_df["start"] + offset_ - 3
            else:
                # BED is end-exclusive, so 3' end index is (end - 1).
                if offset_site == "p":
                    bed_df["P_site"] = bed_df["end"] - offset_ - 1
                else:
                    # Offset points to A-site in this mode; convert to P-site.
                    bed_df["P_site"] = bed_df["end"] - offset_ - 4

            bed_df["A_site"] = bed_df["P_site"] + 3
            bed_df["E_site"] = bed_df["P_site"] - 3

            # Update read_cov_map
            for t_ in bed_df["chrom"].unique():
                if t_ not in fasta_dict:
                    continue
                if t_ not in annotation_df["transcript"].values:
                    continue

                if t_ not in read_cov_map:
                    seq_len = len(fasta_dict[t_].seq)
                    read_cov_map[t_] = {
                        "seq": str(fasta_dict[t_].seq).upper(),
                        "A_site": [0] * seq_len,
                        "P_site": [0] * seq_len,
                        "E_site": [0] * seq_len
                    }

                local_df = bed_df[bed_df["chrom"] == t_]
                for site_ in ["A_site", "P_site", "E_site"]:
                    c_ = local_df[site_].value_counts()
                    for px, cnt_ in c_.items():
                        if 0 <= px < len(read_cov_map[t_][site_]):
                            read_cov_map[t_][site_][px] += cnt_

        # Save read_cov_map to CSV for debugging
        debug_rows = []
        for t_ in read_cov_map:
            seq_len = len(read_cov_map[t_]["seq"])
            A_cov = read_cov_map[t_]["A_site"]
            P_cov = read_cov_map[t_]["P_site"]
            E_cov = read_cov_map[t_]["E_site"]
            for i in range(seq_len):
                debug_rows.append({
                    "Transcript": t_,
                    "Position": i + 1,
                    "A_site": A_cov[i],
                    "P_site": P_cov[i],
                    "E_site": E_cov[i]
                })
        debug_df = pd.DataFrame(debug_rows)
        debug_df.to_csv(os.path.join(debug_dir, f"read_cov_map_{sample_name}.csv"), index=False)

        # -----------
        # Frame usage + primary-site coverage usage
        # -----------
        coverage_usage_per_sample[sample_name] = {}

        for t_ in read_cov_map:
            seq_len = len(read_cov_map[t_]["seq"])
            primary_depth = np.array(read_cov_map[t_][primary_site_key])
            nonzero = primary_depth[primary_depth > 0]
            if len(nonzero) > 0:
                cap_val = np.percentile(nonzero, args.cap_percentile * 100)
            else:
                cap_val = 1
            primary_clipped = np.clip(primary_depth, None, cap_val).astype(int)

            # Save footprint coverage
            p_depth = np.array(read_cov_map[t_]["P_site"])
            a_depth = np.array(read_cov_map[t_]["A_site"])
            e_depth = np.array(read_cov_map[t_]["E_site"])
            foot_df = pd.DataFrame({
                "Position": range(1, seq_len + 1),
                "Nucleotide": list(read_cov_map[t_]["seq"]),
                "P_site": p_depth,
                "A_site": a_depth,
                "E_site": e_depth,
                f"{primary_site_key}_selected_depth": primary_clipped,
            })
            foot_csv_path = os.path.join(foot_dir, f"{t_}_footprint_density.csv")
            foot_df.to_csv(foot_csv_path, index=False)

            plt.figure(figsize=(12, 3))
            plt.bar(
                foot_df["Position"],
                foot_df[f"{primary_site_key}_selected_depth"],
                color=primary_site_plot_color,
                width=3.0,
            )
            plt.title(f"{t_} - {primary_site_label} ({sample_name})")
            plt.xlabel("Position")
            plt.ylabel(f"{primary_site_label} Depth")
            plt.tight_layout()
            plt.savefig(os.path.join(foot_dir, f"{t_}_{primary_site_key}_depth.png"))
            plt.close()

            # Frame usage
            if t_ not in annotation_df["transcript"].values:
                continue

            ann_ = annotation_df[annotation_df["transcript"] == t_].iloc[0]
            start_cdn = ann_["start_codon"]
            # identify frame0/1/2 positions relative to start_codon
            frame0_positions = []
            frame1_positions = []
            frame2_positions = []
            for i in range(seq_len):
                rel = (i - start_cdn) % 3
                if rel == 0:
                    frame0_positions.append(i)
                elif rel == 1:
                    frame1_positions.append(i)
                else:
                    frame2_positions.append(i)

            f0_count = primary_depth[frame0_positions].sum() if frame0_positions else 0
            f1_count = primary_depth[frame1_positions].sum() if frame1_positions else 0
            f2_count = primary_depth[frame2_positions].sum() if frame2_positions else 0
            frame_usage_map[sample_name].setdefault(t_, (f0_count, f1_count, f2_count))

            # Build codon coverage for the selected primary site.
            coverage_usage_per_sample[sample_name].setdefault(t_, {})
            seq_str = read_cov_map[t_]["seq"]

            transcript_stop_codon = ann_["stop_codon"] + 3
            upper_bound = min(transcript_stop_codon, seq_len - 2)

            # step in multiples of 3
            for i in range(start_cdn, upper_bound, 3):
                middle_pos = i + 1
                if (middle_pos + 1) >= seq_len:
                    continue
                codon_seq = seq_str[middle_pos - 1: middle_pos + 2]
                aa_, cat_ = classify_codon(codon_seq)

                cov_here = primary_depth[middle_pos]
                if args.merge_density:
                    if (middle_pos - 1) >= 0:
                        cov_here += primary_depth[middle_pos - 1]
                    if (middle_pos + 1) < seq_len:
                        cov_here += primary_depth[middle_pos + 1]

                coverage_usage_per_sample[sample_name][t_].setdefault((codon_seq, aa_, cat_), 0)
                coverage_usage_per_sample[sample_name][t_][(codon_seq, aa_, cat_)] += cov_here

        # Debug: coverage_usage_per_sample -> flatten and save
        debug_cov_rows = []
        for t_ in coverage_usage_per_sample[sample_name]:
            for (codon, aa_, cat_), c_val in coverage_usage_per_sample[sample_name][t_].items():
                debug_cov_rows.append({
                    "Transcript": t_,
                    "Codon": codon,
                    "AA": aa_,
                    "Category": cat_,
                    "Coverage": c_val
                })
        debug_cov_df = pd.DataFrame(debug_cov_rows)
        debug_cov_df.to_csv(os.path.join(debug_dir, f"coverage_usage_per_sample_{sample_name}.csv"), index=False)

        # -------------------------------------
        # Frame Usage Summaries & Plots
        # -------------------------------------
        frame_t_path = os.path.join(frame_dir, "frame_usage_by_transcript.csv")
        total_f0, total_f1, total_f2 = 0, 0, 0
        with open(frame_t_path, "w") as fw_:
            fw_.write("Transcript,Frame0_Count,Frame1_Count,Frame2_Count,Frame0_Prop,Frame1_Prop,Frame2_Prop\n")
            for t_name, (f0c, f1c, f2c) in frame_usage_map[sample_name].items():
                total_ = f0c + f1c + f2c
                f0p = f0c / total_ if total_ > 0 else 0
                f1p = f1c / total_ if total_ > 0 else 0
                f2p = f2c / total_ if total_ > 0 else 0
                fw_.write(f"{t_name},{f0c},{f1c},{f2c},{f0p},{f1p},{f2p}\n")
                total_f0 += f0c
                total_f1 += f1c
                total_f2 += f2c

        total_all = total_f0 + total_f1 + total_f2
        f0_prop = total_f0 / total_all if total_all > 0 else 0
        f1_prop = total_f1 / total_all if total_all > 0 else 0
        f2_prop = total_f2 / total_all if total_all > 0 else 0

        sample_frame_out = os.path.join(frame_dir, "frame_usage_total.csv")
        with open(sample_frame_out, "w") as fw_:
            fw_.write("Frame,Count,Proportion\n")
            fw_.write(f"Frame0,{total_f0},{f0_prop}\n")
            fw_.write(f"Frame1,{total_f1},{f1_prop}\n")
            fw_.write(f"Frame2,{total_f2},{f2_prop}\n")

        # Frame usage total plot
        frame_df = pd.DataFrame({
            "Frame": ["Frame0", "Frame1", "Frame2"],
            "Count": [total_f0, total_f1, total_f2],
            "Proportion": [f0_prop, f1_prop, f2_prop]
        })
        plt.figure(figsize=(6, 6))
        sns.barplot(x="Frame", y="Proportion", data=frame_df, hue="Frame", dodge=False, palette="Set2", legend=False)
        plt.title(f"Total Frame Usage Proportion ({primary_site_label}) - {sample_name}")
        plt.tight_layout()
        plt.savefig(os.path.join(frame_dir, "frame_usage_total_plot.png"))
        plt.close()

        # By-transcript frame usage
        frame_t_data = []
        with open(frame_t_path) as fin:
            next(fin)  # skip header
            for line in fin:
                t_name, f0c, f1c, f2c, f0p, f1p, f2p = line.strip().split(",")
                f0p, f1p, f2p = float(f0p), float(f1p), float(f2p)
                frame_t_data.append((t_name, "Frame0", f0p))
                frame_t_data.append((t_name, "Frame1", f1p))
                frame_t_data.append((t_name, "Frame2", f2p))

        frame_t_df = pd.DataFrame(frame_t_data, columns=["Transcript", "Frame", "Proportion"])
        plt.figure(figsize=(12, 6))
        sns.barplot(x="Transcript", y="Proportion", hue="Frame", data=frame_t_df, palette="Set3")
        plt.title(f"Frame Usage Proportion by Transcript ({primary_site_label}) - {sample_name}")
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(os.path.join(frame_dir, "frame_usage_by_transcript_plot.png"))
        plt.close()

        # -------------------------------------
        # Process and plot primary-site codon usage (Per-transcript & Overall)
        # -------------------------------------
        codon_usage_dir = os.path.join(sample_out, "codon_usage")
        sample_ov_cov = {}
        sample_ov_freq = {}
        usage_dict_for_sample = coverage_usage_per_sample[sample_name]

        # For each transcript, mask (i.e. completely remove) the final translating codon 
        # so that its P-site coverage is not included in the analysis.
        for t, usage_dict in usage_dict_for_sample.items():
            if (t in annotation_df["transcript"].values) and (t in fasta_dict):
                ann_row = annotation_df[annotation_df["transcript"] == t].iloc[0]
                transcript_stop_codon = ann_row["stop_codon"] + 3
                seq_str = str(fasta_dict[t].seq).upper()
                if transcript_stop_codon > len(seq_str):
                    transcript_stop_codon = len(seq_str)
                last_codon_seq = seq_str[transcript_stop_codon - 3: transcript_stop_codon]
                aa, cat = classify_codon(last_codon_seq)
                masked_key = (last_codon_seq, aa, cat)
                # Remove the masked codon from the usage dictionary.
                if masked_key in usage_dict:
                    usage_dict.pop(masked_key)

        # Now process each transcript for primary-site codon usage analysis.
        rows_list = []
        for t, usage_dict in usage_dict_for_sample.items():
            # Get the transcript frequency dictionary (as computed from the full transcript)
            freq_dict = transcript_codon_counts.get(t, {})
            # Create a local copy so that we can remove the masked codon frequency as well.
            local_freq_dict = freq_dict.copy()
            
            # Compute the masked key again for this transcript.
            masked_key = None
            if (t in annotation_df["transcript"].values) and (t in fasta_dict):
                ann_row = annotation_df[annotation_df["transcript"] == t].iloc[0]
                transcript_stop_codon = ann_row["stop_codon"] + 3
                seq_str = str(fasta_dict[t].seq).upper()
                if transcript_stop_codon > len(seq_str):
                    transcript_stop_codon = len(seq_str)
                last_codon_seq = seq_str[transcript_stop_codon - 3: transcript_stop_codon]
                aa, cat = classify_codon(last_codon_seq)
                masked_key = (last_codon_seq, aa, cat)
                # Remove the masked key from the frequency dictionary so that the codon is not counted.
                if masked_key in local_freq_dict:
                    local_freq_dict.pop(masked_key)
            
            rows = []
            cdf_map = {}
            # Loop over all codon keys. For the masked codon, keep it in output
            # but force zero values so plots still show a complete codon axis.
            for key in all_codon_keys:
                if masked_key is not None and key == masked_key:
                    c_val = 0
                    f_val = 0
                else:
                    c_val = usage_dict.get(key, 0)
                    f_val = local_freq_dict.get(key, 0)
                # Aggregate overall values.
                sample_ov_cov[key] = sample_ov_cov.get(key, 0) + c_val
                sample_ov_freq[key] = sample_ov_freq.get(key, 0) + f_val
                cdf_map[key] = c_val / f_val if f_val > 0 else 0

            sum_cdf = sum(cdf_map.values())
            for key in all_codon_keys:
                if masked_key is not None and key == masked_key:
                    c_val = 0
                    f_val = 0
                else:
                    c_val = usage_dict.get(key, 0)
                    f_val = local_freq_dict.get(key, 0)
                cdf_val = cdf_map.get(key, 0)
                normp = cdf_val / sum_cdf if sum_cdf > 0 else 0
                rows.append({
                    "Transcript": t,
                    "Codon": key[0],
                    "AA": key[1],
                    "Category": key[2],
                    "Coverage": c_val,
                    "Frequency": f_val,
                    "CoverageDivFreq": cdf_val,
                    "NormalizedProp": normp
                })
            
            out_csv = os.path.join(codon_usage_dir, f"codon_usage_{t}.csv")
            df_rows = pd.DataFrame(rows)
            df_rows.to_csv(out_csv, index=False)
            df_rows["Label"] = df_rows["Codon"] + "-" + df_rows["AA"]
            plt.figure(figsize=(20, 6))
            sns.barplot(x="Label", y="CoverageDivFreq", hue="Category",
                        data=df_rows, dodge=False, order=codon_label_order)
            plt.xticks(rotation=90)
            y_max = df_rows["CoverageDivFreq"].max() * 1.2 if not df_rows.empty else 1
            plt.ylim(0, y_max)            
            plt.title(
                f"{t} - {primary_site_label} Codon Usage "
                f"(Frame0{', Merged' if args.merge_density else ''}) - {sample_name}"
            )
            plt.tight_layout()
            out_fig = os.path.join(codon_usage_dir, f"codon_usage_{t}_plot.png")
            plt.savefig(out_fig)
            plt.close()
            rows_list.append(df_rows)

        # Overall A-site usage aggregation and plotting.
        ov_rows = []
        cov_div_freq_map = {}
        PSEUDO = 2  
        # Trying to use empirical-Bayes (pseudocount) shrinkage
        # The pseudocount makes the denominator non-zero even when freq_val == 0, so the if … else 0 guard is no longer needed
        # Adding PSEUDO to both numerator and denominator shrinks extreme values back toward the grand mean.
        for key in all_codon_keys:
            cov_val = sample_ov_cov.get(key, 0)
            freq_val = sample_ov_freq.get(key, 0)
            cdf_val = (cov_val + PSEUDO) / (freq_val + PSEUDO)
            cov_div_freq_map[key] = cdf_val
        S_ov = sum(cov_div_freq_map.values())
        for key in all_codon_keys:
            cov_val = sample_ov_cov.get(key, 0)
            freq_val = sample_ov_freq.get(key, 0)
            cdf_val = cov_div_freq_map[key]
            normp = cdf_val / S_ov if S_ov > 0 else 0
            ov_rows.append({
                "Codon": key[0],
                "AA": key[1],
                "Category": key[2],
                "Coverage": cov_val,
                "Frequency": freq_val,
                "CoverageDivFreq": cdf_val,
                "NormalizedProp": normp
            })
        ov_df = pd.DataFrame(ov_rows)
        ov_csv = os.path.join(codon_usage_dir, "codon_usage_total.csv")
        ov_df.to_csv(ov_csv, index=False)
        ov_df["Label"] = ov_df["Codon"] + "-" + ov_df["AA"]
        plt.figure(figsize=(20, 6))
        sns.barplot(x="Label", y="CoverageDivFreq", hue="Category",
                    data=ov_df, dodge=False, order=codon_label_order)
        plt.xticks(rotation=90)
        overall_y_max = ov_df["CoverageDivFreq"].max() * 1.2 if not ov_df.empty else 1
        plt.ylim(0, overall_y_max)        
        plt.title(
            f"Overall {primary_site_label} Codon Usage "
            f"(Frame0{', Merged' if args.merge_density else ''}) - {sample_name}"
        )
        plt.tight_layout()
        ov_plot = os.path.join(codon_usage_dir, "codon_usage_total_plot.png")
        plt.savefig(ov_plot)
        plt.close()

        # ----------------------------------------
        # A-site Codon Usage (Per-transcript & Overall)
        # Note: For A-site usage, we do NOT mask the stop codon.
        # ----------------------------------------
        coverage_usage_a_site = {}
        for t in read_cov_map:
            seq_len = len(read_cov_map[t]["seq"])
            a_depth = np.array(read_cov_map[t]["A_site"])
            coverage_usage_a_site[t] = {}
            if t not in annotation_df["transcript"].values:
                continue
            ann_ = annotation_df[annotation_df["transcript"] == t].iloc[0]
            start_cdn = ann_["start_codon"]
            transcript_stop_codon = ann_["stop_codon"] + 3
            if transcript_stop_codon > len(read_cov_map[t]["seq"]):
                transcript_stop_codon = len(read_cov_map[t]["seq"])
            upper_bound = min(transcript_stop_codon, seq_len - 2)
            seq_str = read_cov_map[t]["seq"]
            for i in range(start_cdn, upper_bound, 3):
                mid = i + 1
                if (mid + 1) >= seq_len:
                    continue
                codon_seq = seq_str[mid - 1: mid + 2]
                aa_, cat_ = classify_codon(codon_seq)
                cov_here = a_depth[mid]
                if args.merge_density:
                    if (mid - 1) >= 0:
                        cov_here += a_depth[mid - 1]
                    if (mid + 1) < seq_len:
                        cov_here += a_depth[mid + 1]
                coverage_usage_a_site[t].setdefault((codon_seq, aa_, cat_), 0)
                coverage_usage_a_site[t][(codon_seq, aa_, cat_)] += cov_here

        # Save the A-site coverage usage for debugging/analysis.
        debug_a_rows = []
        for t in coverage_usage_a_site:
            for (codon, aa_, cat_), c_val in coverage_usage_a_site[t].items():
                debug_a_rows.append({
                    "Transcript": t,
                    "Codon": codon,
                    "AA": aa_,
                    "Category": cat_,
                    "Coverage": c_val
                })
        dbg_a_df = pd.DataFrame(debug_a_rows)
        dbg_a_df.to_csv(os.path.join(debug_dir, f"coverage_usage_a_site_{sample_name}.csv"), index=False)

        # Build A-site codon summary files used by correlation analysis.
        a_site_ov_cov = {}
        a_site_ov_freq = {}
        for t, usage_dict in coverage_usage_a_site.items():
            freq_dict = transcript_codon_counts.get(t, {})
            for key in all_codon_keys:
                a_site_ov_cov[key] = a_site_ov_cov.get(key, 0) + usage_dict.get(key, 0)
                a_site_ov_freq[key] = a_site_ov_freq.get(key, 0) + freq_dict.get(key, 0)

        a_site_rows = []
        a_site_cdf_map = {}
        PSEUDO_A = 2
        for key in all_codon_keys:
            cov_val = a_site_ov_cov.get(key, 0)
            freq_val = a_site_ov_freq.get(key, 0)
            a_site_cdf_map[key] = (cov_val + PSEUDO_A) / (freq_val + PSEUDO_A)
        a_site_sum = sum(a_site_cdf_map.values())

        for key in all_codon_keys:
            cov_val = a_site_ov_cov.get(key, 0)
            freq_val = a_site_ov_freq.get(key, 0)
            cdf_val = a_site_cdf_map[key]
            normp = cdf_val / a_site_sum if a_site_sum > 0 else 0
            a_site_rows.append({
                "Codon": key[0],
                "AA": key[1],
                "Category": key[2],
                "Coverage": cov_val,
                "Frequency": freq_val,
                "CoverageDivFreq": cdf_val,
                "NormalizedProp": normp
            })

        a_site_df = pd.DataFrame(a_site_rows)
        a_site_csv = os.path.join(codon_usage_dir, "a_site_codon_usage_total.csv")
        a_site_df.to_csv(a_site_csv, index=False)
        a_site_df["Label"] = a_site_df["Codon"] + "-" + a_site_df["AA"]
        plt.figure(figsize=(20, 6))
        sns.barplot(
            x="Label",
            y="CoverageDivFreq",
            hue="Category",
            data=a_site_df,
            dodge=False,
            order=codon_label_order,
        )
        plt.xticks(rotation=90)
        a_site_y_max = a_site_df["CoverageDivFreq"].max() * 1.2 if not a_site_df.empty else 1
        plt.ylim(0, a_site_y_max)
        plt.title(f"Overall A-site Codon Usage (Frame0{', Merged' if args.merge_density else ''}) - {sample_name}")
        plt.tight_layout()
        a_site_plot = os.path.join(codon_usage_dir, "a_site_codon_usage_total_plot.png")
        plt.savefig(a_site_plot)
        plt.close()

        print("In-frame codon analysis complete. merge_density =", args.merge_density)
