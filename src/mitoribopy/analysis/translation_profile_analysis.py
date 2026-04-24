"""Translation-profile analysis for footprint density, frame usage, and codon usage."""

import os
import re

import numpy as np
import pandas as pd
from ..console import iter_with_progress, log_info
from ..plotting.translation_profile_plots import (
    plot_codon_usage_dataframe,
    plot_frame_usage_by_transcript,
    plot_frame_usage_total,
    plot_site_depth_profile,
)

from ..data import (
    build_sequence_display_map,
    load_codon_table,
    resolve_sequence_name,
    resolve_start_codons,
    transcript_display_title,
)
from Bio import SeqIO

def run_translation_profile_analysis(
    sample_dirs,
    selected_offsets_df,
    offset_type,
    fasta_file,
    output_dir,
    args,
    annotation_df,
    filtered_bed_df=None,
    resolved_codon_table=None,
    resolved_start_codons=None,
    selected_offsets_by_sample: dict | None = None,
    site_override: str | None = None,
):
    """Compute translation-profile summaries and write the paired plotting outputs."""

    log_info(
        "PROFILE",
        f"Starting translation-profile analysis (merge_density={args.merge_density}).",
    )
    os.makedirs(output_dir, exist_ok=True)

    # 1) Load FASTA
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    # 2) Pick codon table and start codons
    codon_table = resolved_codon_table
    if codon_table is None:
        codon_table = load_codon_table(preset=getattr(args, "strain", None))
    start_codons = set(
        resolved_start_codons
        if resolved_start_codons is not None
        else resolve_start_codons(
            getattr(args, "strain", "custom"),
            getattr(args, "start_codons", None),
        )
    )
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

    def build_codon_usage_dataframe(
        transcript_name,
        usage_dict,
        freq_dict,
        *,
        masked_key=None,
        aggregate_cov=None,
        aggregate_freq=None,
    ):
        rows = []
        cdf_map = {}

        for key in all_codon_keys:
            if masked_key is not None and key == masked_key:
                c_val = 0
                f_val = 0
            else:
                c_val = usage_dict.get(key, 0)
                f_val = freq_dict.get(key, 0)

            if aggregate_cov is not None:
                aggregate_cov[key] = aggregate_cov.get(key, 0) + c_val
            if aggregate_freq is not None:
                aggregate_freq[key] = aggregate_freq.get(key, 0) + f_val

            cdf_map[key] = c_val / f_val if f_val > 0 else 0

        sum_cdf = sum(cdf_map.values())
        for key in all_codon_keys:
            if masked_key is not None and key == masked_key:
                c_val = 0
                f_val = 0
            else:
                c_val = usage_dict.get(key, 0)
                f_val = freq_dict.get(key, 0)

            cdf_val = cdf_map.get(key, 0)
            normp = cdf_val / sum_cdf if sum_cdf > 0 else 0
            rows.append({
                "Transcript": transcript_name,
                "Codon": key[0],
                "AA": key[1],
                "Category": key[2],
                "Coverage": c_val,
                "Frequency": f_val,
                "CoverageDivFreq": cdf_val,
                "NormalizedProp": normp,
            })

        return pd.DataFrame(rows)

    annotation_df = annotation_df.copy()
    if "start_codon" not in annotation_df.columns and {"l_utr5"}.issubset(annotation_df.columns):
        annotation_df["start_codon"] = annotation_df["l_utr5"]
    if "stop_codon" not in annotation_df.columns and {"l_tr", "l_utr3"}.issubset(annotation_df.columns):
        annotation_df["stop_codon"] = annotation_df["l_tr"] - annotation_df["l_utr3"] - 3

    available_sequence_names = set(fasta_dict.keys())
    annotation_df["resolved_sequence_name"] = annotation_df.apply(
        lambda row: resolve_sequence_name(row, available_sequence_names),
        axis=1,
    )
    resolved_annotation_df = annotation_df[annotation_df["resolved_sequence_name"].notna()].copy()
    annotation_lookup = resolved_annotation_df.drop_duplicates("transcript").set_index("transcript")
    annotation_transcripts = set(annotation_lookup.index)
    annotated_sequence_ids = set(resolved_annotation_df["resolved_sequence_name"])
    sequence_display_map = build_sequence_display_map(resolved_annotation_df, available_sequence_names)

    # 5) Precompute reference codon frequency in frame0
    transcript_codon_counts = {}
    for t_name in iter_with_progress(
        list(annotation_lookup.index),
        component="PROFILE",
        noun="reference transcript",
        labeler=str,
    ):
        ann_row = annotation_lookup.loc[t_name]
        resolved_sequence_name = ann_row["resolved_sequence_name"]
        if resolved_sequence_name not in fasta_dict:
            continue
        seq_obj = fasta_dict[resolved_sequence_name]
        seq_str = str(seq_obj.seq).upper()
        start_pos = ann_row["start_codon"]
        stop_pos = ann_row["stop_codon"] + 3  # typically +3 to include last codon
        start_pos = max(start_pos, 0)
        stop_pos = min(stop_pos, len(seq_str))
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

    if site_override is not None:
        offset_site = str(site_override).lower()
    else:
        offset_site = str(getattr(args, "offset_site", "p")).lower()
    if offset_site not in {"p", "a"}:
        offset_site = "p"
    primary_site_key = "P_site" if offset_site == "p" else "A_site"
    primary_site_label = "P-site" if offset_site == "p" else "A-site"
    primary_site_plot_color = "forestgreen" if offset_site == "p" else "darkorange"

    # -----------
    # Process Each Sample
    # -----------
    for sample_dir in iter_with_progress(
        list(sample_dirs),
        component="PROFILE",
        noun="sample",
        labeler=os.path.basename,
    ):
        sample_name = os.path.basename(sample_dir)
        log_info("PROFILE", f"Processing sample => {sample_name}")

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

        coverage_usage_per_sample[sample_name] = {}
        frame_usage_map[sample_name] = {}

        # ----------
        # read_cov_map:
        #   { transcript => { "A_site": [...], "P_site": [...], "E_site": [...], "seq": <string> } }
        # ----------
        read_cov_map = {}
        if filtered_bed_df is not None:
            sample_bed_df = filtered_bed_df[filtered_bed_df["sample_name"] == sample_name].copy()
            read_length_blocks = [
                (int(read_length), sample_bed_df[sample_bed_df["read_length"] == read_length].copy())
                for read_length in sorted(sample_bed_df["read_length"].dropna().unique())
            ]
        else:
            bed_files = [f for f in os.listdir(sample_dir) if f.endswith(".bed")]
            read_length_blocks = []
            for bf in bed_files:
                rl_match = re.search(r"_(\d+)nt\.bed$", bf)
                if not rl_match:
                    continue
                rl = int(rl_match.group(1))
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
                read_length_blocks.append((rl, bed_df))

        sample_offsets_df = (
            selected_offsets_by_sample.get(sample_name, selected_offsets_df)
            if selected_offsets_by_sample
            else selected_offsets_df
        )
        if sample_offsets_df is None:
            log_info(
                "PROFILE",
                f"{sample_name}: no selected offsets table available; skipping.",
            )
            continue
        for rl, bed_df in read_length_blocks:
            offset_row = sample_offsets_df[sample_offsets_df["Read Length"] == rl]
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
                if t_ not in annotated_sequence_ids:
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
        for t_ in iter_with_progress(
            sorted(read_cov_map),
            component="PROFILE",
            noun="transcript",
            labeler=str,
        ):
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
        transcript_label_map: dict[str, str] = {}

        for sequence_id, coverage_row in read_cov_map.items():
            seq_len = len(coverage_row["seq"])
            primary_depth = np.array(coverage_row[primary_site_key])
            nonzero = primary_depth[primary_depth > 0]
            if len(nonzero) > 0:
                cap_val = np.percentile(nonzero, args.cap_percentile * 100)
            else:
                cap_val = 1
            primary_clipped = np.clip(primary_depth, None, cap_val).astype(int)

            foot_df = pd.DataFrame({
                "Position": range(1, seq_len + 1),
                "Nucleotide": list(coverage_row["seq"]),
                "P_site": np.array(coverage_row["P_site"]),
                "A_site": np.array(coverage_row["A_site"]),
                "E_site": np.array(coverage_row["E_site"]),
                f"{primary_site_key}_selected_depth": primary_clipped,
            })
            foot_csv_path = os.path.join(foot_dir, f"{sequence_id}_footprint_density.csv")
            foot_df.to_csv(foot_csv_path, index=False)
            plot_site_depth_profile(
                foot_df,
                os.path.join(foot_dir, f"{sequence_id}_{primary_site_key}_depth.png"),
                site_column=f"{primary_site_key}_selected_depth",
                site_label=primary_site_label,
                sample_name=sample_name,
                transcript_name=sequence_display_map.get(sequence_id, sequence_id),
                color=primary_site_plot_color,
            )

        for _, ann_ in resolved_annotation_df.iterrows():
            sequence_id = ann_["resolved_sequence_name"]
            transcript_name = ann_["transcript"]
            if sequence_id not in read_cov_map:
                continue

            transcript_label_map[transcript_name] = transcript_display_title(
                ann_,
                include_transcript=True,
            )

            seq_len = len(read_cov_map[sequence_id]["seq"])
            primary_depth = np.array(read_cov_map[sequence_id][primary_site_key])
            start_cdn = int(ann_["start_codon"])
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
            frame_usage_map[sample_name][transcript_name] = (f0_count, f1_count, f2_count)

            coverage_usage_per_sample[sample_name].setdefault(transcript_name, {})
            seq_str = read_cov_map[sequence_id]["seq"]
            transcript_stop_codon = int(ann_["stop_codon"]) + 3
            upper_bound = min(transcript_stop_codon, seq_len - 2)

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

                coverage_usage_per_sample[sample_name][transcript_name].setdefault(
                    (codon_seq, aa_, cat_),
                    0,
                )
                coverage_usage_per_sample[sample_name][transcript_name][(codon_seq, aa_, cat_)] += cov_here

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
        frame_usage_rows = []
        total_f0, total_f1, total_f2 = 0, 0, 0
        for t_name, (f0c, f1c, f2c) in frame_usage_map[sample_name].items():
            total_ = f0c + f1c + f2c
            f0p = f0c / total_ if total_ > 0 else 0
            f1p = f1c / total_ if total_ > 0 else 0
            f2p = f2c / total_ if total_ > 0 else 0
            frame_usage_rows.append(
                {
                    "Transcript": t_name,
                    "DisplayName": transcript_label_map.get(t_name, t_name),
                    "Frame0_Count": f0c,
                    "Frame1_Count": f1c,
                    "Frame2_Count": f2c,
                    "Frame0_Prop": f0p,
                    "Frame1_Prop": f1p,
                    "Frame2_Prop": f2p,
                }
            )
            total_f0 += f0c
            total_f1 += f1c
            total_f2 += f2c

        frame_usage_by_transcript_df = pd.DataFrame(
            frame_usage_rows,
            columns=[
                "Transcript",
                "DisplayName",
                "Frame0_Count",
                "Frame1_Count",
                "Frame2_Count",
                "Frame0_Prop",
                "Frame1_Prop",
                "Frame2_Prop",
            ],
        )
        frame_usage_by_transcript_df.to_csv(frame_t_path, index=False)

        total_all = total_f0 + total_f1 + total_f2
        f0_prop = total_f0 / total_all if total_all > 0 else 0
        f1_prop = total_f1 / total_all if total_all > 0 else 0
        f2_prop = total_f2 / total_all if total_all > 0 else 0

        sample_frame_out = os.path.join(frame_dir, "frame_usage_total.csv")
        frame_df = pd.DataFrame({
            "Frame": ["Frame0", "Frame1", "Frame2"],
            "Count": [total_f0, total_f1, total_f2],
            "Proportion": [f0_prop, f1_prop, f2_prop]
        })
        frame_df.to_csv(sample_frame_out, index=False)
        plot_frame_usage_total(
            frame_df,
            os.path.join(frame_dir, "frame_usage_total_plot.png"),
            primary_site_label=primary_site_label,
            sample_name=sample_name,
        )

        frame_t_df = frame_usage_by_transcript_df.melt(
            id_vars=["Transcript", "DisplayName"],
            value_vars=["Frame0_Prop", "Frame1_Prop", "Frame2_Prop"],
            var_name="Frame",
            value_name="Proportion",
        )
        frame_t_df["Frame"] = frame_t_df["Frame"].str.replace("_Prop", "", regex=False)
        frame_t_df["Transcript"] = frame_t_df["DisplayName"]
        plot_frame_usage_by_transcript(
            frame_t_df,
            os.path.join(frame_dir, "frame_usage_by_transcript_plot.png"),
            primary_site_label=primary_site_label,
            sample_name=sample_name,
        )

        # -------------------------------------
        # Process and plot primary-site codon usage (Per-transcript & Overall)
        # -------------------------------------
        codon_usage_dir = os.path.join(sample_out, "codon_usage")
        sample_ov_cov = {}
        sample_ov_freq = {}
        usage_dict_for_sample = coverage_usage_per_sample[sample_name]
        mask_primary_stop_codon = primary_site_key == "P_site"

        # Now process each transcript for primary-site codon usage analysis.
        for t, usage_dict in usage_dict_for_sample.items():
            freq_dict = transcript_codon_counts.get(t, {})
            masked_key = None
            if mask_primary_stop_codon and (t in annotation_transcripts):
                ann_row = annotation_lookup.loc[t]
                resolved_sequence_name = ann_row["resolved_sequence_name"]
                if resolved_sequence_name not in fasta_dict:
                    continue
                transcript_stop_codon = ann_row["stop_codon"] + 3
                seq_str = str(fasta_dict[resolved_sequence_name].seq).upper()
                if transcript_stop_codon > len(seq_str):
                    transcript_stop_codon = len(seq_str)
                last_codon_seq = seq_str[transcript_stop_codon - 3: transcript_stop_codon]
                aa, cat = classify_codon(last_codon_seq)
                masked_key = (last_codon_seq, aa, cat)

            df_rows = build_codon_usage_dataframe(
                t,
                usage_dict,
                freq_dict,
                masked_key=masked_key,
                aggregate_cov=sample_ov_cov,
                aggregate_freq=sample_ov_freq,
            )
            out_csv = os.path.join(codon_usage_dir, f"codon_usage_{t}.csv")
            df_rows.to_csv(out_csv, index=False)
            plot_codon_usage_dataframe(
                df_rows,
                os.path.join(codon_usage_dir, f"codon_usage_{t}_plot.png"),
                f"{transcript_label_map.get(t, t)} - {primary_site_label} Codon Usage "
                f"(Frame0{', Merged' if args.merge_density else ''}) - {sample_name}",
                codon_label_order,
            )

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
        plot_codon_usage_dataframe(
            ov_df,
            os.path.join(codon_usage_dir, "codon_usage_total_plot.png"),
            f"Overall {primary_site_label} Codon Usage "
            f"(Frame0{', Merged' if args.merge_density else ''}) - {sample_name}",
            codon_label_order,
        )

        # ----------------------------------------
        # A-site Codon Usage (Per-transcript & Overall)
        # Note: For A-site usage, we do NOT mask the stop codon.
        # ----------------------------------------
        coverage_usage_a_site = {}
        for _, ann_ in resolved_annotation_df.iterrows():
            t = ann_["transcript"]
            sequence_id = ann_["resolved_sequence_name"]
            if sequence_id not in read_cov_map:
                continue
            seq_len = len(read_cov_map[sequence_id]["seq"])
            a_depth = np.array(read_cov_map[sequence_id]["A_site"])
            coverage_usage_a_site[t] = {}
            start_cdn = int(ann_["start_codon"])
            transcript_stop_codon = int(ann_["stop_codon"]) + 3
            if transcript_stop_codon > len(read_cov_map[sequence_id]["seq"]):
                transcript_stop_codon = len(read_cov_map[sequence_id]["seq"])
            upper_bound = min(transcript_stop_codon, seq_len - 2)
            seq_str = read_cov_map[sequence_id]["seq"]
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

        a_site_ov_cov = {}
        a_site_ov_freq = {}
        for t, usage_dict in coverage_usage_a_site.items():
            freq_dict = transcript_codon_counts.get(t, {})
            a_site_df_rows = build_codon_usage_dataframe(
                t,
                usage_dict,
                freq_dict,
                aggregate_cov=a_site_ov_cov,
                aggregate_freq=a_site_ov_freq,
            )
            a_site_transcript_csv = os.path.join(codon_usage_dir, f"a_site_codon_usage_{t}.csv")
            a_site_df_rows.to_csv(a_site_transcript_csv, index=False)
            plot_codon_usage_dataframe(
                a_site_df_rows,
                os.path.join(codon_usage_dir, f"a_site_codon_usage_{t}_plot.png"),
                f"{transcript_label_map.get(t, t)} - A-site Codon Usage "
                f"(Frame0{', Merged' if args.merge_density else ''}) - {sample_name}",
                codon_label_order,
            )

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
        plot_codon_usage_dataframe(
            a_site_df,
            os.path.join(codon_usage_dir, "a_site_codon_usage_total_plot.png"),
            f"Overall A-site Codon Usage (Frame0{', Merged' if args.merge_density else ''}) - {sample_name}",
            codon_label_order,
        )

    log_info(
        "PROFILE",
        f"Translation-profile analysis complete (merge_density={args.merge_density}).",
    )
