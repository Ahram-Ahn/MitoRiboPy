#!/usr/bin/env python3
"""Entry point for the mitochondrial ribosome profiling pipeline."""

from __future__ import annotations

import argparse
import os

from mitoribopy.io import (
    compute_total_counts,
    compute_unfiltered_read_length_summary,
    plot_unfiltered_read_length_heatmap,
    process_bed_files,
)
from mitoribopy.analysis import (
    create_csv_for_offset_enrichment,
    determine_p_site_offsets,
    plot_offset_enrichment,
    run_codon_correlation,
    run_inframe_codon_analysis,
    run_rna_seq_analysis,
)
from mitoribopy.data import human_annotation_df, yeast_annotation_df
from mitoribopy.config import DEFAULT_CONFIG, load_user_config, resolve_rpf_range
from mitoribopy.plotting import (
    apply_publication_style,
    run_igv_style_plot,
    run_varna_plot,
)


def configure_plot_defaults() -> None:
    apply_publication_style()


def build_parser(defaults: dict) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Ribo-seq pipeline with optional RNA-seq and VARNA output."
    )

    parser.add_argument("--config", default=None, help="Path to JSON config file.")
    parser.add_argument("-f", "--fasta", required=True, help="Path to FASTA file.")
    parser.add_argument(
        "-s",
        "--strain",
        choices=["y", "h"],
        default=defaults["strain"],
        help="Strain: y=yeast, h=human.",
    )
    parser.add_argument(
        "-d",
        "--directory",
        default=defaults["directory"],
        help="Directory with Ribo-seq BED files.",
    )
    parser.add_argument(
        "-rpf",
        nargs=2,
        type=int,
        default=defaults["rpf"],
        help="RPF length range. Example: -rpf 28 35",
    )
    parser.add_argument(
        "-a",
        "--align",
        choices=["start", "stop"],
        default=defaults["align"],
        help="Align to start or stop codon.",
    )
    parser.add_argument(
        "-r",
        "--range",
        type=int,
        default=defaults["range"],
        help="Offset plotting range (uses -range..+range).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=defaults["output"],
        help="Base output directory.",
    )
    parser.add_argument(
        "--downstream_dir",
        default=defaults["downstream_dir"],
        help="Subdirectory name for in-frame outputs.",
    )
    parser.add_argument(
        "--min_offset",
        type=int,
        default=defaults["min_offset"],
        help="Minimum absolute offset used to select P-site offsets.",
    )
    parser.add_argument(
        "--max_offset",
        type=int,
        default=defaults["max_offset"],
        help="Maximum absolute offset used to select P-site offsets.",
    )
    parser.add_argument(
        "--offset_pick_reference",
        choices=["selected_site", "p_site"],
        default=defaults["offset_pick_reference"],
        help=(
            "Reference space used when selecting most-enriched offsets: "
            "'selected_site' picks directly in the chosen offset_site space "
            "(legacy behavior), 'p_site' picks in canonical P-site space and "
            "then transforms to A-site if needed (recommended for stable P/A coupling)."
        ),
    )
    parser.add_argument(
        "--offset_type",
        choices=["5", "3"],
        default=defaults["offset_type"],
        help="Offset type for in-frame and IGV analysis.",
    )
    parser.add_argument(
        "--offset_site",
        choices=["p", "a"],
        default=defaults["offset_site"],
        help="Whether selected offsets represent P-site or A-site positions.",
    )
    parser.add_argument(
        "--codon_overlap_mode",
        choices=["full", "any"],
        default=defaults["codon_overlap_mode"],
        help="Use only full codon overlap (recommended) or any partial overlap for offset enrichment.",
    )
    parser.add_argument(
        "--plot_dir",
        default=defaults["plot_dir"],
        help="Subdirectory for offset plots and CSV files.",
    )
    parser.add_argument(
        "-fmt",
        "--plot_format",
        choices=["png", "pdf", "svg"],
        default=defaults["plot_format"],
        help="Plot format.",
    )
    parser.add_argument(
        "--x_breaks",
        nargs="+",
        type=int,
        default=defaults["x_breaks"],
        help="X-axis tick breaks for offset line plots.",
    )
    parser.add_argument(
        "--line_plot_style",
        choices=["combined", "separate"],
        default=defaults["line_plot_style"],
        help="Offset line plot style.",
    )
    parser.add_argument(
        "--cap_percentile",
        type=float,
        default=defaults["cap_percentile"],
        help="Coverage capping percentile.",
    )
    parser.add_argument(
        "-m",
        "--merge_density",
        action="store_true",
        default=defaults["merge_density"],
        help="Merge frame1 and frame2 into frame0 for codon usage.",
    )
    parser.add_argument(
        "-p",
        "--psite_offset",
        type=int,
        default=defaults["psite_offset"],
        help="Fixed offset for all reads (skips overlap-based offset logic).",
    )
    parser.add_argument(
        "--read_counts_file",
        default=defaults["read_counts_file"],
        help="Path to read counts summary CSV.",
    )
    parser.add_argument(
        "--read_counts_sample_col",
        default=defaults["read_counts_sample_col"],
        help="Sample column name in read-counts table (auto-detected when omitted).",
    )
    parser.add_argument(
        "--read_counts_reads_col",
        default=defaults["read_counts_reads_col"],
        help="Read-count column name in read-counts table (auto-detected when omitted).",
    )
    parser.add_argument(
        "--rpm_norm_mode",
        choices=["total", "mt_mrna"],
        default=defaults["rpm_norm_mode"],
        help=(
            "RPM denominator mode for IGV/unfiltered-RPM normalization: "
            "'total' uses all reads in count file; 'mt_mrna' uses only mt-mRNA references."
        ),
    )
    parser.add_argument(
        "--read_counts_reference_col",
        default=defaults["read_counts_reference_col"],
        help=(
            "Reference column name in read-counts table. Needed for "
            "--rpm_norm_mode mt_mrna when auto-detection is not sufficient."
        ),
    )
    parser.add_argument(
        "--mrna_ref_patterns",
        nargs="+",
        default=defaults["mrna_ref_patterns"],
        help=(
            "Substring pattern(s) used to identify mt-mRNA rows in read-counts table "
            "for --rpm_norm_mode mt_mrna."
        ),
    )
    parser.add_argument(
        "-v",
        "--varna",
        action="store_true",
        default=defaults["varna"],
        help="Generate VARNA color output.",
    )
    parser.add_argument(
        "--varna_norm_perc",
        type=float,
        default=defaults["varna_norm_perc"],
        help="Percentile for VARNA scaling.",
    )
    parser.add_argument(
        "--order_samples",
        nargs="+",
        default=defaults["order_samples"],
        help="Optional sample order for Ribo-seq outputs.",
    )
    parser.add_argument(
        "--cor_plot",
        action="store_true",
        default=defaults["cor_plot"],
        help="Generate codon-correlation plots.",
    )
    parser.add_argument(
        "--base_sample",
        default=defaults["base_sample"],
        help="Reference sample for codon-correlation plotting.",
    )
    parser.add_argument(
        "--cor_mask_method",
        choices=["percentile", "fixed", "none"],
        default="percentile",
        help="Outlier masking strategy for codon correlation.",
    )
    parser.add_argument(
        "--cor_mask_percentile",
        type=float,
        default=0.99,
        help="Percentile cutoff used when --cor_mask_method percentile.",
    )
    parser.add_argument(
        "--cor_mask_threshold",
        type=float,
        default=None,
        help="Fixed cutoff used when --cor_mask_method fixed.",
    )

    parser.add_argument(
        "--use_rna_seq",
        action="store_true",
        default=defaults["use_rna_seq"],
        help="Run RNA-seq analysis.",
    )
    parser.add_argument(
        "--rna_seq_dir",
        default=defaults["rna_seq_dir"],
        help="Directory containing RNA-seq BED files.",
    )
    parser.add_argument(
        "--rna_order",
        nargs="+",
        default=defaults["rna_order"],
        help="RNA sample order.",
    )
    parser.add_argument(
        "--rna_out_dir",
        default=defaults["rna_out_dir"],
        help="RNA-seq output directory name.",
    )
    parser.add_argument(
        "--do_rna_ribo_ratio",
        action="store_true",
        default=defaults["do_rna_ribo_ratio"],
        help="Merge Ribo-seq and RNA-seq outputs into ratios.",
    )

    return parser


def get_args() -> argparse.Namespace:
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", default=None)
    pre_args, _ = pre.parse_known_args()

    defaults = dict(DEFAULT_CONFIG)
    defaults["directory"] = os.getcwd()
    defaults.update(load_user_config(pre_args.config))

    parser = build_parser(defaults)
    args = parser.parse_args()

    if args.min_offset > args.max_offset:
        parser.error("--min_offset cannot be greater than --max_offset")
    if args.range <= 0:
        parser.error("--range must be a positive integer")
    if args.cap_percentile <= 0 or args.cap_percentile > 1:
        parser.error("--cap_percentile must be in (0, 1]")
    if args.cor_mask_method == "percentile" and (args.cor_mask_percentile <= 0 or args.cor_mask_percentile >= 1):
        parser.error("--cor_mask_percentile must be in (0, 1)")
    if args.use_rna_seq and (not args.rna_seq_dir or not args.rna_order):
        parser.error("--use_rna_seq requires both --rna_seq_dir and --rna_order")

    return args


def run_pipeline(args: argparse.Namespace) -> None:
    base_out = os.path.abspath(args.output)
    os.makedirs(base_out, exist_ok=True)

    plot_out_dir = os.path.join(base_out, args.plot_dir)
    os.makedirs(plot_out_dir, exist_ok=True)

    if args.strain == "y":
        annotation_df = yeast_annotation_df.copy()
    else:
        annotation_df = human_annotation_df.copy()
    rpf_range = resolve_rpf_range(args.strain, args.rpf)

    annotation_df["start_codon"] = annotation_df["l_utr5"]
    # 0-based first base of stop codon: transcript_length - utr3 - 3
    annotation_df["stop_codon"] = annotation_df["l_tr"] - annotation_df["l_utr3"] - 3

    read_counts_path = args.read_counts_file
    if not os.path.exists(read_counts_path):
        print(f"[MAIN] File not found: {read_counts_path}. Cannot compute total read counts.")
        total_counts_map = {}
    else:
        try:
            total_counts_map, _ = compute_total_counts(
                read_counts_path,
                sample_col=args.read_counts_sample_col,
                reads_col=args.read_counts_reads_col,
                normalization_mode=args.rpm_norm_mode,
                reference_col=args.read_counts_reference_col,
                mrna_ref_patterns=args.mrna_ref_patterns,
            )
        except Exception as exc:
            print(f"[MAIN] Failed to parse read-count file '{read_counts_path}': {exc}")
            total_counts_map = {}

    unfiltered_summary_csv = os.path.join(
        plot_out_dir, "unfiltered_read_length_summary_15_50.csv"
    )
    compute_unfiltered_read_length_summary(
        input_dir=args.directory,
        output_csv=unfiltered_summary_csv,
        total_counts_map=total_counts_map,
        read_length_range=(15, 50),
    )
    heatmap_base = os.path.join(plot_out_dir, "unfiltered_heatmap_15_50")
    plot_unfiltered_read_length_heatmap(
        summary_csv_path=unfiltered_summary_csv,
        output_png_base=heatmap_base,
        value_col="count",
    )
    plot_unfiltered_read_length_heatmap(
        summary_csv_path=unfiltered_summary_csv,
        output_png_base=heatmap_base,
        value_col="RPM",
    )

    filtered_bed_df, sample_dirs = process_bed_files(
        input_dir=args.directory,
        output_dir=plot_out_dir,
        organism=args.strain,
        annotation_df=annotation_df,
        rpf_range=rpf_range,
    )
    if filtered_bed_df.empty:
        print("[MAIN] No data after filtering. Pipeline finished early.")
        return

    if args.order_samples:
        name_to_dir = {os.path.basename(sd): sd for sd in sample_dirs}
        ordered_dirs = []
        missing_requested = []
        for requested in args.order_samples:
            if requested in name_to_dir:
                ordered_dirs.append(name_to_dir.pop(requested))
            else:
                missing_requested.append(requested)
        # Append any remaining samples deterministically.
        ordered_dirs.extend(
            sorted(name_to_dir.values(), key=lambda p: os.path.basename(p))
        )
        sample_dirs = ordered_dirs
        if missing_requested:
            print(
                "[MAIN] WARNING: Requested sample(s) not found in filtered BEDs: "
                + ", ".join(missing_requested)
            )

    offset_csv_path = os.path.join(plot_out_dir, f"offset_{args.align}.csv")
    summary_df, offsets_df = create_csv_for_offset_enrichment(
        bed_df=filtered_bed_df,
        annotation_df=annotation_df,
        align_to=args.align,
        rpf_range=rpf_range,
        output_csv=offset_csv_path,
        offset_limit=args.range,
        manual_offset=args.psite_offset,
        offset_site=args.offset_site,
        codon_overlap_mode=args.codon_overlap_mode,
        strain=args.strain,
    )
    if summary_df is None:
        print("[MAIN] No offsets generated. Pipeline finished early.")
        return

    p_site_offset_file = os.path.join(plot_out_dir, f"p_site_offsets_{args.align}.csv")
    p_site_offsets = determine_p_site_offsets(
        offsets_df=offsets_df,
        align_to=args.align,
        out_file=p_site_offset_file,
        offset_min=args.min_offset,
        offset_max=args.max_offset,
        offset_site=args.offset_site,
        selection_reference=args.offset_pick_reference,
    )

    plot_offset_enrichment(
        summary_df=summary_df,
        align_to=args.align,
        plot_dir=plot_out_dir,
        plot_format=args.plot_format,
        x_breaks=args.x_breaks,
        line_plot_style=args.line_plot_style,
        offset_limit=args.range,
        p_site_offsets=p_site_offsets,
        offset_min=args.min_offset,
        offset_max=args.max_offset,
    )

    if p_site_offsets is None:
        print("[MAIN] No P-site offsets. Skipping in-frame, IGV, and VARNA.")
    else:
        inframe_out_dir = base_out
        os.makedirs(inframe_out_dir, exist_ok=True)
        run_inframe_codon_analysis(
            sample_dirs=sample_dirs,
            a_site_offsets_df=p_site_offsets,
            offset_type=args.offset_type,
            fasta_file=args.fasta,
            output_dir=inframe_out_dir,
            args=args,
            annotation_df=annotation_df,
        )

        igv_out_dir = os.path.join(base_out, "igv_style_plots")
        os.makedirs(igv_out_dir, exist_ok=True)
        args.total_mrna_map = total_counts_map
        run_igv_style_plot(
            sample_dirs=sample_dirs,
            p_site_offsets_df=p_site_offsets,
            offset_type=args.offset_type,
            fasta_file=args.fasta,
            output_dir=igv_out_dir,
            args=args,
            annotation_df=annotation_df,
        )

        if args.varna:
            varna_out_dir = os.path.join(base_out, "varna")
            varna_site_col = "P_site" if args.offset_site == "p" else "A_site"
            run_varna_plot(
                inframe_out_dir=inframe_out_dir,
                varna_norm_perc=args.varna_norm_perc,
                varna_out_dir=varna_out_dir,
                varna_site_col=varna_site_col,
            )
        else:
            print("[MAIN] VARNA skipped. Use --varna to enable it.")

    if args.cor_plot and args.base_sample:
        sample_names = [os.path.basename(sd) for sd in sample_dirs]
        if args.base_sample not in sample_names:
            print(f"[MAIN] base_sample '{args.base_sample}' not present. Skipping correlation.")
        else:
            cor_out_dir = os.path.join(base_out, "codon_correlation")
            run_codon_correlation(
                inframe_analysis_dir=base_out,
                samples=sample_names,
                base_sample=args.base_sample,
                column="CoverageDivFreq",
                output_dir=cor_out_dir,
                mask_method=args.cor_mask_method,
                mask_percentile=args.cor_mask_percentile,
                mask_threshold=(args.cor_mask_threshold if args.cor_mask_method == "fixed" else None),
            )

    if args.use_rna_seq:
        rna_out_dir = os.path.join(base_out, args.rna_out_dir)
        os.makedirs(rna_out_dir, exist_ok=True)
        run_rna_seq_analysis(
            rna_seq_dir=args.rna_seq_dir,
            rna_order=args.rna_order,
            annotation_df=annotation_df,
            fasta_file=args.fasta,
            output_dir=rna_out_dir,
            inframe_out_dir=base_out,
            do_merge_with_ribo=args.do_rna_ribo_ratio,
            ribo_order=args.order_samples,
        )
    else:
        print("[MAIN] RNA-seq skipped.")


def main() -> None:
    configure_plot_defaults()
    args = get_args()
    run_pipeline(args)


if __name__ == "__main__":
    main()
