"""Standalone package-native pipeline runner."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Callable, Iterable

from ..cli.common import MitoRiboPyHelpFormatter
from ..console import configure_file_logging, log_message
from ..config import DEFAULT_CONFIG, load_user_config
from ..data import available_codon_table_names

StatusWriter = Callable[[str], None]


def configure_plot_defaults() -> None:
    """Apply global plotting defaults before pipeline execution."""
    from ..plotting import apply_publication_style

    apply_publication_style()


def build_parser(defaults: dict) -> argparse.ArgumentParser:
    """Build the standalone package CLI parser."""
    parser = argparse.ArgumentParser(
        prog="mitoribopy",
        description=(
            "Run the standalone MitoRiboPy Ribo-seq pipeline.\n"
            "This CLI filters BED inputs, estimates offsets, and then generates\n"
            "translation-profile, codon, coverage-profile, and optional\n"
            "RNA-seq/structure-density outputs."
        ),
        epilog=(
            "Examples:\n"
            "  mitoribopy -s h -f ref.fa --directory beds -rpf 29 34 --align stop \\\n"
            "    --offset_type 5 --offset_site p --offset_mask_nt 5 --output output -m\n"
            "  mitoribopy ... --min_5_offset 10 --max_5_offset 22 \\\n"
            "    --min_3_offset 12 --max_3_offset 30\n"
            "  mitoribopy -s custom -f ref.fa --directory beds --annotation_file ann.csv \\\n"
            "    --codon_tables_file codon_tables.json --codon_table_name standard\n"
            "\n"
            "Tip: prefer the end-specific 5'/3' offset bounds.\n"
            "Use --min_offset/--max_offset only as shared fallback bounds."
        ),
        formatter_class=MitoRiboPyHelpFormatter,
    )

    core_group = parser.add_argument_group("Core Inputs")
    offset_group = parser.add_argument_group("Offset Enrichment and Selection")
    normalization_group = parser.add_argument_group("Read-Count Normalization")
    output_group = parser.add_argument_group("Outputs and Plotting")
    optional_group = parser.add_argument_group("Optional Modules")

    core_group.add_argument(
        "--config",
        default=None,
        metavar="CONFIG.json",
        help="Optional JSON config file. CLI arguments override values from this file.",
    )
    core_group.add_argument(
        "-f",
        "--fasta",
        required=True,
        metavar="REF_FASTA",
        help="Reference FASTA file used to build the transcript/annotation context.",
    )
    core_group.add_argument(
        "-s",
        "--strain",
        choices=["y", "h", "vm", "ym", "custom"],
        default=defaults["strain"],
        help=(
            "Reference preset:\n"
            "  h      human mt (ships annotation + vertebrate_mitochondrial codons)\n"
            "  y      yeast mt (ships annotation + yeast_mitochondrial codons)\n"
            "  vm     any vertebrate mt (uses vertebrate_mitochondrial codons;\n"
            "         requires --annotation_file and an explicit -rpf range)\n"
            "  ym     any fungus with yeast-mito code (uses yeast_mitochondrial\n"
            "         codons; requires --annotation_file and explicit -rpf)\n"
            "  custom fully user-specified (annotation + codon table + -rpf)"
        ),
    )
    directory_action = core_group.add_argument(
        "-d",
        "--directory",
        default=defaults["directory"],
        metavar="BED_DIR",
        help=(
            "Directory containing Ribo-seq input files.\n"
            "Both .bed and .bam are accepted; BAM files are auto-converted\n"
            "to BED6 under <output>/bam_converted/ via pysam."
        ),
    )
    directory_action.default_display = "current working directory"
    core_group.add_argument(
        "--bam_mapq",
        type=int,
        default=defaults["bam_mapq"],
        metavar="Q",
        help=(
            "MAPQ threshold applied to BAM inputs before BAM->BED6 conversion.\n"
            "Set to 0 to disable. Default 10 is the same NUMT-suppression\n"
            "default used by 'mitoribopy align'."
        ),
    )
    rpf_action = core_group.add_argument(
        "-rpf",
        nargs=2,
        type=int,
        metavar=("MIN_LEN", "MAX_LEN"),
        default=defaults["rpf"],
        help="Inclusive read-length filter range, for example: -rpf 29 34",
    )
    rpf_action.default_display = (
        "monosome h/vm: 28-34, y/ym: 37-41; "
        "disome h/vm: 60-90, y/ym: 65-95"
    )
    footprint_action = core_group.add_argument(
        "--footprint_class",
        choices=["monosome", "disome", "custom"],
        default=defaults["footprint_class"],
        help=(
            "Expected RPF class, selects sensible --rpf and\n"
            "--unfiltered_read_length_range defaults when the user does\n"
            "not pass them explicitly:\n"
            "  monosome  single-ribosome footprint (default)\n"
            "  disome    collided-ribosome footprint (~60-90 nt, e.g.\n"
            "            eIF5A-depletion / queueing studies)\n"
            "  custom    no biological default; pass -rpf / \n"
            "            --unfiltered_read_length_range yourself"
        ),
    )
    footprint_action.default_display = "monosome"
    core_group.add_argument(
        "--annotation_file",
        default=defaults["annotation_file"],
        metavar="ANNOTATION.csv",
        help=(
            "Optional annotation CSV override.\n"
            "Required for --strain custom / vm / ym; built-in yeast (y)\n"
            "and human (h) tables are used otherwise."
        ),
    )
    core_group.add_argument(
        "--codon_tables_file",
        default=defaults["codon_tables_file"],
        metavar="CODON_TABLES.json",
        help=(
            "Optional codon-table JSON override.\n"
            "Supports either one flat 64-codon table or multiple named tables."
        ),
    )
    codon_table_action = core_group.add_argument(
        "--codon_table_name",
        default=defaults["codon_table_name"],
        metavar="TABLE_NAME",
        help=(
            "Codon-table name to load from built-ins or from --codon_tables_file.\n"
            f"Built-in names: {', '.join(available_codon_table_names())}"
        ),
    )
    codon_table_action.default_display = "y: yeast_mitochondrial, h: vertebrate_mitochondrial"
    start_codons_action = core_group.add_argument(
        "--start_codons",
        nargs="+",
        default=defaults["start_codons"],
        metavar="CODON",
        help=(
            "Allowed translation start codons used for codon classification.\n"
            "Defaults to strain presets; custom organisms can override them here."
        ),
    )
    start_codons_action.default_display = "y: ATG, h: ATG ATA, custom: ATG"
    core_group.add_argument(
        "--atp8_atp6_baseline",
        choices=["ATP6", "ATP8"],
        default=defaults["atp8_atp6_baseline"],
        help=(
            "Preferred baseline name when resolving the shared ATP8/ATP6 bicistronic transcript.\n"
            "Titles remain ATP8/ATP6."
        ),
    )
    core_group.add_argument(
        "--nd4l_nd4_baseline",
        choices=["ND4", "ND4L"],
        default=defaults["nd4l_nd4_baseline"],
        help=(
            "Preferred baseline name when resolving the shared ND4L/ND4 bicistronic transcript.\n"
            "Titles remain ND4L/ND4."
        ),
    )
    offset_group.add_argument(
        "-a",
        "--align",
        choices=["start", "stop"],
        default=defaults["align"],
        help="Anchor offset enrichment around the start or stop codon.",
    )
    offset_group.add_argument(
        "-r",
        "--range",
        type=int,
        default=defaults["range"],
        metavar="NT",
        help="Plot offsets from -range to +range around the chosen anchor codon.",
    )
    offset_group.add_argument(
        "--min_offset",
        type=int,
        default=defaults["min_offset"],
        metavar="NT",
        help=(
            "Backward-compatible shared minimum absolute offset used only when "
            "the 5'/3' end-specific bounds are not provided."
        ),
    )
    offset_group.add_argument(
        "--max_offset",
        type=int,
        default=defaults["max_offset"],
        metavar="NT",
        help=(
            "Backward-compatible shared maximum absolute offset used only when "
            "the 5'/3' end-specific bounds are not provided."
        ),
    )
    min_5_action = offset_group.add_argument(
        "--min_5_offset",
        type=int,
        default=defaults["min_5_offset"],
        metavar="NT",
        help="Recommended minimum absolute 5' offset considered during offset selection.",
    )
    min_5_action.default_display = "same as --min_offset"
    max_5_action = offset_group.add_argument(
        "--max_5_offset",
        type=int,
        default=defaults["max_5_offset"],
        metavar="NT",
        help="Recommended maximum absolute 5' offset considered during offset selection.",
    )
    max_5_action.default_display = "same as --max_offset"
    min_3_action = offset_group.add_argument(
        "--min_3_offset",
        type=int,
        default=defaults["min_3_offset"],
        metavar="NT",
        help="Recommended minimum absolute 3' offset considered during offset selection.",
    )
    min_3_action.default_display = "same as --min_offset"
    max_3_action = offset_group.add_argument(
        "--max_3_offset",
        type=int,
        default=defaults["max_3_offset"],
        metavar="NT",
        help="Recommended maximum absolute 3' offset considered during offset selection.",
    )
    max_3_action.default_display = "same as --max_offset"
    offset_group.add_argument(
        "--offset_mask_nt",
        type=int,
        default=defaults["offset_mask_nt"],
        metavar="NT",
        help=(
            "Mask near-anchor bins from -N..-1 and +1..+N in offset summaries,\n"
            "line plots, and heatmaps. The first visible bins are -(N+1) and +(N+1)."
        ),
    )
    offset_group.add_argument(
        "--offset_pick_reference",
        choices=["selected_site", "p_site"],
        default=defaults["offset_pick_reference"],
        help=(
            "How the best offset is chosen from enrichment tables:\n"
            "  p_site         choose in canonical P-site space first, then convert if needed\n"
            "  selected_site  choose directly in the final --offset_site space (legacy)"
        ),
    )
    offset_group.add_argument(
        "--offset_type",
        choices=["5", "3"],
        default=defaults["offset_type"],
        help=(
            "Which read end defines downstream site placement after offsets are selected:\n"
            "  5  measure offsets from the read 5' end\n"
            "  3  measure offsets from the read 3' end"
        ),
    )
    offset_group.add_argument(
        "--offset_site",
        choices=["p", "a"],
        default=defaults["offset_site"],
        help=(
            "Coordinate space for the SELECTED OFFSETS table. p = P-site,\n"
            "a = A-site. This controls the column values in\n"
            "p_site_offsets_<align>.csv only. To control which downstream\n"
            "outputs (codon usage, coverage plots) are generated, use\n"
            "--analysis_sites."
        ),
    )
    offset_group.add_argument(
        "--analysis_sites",
        choices=["p", "a", "both"],
        default=defaults.get("analysis_sites", "both"),
        help=(
            "Which downstream outputs to generate per sample.\n"
            "  both  (default) write both P-site and A-site codon usage and\n"
            "        coverage plots, side by side under per-site subdirs.\n"
            "  p     P-site outputs only.\n"
            "  a     A-site outputs only.\n"
            "Independent of --offset_site, which only controls the offset\n"
            "selection coordinate space."
        ),
    )
    offset_group.add_argument(
        "--codon_overlap_mode",
        choices=["full", "any"],
        default=defaults["codon_overlap_mode"],
        help=(
            "How reads count toward offset enrichment:\n"
            "  full  the read must span all 3 nt of the anchor codon\n"
            "  any   any partial codon overlap is enough"
        ),
    )
    offset_group.add_argument(
        "-p",
        "--psite_offset",
        type=int,
        default=defaults["psite_offset"],
        metavar="NT",
        help=(
            "Use one fixed offset for every read length and sample.\n"
            "This bypasses enrichment-based offset selection."
        ),
    )
    offset_group.add_argument(
        "--offset_mode",
        choices=["per_sample", "combined"],
        default=defaults.get("offset_mode", "per_sample"),
        help=(
            "How offsets drive the downstream translation-profile and "
            "coverage plots.\n"
            "  per_sample  (default) run enrichment + offset selection per\n"
            "              sample; each sample uses its own offsets for\n"
            "              downstream outputs. Offset drift across samples\n"
            "              is surfaced in offset_drift_<align>.svg.\n"
            "  combined    select one offset table from every sample pooled\n"
            "              together and apply it uniformly; matches the\n"
            "              v0.3.x behavior. Use when you have very low\n"
            "              coverage per sample and need the pooled signal."
        ),
    )
    output_group.add_argument(
        "-o",
        "--output",
        default=defaults["output"],
        metavar="OUTPUT_DIR",
        help="Base output directory for all pipeline results, including mitoribopy.log.",
    )
    output_group.add_argument(
        "--downstream_dir",
        default=defaults["downstream_dir"],
        help="Per-sample subdirectory name for frame and codon analyses.",
    )
    output_group.add_argument(
        "--plot_dir",
        default=defaults["plot_dir"],
        help="Shared subdirectory name for offset CSV files and plots.",
    )
    output_group.add_argument(
        "-fmt",
        "--plot_format",
        choices=["png", "pdf", "svg"],
        default=defaults["plot_format"],
        help="File format used for saved plots.",
    )
    output_group.add_argument(
        "--x_breaks",
        nargs="+",
        type=int,
        metavar="NT",
        default=defaults["x_breaks"],
        help="Optional custom x-axis tick marks for offset line plots.",
    )
    output_group.add_argument(
        "--line_plot_style",
        choices=["combined", "separate"],
        default=defaults["line_plot_style"],
        help="Draw offset line plots in one combined panel or separate 5'/3' panels.",
    )
    output_group.add_argument(
        "--cap_percentile",
        type=float,
        default=defaults["cap_percentile"],
        help="Upper percentile cap for coverage-style plots (0 < value <= 1).",
    )
    output_group.add_argument(
        "-m",
        "--merge_density",
        action="store_true",
        default=defaults["merge_density"],
        help="Collapse frame 1 and frame 2 into frame 0 for codon-density summaries.",
    )
    output_group.add_argument(
        "--order_samples",
        nargs="+",
        default=defaults["order_samples"],
        help="Optional sample order for Ribo-seq plots and aggregated outputs.",
    )
    normalization_group.add_argument(
        "--read_counts_file",
        default=defaults["read_counts_file"],
        metavar="COUNTS_TABLE",
        help="Read-count table for RPM normalization (.csv, .tsv, or .txt; delimiter auto-detected).",
    )
    normalization_group.add_argument(
        "--read_counts_sample_col",
        default=defaults["read_counts_sample_col"],
        help=(
            "Optional sample column override.\n"
            "If omitted, matching is case-insensitive and then falls back to column 1."
        ),
    )
    normalization_group.add_argument(
        "--read_counts_reads_col",
        default=defaults["read_counts_reads_col"],
        help=(
            "Optional read-count column override.\n"
            "If omitted, names like reads/read_count/counts are matched first, then column 3."
        ),
    )
    normalization_group.add_argument(
        "--unfiltered_read_length_range",
        nargs=2,
        type=int,
        metavar=("MIN_LEN", "MAX_LEN"),
        default=defaults["unfiltered_read_length_range"],
        help=(
            "Inclusive read-length range used for the unfiltered QC summary and heatmaps.\n"
            "Broaden this when you want to inspect longer footprints such as disomes."
        ),
    )
    normalization_group.add_argument(
        "--rpm_norm_mode",
        choices=["total", "mt_mrna"],
        default=defaults["rpm_norm_mode"],
        help=(
            "RPM denominator used for unfiltered and coverage-profile normalization:\n"
            "  total    sum all rows in the read-count table\n"
            "  mt_mrna  sum only rows whose reference matches --mrna_ref_patterns"
        ),
    )
    normalization_group.add_argument(
        "--read_counts_reference_col",
        default=defaults["read_counts_reference_col"],
        help=(
            "Optional reference column override.\n"
            "Needed for --rpm_norm_mode mt_mrna if auto-detection fails; otherwise falls back to column 2."
        ),
    )
    normalization_group.add_argument(
        "--mrna_ref_patterns",
        nargs="+",
        metavar="PATTERN",
        default=defaults["mrna_ref_patterns"],
        help=(
            "Substring pattern(s) used to identify mt-mRNA rows when\n"
            "--rpm_norm_mode mt_mrna is selected."
        ),
    )
    optional_group.add_argument(
        "--structure_density",
        action="store_true",
        default=defaults["structure_density"],
        help="Generate structure-density exports from footprint-density tables.",
    )
    optional_group.add_argument(
        "--structure_density_norm_perc",
        type=float,
        default=defaults["structure_density_norm_perc"],
        help="Upper percentile used to cap and scale structure-density values.",
    )
    optional_group.add_argument(
        "--cor_plot",
        action="store_true",
        default=defaults["cor_plot"],
        help="Generate codon-correlation plots.",
    )
    optional_group.add_argument(
        "--base_sample",
        default=defaults["base_sample"],
        help="Reference sample for codon-correlation comparisons.",
    )
    optional_group.add_argument(
        "--cor_mask_method",
        choices=["percentile", "fixed", "none"],
        default=defaults["cor_mask_method"],
        help="Masking rule for extreme codon-correlation outliers.",
    )
    optional_group.add_argument(
        "--cor_mask_percentile",
        type=float,
        default=defaults["cor_mask_percentile"],
        help="Percentile cutoff used when --cor_mask_method percentile.",
    )
    optional_group.add_argument(
        "--cor_mask_threshold",
        type=float,
        default=defaults["cor_mask_threshold"],
        help="Fixed absolute cutoff used when --cor_mask_method fixed.",
    )
    optional_group.add_argument(
        "--use_rna_seq",
        action="store_true",
        default=defaults["use_rna_seq"],
        help=(
            "[DEPRECATED in v0.3.0; removed in v0.4.0] Run the legacy\n"
            "RNA-seq ratio module. Use the dedicated 'mitoribopy rnaseq'\n"
            "subcommand instead for DE-based TE and delta-TE analysis\n"
            "with a SHA256 reference-consistency gate."
        ),
    )
    optional_group.add_argument(
        "--rna_seq_dir",
        default=defaults["rna_seq_dir"],
        help="Directory containing RNA-seq BED files.",
    )
    optional_group.add_argument(
        "--rna_order",
        nargs="+",
        default=defaults["rna_order"],
        help="RNA-seq sample order for paired outputs and plots.",
    )
    optional_group.add_argument(
        "--rna_out_dir",
        default=defaults["rna_out_dir"],
        help="RNA-seq output directory name.",
    )
    optional_group.add_argument(
        "--do_rna_ribo_ratio",
        action="store_true",
        default=defaults["do_rna_ribo_ratio"],
        help="Merge Ribo-seq and RNA-seq outputs into ratios.",
    )

    return parser


def parse_pipeline_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    """Parse and validate pipeline CLI arguments."""
    argv_list = None if argv is None else list(argv)

    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument("--config", default=None)
    pre_args, _ = pre_parser.parse_known_args(argv_list)

    defaults = dict(DEFAULT_CONFIG)
    defaults["directory"] = os.getcwd()
    defaults.update(load_user_config(pre_args.config))

    parser = build_parser(defaults)
    args = parser.parse_args(argv_list)

    if args.min_5_offset is None:
        args.min_5_offset = args.min_offset
    if args.max_5_offset is None:
        args.max_5_offset = args.max_offset
    if args.min_3_offset is None:
        args.min_3_offset = args.min_offset
    if args.max_3_offset is None:
        args.max_3_offset = args.max_offset

    if args.min_offset > args.max_offset:
        parser.error("--min_offset cannot be greater than --max_offset")
    if args.min_5_offset > args.max_5_offset:
        parser.error("--min_5_offset cannot be greater than --max_5_offset")
    if args.min_3_offset > args.max_3_offset:
        parser.error("--min_3_offset cannot be greater than --max_3_offset")
    if args.offset_mask_nt < 0:
        parser.error("--offset_mask_nt must be zero or a positive integer")
    if args.range <= 0:
        parser.error("--range must be a positive integer")
    if args.cap_percentile <= 0 or args.cap_percentile > 1:
        parser.error("--cap_percentile must be in (0, 1]")
    if len(args.unfiltered_read_length_range) != 2:
        parser.error("--unfiltered_read_length_range requires MIN_LEN MAX_LEN")
    if int(args.unfiltered_read_length_range[0]) > int(args.unfiltered_read_length_range[1]):
        parser.error(
            "--unfiltered_read_length_range MIN_LEN cannot be greater than MAX_LEN"
        )
    if args.structure_density_norm_perc <= 0 or args.structure_density_norm_perc > 1:
        parser.error("--structure_density_norm_perc must be in (0, 1]")
    if args.cor_mask_method == "percentile" and (
        args.cor_mask_percentile <= 0 or args.cor_mask_percentile >= 1
    ):
        parser.error("--cor_mask_percentile must be in (0, 1)")
    if args.use_rna_seq and (not args.rna_seq_dir or not args.rna_order):
        parser.error("--use_rna_seq requires both --rna_seq_dir and --rna_order")
    # Strain-preset requirements. Order of checks matches how a user
    # would experience the failures.
    if args.strain == "custom":
        if not args.annotation_file:
            parser.error("--strain custom requires --annotation_file")
        if args.rpf is None:
            parser.error(
                "--strain custom requires an explicit -rpf MIN_LEN MAX_LEN range"
            )
        if not (args.codon_tables_file or args.codon_table_name):
            parser.error(
                "--strain custom requires --codon_tables_file or --codon_table_name"
            )
    if args.strain in {"vm", "ym"}:
        if not args.annotation_file:
            parser.error(
                f"--strain {args.strain} requires --annotation_file "
                "(only h and y ship a built-in annotation)"
            )
        if args.rpf is None:
            parser.error(
                f"--strain {args.strain} requires an explicit -rpf MIN_LEN MAX_LEN "
                "range (only h and y have a built-in default)"
            )

    # --footprint_class=custom requires an explicit -rpf even for h/y,
    # because the whole point of 'custom' is "I know my footprint class,
    # don't pick one for me".
    if args.footprint_class == "custom" and args.rpf is None and args.strain not in {"h", "y"}:
        parser.error(
            "--footprint_class custom requires an explicit -rpf MIN_LEN MAX_LEN range"
        )

    # Inject footprint-class unfiltered-length defaults when the user did
    # not override --unfiltered_read_length_range. We compare to the
    # monosome default because that is what DEFAULT_CONFIG ships; if the
    # user passed their own range, it is already different and we leave
    # it alone.
    from ..config.runtime import FOOTPRINT_CLASS_DEFAULTS  # local import to avoid cycle
    class_defaults = FOOTPRINT_CLASS_DEFAULTS.get(args.footprint_class, {})
    user_kept_default = (
        list(args.unfiltered_read_length_range)
        == list(DEFAULT_CONFIG["unfiltered_read_length_range"])
    )
    if user_kept_default and "unfiltered_read_length_range" in class_defaults:
        args.unfiltered_read_length_range = list(
            class_defaults["unfiltered_read_length_range"]
        )

    return args


def run_pipeline(
    args: argparse.Namespace,
    *,
    emit_status: StatusWriter | None = None,
) -> int:
    """Execute the standalone package-native pipeline."""
    from .steps import (
        TOTAL_PIPELINE_STEPS,
        build_pipeline_context,
        compute_offset_enrichment_step,
        filter_bed_inputs,
        load_total_read_counts,
        run_downstream_modules,
        run_unfiltered_read_length_qc,
        select_offsets_and_plot,
    )

    status_writer = log_message if emit_status is None else emit_status

    context = build_pipeline_context(args)
    status_writer(
        "[PIPELINE] Step 1/"
        f"{TOTAL_PIPELINE_STEPS} OK: initialized outputs, annotations, and RPF range."
    )

    load_total_read_counts(context, status_writer)
    run_unfiltered_read_length_qc(context, status_writer)

    if not filter_bed_inputs(context, status_writer):
        return 0
    if not compute_offset_enrichment_step(context, status_writer):
        return 0

    select_offsets_and_plot(context, status_writer)
    run_downstream_modules(context, status_writer)
    return 0


def run_pipeline_cli(argv: Iterable[str] | None = None) -> int:
    """Parse arguments and run the standalone package-native pipeline."""
    import sys as _sys

    args = parse_pipeline_args(argv)
    if getattr(args, "use_rna_seq", False):
        _sys.stderr.write(
            "[mitoribopy rpf] DEPRECATION: --use_rna_seq is deprecated in "
            "v0.3.0 and will be removed in v0.4.0. Use 'mitoribopy rnaseq' "
            "instead - it enforces a SHA256 reference-consistency gate and "
            "runs on DE output (DESeq2 / Xtail / Anota2Seq).\n"
        )
    log_path = configure_file_logging(Path(args.output) / "mitoribopy.log")
    log_message(f"[PIPELINE] Log file: {log_path}")
    configure_plot_defaults()
    return run_pipeline(args)


def main(argv: Iterable[str] | None = None) -> int:
    """Compatibility entrypoint used by both package and wrapper scripts."""
    return run_pipeline_cli(argv)
