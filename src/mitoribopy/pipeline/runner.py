"""Standalone package-native pipeline runner."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Callable, Iterable

from ..cli.common import MitoRiboPyHelpFormatter
from ..console import configure_file_logging, log_message
from ..config import DEFAULT_CONFIG, load_user_config
from ..data import available_codon_table_names
from ..progress import Stopwatch, format_duration, render_step_timeline

StatusWriter = Callable[[str], None]


def configure_plot_defaults() -> None:
    """Apply global plotting defaults before pipeline execution."""
    from ..plotting import apply_publication_style

    apply_publication_style()


def build_parser(defaults: dict) -> argparse.ArgumentParser:
    """Build the standalone package CLI parser."""
    parser = argparse.ArgumentParser(
        prog="mitoribopy rpf",
        description=(
            "Run the MitoRiboPy Ribo-seq analysis stage on BED / BAM inputs.\n"
            "This subcommand filters reads, estimates P-site / A-site offsets,\n"
            "and then generates translation-profile, codon, coverage-profile,\n"
            "and optional RNA-seq / structure-density outputs."
        ),
        epilog=(
            "Examples:\n"
            "  mitoribopy rpf --strain h.sapiens --fasta ref.fa --directory beds \\\n"
            "                 -rpf 29 34 --align stop --offset-type 5 \\\n"
            "                 --offset-site p --offset-mask-nt 5 --output out \\\n"
            "                 --codon-density-window\n"
            "  mitoribopy rpf ... --min-5-offset 10 --max-5-offset 22 \\\n"
            "                     --min-3-offset 12 --max-3-offset 30\n"
            "  mitoribopy rpf --strain custom --fasta ref.fa --directory beds \\\n"
            "                 --annotation-file ann.csv \\\n"
            "                 --codon-tables-file codon_tables.json \\\n"
            "                 --codon-table-name standard\n"
            "\n"
            "Tip: prefer the end-specific 5'/3' offset bounds.\n"
            "Use --min-offset / --max-offset only as shared fallback bounds.\n"
            "Underscore-style flags (--offset_type, --min_5_offset, ...) are\n"
            "still accepted for one transition cycle but no longer shown here."
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
        metavar="CONFIG",
        help=(
            "Optional config file (.yaml, .yml, .json, or .toml; format "
            "auto-detected by extension). CLI arguments override values "
            "from this file."
        ),
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
        choices=["h.sapiens", "s.cerevisiae", "custom", "h", "y"],
        default=defaults["strain"],
        metavar="STRAIN",
        help=(
            "Reference preset (organism). Built-in strains ship a complete\n"
            "annotation table and codon table, so the user only needs to\n"
            "supply the FASTA reference. Other organisms must use\n"
            "'custom' and supply their own annotation + codon-table\n"
            "choice (see --annotation_file / --codon_table_name).\n"
            "  h.sapiens     human mt (default; vertebrate_mitochondrial\n"
            "                codon table; ATG/ATA start codons)\n"
            "  s.cerevisiae  budding yeast mt (yeast_mitochondrial codon\n"
            "                table; ATG start codon)\n"
            "  custom        any other organism (mouse, fly, A. thaliana,\n"
            "                fungi, ...); requires --annotation_file plus\n"
            "                a --codon_table_name picked from the\n"
            "                built-in NCBI Genetic Codes list (see\n"
            "                --codon_table_name help) or a\n"
            "                --codon_tables_file you supply.\n"
            "  h, y          deprecated short aliases for h.sapiens and\n"
            "                s.cerevisiae; still accepted for one cycle."
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
        "--bam-mapq",
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
        "short    h.sapiens / s.cerevisiae: 16-24; "
        "monosome h.sapiens: 28-34, s.cerevisiae: 37-41; "
        "disome   h.sapiens: 50-70, s.cerevisiae: 60-90"
    )
    footprint_action = core_group.add_argument(
        "--footprint-class",
        "--footprint_class",
        choices=["short", "monosome", "disome", "custom"],
        default=defaults["footprint_class"],
        help=(
            "Expected ribosome-protected fragment class. Selects the\n"
            "default --rpf and --unfiltered_read_length_range windows for\n"
            "the chosen --strain when the user does not pass them\n"
            "explicitly. Built-in defaults exist for h.sapiens and\n"
            "s.cerevisiae; for --strain custom you must also pass -rpf.\n"
            "  short     truncated RNase products (~16-24 nt). Sit just\n"
            "            below the canonical monosome window; useful for\n"
            "            mapping context-dependent pausing and as a QC\n"
            "            indicator of digest aggressiveness.\n"
            "  monosome  single-ribosome footprint (default).\n"
            "            h.sapiens 28-34 nt, s.cerevisiae 37-41 nt.\n"
            "  disome    collided-ribosome footprint. h.sapiens 50-70 nt,\n"
            "            s.cerevisiae 60-90 nt. Used for ribosome-stalling /\n"
            "            queueing studies (e.g. eIF5A-depletion).\n"
            "  custom    no biological default; supply -rpf and\n"
            "            --unfiltered_read_length_range yourself."
        ),
    )
    footprint_action.default_display = "monosome"
    core_group.add_argument(
        "--annotation-file",
        "--annotation_file",
        default=defaults["annotation_file"],
        metavar="ANNOTATION.csv",
        help=(
            "Per-transcript annotation CSV. Required when --strain\n"
            "custom; the built-in human (h.sapiens) and yeast\n"
            "(s.cerevisiae) tables are used otherwise. The CSV must have\n"
            "one row per transcript with the columns: transcript,\n"
            "sequence_name, display_name, sequence_aliases, l_utr5,\n"
            "l_cds, l_utr3, start_codon, stop_codon. See the 'Custom\n"
            "organisms' section of the README for a worked example."
        ),
    )
    core_group.add_argument(
        "--codon-tables-file",
        "--codon_tables_file",
        default=defaults["codon_tables_file"],
        metavar="CODON_TABLES.json",
        help=(
            "Optional path to a codon-table JSON file. Supports either\n"
            "(a) ONE flat 64-codon mapping {codon: amino_acid}, or\n"
            "(b) multiple named tables {table_name: {codon: amino_acid}}.\n"
            "Combine with --codon_table_name when (b). Use this only when\n"
            "your organism's genetic code is not already in the built-in\n"
            "list — see --codon_table_name help."
        ),
    )
    codon_table_action = core_group.add_argument(
        "--codon-table-name",
        "--codon_table_name",
        default=defaults["codon_table_name"],
        metavar="TABLE_NAME",
        help=(
            "Codon-table name to load from built-ins or from\n"
            "--codon_tables_file. The built-in tables are the NCBI\n"
            "Genetic Codes (ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)\n"
            "with NCBI numbers converted to descriptive names. Pick the\n"
            "one matching your organism's MITOCHONDRIAL code (or NUCLEAR\n"
            "code for ciliate / candida / etc.):\n"
            "  vertebrate_mitochondrial   (NCBI #2)  mouse, fish, frog,\n"
            "                                       any non-human vertebrate\n"
            "  yeast_mitochondrial        (NCBI #3)  S. cerevisiae and\n"
            "                                       close relatives\n"
            "  mold_mitochondrial         (NCBI #4)  Neurospora, Aspergillus,\n"
            "                                       Trichoderma, mycoplasmas\n"
            "  invertebrate_mitochondrial (NCBI #5)  fruit fly, mosquito,\n"
            "                                       nematodes\n"
            "  echinoderm_mitochondrial   (NCBI #9)  sea urchin, starfish\n"
            "  ascidian_mitochondrial     (NCBI #13) sea squirts\n"
            "  alternative_flatworm_mitochondrial (NCBI #14) flatworms\n"
            "  trematode_mitochondrial    (NCBI #21) trematodes\n"
            "  standard                   (NCBI #1)  use this when the\n"
            "                                       organism uses the\n"
            "                                       standard genetic code\n"
            "                                       (most plants, including\n"
            "                                       A. thaliana mt + plastid)\n"
            f"Full list of {len(available_codon_table_names())} built-in tables: "
            f"{', '.join(available_codon_table_names())}."
        ),
    )
    codon_table_action.default_display = (
        "h.sapiens: vertebrate_mitochondrial, "
        "s.cerevisiae: yeast_mitochondrial"
    )
    start_codons_action = core_group.add_argument(
        "--start-codons",
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
        "--atp8-atp6-baseline",
        "--atp8_atp6_baseline",
        choices=["ATP6", "ATP8"],
        default=defaults["atp8_atp6_baseline"],
        help=(
            "Preferred baseline name when resolving the shared ATP8/ATP6 bicistronic transcript.\n"
            "Titles remain ATP8/ATP6."
        ),
    )
    core_group.add_argument(
        "--nd4l-nd4-baseline",
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
        "--min-offset",
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
        "--max-offset",
        "--max_offset",
        type=int,
        default=defaults["max_offset"],
        metavar="NT",
        help=(
            "Backward-compatible shared maximum absolute offset used only when "
            "the 5'/3' end-specific bounds are not provided."
        ),
    )
    offset_group.add_argument(
        "--rpf-min-count-frac",
        "--rpf_min_count_frac",
        type=float,
        default=defaults["rpf_min_count_frac"],
        metavar="FRAC",
        help=(
            "Drop read-length bins from the RPF window whose total count\n"
            "across all samples is below FRAC x the most-enriched length.\n"
            "FRAC=0 disables the filter. Default 0.20 keeps only read\n"
            "lengths with >=20%% of the dominant-length count, so noisy\n"
            "low-count bins do not pollute offset selection. Pair with a\n"
            "wide --rpf range (e.g. 27 36 for human) to let the data\n"
            "decide which lengths actually carry signal."
        ),
    )
    min_5_action = offset_group.add_argument(
        "--min-5-offset",
        "--min_5_offset",
        type=int,
        default=defaults["min_5_offset"],
        metavar="NT",
        help="Recommended minimum absolute 5' offset considered during offset selection.",
    )
    min_5_action.default_display = "same as --min_offset"
    max_5_action = offset_group.add_argument(
        "--max-5-offset",
        "--max_5_offset",
        type=int,
        default=defaults["max_5_offset"],
        metavar="NT",
        help="Recommended maximum absolute 5' offset considered during offset selection.",
    )
    max_5_action.default_display = "same as --max_offset"
    min_3_action = offset_group.add_argument(
        "--min-3-offset",
        "--min_3_offset",
        type=int,
        default=defaults["min_3_offset"],
        metavar="NT",
        help="Recommended minimum absolute 3' offset considered during offset selection.",
    )
    min_3_action.default_display = "same as --min_offset"
    max_3_action = offset_group.add_argument(
        "--max-3-offset",
        "--max_3_offset",
        type=int,
        default=defaults["max_3_offset"],
        metavar="NT",
        help="Recommended maximum absolute 3' offset considered during offset selection.",
    )
    max_3_action.default_display = "same as --max_offset"
    offset_group.add_argument(
        "--offset-mask-nt",
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
        "--offset-pick-reference",
        "--offset_pick_reference",
        choices=["reported_site", "p_site", "selected_site"],
        default=defaults["offset_pick_reference"],
        help=(
            "Which coordinate space the best offset is chosen in.\n"
            "  p_site         (default) pick the offset in canonical P-site\n"
            "                 space first, then convert into the space\n"
            "                 named by --offset_site if that is 'a'.\n"
            "                 Use this when comparing across samples with\n"
            "                 different offset_site choices.\n"
            "  reported_site  pick the offset directly in the same space\n"
            "                 named by --offset_site (no P<->A conversion).\n"
            "                 Use this when you want what you see in the\n"
            "                 enrichment table for --offset_site to be the\n"
            "                 exact value picked.\n"
            "  selected_site  DEPRECATED alias for 'reported_site'."
        ),
    )
    offset_group.add_argument(
        "--offset-type",
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
        "--offset-site",
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
        "--analysis-sites",
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
        "--codon-overlap-mode",
        "--codon_overlap_mode",
        choices=["full", "any"],
        default=defaults["codon_overlap_mode"],
        help=(
            "How reads count toward the codon-level offset-enrichment\n"
            "table at each anchor codon (start or stop, per --align).\n"
            "  full  (default) the read must span ALL 3 nt of the anchor\n"
            "        codon. Example: anchor codon at positions 101-103;\n"
            "        a 30-nt read at 100-129 counts (read covers 101, 102,\n"
            "        103); a read at 102-131 does NOT count (does not\n"
            "        cover position 101). Strict and the right default for\n"
            "        ribo-seq footprints (>= ~26 nt).\n"
            "  any   any 1+ nt overlap with the anchor codon counts. The\n"
            "        102-131 read above WOULD count. Use only for very\n"
            "        short reads relative to a codon (rare for ribo-seq)."
        ),
    )
    offset_group.add_argument(
        "-p",
        "--psite-offset",
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
        "--offset-mode",
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
        "--downstream-dir",
        "--downstream_dir",
        default=defaults["downstream_dir"],
        help="Per-sample subdirectory name for frame and codon analyses.",
    )
    output_group.add_argument(
        "--plot-dir",
        "--plot_dir",
        default=defaults["plot_dir"],
        help="Shared subdirectory name for offset CSV files and plots.",
    )
    output_group.add_argument(
        "-fmt",
        "--plot-format",
        "--plot_format",
        choices=["png", "pdf", "svg"],
        default=defaults["plot_format"],
        help="File format used for saved plots.",
    )
    output_group.add_argument(
        "--x-breaks",
        "--x_breaks",
        nargs="+",
        type=int,
        metavar="NT",
        default=defaults["x_breaks"],
        help="Optional custom x-axis tick marks for offset line plots.",
    )
    output_group.add_argument(
        "--line-plot-style",
        "--line_plot_style",
        choices=["combined", "separate"],
        default=defaults["line_plot_style"],
        help="Draw offset line plots in one combined panel or separate 5'/3' panels.",
    )
    output_group.add_argument(
        "--cap-percentile",
        "--cap_percentile",
        type=float,
        default=defaults["cap_percentile"],
        help="Upper percentile cap for coverage-style plots (0 < value <= 1).",
    )
    output_group.add_argument(
        "-m",
        "--codon-density-window",
        "--codon_density_window",
        "--merge_density",
        dest="codon_density_window",
        action="store_true",
        default=defaults["codon_density_window"],
        help=(
            "When set, the codon-level density at each codon centre is "
            "summed with its +/-1 nt neighbours (a 3-nt sliding window) "
            "before being written into the codon coverage / usage tables "
            "and plots. Smooths short-window noise around the codon "
            "centre; does NOT collapse reading frames despite the "
            "historical flag name '--merge_density' (which is kept as a "
            "deprecated alias)."
        ),
    )
    output_group.add_argument(
        "--order-samples",
        "--order_samples",
        nargs="+",
        default=defaults["order_samples"],
        help="Optional sample order for Ribo-seq plots and aggregated outputs.",
    )
    normalization_group.add_argument(
        "--read-counts-file",
        "--read_counts_file",
        default=defaults["read_counts_file"],
        metavar="COUNTS_TABLE",
        help="Read-count table for RPM normalization (.csv, .tsv, or .txt; delimiter auto-detected).",
    )
    normalization_group.add_argument(
        "--read-counts-sample-col",
        "--read_counts_sample_col",
        default=defaults["read_counts_sample_col"],
        help=(
            "Optional sample column override.\n"
            "If omitted, matching is case-insensitive and then falls back to column 1."
        ),
    )
    normalization_group.add_argument(
        "--read-counts-reads-col",
        "--read_counts_reads_col",
        default=defaults["read_counts_reads_col"],
        help=(
            "Optional read-count column override.\n"
            "If omitted, names like reads/read_count/counts are matched first, then column 3."
        ),
    )
    normalization_group.add_argument(
        "--unfiltered-read-length-range",
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
        "--rpm-norm-mode",
        "--rpm_norm_mode",
        choices=["total", "mt_mrna"],
        default=defaults["rpm_norm_mode"],
        help=(
            "RPM denominator used for unfiltered and coverage-profile normalization:\n"
            "  total    sum all rows in the read-count table\n"
            "  mt_mrna  sum only rows whose reference matches\n"
            "           --mt_mrna_substring_patterns"
        ),
    )
    normalization_group.add_argument(
        "--read-counts-reference-col",
        "--read_counts_reference_col",
        default=defaults["read_counts_reference_col"],
        help=(
            "Optional reference column override.\n"
            "Needed for --rpm_norm_mode mt_mrna if auto-detection fails; otherwise falls back to column 2."
        ),
    )
    normalization_group.add_argument(
        "--mt-mrna-substring-patterns",
        "--mt_mrna_substring_patterns",
        "--mrna_ref_patterns",
        dest="mt_mrna_substring_patterns",
        nargs="+",
        metavar="PATTERN",
        default=defaults["mt_mrna_substring_patterns"],
        help=(
            "When --rpm_norm_mode mt_mrna is selected, the RPM denominator\n"
            "uses only rows from the read-count table whose value in the\n"
            "reference column (set via --read_counts_reference_col, or\n"
            "auto-detected) contains ANY of these substrings. Default\n"
            "matches reference names like 'mt_genome', 'mt-mrna', and\n"
            "'mt_mrna'. Example: in a read-count table with rows for\n"
            "'rrna', 'trna', 'mt_mrna_ND1', and 'mt_mrna_COX1', the\n"
            "default pattern keeps the latter two and drops the\n"
            "(r/t)RNA rows. The legacy flag name '--mrna_ref_patterns'\n"
            "is kept as a deprecated alias."
        ),
    )
    optional_group.add_argument(
        "--structure-density",
        "--structure_density",
        action="store_true",
        default=defaults["structure_density"],
        help="Generate structure-density exports from footprint-density tables.",
    )
    optional_group.add_argument(
        "--structure-density-norm-perc",
        "--structure_density_norm_perc",
        type=float,
        default=defaults["structure_density_norm_perc"],
        help="Upper percentile used to cap and scale structure-density values.",
    )
    optional_group.add_argument(
        "--cor-plot",
        "--cor_plot",
        action="store_true",
        default=defaults["cor_plot"],
        help="Generate codon-correlation plots.",
    )
    optional_group.add_argument(
        "--base-sample",
        "--base_sample",
        default=defaults["base_sample"],
        help="Reference sample for codon-correlation comparisons.",
    )
    optional_group.add_argument(
        "--cor-mask-method",
        "--cor_mask_method",
        choices=["percentile", "fixed", "none"],
        default=defaults["cor_mask_method"],
        help="Masking rule for extreme codon-correlation outliers.",
    )
    optional_group.add_argument(
        "--cor-mask-percentile",
        "--cor_mask_percentile",
        type=float,
        default=defaults["cor_mask_percentile"],
        help="Percentile cutoff used when --cor_mask_method percentile.",
    )
    optional_group.add_argument(
        "--cor-mask-threshold",
        "--cor_mask_threshold",
        type=float,
        default=defaults["cor_mask_threshold"],
        help="Fixed absolute cutoff used when --cor_mask_method fixed.",
    )
    optional_group.add_argument(
        "--cor-metric",
        "--cor_metric",
        choices=["log2_density_rpm", "log2_rpm", "linear", "raw_count"],
        default=defaults["cor_metric"],
        help=(
            "Primary codon-correlation metric. log2_density_rpm (default) "
            "log-transforms RPM-normalised codon density; raw_count "
            "reproduces the legacy depth-dominated scatter (and emits "
            "W_CODON_RAW_COUNT_PRIMARY)."
        ),
    )
    optional_group.add_argument(
        "--cor-regression",
        "--cor_regression",
        choices=["theil_sen", "ols", "none"],
        default=defaults["cor_regression"],
        help=(
            "Robust regression method drawn through the codon-correlation "
            "scatter. theil_sen (default) is median-based and tolerant of "
            "outliers; ols reproduces the legacy linear fit."
        ),
    )
    optional_group.add_argument(
        "--cor-support-min-raw",
        "--cor_support_min_raw",
        type=int,
        default=defaults["cor_support_min_raw"],
        help=(
            "Minimum raw value required in BOTH samples for a codon to "
            "be marked include_primary. Low-support codons remain in the "
            "metrics TSV but are dimmed in the figure and excluded from "
            "the residual ranking."
        ),
    )
    optional_group.add_argument(
        "--cor-label-top-n",
        "--cor_label_top_n",
        type=int,
        default=defaults["cor_label_top_n"],
        help="Number of codons to label on the scatter (by support-aware label score).",
    )
    optional_group.add_argument(
        "--cor-pseudocount",
        "--cor_pseudocount",
        default=defaults["cor_pseudocount"],
        help="Additive pseudocount before log2. 'auto' or float.",
    )
    optional_group.add_argument(
        "--cor-raw-panel",
        "--cor_raw_panel",
        choices=["qc_only", "off"],
        default=defaults["cor_raw_panel"],
        help=(
            "Where to write the raw-count panel when metric=raw_count. "
            "'qc_only' (default) routes it to raw_count_qc/; 'off' skips it."
        ),
    )
    optional_group.add_argument(
        "--read-coverage-raw",
        "--read_coverage_raw",
        action=argparse.BooleanOptionalAction,
        default=defaults["read_coverage_raw"],
        help=(
            "Write read-coverage plots in raw counts under "
            "coverage_profile_plots/read_coverage_raw[_codon]/. "
            "Use --no-read_coverage_raw to skip."
        ),
    )
    optional_group.add_argument(
        "--read-coverage-rpm",
        "--read_coverage_rpm",
        action=argparse.BooleanOptionalAction,
        default=defaults["read_coverage_rpm"],
        help=(
            "Write read-coverage plots in RPM under "
            "coverage_profile_plots/read_coverage_rpm[_codon]/. "
            "Use --no-read_coverage_rpm to skip."
        ),
    )
    optional_group.add_argument(
        "--igv-export",
        "--igv_export",
        action=argparse.BooleanOptionalAction,
        default=defaults["igv_export"],
        help=(
            "Export per-sample BedGraph tracks (P-site / A-site) under "
            "<output>/igv_tracks/<sample>/, suitable for opening in IGV."
        ),
    )
    # ----- Periodicity QC knobs (spec-defined, see docs/release-notes/v0.6.2.md)
    # The pipeline ALWAYS emits the periodicity bundle; these flags let
    # users tune thresholds / metric toggles without touching code.
    # `periodicity:` YAML section maps 1:1 to these via the standard
    # snake_case-to-hyphen rule (`periodicity.exclude_start_codons:` ->
    # `--periodicity-exclude-start-codons`).
    optional_group.add_argument(
        "--periodicity-enabled",
        "--periodicity_enabled",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Skip the periodicity QC step entirely when set to false.",
    )
    optional_group.add_argument(
        "--periodicity-exclude-start-codons",
        "--periodicity_exclude_start_codons",
        type=int,
        default=None,
        help=(
            "Number of codons after CDS start to mask out of frame "
            "statistics (initiation pause). Default: 6 (per spec)."
        ),
    )
    optional_group.add_argument(
        "--periodicity-exclude-stop-codons",
        "--periodicity_exclude_stop_codons",
        type=int,
        default=None,
        help=(
            "Number of codons before CDS stop to mask (termination "
            "pause). Default: 3 (per spec)."
        ),
    )
    optional_group.add_argument(
        "--periodicity-phase-score",
        "--periodicity_phase_score",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "Compute a ribotricer-style gene-level phase_score column "
            "in gene_periodicity.tsv. Default: enabled."
        ),
    )
    optional_group.add_argument(
        "--periodicity-fourier-spectrum",
        "--periodicity_fourier_spectrum",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "Compute the Wakigawa-style amplitude spectrum per "
            "(sample, read_length, gene, region). Writes "
            "fourier_spectrum.tsv, fourier_period3_score.tsv, "
            "fourier_period3_summary.tsv plus per-(sample, length) "
            "two-panel overlay plots under fourier_spectrum/. "
            "Default: enabled."
        ),
    )
    optional_group.add_argument(
        "--periodicity-fourier-window-nt",
        "--periodicity_fourier_window_nt",
        type=int,
        default=None,
        help=(
            "Window (nt) on each side of the canonical stop codon for "
            "the Fourier spectrum. Default: 100 (Wakigawa published "
            "value)."
        ),
    )
    optional_group.add_argument(
        "--periodicity-metagene-nt",
        "--periodicity_metagene_nt",
        type=int,
        default=None,
        help=(
            "Window (nt) up/downstream of start/stop codons for the "
            "metagene plots. Default: 300 (legacy MitoRiboPy window). "
            "Spec recommends 90 for tighter publication-ready plots."
        ),
    )
    optional_group.add_argument(
        "--periodicity-min-reads-per-length",
        "--periodicity_min_reads_per_length",
        type=int,
        default=None,
        help=(
            "Minimum CDS sites required to score a read length. "
            "Default: 200 (pipeline). Bump to 1000 to match the spec."
        ),
    )
    return parser


def _warn_deprecated(message: str) -> None:
    """Emit a single-line deprecation notice on stderr."""
    print(f"[mitoribopy] DEPRECATED: {message}", file=sys.stderr)


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

    # Canonicalise the deprecated --offset-pick-reference value at the
    # CLI boundary so every downstream consumer sees the new name and
    # the user gets one (and only one) deprecation line.
    if getattr(args, "offset_pick_reference", None) == "selected_site":
        _warn_deprecated(
            "--offset-pick-reference selected_site -> reported_site "
            "(same behaviour, clearer name)."
        )
        args.offset_pick_reference = "reported_site"

    if args.min_5_offset is None:
        args.min_5_offset = args.min_offset
    if args.max_5_offset is None:
        args.max_5_offset = args.max_offset
    if args.min_3_offset is None:
        args.min_3_offset = args.min_offset
    if args.max_3_offset is None:
        args.max_3_offset = args.max_offset

    if args.min_offset > args.max_offset:
        parser.error("--min-offset cannot be greater than --max-offset")
    if args.min_5_offset > args.max_5_offset:
        parser.error("--min-5-offset cannot be greater than --max-5-offset")
    if args.min_3_offset > args.max_3_offset:
        parser.error("--min-3-offset cannot be greater than --max-3-offset")
    if args.offset_mask_nt < 0:
        parser.error("--offset-mask-nt must be zero or a positive integer")
    if args.range <= 0:
        parser.error("--range must be a positive integer")
    if args.cap_percentile <= 0 or args.cap_percentile > 1:
        parser.error("--cap-percentile must be in (0, 1]")
    if len(args.unfiltered_read_length_range) != 2:
        parser.error("--unfiltered-read-length-range requires MIN_LEN MAX_LEN")
    if int(args.unfiltered_read_length_range[0]) > int(args.unfiltered_read_length_range[1]):
        parser.error(
            "--unfiltered-read-length-range MIN_LEN cannot be greater than MAX_LEN"
        )
    if args.structure_density_norm_perc <= 0 or args.structure_density_norm_perc > 1:
        parser.error("--structure-density-norm-perc must be in (0, 1]")
    if args.cor_mask_method == "percentile" and (
        args.cor_mask_percentile <= 0 or args.cor_mask_percentile >= 1
    ):
        parser.error("--cor-mask-percentile must be in (0, 1)")
    # Canonicalise the deprecated short strain aliases (h, y) to their
    # full species names so the rest of the pipeline only sees the
    # canonical form. Emit one DEPRECATED line so the user knows to
    # update their config.
    from ..data.reference_data import (
        BUILTIN_STRAIN_PRESETS,
        STRAIN_ALIASES,
        canonical_strain,
    )

    if args.strain in STRAIN_ALIASES:
        canonical = STRAIN_ALIASES[args.strain]
        _warn_deprecated(
            f"--strain {args.strain} -> {canonical} (use the full species name)."
        )
        args.strain = canonical

    # Strain-preset requirements. Built-in strains (h.sapiens,
    # s.cerevisiae) ship a complete annotation + codon table and need
    # no extra files. 'custom' MUST supply an annotation, an explicit
    # -rpf range, and a codon-table source.
    if args.strain == "custom":
        if not args.annotation_file:
            parser.error("--strain custom requires --annotation-file")
        if args.rpf is None:
            parser.error(
                "--strain custom requires an explicit -rpf MIN_LEN MAX_LEN range"
            )
        if not (args.codon_tables_file or args.codon_table_name):
            parser.error(
                "--strain custom requires --codon-table-name (pick from "
                "the built-in NCBI Genetic Codes list, see "
                "--codon-table-name help) or --codon-tables-file (your "
                "own codon-table JSON)"
            )

    # --footprint-class=custom requires an explicit -rpf even for the
    # built-in strains, because the whole point of 'custom' is
    # "I know my footprint class, don't pick one for me".
    if (
        args.footprint_class == "custom"
        and args.rpf is None
        and args.strain not in BUILTIN_STRAIN_PRESETS
    ):
        parser.error(
            "--footprint-class custom requires an explicit -rpf MIN_LEN MAX_LEN range"
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

    # Per-step timing: each step is wrapped in a Stopwatch, the duration
    # is emitted immediately after the step's own OK banner, and a final
    # summary table is rendered at the end. Mirrors the per-sample timing
    # already in `mitoribopy align`.
    step_timings: list[tuple[str, float]] = []
    wall_sw = Stopwatch()
    wall_sw.__enter__()

    def _run_step(label: str, fn, *fargs, **fkwargs):
        sw = Stopwatch()
        sw.__enter__()
        try:
            result = fn(*fargs, **fkwargs)
        finally:
            sw.__exit__(None, None, None)
            step_timings.append((label, sw.seconds))
            status_writer(
                f"[PIPELINE] {label}: {format_duration(sw.seconds)}"
            )
        return result

    context = _run_step("step 1 initialize", build_pipeline_context, args)
    status_writer(
        "[PIPELINE] Step 1/"
        f"{TOTAL_PIPELINE_STEPS} OK: initialized outputs, annotations, and RPF range."
    )

    _run_step(
        "step 2 read-counts", load_total_read_counts, context, status_writer
    )
    _run_step(
        "step 3 read-length QC",
        run_unfiltered_read_length_qc,
        context,
        status_writer,
    )

    # Steps 4 + 5 each may bail out (no rows pass the BED filter, or no
    # offsets could be enriched). The summary still renders for the
    # steps that did run, which is exactly what the user wants when
    # diagnosing why a run stopped early.
    if _run_step(
        "step 4 filter BED", filter_bed_inputs, context, status_writer
    ) and _run_step(
        "step 5 offset enrichment",
        compute_offset_enrichment_step,
        context,
        status_writer,
    ):
        _run_step(
            "step 6 offset selection",
            select_offsets_and_plot,
            context,
            status_writer,
        )
        _run_step(
            "step 7 downstream",
            run_downstream_modules,
            context,
            status_writer,
        )

    wall_sw.__exit__(None, None, None)
    for line in render_step_timeline(
        step_timings, wall_seconds=wall_sw.seconds
    ):
        status_writer(f"[PIPELINE] {line}")
    return 0


def run_pipeline_cli(argv: Iterable[str] | None = None) -> int:
    """Parse arguments and run the standalone package-native pipeline."""
    args = parse_pipeline_args(argv)
    log_path = configure_file_logging(Path(args.output) / "mitoribopy.log")
    log_message(f"[PIPELINE] Log file: {log_path}")
    configure_plot_defaults()
    return run_pipeline(args)


def main(argv: Iterable[str] | None = None) -> int:
    """Compatibility entrypoint used by both package and wrapper scripts."""
    return run_pipeline_cli(argv)
