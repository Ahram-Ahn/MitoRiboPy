"""Standalone pipeline step helpers built on package-native modules."""

from __future__ import annotations

from pathlib import Path
from typing import Callable

from ..analysis import (
    create_csv_for_offset_enrichment,
    determine_p_site_offsets,
    run_codon_correlation,
    run_rna_seq_analysis,
    run_translation_profile_analysis,
)
from ..config import resolve_rpf_range
from ..data import load_annotation_table, load_codon_table, resolve_start_codons
from ..io import (
    compute_total_counts,
    compute_unfiltered_read_length_summary,
    process_bed_files,
)
from ..plotting import (
    plot_offset_enrichment,
    plot_unfiltered_read_length_heatmap,
    run_coverage_profile_plots,
    run_structure_density_export,
)
from .context import PipelineContext

StatusWriter = Callable[[str], None]
TOTAL_PIPELINE_STEPS = 7


def _emit_step_ok(step_number: int, message: str, emit_status: StatusWriter) -> None:
    emit_status(f"[PIPELINE] Step {step_number}/{TOTAL_PIPELINE_STEPS} OK: {message}")


def _emit_step_stop(step_number: int, message: str, emit_status: StatusWriter) -> None:
    emit_status(f"[PIPELINE] Step {step_number}/{TOTAL_PIPELINE_STEPS} STOP: {message}")


def _resolved_offset_plot_limit(args) -> int:
    """Expand the plotted offset window so it always covers selectable offsets."""
    return max(
        int(args.range),
        int(args.max_offset),
        int(args.max_5_offset),
        int(args.max_3_offset),
    )


def build_pipeline_context(args) -> PipelineContext:
    """Create shared runtime state for a package-native pipeline run."""
    base_output_dir = Path(args.output).resolve()
    base_output_dir.mkdir(parents=True, exist_ok=True)

    plot_output_dir = base_output_dir / args.plot_dir
    plot_output_dir.mkdir(parents=True, exist_ok=True)

    annotation_df = load_annotation_table(
        preset=(args.strain if args.strain in {"y", "h"} else None),
        annotation_file=args.annotation_file,
        atp8_atp6_baseline=args.atp8_atp6_baseline,
        nd4l_nd4_baseline=args.nd4l_nd4_baseline,
    )
    args.resolved_codon_table = load_codon_table(
        preset=args.strain,
        table_name=args.codon_table_name,
        table_file=args.codon_tables_file,
    )
    args.resolved_start_codons = resolve_start_codons(args.strain, args.start_codons)

    return PipelineContext(
        args=args,
        base_output_dir=base_output_dir,
        plot_output_dir=plot_output_dir,
        annotation_df=annotation_df,
        rpf_range=resolve_rpf_range(args.strain, args.rpf),
    )


def load_total_read_counts(context: PipelineContext, emit_status: StatusWriter) -> None:
    """Load optional read-count totals used for RPM normalization."""
    read_counts_path = Path(context.args.read_counts_file)
    if not read_counts_path.exists():
        emit_status(
            f"[PIPELINE] Read-count file not found: {read_counts_path}. Continuing without RPM totals."
        )
        context.total_counts_map = {}
    else:
        try:
            context.total_counts_map, _ = compute_total_counts(
                str(read_counts_path),
                sample_col=context.args.read_counts_sample_col,
                reads_col=context.args.read_counts_reads_col,
                normalization_mode=context.args.rpm_norm_mode,
                reference_col=context.args.read_counts_reference_col,
                mrna_ref_patterns=context.args.mrna_ref_patterns,
            )
        except Exception as exc:
            emit_status(
                f"[PIPELINE] Failed to parse read-count file '{read_counts_path}': {exc}"
            )
            context.total_counts_map = {}

    context.args.total_mrna_map = context.total_counts_map
    _emit_step_ok(
        2,
        f"loaded read-count totals for {len(context.total_counts_map)} sample(s).",
        emit_status,
    )


def run_unfiltered_read_length_qc(
    context: PipelineContext,
    emit_status: StatusWriter,
) -> None:
    """Generate unfiltered read-length summary tables and heatmaps."""
    unfiltered_summary_csv = context.plot_output_dir / "unfiltered_read_length_summary_15_50.csv"
    compute_unfiltered_read_length_summary(
        input_dir=context.args.directory,
        output_csv=str(unfiltered_summary_csv),
        total_counts_map=context.total_counts_map,
        read_length_range=(15, 50),
    )

    heatmap_base = context.plot_output_dir / "unfiltered_heatmap_15_50"
    plot_unfiltered_read_length_heatmap(
        summary_csv_path=str(unfiltered_summary_csv),
        output_png_base=str(heatmap_base),
        value_col="count",
    )
    plot_unfiltered_read_length_heatmap(
        summary_csv_path=str(unfiltered_summary_csv),
        output_png_base=str(heatmap_base),
        value_col="RPM",
    )

    _emit_step_ok(
        3,
        f"wrote unfiltered read-length QC outputs to {context.plot_output_dir}.",
        emit_status,
    )


def _order_sample_dirs(sample_dirs: list[str], requested_order: list[str] | None, emit_status: StatusWriter) -> list[str]:
    if not requested_order:
        return sample_dirs

    name_to_dir = {Path(sample_dir).name: sample_dir for sample_dir in sample_dirs}
    ordered_dirs: list[str] = []
    missing_requested: list[str] = []

    for requested_sample in requested_order:
        if requested_sample in name_to_dir:
            ordered_dirs.append(name_to_dir.pop(requested_sample))
        else:
            missing_requested.append(requested_sample)

    ordered_dirs.extend(sorted(name_to_dir.values(), key=lambda path: Path(path).name))
    if missing_requested:
        emit_status(
            "[PIPELINE] Requested sample(s) not found in filtered BEDs: "
            + ", ".join(missing_requested)
        )
    return ordered_dirs


def filter_bed_inputs(context: PipelineContext, emit_status: StatusWriter) -> bool:
    """Filter BED inputs and resolve the sample order for downstream steps."""
    filtered_bed_df, sample_dirs = process_bed_files(
        input_dir=context.args.directory,
        output_dir=str(context.plot_output_dir),
        organism=context.args.strain,
        annotation_df=context.annotation_df,
        rpf_range=context.rpf_range,
    )
    context.filtered_bed_df = filtered_bed_df

    if filtered_bed_df.empty:
        context.sample_dirs = []
        _emit_step_stop(4, "no data remained after BED filtering; pipeline finished early.", emit_status)
        return False

    context.sample_dirs = _order_sample_dirs(sample_dirs, context.args.order_samples, emit_status)
    _emit_step_ok(
        4,
        f"retained {len(filtered_bed_df)} filtered reads across {len(context.sample_dirs)} sample(s).",
        emit_status,
    )
    return True


def compute_offset_enrichment_step(context: PipelineContext, emit_status: StatusWriter) -> bool:
    """Build the offset enrichment summary and detailed offset tables."""
    effective_offset_limit = _resolved_offset_plot_limit(context.args)
    context.extra["effective_offset_limit"] = effective_offset_limit
    if effective_offset_limit > int(context.args.range):
        emit_status(
            "[PIPELINE] Expanded offset plotting window from "
            f"{context.args.range} to {effective_offset_limit} nt to cover the selectable offsets."
        )

    offset_csv_path = context.plot_output_dir / f"offset_{context.args.align}.csv"
    summary_df, offsets_df = create_csv_for_offset_enrichment(
        bed_df=context.filtered_bed_df,
        annotation_df=context.annotation_df,
        align_to=context.args.align,
        rpf_range=context.rpf_range,
        output_csv=str(offset_csv_path),
        offset_limit=effective_offset_limit,
        manual_offset=context.args.psite_offset,
        offset_mask_nt=context.args.offset_mask_nt,
        offset_site=context.args.offset_site,
        codon_overlap_mode=context.args.codon_overlap_mode,
        strain=context.args.strain,
    )
    context.offset_summary_df = summary_df
    context.offset_details_df = offsets_df

    if summary_df is None or offsets_df is None:
        _emit_step_stop(5, "no offsets were generated; pipeline finished early.", emit_status)
        return False

    _emit_step_ok(
        5,
        f"computed offset enrichment tables for {len(context.rpf_range)} read length(s).",
        emit_status,
    )
    return True


def select_offsets_and_plot(context: PipelineContext, emit_status: StatusWriter) -> None:
    """Select offsets and render the diagnostic offset plots."""
    p_site_offset_file = context.plot_output_dir / f"p_site_offsets_{context.args.align}.csv"
    context.selected_offsets_df = determine_p_site_offsets(
        offsets_df=context.offset_details_df,
        align_to=context.args.align,
        out_file=str(p_site_offset_file),
        offset_min=context.args.min_offset,
        offset_max=context.args.max_offset,
        five_offset_min=context.args.min_5_offset,
        five_offset_max=context.args.max_5_offset,
        three_offset_min=context.args.min_3_offset,
        three_offset_max=context.args.max_3_offset,
        offset_mask_nt=context.args.offset_mask_nt,
        offset_site=context.args.offset_site,
        selection_reference=context.args.offset_pick_reference,
    )

    plot_offset_enrichment(
        summary_df=context.offset_summary_df,
        align_to=context.args.align,
        plot_dir=str(context.plot_output_dir),
        plot_format=context.args.plot_format,
        x_breaks=context.args.x_breaks,
        line_plot_style=context.args.line_plot_style,
        offset_limit=context.extra.get("effective_offset_limit", context.args.range),
        offset_mask_nt=context.args.offset_mask_nt,
        selected_offsets=context.selected_offsets_df,
        offset_min=context.args.min_offset,
        offset_max=context.args.max_offset,
        five_offset_min=context.args.min_5_offset,
        five_offset_max=context.args.max_5_offset,
        three_offset_min=context.args.min_3_offset,
        three_offset_max=context.args.max_3_offset,
    )

    if context.selected_offsets_df is None:
        _emit_step_ok(
            6,
            "rendered offset diagnostics but no P-site offsets were selected.",
            emit_status,
        )
        return

    _emit_step_ok(
        6,
        f"selected offsets for {len(context.selected_offsets_df)} read length(s) and rendered diagnostics.",
        emit_status,
    )


def run_downstream_modules(context: PipelineContext, emit_status: StatusWriter) -> None:
    """Run downstream analyses that consume the selected offsets."""
    ran_modules: list[str] = []
    skipped_modules: list[str] = []

    if context.selected_offsets_df is None:
        skipped_modules.extend(
            ["translation-profile analysis", "coverage-profile plots", "structure-density export"]
        )
    else:
        run_translation_profile_analysis(
            sample_dirs=context.sample_dirs,
            selected_offsets_df=context.selected_offsets_df,
            offset_type=context.args.offset_type,
            fasta_file=context.args.fasta,
            output_dir=str(context.base_output_dir),
            args=context.args,
            annotation_df=context.annotation_df,
            filtered_bed_df=context.filtered_bed_df,
        )
        ran_modules.append("translation-profile analysis")

        coverage_plot_dir = context.base_output_dir / "coverage_profile_plots"
        coverage_plot_dir.mkdir(parents=True, exist_ok=True)
        run_coverage_profile_plots(
            sample_dirs=context.sample_dirs,
            selected_offsets_df=context.selected_offsets_df,
            offset_type=context.args.offset_type,
            fasta_file=context.args.fasta,
            output_dir=str(coverage_plot_dir),
            args=context.args,
            annotation_df=context.annotation_df,
            filtered_bed_df=context.filtered_bed_df,
        )
        ran_modules.append("coverage-profile plots")

        if context.args.structure_density:
            structure_density_dir = context.base_output_dir / "structure_density"
            site_column = "P_site" if context.args.offset_site == "p" else "A_site"
            run_structure_density_export(
                translation_profile_dir=str(context.base_output_dir),
                structure_density_norm_perc=context.args.structure_density_norm_perc,
                output_dir=str(structure_density_dir),
                site_column=site_column,
            )
            ran_modules.append("structure-density export")
        else:
            skipped_modules.append("structure-density export")

    if context.args.cor_plot and context.args.base_sample:
        if context.selected_offsets_df is None:
            skipped_modules.append("codon correlation")
        else:
            sample_names = [Path(sample_dir).name for sample_dir in context.sample_dirs]
            if context.args.base_sample not in sample_names:
                emit_status(
                    f"[PIPELINE] base_sample '{context.args.base_sample}' not present. Skipping codon correlation."
                )
                skipped_modules.append("codon correlation")
            else:
                cor_out_dir = context.base_output_dir / "codon_correlation"
                run_codon_correlation(
                    translation_profile_dir=str(context.base_output_dir),
                    samples=sample_names,
                    base_sample=context.args.base_sample,
                    column="CoverageDivFreq",
                    output_dir=str(cor_out_dir),
                    mask_method=context.args.cor_mask_method,
                    mask_percentile=context.args.cor_mask_percentile,
                    mask_threshold=(
                        context.args.cor_mask_threshold
                        if context.args.cor_mask_method == "fixed"
                        else None
                    ),
                )
                ran_modules.append("codon correlation")
    else:
        skipped_modules.append("codon correlation")

    if context.args.use_rna_seq:
        rna_out_dir = context.base_output_dir / context.args.rna_out_dir
        rna_out_dir.mkdir(parents=True, exist_ok=True)
        run_rna_seq_analysis(
            rna_seq_dir=context.args.rna_seq_dir,
            rna_order=context.args.rna_order,
            annotation_df=context.annotation_df,
            fasta_file=context.args.fasta,
            output_dir=str(rna_out_dir),
            translation_profile_dir=str(context.base_output_dir),
            do_merge_with_ribo=context.args.do_rna_ribo_ratio,
            ribo_order=context.args.order_samples,
        )
        ran_modules.append("RNA-seq integration")
    else:
        skipped_modules.append("RNA-seq integration")

    summary_parts: list[str] = []
    if ran_modules:
        summary_parts.append("ran " + ", ".join(ran_modules))
    if skipped_modules:
        summary_parts.append("skipped " + ", ".join(skipped_modules))

    _emit_step_ok(7, "; ".join(summary_parts) + ".", emit_status)
