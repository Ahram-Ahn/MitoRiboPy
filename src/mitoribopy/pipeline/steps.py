"""Standalone pipeline step helpers built on package-native modules."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Callable

import pandas as pd

from ..analysis import (
    build_per_sample_summaries,
    create_csv_for_offset_enrichment,
    determine_p_site_offsets,
    run_codon_correlation,
    run_rna_seq_analysis,
    run_translation_profile_analysis,
)
from ..config import resolve_rpf_range
from ..data import load_annotation_table, load_codon_table, resolve_start_codons
from ..data.reference_data import BUILTIN_ANNOTATION_PRESETS
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


def build_pipeline_context(args: argparse.Namespace) -> PipelineContext:
    """Create shared runtime state for a package-native pipeline run."""
    base_output_dir = Path(args.output).resolve()
    base_output_dir.mkdir(parents=True, exist_ok=True)

    plot_output_dir = base_output_dir / args.plot_dir
    plot_output_dir.mkdir(parents=True, exist_ok=True)

    # Only h and y ship a built-in annotation. vm / ym / custom require
    # the user to pass --annotation_file (the parser enforces this).
    annotation_df = load_annotation_table(
        preset=(args.strain if args.strain in BUILTIN_ANNOTATION_PRESETS else None),
        annotation_file=args.annotation_file,
        atp8_atp6_baseline=args.atp8_atp6_baseline,
        nd4l_nd4_baseline=args.nd4l_nd4_baseline,
    )
    resolved_codon_table = load_codon_table(
        preset=args.strain,
        table_name=args.codon_table_name,
        table_file=args.codon_tables_file,
    )
    resolved_start_codons = resolve_start_codons(args.strain, args.start_codons)
    unfiltered_read_length_range = tuple(
        int(value) for value in args.unfiltered_read_length_range
    )
    footprint_class = getattr(args, "footprint_class", "monosome")

    return PipelineContext(
        args=args,
        base_output_dir=base_output_dir,
        plot_output_dir=plot_output_dir,
        annotation_df=annotation_df,
        resolved_codon_table=resolved_codon_table,
        resolved_start_codons=resolved_start_codons,
        rpf_range=resolve_rpf_range(
            args.strain, args.rpf, footprint_class=footprint_class
        ),
        unfiltered_read_length_range=unfiltered_read_length_range,
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
                mt_mrna_substring_patterns=context.args.mt_mrna_substring_patterns,
            )
        except Exception as exc:
            emit_status(
                f"[PIPELINE] Failed to parse read-count file '{read_counts_path}': {exc}"
            )
            context.total_counts_map = {}

    # Task 5c: derived state lives on the context. The
    # context.args.total_mrna_map assignment below is a compatibility
    # shim for external code that still reads from the namespace; it
    # emits a deprecation warning on access and will be removed in
    # v0.4.0. Internal consumers (run_coverage_profile_plots) now take
    # total_mrna_map as an explicit kwarg.
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
    min_len, max_len = context.unfiltered_read_length_range
    unfiltered_summary_csv = (
        context.plot_output_dir
        / f"unfiltered_read_length_summary_{min_len}_{max_len}.csv"
    )
    compute_unfiltered_read_length_summary(
        input_dir=context.args.directory,
        output_csv=str(unfiltered_summary_csv),
        total_counts_map=context.total_counts_map,
        read_length_range=context.unfiltered_read_length_range,
    )

    heatmap_base = context.plot_output_dir / f"unfiltered_heatmap_{min_len}_{max_len}"
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
    """Filter BED/BAM inputs and resolve the sample order for downstream steps."""
    filtered_bed_df, sample_dirs = process_bed_files(
        input_dir=context.args.directory,
        output_dir=str(context.plot_output_dir),
        organism=context.args.strain,
        annotation_df=context.annotation_df,
        rpf_range=context.rpf_range,
        bam_mapq=int(getattr(context.args, "bam_mapq", 0) or 0),
    )
    context.filtered_bed_df = filtered_bed_df

    if filtered_bed_df.empty:
        context.sample_dirs = []
        _emit_step_stop(4, "no data remained after BED filtering; pipeline finished early.", emit_status)
        return False

    context.sample_dirs = _order_sample_dirs(sample_dirs, context.args.order_samples, emit_status)

    # Emit the per-sample per-gene RPF count table consumed by
    # 'mitoribopy rnaseq'. This is the provenance spine that ties rpf
    # outputs to the downstream TE / delta-TE step.
    write_rpf_counts_table(
        filtered_bed_df=filtered_bed_df,
        output_path=context.base_output_dir / "rpf_counts.tsv",
    )

    # Emit run_settings.json with reference_checksum so the Phase 5
    # rnaseq subcommand's reference-consistency gate can verify that
    # the Ribo-seq and RNA-seq sides used the same transcript set.
    write_rpf_run_settings(context)

    _emit_step_ok(
        4,
        f"retained {len(filtered_bed_df)} filtered reads across {len(context.sample_dirs)} sample(s).",
        emit_status,
    )
    return True


def write_rpf_run_settings(context: PipelineContext) -> Path | None:
    """Write ``<output>/run_settings.json`` with the rpf-side provenance.

    Key field for downstream rnaseq: ``reference_checksum`` is the
    SHA-256 of the FASTA passed via ``-f / --fasta``. The rnaseq
    subcommand's reference-consistency gate reads this and refuses to
    proceed unless the rnaseq side matches.
    """
    import hashlib
    import json

    from .. import __version__

    fasta_path = Path(context.args.fasta)
    if not fasta_path.is_file():
        # Do not emit a manifest with a bogus hash; rnaseq will fail
        # loudly on missing checksum, which is the correct behavior.
        return None

    digest = hashlib.sha256()
    with fasta_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)

    settings = {
        "subcommand": "rpf",
        "mitoribopy_version": __version__,
        "reference_fasta": str(fasta_path),
        "reference_checksum": digest.hexdigest(),
        "strain": context.args.strain,
        "rpf_range": list(context.rpf_range),
        "align": context.args.align,
        "offset_type": context.args.offset_type,
        "offset_site": context.args.offset_site,
        "codon_overlap_mode": context.args.codon_overlap_mode,
        "bam_mapq": int(getattr(context.args, "bam_mapq", 0) or 0),
    }

    out_path = context.base_output_dir / "run_settings.json"
    out_path.write_text(
        json.dumps(settings, indent=2, sort_keys=True), encoding="utf-8"
    )
    return out_path


def write_rpf_counts_table(
    *,
    filtered_bed_df,
    output_path,
) -> None:
    """Write the per-sample per-gene RPF read-count TSV.

    Columns: ``sample``, ``gene``, ``count``. One row per unique
    (sample_name, chrom) pair. Used by ``mitoribopy rnaseq`` to load
    Ribo-seq counts for TE / delta-TE computation.
    """
    from pathlib import Path

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if filtered_bed_df is None or filtered_bed_df.empty:
        output_path.write_text("sample\tgene\tcount\n", encoding="utf-8")
        return

    grouped = (
        filtered_bed_df.groupby(["sample_name", "chrom"]).size().reset_index(name="count")
    )
    grouped = grouped.rename(columns={"sample_name": "sample", "chrom": "gene"})
    grouped = grouped.sort_values(["sample", "gene"]).reset_index(drop=True)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("sample\tgene\tcount\n")
        for _, row in grouped.iterrows():
            handle.write(f"{row['sample']}\t{row['gene']}\t{int(row['count'])}\n")


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

    # Build per-sample summaries (used when --offset_mode=per_sample).
    context.per_sample_offset_summaries = build_per_sample_summaries(
        offsets_df=offsets_df,
        rpf_range=context.rpf_range,
        offset_limit=effective_offset_limit,
        offset_mask_nt=context.args.offset_mask_nt,
    )

    # Persist per-sample enrichment tables so the user can compare them
    # against the combined diagnostic.
    if context.per_sample_offset_summaries:
        per_sample_dir = context.plot_output_dir / "per_sample"
        per_sample_dir.mkdir(parents=True, exist_ok=True)
        for sample_name, per_summary in context.per_sample_offset_summaries.items():
            sample_subdir = per_sample_dir / sample_name
            sample_subdir.mkdir(parents=True, exist_ok=True)
            per_summary.to_csv(
                sample_subdir / f"offset_{context.args.align}.csv", index=False
            )

    _emit_step_ok(
        5,
        f"computed offset enrichment tables for {len(context.rpf_range)} read length(s) "
        f"(combined + {len(context.per_sample_offset_summaries)} per-sample).",
        emit_status,
    )
    return True


def _plot_offset_drift(
    *,
    selected_offsets_by_sample: dict[str, pd.DataFrame],
    combined: pd.DataFrame | None,
    align_to: str,
    output_path: Path,
    plot_format: str,
) -> Path | None:
    """Render a per-sample 5' / 3' offset comparison by read length.

    Returns ``None`` when there is nothing to plot (fewer than two
    samples and no combined reference).
    """
    if not selected_offsets_by_sample:
        return None

    rows: list[dict] = []
    for sample, df in selected_offsets_by_sample.items():
        if df is None:
            continue
        for _, row in df.iterrows():
            rows.append(
                {
                    "sample": sample,
                    "Read Length": int(row["Read Length"]),
                    "5' Offset": row.get("Most Enriched 5' Offset"),
                    "3' Offset": row.get("Most Enriched 3' Offset"),
                }
            )
    if combined is not None:
        for _, row in combined.iterrows():
            rows.append(
                {
                    "sample": "[combined]",
                    "Read Length": int(row["Read Length"]),
                    "5' Offset": row.get("Most Enriched 5' Offset"),
                    "3' Offset": row.get("Most Enriched 3' Offset"),
                }
            )
    if not rows:
        return None

    drift_df = pd.DataFrame(rows)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path = output_path.with_suffix(".csv")
    drift_df.to_csv(csv_path, index=False)

    import matplotlib.pyplot as plt

    samples = [s for s in drift_df["sample"].unique() if s != "[combined]"]
    if "[combined]" in drift_df["sample"].unique():
        samples.append("[combined]")
    if not samples:
        return csv_path

    fig, (ax5, ax3) = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for sample in samples:
        sub = drift_df[drift_df["sample"] == sample].sort_values("Read Length")
        linestyle = "--" if sample == "[combined]" else "-"
        linewidth = 2.2 if sample == "[combined]" else 1.4
        ax5.plot(
            sub["Read Length"],
            sub["5' Offset"],
            marker="o",
            linestyle=linestyle,
            linewidth=linewidth,
            label=sample,
        )
        ax3.plot(
            sub["Read Length"],
            sub["3' Offset"],
            marker="o",
            linestyle=linestyle,
            linewidth=linewidth,
            label=sample,
        )
    ax5.set_title(f"5' offset by sample (align to {align_to})")
    ax3.set_title(f"3' offset by sample (align to {align_to})")
    for ax in (ax5, ax3):
        ax.set_xlabel("Read length (nt)")
        ax.set_ylabel("Offset (nt)")
        ax.grid(True, alpha=0.2)
    ax3.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(str(output_path), format=plot_format)
    plt.close(fig)
    return output_path


def select_offsets_and_plot(context: PipelineContext, emit_status: StatusWriter) -> None:
    """Select offsets and render the diagnostic offset plots.

    Always produces a combined-across-samples selection as the
    diagnostic spine (written to ``p_site_offsets_<align>.csv`` and
    recorded on ``context.selected_offsets_df``). When
    ``--offset_mode=per_sample`` is active (the default) it additionally
    runs per-sample selection and records the per-sample picks on
    ``context.selected_offsets_by_sample``. Downstream modules honor
    ``--offset_mode`` when choosing which set of offsets to apply.
    """
    combined_offset_file = context.plot_output_dir / f"p_site_offsets_{context.args.align}.csv"
    context.selected_offsets_df = determine_p_site_offsets(
        offsets_df=context.offset_details_df,
        align_to=context.args.align,
        out_file=str(combined_offset_file),
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

    # Per-sample selection. Skipped when the user pinned a fixed
    # manual offset (--psite_offset) since per-sample selection would
    # just echo the global choice.
    offset_mode = getattr(context.args, "offset_mode", "per_sample")
    manual_offset = getattr(context.args, "psite_offset", None)
    if offset_mode == "per_sample" and manual_offset is None and context.offset_details_df is not None:
        per_sample_selections: dict[str, pd.DataFrame] = {}
        sample_dir = context.plot_output_dir / "per_sample"
        sample_dir.mkdir(parents=True, exist_ok=True)
        if "sample_name" in context.offset_details_df.columns:
            for sample_name, sub in context.offset_details_df.groupby(
                "sample_name", sort=True
            ):
                out_csv = sample_dir / sample_name / f"p_site_offsets_{context.args.align}.csv"
                out_csv.parent.mkdir(parents=True, exist_ok=True)
                selection = determine_p_site_offsets(
                    offsets_df=sub,
                    align_to=context.args.align,
                    out_file=str(out_csv),
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
                if selection is not None:
                    per_sample_selections[str(sample_name)] = selection
                else:
                    emit_status(
                        "[PIPELINE] No per-sample offsets selected for "
                        f"{sample_name}; will fall back to the combined "
                        "offsets table for that sample."
                    )
        context.selected_offsets_by_sample = per_sample_selections
    else:
        context.selected_offsets_by_sample = {}

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

    # Per-sample offset drift plot — makes cross-sample offset variance
    # visible before the user reads any downstream output.
    if context.selected_offsets_by_sample:
        drift_plot = context.plot_output_dir / f"offset_drift_{context.args.align}.{context.args.plot_format}"
        _plot_offset_drift(
            selected_offsets_by_sample=context.selected_offsets_by_sample,
            combined=context.selected_offsets_df,
            align_to=context.args.align,
            output_path=drift_plot,
            plot_format=context.args.plot_format,
        )

    if context.selected_offsets_df is None and not context.selected_offsets_by_sample:
        _emit_step_ok(
            6,
            "rendered offset diagnostics but no P-site offsets were selected.",
            emit_status,
        )
        return

    combined_n = (
        len(context.selected_offsets_df) if context.selected_offsets_df is not None else 0
    )
    _emit_step_ok(
        6,
        f"selected offsets for {combined_n} read length(s) (combined) "
        f"+ {len(context.selected_offsets_by_sample)} per-sample; "
        "rendered offset diagnostics.",
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
        analysis_sites = getattr(context.args, "analysis_sites", "both")
        if analysis_sites == "both":
            requested_sites = ["p", "a"]
        else:
            requested_sites = [analysis_sites]

        for site in requested_sites:
            site_label = "P-site" if site == "p" else "A-site"
            if len(requested_sites) == 1:
                # Single-site run: write to the legacy top-level layout
                # so v0.3 consumers see no surprise rename.
                profile_out = context.base_output_dir
                coverage_out = context.base_output_dir / "coverage_profile_plots"
            else:
                # Multi-site run: per-site subdirectory keeps the two
                # outputs side by side without overwriting each other.
                profile_out = context.base_output_dir / f"translation_profile_{site}"
                coverage_out = context.base_output_dir / f"coverage_profile_plots_{site}"
            profile_out.mkdir(parents=True, exist_ok=True)
            coverage_out.mkdir(parents=True, exist_ok=True)

            run_translation_profile_analysis(
                sample_dirs=context.sample_dirs,
                selected_offsets_df=context.selected_offsets_df,
                offset_type=context.args.offset_type,
                fasta_file=context.args.fasta,
                output_dir=str(profile_out),
                args=context.args,
                annotation_df=context.annotation_df,
                filtered_bed_df=context.filtered_bed_df,
                resolved_codon_table=context.resolved_codon_table,
                resolved_start_codons=context.resolved_start_codons,
                selected_offsets_by_sample=context.selected_offsets_by_sample or None,
                site_override=site,
            )
            run_coverage_profile_plots(
                sample_dirs=context.sample_dirs,
                selected_offsets_df=context.selected_offsets_df,
                offset_type=context.args.offset_type,
                fasta_file=context.args.fasta,
                output_dir=str(coverage_out),
                args=context.args,
                annotation_df=context.annotation_df,
                filtered_bed_df=context.filtered_bed_df,
                total_mrna_map=context.total_counts_map,
                selected_offsets_by_sample=context.selected_offsets_by_sample or None,
                site_override=site,
            )
            ran_modules.append(f"translation-profile analysis ({site_label})")
            ran_modules.append(f"coverage-profile plots ({site_label})")

        if context.args.structure_density:
            structure_density_dir = context.base_output_dir / "structure_density"
            # Pick the structure-density site to match the user's
            # analysis_sites preference. When both sites were generated,
            # default to P-site (matches v0.3.x behaviour).
            structure_site = "p" if "p" in requested_sites else "a"
            site_column = "P_site" if structure_site == "p" else "A_site"
            if len(requested_sites) > 1:
                tp_source = context.base_output_dir / f"translation_profile_{structure_site}"
            else:
                tp_source = context.base_output_dir
            run_structure_density_export(
                translation_profile_dir=str(tp_source),
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
                # Codon correlation reads from a translation_profile dir.
                # When both sites were generated, default to the P-site
                # subdir (matches v0.3.x).
                cor_site = "p" if "p" in requested_sites else "a"
                if len(requested_sites) > 1:
                    cor_source = context.base_output_dir / f"translation_profile_{cor_site}"
                else:
                    cor_source = context.base_output_dir
                cor_out_dir = context.base_output_dir / "codon_correlation"
                run_codon_correlation(
                    translation_profile_dir=str(cor_source),
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
