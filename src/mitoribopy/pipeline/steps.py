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
    run_translation_profile_analysis,
)
from ..config import resolve_rpf_range
from ..data import load_annotation_table, load_codon_table, resolve_start_codons
from ..data.reference_data import BUILTIN_ANNOTATION_PRESETS
from ..io import (
    compute_total_counts,
    compute_unfiltered_read_length_summary,
    prepare_bam_inputs,
    process_bed_files,
)
from ..plotting import (
    plot_offset_enrichment,
    plot_unfiltered_read_length_heatmap,
    run_coverage_profile_plots,
    run_structure_density_export,
)
from ..plotting.igv_export import run_igv_export
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
    # Separate CSV and plot outputs under <plot_output_dir>/{csv,plots}/
    # so the user does not have to wade through interleaved file types.
    csv_dir = plot_output_dir / "csv"
    plot_subdir = plot_output_dir / "plots"
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_subdir.mkdir(parents=True, exist_ok=True)

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
        csv_dir=csv_dir,
        plot_subdir=plot_subdir,
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

    # Derived state lives on the context. The context.args.total_mrna_map
    # assignment below is a compatibility shim for external code that
    # still reads from the namespace; internal consumers
    # (run_coverage_profile_plots) take total_mrna_map as an explicit kwarg.
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
        context.csv_dir
        / f"unfiltered_read_length_summary_{min_len}_{max_len}.csv"
    )
    converted_bed_paths = prepare_bam_inputs(
        input_dir=Path(context.args.directory),
        converted_dir=context.plot_output_dir / "bam_converted",
        mapq_threshold=int(getattr(context.args, "bam_mapq", 0) or 0),
    )
    context.extra["converted_bed_paths"] = [str(p) for p in converted_bed_paths]
    compute_unfiltered_read_length_summary(
        input_dir=context.args.directory,
        output_csv=str(unfiltered_summary_csv),
        total_counts_map=context.total_counts_map,
        read_length_range=context.unfiltered_read_length_range,
        converted_bed_paths=converted_bed_paths,
    )

    heatmap_base = context.plot_subdir / f"unfiltered_heatmap_{min_len}_{max_len}"
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


def _apply_rpf_count_filter(
    context: PipelineContext, emit_status: StatusWriter
) -> None:
    """Drop read-length bins whose total count is below FRAC x the dominant length.

    Pure data-driven prune: keep the read lengths that carry the bulk of
    the signal, drop the noisy long tail. The user controls FRAC via
    ``--rpf_min_count_frac`` (default 0.20). Setting FRAC <= 0 disables.
    Filters on the COMBINED-across-samples count so the surviving
    read-length set is consistent between samples (per-sample offset
    selection still varies within that shared set).
    """
    # Snapshot the user-requested range up-front so the metadata
    # sidecar can record both "what they asked for" and "what survived".
    context.requested_rpf_range = sorted(int(x) for x in context.rpf_range)

    frac = float(getattr(context.args, "rpf_min_count_frac", 0.0) or 0.0)
    df = context.filtered_bed_df
    if df is not None and not df.empty and "read_length" in df.columns:
        counts = df["read_length"].value_counts()
        # Always record observed counts so the metadata sidecar carries
        # them whether or not the auto-filter kicks in.
        context.observed_counts_by_length = {
            int(k): int(v) for k, v in counts.items()
        }

    if frac <= 0:
        context.read_length_filter_threshold_rule = "disabled"
        context.read_length_filter_threshold_value = 0.0
        return  # disabled
    if df is None or df.empty or "read_length" not in df.columns:
        return
    counts = df["read_length"].value_counts()
    if counts.empty:
        return
    max_count = int(counts.max())
    threshold = max_count * frac
    context.read_length_filter_threshold_rule = "dominant_fraction"
    context.read_length_filter_threshold_value = float(frac)
    kept = sorted(int(rl) for rl, c in counts.items() if c >= threshold)
    if not kept:
        return  # belt-and-braces; should never trigger because the max bin clears its own threshold
    dropped = sorted(set(int(x) for x in context.rpf_range) - set(kept))
    context.dropped_lengths = dropped
    if not dropped:
        return  # nothing pruned; current range is already tight
    emit_status(
        f"[BED] read-length auto-filter: kept {kept} (>={frac:.0%} of "
        f"max={max_count:,} reads); dropped {dropped}."
    )
    context.rpf_range = kept
    context.filtered_bed_df = df[df["read_length"].isin(kept)].copy()


def filter_bed_inputs(context: PipelineContext, emit_status: StatusWriter) -> bool:
    """Filter BED/BAM inputs and resolve the sample order for downstream steps."""
    filtered_bed_df, sample_dirs = process_bed_files(
        input_dir=context.args.directory,
        output_dir=str(context.plot_output_dir),
        organism=context.args.strain,
        annotation_df=context.annotation_df,
        rpf_range=context.rpf_range,
        bam_mapq=int(getattr(context.args, "bam_mapq", 0) or 0),
        plot_dir=str(context.plot_subdir),
        csv_dir=str(context.csv_dir),
        converted_bed_paths=[
            Path(p) for p in context.extra.get("converted_bed_paths", [])
        ] or None,
    )
    context.filtered_bed_df = filtered_bed_df

    if filtered_bed_df.empty:
        context.sample_dirs = []
        _emit_step_stop(4, "no data remained after BED filtering; pipeline finished early.", emit_status)
        return False

    context.sample_dirs = _order_sample_dirs(sample_dirs, context.args.order_samples, emit_status)

    # Auto-prune read-length bins whose total count is far below the
    # dominant length. This lets users pass a wide --rpf range (e.g.
    # 27 36 for human) and rely on the data to decide which lengths
    # carry real signal. Quiet when nothing is pruned.
    _apply_rpf_count_filter(context, emit_status)

    # Emit the per-sample per-gene RPF count table consumed by
    # 'mitoribopy rnaseq'. This is the provenance spine that ties rpf
    # outputs to the downstream TE / delta-TE step.
    write_rpf_counts_table(
        filtered_bed_df=context.filtered_bed_df,
        output_path=context.base_output_dir / "rpf_counts.tsv",
    )
    write_rpf_counts_metadata(
        context=context,
        output_path=context.base_output_dir / "rpf_counts.metadata.json",
    )

    # Emit run_settings.json with reference_checksum so the rnaseq
    # subcommand's reference-consistency gate can verify that the
    # Ribo-seq and RNA-seq sides used the same transcript set.
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


def write_rpf_counts_metadata(
    *,
    context: PipelineContext,
    output_path,
) -> Path:
    """Write the rpf_counts.tsv provenance sidecar.

    Documents the read-length filtering decision, the offset mode the
    counts were produced under, and the reference state — every field a
    reviewer of a TE / ΔTE table needs to defend the upstream counts.

    Schema (v1.0.0)::

        {
          "schema_version": "1.0.0",
          "mitoribopy_version": "...",
          "subcommand": "rpf",
          "counts_path": "rpf_counts.tsv",
          "n_rows": int,            # rows in the counts table
          "n_samples": int,
          "n_genes": int,
          "total_reads": int,       # sum of count column
          "read_length_filter": {
            "requested_rpf_range": [29, 30, 31, 32, 33, 34],
            "observed_counts_by_length": {"30": 12345, ...},
            "retained_lengths":  [29, 30, 31, 32],
            "dropped_lengths":   [33, 34],
            "threshold_rule":    "dominant_fraction" | "disabled",
            "threshold_value":   0.20
          },
          "offset_mode": "per_sample" | "combined",
          "offset_type": "5",
          "offset_site": "p",
          "footprint_class": "monosome",
          "strain": "h.sapiens",
          "reference_fasta": "...",
          "reference_checksum": "<sha256>" | None
        }

    Always written regardless of whether the auto-filter pruned
    anything, so a downstream tool can read the sidecar and trust that
    every field is present.
    """
    import hashlib
    import json

    from .. import __version__

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = context.filtered_bed_df
    if df is None or df.empty:
        n_rows = 0
        n_samples = 0
        n_genes = 0
        total_reads = 0
    else:
        grouped = df.groupby(["sample_name", "chrom"]).size()
        n_rows = int(len(grouped))
        n_samples = int(df["sample_name"].nunique())
        n_genes = int(df["chrom"].nunique())
        total_reads = int(len(df))

    fasta_path = Path(context.args.fasta) if getattr(context.args, "fasta", None) else None
    reference_checksum: str | None = None
    if fasta_path and fasta_path.is_file():
        digest = hashlib.sha256()
        with fasta_path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(65536), b""):
                digest.update(chunk)
        reference_checksum = digest.hexdigest()

    offset_mode = str(getattr(context.args, "offset_mode", "per_sample"))

    metadata = {
        "schema_version": "1.0.0",
        "mitoribopy_version": __version__,
        "subcommand": "rpf",
        "counts_path": "rpf_counts.tsv",
        "n_rows": n_rows,
        "n_samples": n_samples,
        "n_genes": n_genes,
        "total_reads": total_reads,
        "read_length_filter": {
            "requested_rpf_range": list(context.requested_rpf_range)
            or list(context.rpf_range),
            "observed_counts_by_length": {
                str(k): int(v) for k, v in sorted(
                    context.observed_counts_by_length.items()
                )
            },
            "retained_lengths": sorted(int(x) for x in context.rpf_range),
            "dropped_lengths": list(context.dropped_lengths),
            "threshold_rule": context.read_length_filter_threshold_rule,
            "threshold_value": context.read_length_filter_threshold_value,
        },
        "offset_mode": offset_mode,
        "offset_type": str(getattr(context.args, "offset_type", "5")),
        "offset_site": str(getattr(context.args, "offset_site", "p")),
        "footprint_class": str(getattr(context.args, "footprint_class", "")),
        "strain": str(getattr(context.args, "strain", "")),
        "reference_fasta": str(fasta_path) if fasta_path else None,
        "reference_checksum": reference_checksum,
    }

    output_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8"
    )
    return output_path


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

    offset_csv_path = context.csv_dir / f"offset_{context.args.align}.csv"
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
        per_sample_dir = context.csv_dir / "per_sample_offset"
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
    # Anchor the legend outside the plot area so per-sample lines stay
    # readable even when many samples crowd the curve.
    ax3.legend(
        loc="upper left", bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0, fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(str(output_path), format=plot_format, bbox_inches="tight")
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
    combined_offset_file = context.csv_dir / f"p_site_offsets_{context.args.align}.csv"
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
        sample_dir = context.csv_dir / "per_sample_offset"
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
        plot_dir=str(context.plot_subdir),
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
    # visible before the user reads any downstream output. The PNG/SVG
    # goes under plots/ but ``_plot_offset_drift`` also writes a sister
    # CSV next to it; we redirect that to csv_dir afterwards so the two
    # file types stay separated.
    if context.selected_offsets_by_sample:
        drift_plot = context.plot_subdir / f"offset_drift_{context.args.align}.{context.args.plot_format}"
        _plot_offset_drift(
            selected_offsets_by_sample=context.selected_offsets_by_sample,
            combined=context.selected_offsets_df,
            align_to=context.args.align,
            output_path=drift_plot,
            plot_format=context.args.plot_format,
        )
        sister_csv = drift_plot.with_suffix(".csv")
        if sister_csv.exists():
            target_csv = context.csv_dir / sister_csv.name
            sister_csv.replace(target_csv)

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


def _write_per_sample_offset_audit(context: PipelineContext) -> None:
    """Snapshot the offset row each sample will use to disk.

    Pure diagnostic so a reviewer can confirm from the filesystem that
    per-sample offset selection was honoured downstream. No effect on
    behaviour. The CSVs land at
    ``<plot_dir>/csv/per_sample_offset/<sample>/offset_applied.csv``.
    """
    if context.selected_offsets_df is None and not context.selected_offsets_by_sample:
        return
    audit_root = context.csv_dir / "per_sample_offset"
    audit_root.mkdir(parents=True, exist_ok=True)
    for sample_dir in context.sample_dirs:
        sample_name = Path(sample_dir).name
        applied = (
            context.selected_offsets_by_sample.get(sample_name)
            if context.selected_offsets_by_sample
            else None
        )
        if applied is None:
            applied = context.selected_offsets_df
        if applied is None:
            continue
        target = audit_root / sample_name / "offset_applied.csv"
        target.parent.mkdir(parents=True, exist_ok=True)
        applied.to_csv(target, index=False)


def run_downstream_modules(context: PipelineContext, emit_status: StatusWriter) -> None:
    """Run downstream analyses that consume the selected offsets."""
    ran_modules: list[str] = []
    skipped_modules: list[str] = []

    requested_sites: list[str] = []
    translation_profile_root = context.base_output_dir / "translation_profile"
    coverage_root = context.base_output_dir / "coverage_profile_plots"

    if context.selected_offsets_df is None:
        skipped_modules.extend(
            [
                "translation-profile analysis",
                "coverage-profile plots",
                "structure-density export",
            ]
        )
    else:
        analysis_sites = getattr(context.args, "analysis_sites", "both")
        if analysis_sites == "both":
            requested_sites = ["p", "a"]
        else:
            requested_sites = [analysis_sites]

        # Flat layout (v0.4.x):
        #   <output>/translation_profile/<sample>/{codon_usage,footprint_density,translating_frame}/
        #   <output>/coverage_profile_plots/{p_site,a_site}_density_*/, read_coverage_*/
        # Filenames carry the site prefix (``p_site_*`` / ``a_site_*``)
        # so a single per-sample folder can host both sites without
        # ambiguity.
        translation_profile_root.mkdir(parents=True, exist_ok=True)
        coverage_root.mkdir(parents=True, exist_ok=True)

        # Audit the offset row applied to each sample for reviewers.
        _write_per_sample_offset_audit(context)

        # Periodicity + frame + strand QC. Cheap (one pass over the
        # filtered BED with the per-sample offsets we just selected) and
        # the only routine output that lets a reviewer judge whether the
        # offset choices yielded a periodic profile, which is the
        # primary defensibility argument for any Ribo-seq run.
        from ..analysis.periodicity import run_periodicity_qc

        sample_names = [Path(s).name for s in context.sample_dirs]
        qc_dir = context.base_output_dir / "qc"
        # Periodicity-section knobs (when present in the YAML/CLI),
        # falling back to the spec-aligned defaults baked into
        # ``run_periodicity_qc``. Allows users to set
        # ``periodicity.enabled: false`` to skip the step cleanly,
        # or to tune thresholds / exclude windows / metric toggles
        # without touching code. Each knob is read with getattr so
        # existing call sites that don't supply it still work.
        period_enabled = bool(getattr(context.args, "periodicity_enabled", True))
        if period_enabled:
            period_kwargs: dict = {}
            for cli_attr, kw in (
                ("periodicity_fourier_window_nt", "fourier_window_nt"),
                ("periodicity_metagene_nt", "metagene_nt"),
                ("periodicity_metagene_normalize", "metagene_normalize"),
                ("periodicity_fourier_bootstrap_n", "fourier_n_bootstrap"),
                ("periodicity_fourier_permutations_n", "fourier_n_permutations"),
                ("periodicity_fourier_ci_alpha", "fourier_ci_alpha"),
                ("periodicity_fourier_random_seed", "fourier_random_seed"),
            ):
                val = getattr(context.args, cli_attr, None)
                if val is not None:
                    period_kwargs[kw] = val
            # The --periodicity-no-fourier-stats flag is the inverse of
            # `fourier_compute_stats`; map it explicitly.
            no_stats = bool(getattr(context.args, "periodicity_no_fourier_stats", False))
            if no_stats:
                period_kwargs["fourier_compute_stats"] = False
            run_periodicity_qc(
                bed_df=context.filtered_bed_df,
                annotation_df=context.annotation_df,
                samples=sample_names,
                selected_offsets_by_sample=context.selected_offsets_by_sample or None,
                selected_offsets_combined=context.selected_offsets_df,
                offset_type=str(context.args.offset_type),
                offset_site=str(context.args.offset_site),
                output_dir=qc_dir,
                **period_kwargs,
            )

        run_translation_profile_analysis(
            sample_dirs=context.sample_dirs,
            selected_offsets_df=context.selected_offsets_df,
            offset_type=context.args.offset_type,
            fasta_file=context.args.fasta,
            output_dir=str(translation_profile_root),
            args=context.args,
            annotation_df=context.annotation_df,
            filtered_bed_df=context.filtered_bed_df,
            resolved_codon_table=context.resolved_codon_table,
            resolved_start_codons=context.resolved_start_codons,
            selected_offsets_by_sample=context.selected_offsets_by_sample or None,
            requested_sites=requested_sites,
        )
        run_coverage_profile_plots(
            sample_dirs=context.sample_dirs,
            selected_offsets_df=context.selected_offsets_df,
            offset_type=context.args.offset_type,
            fasta_file=context.args.fasta,
            output_dir=str(coverage_root),
            read_coverage_dir=str(coverage_root),
            write_read_coverage_raw=getattr(context.args, "read_coverage_raw", True),
            write_read_coverage_rpm=getattr(context.args, "read_coverage_rpm", True),
            args=context.args,
            annotation_df=context.annotation_df,
            filtered_bed_df=context.filtered_bed_df,
            total_mrna_map=context.total_counts_map,
            selected_offsets_by_sample=context.selected_offsets_by_sample or None,
            requested_sites=requested_sites,
        )
        ran_modules.append(
            "translation-profile analysis ("
            + ", ".join("P-site" if s == "p" else "A-site" for s in requested_sites)
            + ")"
        )
        ran_modules.append(
            "coverage-profile plots ("
            + ", ".join("P-site" if s == "p" else "A-site" for s in requested_sites)
            + ")"
        )

        if context.args.structure_density:
            structure_density_dir = context.base_output_dir / "structure_density"
            structure_site = "p" if "p" in requested_sites else "a"
            site_column = "P_site" if structure_site == "p" else "A_site"
            run_structure_density_export(
                translation_profile_dir=str(translation_profile_root),
                structure_density_norm_perc=context.args.structure_density_norm_perc,
                output_dir=str(structure_density_dir),
                site_column=site_column,
            )
            ran_modules.append("structure-density export")
        else:
            skipped_modules.append("structure-density export")

        if getattr(context.args, "igv_export", False):
            igv_dir = context.base_output_dir / "igv_tracks"
            run_igv_export(
                translation_profile_dir=str(translation_profile_root),
                output_dir=str(igv_dir),
                sites=requested_sites,
            )
            ran_modules.append("IGV BedGraph export")
        else:
            skipped_modules.append("IGV BedGraph export")

    # Codon correlation: per-site, with explicit logging of every skip
    # path. Three of four skip branches were silent in v0.3.x — that is
    # the bug behind 'cor_plot: true did nothing'. Each branch now
    # emit_status's its reason so the user can diagnose from the log.
    if not context.args.cor_plot:
        emit_status("[PIPELINE] cor_plot is false; skipping codon correlation.")
        skipped_modules.append("codon correlation")
    elif not context.args.base_sample:
        emit_status(
            "[PIPELINE] cor_plot is true but no --base_sample was set; "
            "skipping codon correlation."
        )
        skipped_modules.append("codon correlation")
    elif context.selected_offsets_df is None:
        emit_status(
            "[PIPELINE] cor_plot requested but no offsets were selected; "
            "skipping codon correlation."
        )
        skipped_modules.append("codon correlation")
    else:
        sample_names = [Path(sample_dir).name for sample_dir in context.sample_dirs]
        if context.args.base_sample not in sample_names:
            emit_status(
                f"[PIPELINE] base_sample '{context.args.base_sample}' is not in "
                f"the sample list {sample_names}. Skipping codon correlation."
            )
            skipped_modules.append("codon correlation")
        else:
            cor_root = context.base_output_dir / "codon_correlation"
            ran_sites: list[str] = []
            for site in requested_sites:
                site_dir = "p_site" if site == "p" else "a_site"
                cor_out_dir = cor_root / site_dir
                run_codon_correlation(
                    translation_profile_dir=str(translation_profile_root),
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
                    site=site,
                    metric=getattr(context.args, "cor_metric", "log2_density_rpm"),
                    regression=getattr(context.args, "cor_regression", "theil_sen"),
                    support_min_raw=int(
                        getattr(context.args, "cor_support_min_raw", 10)
                    ),
                    label_top_n=int(
                        getattr(context.args, "cor_label_top_n", 10)
                    ),
                    pseudocount=getattr(context.args, "cor_pseudocount", "auto"),
                    raw_panel=getattr(context.args, "cor_raw_panel", "qc_only"),
                )
                ran_sites.append("P-site" if site == "p" else "A-site")
            if ran_sites:
                ran_modules.append(
                    "codon correlation (" + ", ".join(ran_sites) + ")"
                )
            else:
                emit_status(
                    "[PIPELINE] cor_plot requested but no analysis sites were "
                    "generated upstream; skipping codon correlation."
                )
                skipped_modules.append("codon correlation")

    summary_parts: list[str] = []
    if ran_modules:
        summary_parts.append("ran " + ", ".join(ran_modules))
    if skipped_modules:
        summary_parts.append("skipped " + ", ".join(skipped_modules))

    _emit_step_ok(7, "; ".join(summary_parts) + ".", emit_status)
