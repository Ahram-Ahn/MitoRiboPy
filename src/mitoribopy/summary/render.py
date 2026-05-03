"""Human-readable run summary renderer (P1.6).

Builds a Markdown document from the manifest + summary_qc rows that
gives a reviewer everything they need to judge a finished run at a
glance, with deep-link pointers to the underlying files.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Iterable


__all__ = ["render_summary_md"]


def _section(title: str, body: str) -> str:
    return f"## {title}\n\n{body}\n\n"


def _stage_status_line(stage: str, info: dict | None) -> str:
    if not info:
        return f"- **{stage}**: not configured"
    status = info.get("status", "unknown")
    runtime = info.get("runtime_seconds")
    reason = info.get("reason")
    bits = [f"`{status}`"]
    if runtime is not None:
        bits.append(f"runtime={runtime:.1f}s")
    if reason:
        bits.append(reason)
    return f"- **{stage}**: " + ", ".join(bits)


def _format_row(
    columns: Iterable[str], row: dict[str, str], widths: dict[str, int]
) -> str:
    cells = [str(row.get(col, "") or "").ljust(widths[col]) for col in columns]
    return "| " + " | ".join(cells) + " |"


def _markdown_table(
    columns: Iterable[str], rows: list[dict[str, str]]
) -> str:
    cols = list(columns)
    widths = {
        col: max(len(col), max((len(str(r.get(col) or "")) for r in rows), default=0))
        for col in cols
    }
    header = "| " + " | ".join(col.ljust(widths[col]) for col in cols) + " |"
    sep = "| " + " | ".join("-" * widths[col] for col in cols) + " |"
    body_lines = [_format_row(cols, r, widths) for r in rows]
    return "\n".join([header, sep, *body_lines])


def _format_periodicity_confidence_section(
    score_table_path: Path,
) -> str | None:
    """Read fourier_period3_score_combined.tsv and render a per-sample
    statistical-confidence table for the SUMMARY.md.

    Picks the headline row per sample — the ``combined`` gene_set,
    ``orf_start`` region, highest ``n_sites_total`` read length — and
    reports the spectral ratio with its bootstrap CI and permutation
    p-value. Returns ``None`` when the score TSV is missing or the
    headline row cannot be located.
    """
    if not score_table_path.is_file():
        return None
    try:
        import pandas as pd
    except ImportError:  # pragma: no cover - pandas is a hard dep
        return None
    try:
        score = pd.read_csv(score_table_path, sep="\t", comment="#")
    except (OSError, pd.errors.EmptyDataError, pd.errors.ParserError):
        return None
    if score.empty:
        return None

    sub = score[
        (score.get("gene_set") == "combined")
        & (score.get("region") == "orf_start")
    ].copy()
    if sub.empty:
        return None
    # Pick the read length with the most sites per sample.
    sub_sorted = sub.sort_values(
        ["sample", "n_sites_total"], ascending=[True, False],
    )
    headline = sub_sorted.drop_duplicates(subset=["sample"], keep="first")

    rows: list[dict[str, str]] = []
    for _, r in headline.iterrows():
        ratio = r.get("spectral_ratio_3nt")
        ci_lo = r.get("spectral_ratio_3nt_ci_low")
        ci_hi = r.get("spectral_ratio_3nt_ci_high")
        perm_p = r.get("permutation_p")
        snr = r.get("snr_call", "")
        if pd.notna(ci_lo) and pd.notna(ci_hi):
            ci_str = f"[{float(ci_lo):.2f}, {float(ci_hi):.2f}]"
        else:
            ci_str = "—"
        if pd.notna(perm_p):
            p_val = float(perm_p)
            p_str = "<0.001" if p_val < 0.001 else f"{p_val:.3f}"
        else:
            p_str = "—"
        rows.append({
            "sample": str(r.get("sample", "")),
            "read_length": str(int(r.get("read_length", 0))),
            "n_genes": str(int(r.get("n_genes", 0))),
            "spectral_ratio_3nt": (
                f"{float(ratio):.2f}x" if pd.notna(ratio) else "—"
            ),
            "ci_90pct": ci_str,
            "permutation_p": p_str,
            "snr_call": str(snr or "—"),
        })

    if not rows:
        return None
    columns = (
        "sample", "read_length", "n_genes",
        "spectral_ratio_3nt", "ci_90pct", "permutation_p", "snr_call",
    )
    return _markdown_table(columns, rows)


def render_summary_md(
    *,
    run_root: Path,
    manifest: dict,
    summary_qc_rows: list[dict[str, str]],
    warnings: list[dict] | None = None,
) -> str:
    """Return the rendered SUMMARY.md content for one finished run."""
    warnings = warnings or []

    # Header
    body = "# MitoRiboPy run summary\n\n"
    body += (
        f"- mitoribopy version: `{manifest.get('mitoribopy_version', '?')}`\n"
        f"- manifest schema: `{manifest.get('schema_version', '?')}`\n"
        f"- git commit: `{manifest.get('git_commit') or 'n/a'}`\n"
        f"- command: `{manifest.get('command', '?')}`\n"
        f"- config source: `{manifest.get('config_source', '?')}`\n"
    )
    sheet = manifest.get("sample_sheet")
    if sheet:
        body += f"- sample sheet: `{sheet}`\n"
    ref = manifest.get("reference_checksum")
    if ref:
        body += f"- reference checksum: `{ref[:16]}…`\n"
    body += "\n"

    # Stages
    stages = manifest.get("stages") or {}
    stage_lines = "\n".join(
        _stage_status_line(name, stages.get(name))
        for name in ("align", "rpf", "rnaseq")
    )
    body += _section("Stages", stage_lines)

    # Sample summary
    if summary_qc_rows:
        by_assay = Counter(row.get("assay", "") for row in summary_qc_rows)
        by_condition = Counter(row.get("condition", "") for row in summary_qc_rows)
        sample_lines = [
            f"- total active samples: {len(summary_qc_rows)}",
        ]
        for assay, n in sorted(by_assay.items()):
            sample_lines.append(f"  - assay `{assay}`: {n}")
        cond_lines = [
            f"  - condition `{cond or '(empty)'}`: {n}"
            for cond, n in sorted(by_condition.items())
        ]
        body += _section(
            "Samples",
            "\n".join(sample_lines + cond_lines),
        )

        # Per-sample QC table.
        qc_columns = (
            "sample_id", "assay", "condition", "kit_applied",
            "post_trim_reads", "mt_mrna_fraction", "qc_status", "qc_notes",
        )
        body += _section(
            "Per-sample QC (subset)",
            _markdown_table(qc_columns, summary_qc_rows),
        )

    # v0.9.0: Statistical confidence on the metagene Fourier QC.
    # Pulls the headline (sample, gene_set=combined, region=orf_start,
    # max-sites read_length) row out of fourier_period3_score_combined.tsv
    # and renders the spectral ratio with its bootstrap CI and the
    # circular-shift permutation p-value. Skipped silently when the
    # score TSV is missing or the headline row cannot be located (e.g.
    # very small fixtures).
    score_table = run_root / "rpf" / "qc" / "fourier_period3_score_combined.tsv"
    confidence_table = _format_periodicity_confidence_section(score_table)
    if confidence_table is not None:
        body += _section(
            "Periodicity statistical confidence",
            confidence_table
            + "\n\nColumns: spectral_ratio_3nt = amp(3) / median(amp at "
            "non-period-3 background); ci_90pct = 90 % bootstrap CI over "
            "genes; permutation_p = Laplace-smoothed circular-shift null "
            "p; snr_call = four-tier verdict (excellent ≥ 10×, "
            "healthy ≥ 5×, modest ≥ 2×, broken < 2×).",
        )

    # Warnings
    if warnings:
        warn_lines = []
        for w in warnings[:25]:
            sample = w.get("sample") or "-"
            code = w.get("code") or "-"
            warn_lines.append(
                f"- `{w.get('component', '?')}` `{code}` (sample `{sample}`): "
                f"{w.get('message', '')}"
            )
        if len(warnings) > 25:
            warn_lines.append(f"- … and {len(warnings) - 25} more (see warnings.tsv)")
        body += _section("Warnings", "\n".join(warn_lines))
    else:
        body += _section("Warnings", "_no structured warnings recorded_")

    # Output index
    output_lines = [
        f"- `{run_root.name}/run_manifest.json` — provenance + recorded hashes",
        f"- `{run_root.name}/canonical_config.yaml` — fully-resolved config",
        f"- `{run_root.name}/resource_plan.json` — CPU / memory / parallelism plan",
        f"- `{run_root.name}/outputs_index.tsv` — index of every produced TSV / JSON / plot",
        f"- `{run_root.name}/warnings.tsv` — structured warnings collected during the run",
        f"- `{run_root.name}/summary_qc.tsv` — per-sample QC roll-up",
        f"- `{run_root.name}/figure_qc.tsv` — per-plot mechanical QC (run by validate-figures)",
        f"- `{run_root.name}/align/` — per-sample read counts, kit resolution, BEDs",
        f"- `{run_root.name}/align/umi_qc.tsv` — UMI / dedup audit (when align ran)",
        f"- `{run_root.name}/rpf/` — RPF counts, offsets, P-site / A-site outputs",
        f"- `{run_root.name}/rpf/qc/` — metagene Fourier QC "
        "(fourier_spectrum_combined.tsv, fourier_period3_score_combined.tsv, "
        "metagene_start.tsv, metagene_stop.tsv, fourier_spectrum/<sample>/)",
        f"- `{run_root.name}/rpf/codon_correlation/{{p_site,a_site}}/` — "
        "codon-correlation metrics + sidecar metadata",
        f"- `{run_root.name}/rnaseq/` — TE / ΔTE tables and plots (when configured)",
    ]
    body += _section("Outputs", "\n".join(output_lines))

    schemas = manifest.get("output_schemas")
    if isinstance(schemas, dict):
        schema_lines = "\n".join(
            f"- `{k}`: {v}" for k, v in sorted(schemas.items())
        )
        body += _section("Output schema versions", schema_lines)

    return body
