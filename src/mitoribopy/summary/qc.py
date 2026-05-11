"""Per-sample QC roll-up writer (P1.7).

Reads the existing align / rpf / rnaseq stage outputs from a finished
``mitoribopy all`` run root and aggregates them into a single
``summary_qc.tsv`` whose columns are:

* identity / metadata: ``sample_id``, ``assay``, ``condition``
* align: ``kit_applied``, ``adapter_match_rate``, ``post_trim_reads``,
  ``contam_fraction``, ``mt_mrna_fraction``, ``mapq_kept_fraction``,
  ``dedup_removed_fraction``, ``umi_source``
* rpf: ``dominant_rpf_length``, ``frame0_fraction``, ``offset_mode``,
  ``offset_confidence``
* qc: ``qc_status`` (``pass`` / ``warn``), ``qc_notes`` (comma-joined
  reasons for warnings)

Every metric is best-effort: when a source file is missing or a column
is absent, the corresponding cell is empty. This keeps the writer
robust on partial runs and on tiny test fixtures that don't carry the
full output set.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable


__all__ = [
    "QC_THRESHOLDS",
    "SUMMARY_QC_COLUMNS",
    "build_summary_qc",
    "write_summary_qc",
]


SUMMARY_QC_COLUMNS: tuple[str, ...] = (
    "sample_id",
    "assay",
    "condition",
    "kit_applied",
    "adapter_match_rate",
    "post_trim_reads",
    "contam_fraction",
    "mt_mrna_fraction",
    "mapq_kept_fraction",
    "dedup_removed_fraction",
    "umi_source",
    "dominant_rpf_length",
    "frame0_fraction",
    "offset_mode",
    "offset_confidence",
    "qc_status",
    "qc_notes",
)


# Conservative defaults; the orchestrator can override per-run via the
# config later. The thresholds are directional only — falling below
# `mt_mrna_fraction_min` flags the sample, falling above
# `dedup_removed_fraction_max` flags it too (both are unusual).
QC_THRESHOLDS: dict[str, float] = {
    "post_trim_reads_min": 100_000,
    "mt_mrna_fraction_min": 0.05,
    "frame0_fraction_min": 0.45,
    "dedup_removed_fraction_max": 0.95,
}


# ---------- helpers ---------------------------------------------------------


def _read_tsv_with_comment_skip(path: Path) -> list[dict[str, str]]:
    """Read a tab-delimited file, skipping `#`-prefixed comment lines."""
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8") as handle:
        non_comment = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(non_comment, delimiter="\t")
        return list(reader)


def _safe_float(value: object) -> float | None:
    if value in (None, "", "NA"):
        return None
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def _safe_int(value: object) -> int | None:
    if value in (None, "", "NA"):
        return None
    try:
        return int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def _ratio(num: float | None, denom: float | None) -> float | None:
    """Return ``num/denom`` clamped to [0,1]; ``None`` when undefined."""
    if num is None or denom is None or denom == 0:
        return None
    if denom < 0 or num < 0:
        return None
    return max(0.0, min(1.0, num / denom))


# ---------- per-stage extractors --------------------------------------------


def _load_sample_metadata(run_root: Path, manifest: dict) -> list[dict[str, str]]:
    """Return per-sample identity rows: sample_id, assay, condition.

    Source: the unified sample sheet recorded in the manifest. Falls
    back to deriving identities from align/read_counts.tsv when no
    sample sheet is set.
    """
    sheet_path = manifest.get("sample_sheet")
    if sheet_path:
        try:
            from ..sample_sheet import load_sample_sheet

            sheet = load_sample_sheet(sheet_path)
        except Exception:
            sheet = None
        if sheet is not None:
            return [
                {
                    "sample_id": r.sample_id,
                    "assay": r.assay,
                    "condition": r.condition,
                }
                for r in sheet.active()
            ]

    # Fallback: derive from align/read_counts.tsv
    rows = _read_tsv_with_comment_skip(run_root / "align" / "read_counts.tsv")
    return [
        {
            "sample_id": r.get("sample", ""),
            "assay": "ribo",
            "condition": "",
        }
        for r in rows
    ]


def _load_align_metrics(run_root: Path) -> dict[str, dict[str, str]]:
    """Return ``{sample_id: {metric: value_str}}`` from align outputs."""
    counts_rows = _read_tsv_with_comment_skip(
        run_root / "align" / "read_counts.tsv"
    )
    kit_rows = _read_tsv_with_comment_skip(
        run_root / "align" / "kit_resolution.tsv"
    )
    by_sample: dict[str, dict[str, str]] = {}

    for row in counts_rows:
        sample = row.get("sample") or ""
        if not sample:
            continue
        post_trim = _safe_int(row.get("post_trim"))
        post_rrna = _safe_int(row.get("post_rrna_filter"))
        mt_aligned = _safe_int(row.get("mt_aligned"))
        mt_after_mapq = _safe_int(row.get("mt_aligned_after_mapq"))
        mt_after_dedup = _safe_int(row.get("mt_aligned_after_dedup"))

        contam_frac = _ratio(
            (post_trim or 0) - (post_rrna or 0)
            if post_trim is not None and post_rrna is not None
            else None,
            post_trim,
        )
        mt_frac = _ratio(mt_aligned, post_rrna)
        mapq_frac = _ratio(mt_after_mapq, mt_aligned)
        dedup_removed = (
            None
            if (mt_after_mapq is None or mt_after_dedup is None or mt_after_mapq == 0)
            else max(0.0, 1.0 - mt_after_dedup / mt_after_mapq)
        )

        by_sample[sample] = {
            "post_trim_reads": str(post_trim) if post_trim is not None else "",
            "contam_fraction": _fmt_frac(contam_frac),
            "mt_mrna_fraction": _fmt_frac(mt_frac),
            "mapq_kept_fraction": _fmt_frac(mapq_frac),
            "dedup_removed_fraction": _fmt_frac(dedup_removed),
        }

    for row in kit_rows:
        sample = row.get("sample") or ""
        entry = by_sample.setdefault(sample, {})
        entry["kit_applied"] = row.get("applied_kit") or ""
        entry["adapter_match_rate"] = row.get("detection_match_rate") or ""
        entry["umi_source"] = row.get("umi_source") or ""

    return by_sample


def _fmt_frac(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.4f}"


def _load_rpf_metrics(run_root: Path) -> dict[str, dict[str, str]]:
    """Return ``{sample_id: {metric: value_str}}`` from rpf outputs.

    Best-effort: rpf does not currently expose a single per-sample QC
    table, so we surface what is reliably available — the ``offset_mode``
    from rpf's ``run_settings.json`` and (when present) the per-sample
    offset table's ``confidence`` column.
    """
    settings_path = run_root / "rpf" / "run_settings.json"
    if not settings_path.is_file():
        return {}
    try:
        import json

        settings = json.loads(settings_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    offset_mode = settings.get("offset_mode") or settings.get("offset_pick_reference") or ""

    by_sample: dict[str, dict[str, str]] = {}
    # Per-sample offset confidence: probe the standard offset diagnostics
    # CSV when it exists. Layout is open-ended; we only consume two
    # columns (``sample`` and ``confidence``).
    diag_dir = run_root / "rpf" / "offset_diagnostics" / "csv"
    confidence_csv = diag_dir / "offset_confidence.csv"
    if confidence_csv.is_file():
        with confidence_csv.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                sample = row.get("sample") or ""
                if not sample:
                    continue
                by_sample.setdefault(sample, {})["offset_confidence"] = (
                    row.get("confidence") or ""
                )
    if offset_mode:
        for entry in by_sample.values():
            entry.setdefault("offset_mode", offset_mode)
        # Make sure samples that had no offset row still show the mode.
        # Caller fills missing entries downstream.
    return by_sample


# ---------- top-level builder + writer --------------------------------------


def _classify_qc(row: dict[str, str]) -> tuple[str, str]:
    """Return ``(qc_status, qc_notes)`` based on QC_THRESHOLDS."""
    notes: list[str] = []
    post_trim = _safe_int(row.get("post_trim_reads"))
    if post_trim is not None and post_trim < QC_THRESHOLDS["post_trim_reads_min"]:
        notes.append(f"post_trim_reads<{QC_THRESHOLDS['post_trim_reads_min']}")
    mt_frac = _safe_float(row.get("mt_mrna_fraction"))
    if mt_frac is not None and mt_frac < QC_THRESHOLDS["mt_mrna_fraction_min"]:
        notes.append(f"mt_mrna_fraction<{QC_THRESHOLDS['mt_mrna_fraction_min']}")
    dedup_removed = _safe_float(row.get("dedup_removed_fraction"))
    if (
        dedup_removed is not None
        and dedup_removed > QC_THRESHOLDS["dedup_removed_fraction_max"]
    ):
        notes.append(
            f"dedup_removed_fraction>{QC_THRESHOLDS['dedup_removed_fraction_max']}"
        )
    return ("warn" if notes else "pass", ",".join(notes))


def build_summary_qc(
    run_root: Path,
    manifest: dict,
) -> list[dict[str, str]]:
    """Return a list of per-sample QC rows for *run_root*.

    The manifest is the live ``run_manifest.json`` dict (after
    serialisation by ``mitoribopy all``). Used here to discover the
    sample sheet path; the remaining metrics come from the per-stage
    output TSVs on disk.
    """
    run_root = Path(run_root)
    sample_metadata = _load_sample_metadata(run_root, manifest)
    align_metrics = _load_align_metrics(run_root)
    rpf_metrics = _load_rpf_metrics(run_root)

    rows: list[dict[str, str]] = []
    for meta in sample_metadata:
        sid = meta["sample_id"]
        align = align_metrics.get(sid, {})
        rpf = rpf_metrics.get(sid, {})
        merged = {**meta, **align, **rpf}
        # Fill missing columns with empty strings so the TSV layout is stable.
        for col in SUMMARY_QC_COLUMNS:
            merged.setdefault(col, "")
        merged["qc_status"], merged["qc_notes"] = _classify_qc(merged)
        rows.append(merged)
    return rows


def write_summary_qc(
    rows: Iterable[dict[str, str]],
    output: Path,
) -> Path:
    """Write *rows* to ``output`` with the schema-version header line."""
    from ..io.schema_versions import schema_header_line

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    with output.open("w", encoding="utf-8") as handle:
        handle.write(schema_header_line("summary_qc.tsv"))
        handle.write("\t".join(SUMMARY_QC_COLUMNS) + "\n")
        for row in rows:
            handle.write(
                "\t".join(str(row.get(col, "") or "") for col in SUMMARY_QC_COLUMNS)
                + "\n"
            )
    return output
