"""Mechanical QC for MitoRiboPy plot outputs.

Powers ``mitoribopy validate-figures`` (assessment §4 / §6). A finished
run root is walked for every PNG / SVG produced by the rpf and rnaseq
stages; each plot is scored for a small set of objective failure modes:

* point-count drift between the plot's metadata sidecar and the
  source TSV that fed it (off-by-one, accidentally dropped genes);
* SVG text editability (the publication promise of editable Illustrator
  outputs is null if ``svg.fonttype`` regressed to ``path``);
* PNG dpi ≥ 300 (hard publication requirement).

The validator is intentionally conservative — it only flags conditions
that *can* be checked from disk without re-running matplotlib. Things
like label-label overlap that need a live ``Figure`` object are covered
in unit tests under :mod:`tests.test_figure_validator`; when a plot
records ``label_overlap_count`` in its sidecar, the validator promotes
that to a fail/warn here.

Output:

* ``figure_qc.tsv`` at the run root (schema 1.0); columns documented in
  :data:`_FIGURE_QC_COLUMNS`.
* Each problem also surfaces through
  :func:`mitoribopy.io.warnings_log.record` so it shows up in
  ``warnings.tsv`` and the manifest's ``warnings`` array.

Sidecar shape (written by :func:`write_plot_metadata`)::

    {
      "plot_type": "delta_te_volcano",
      "stage": "rnaseq",
      "source_data": "rnaseq/delta_te.tsv",
      "n_points_expected": 13,
      "n_points_drawn": 13,
      "n_labels": 13,
      "labels_drawn": 13,
      "label_policy": "all_mt_genes",
      "label_overlap_count": 0,
      "label_outside_axes_count": 0,
      "palette": "Okabe-Ito",
      "formats": ["png", "svg"],
      "dpi": 300,
      "min_font_size": 8
    }
"""

from __future__ import annotations

import json
import struct
import zlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

from ..io.schema_versions import schema_header_line


__all__ = [
    "FigureRecord",
    "build_figure_qc_rows",
    "discover_plots",
    "load_plot_metadata",
    "metadata_sidecar_path",
    "validate_figures",
    "write_figure_qc",
    "write_plot_metadata",
]


_FIGURE_QC_COLUMNS: tuple[str, ...] = (
    "plot_path",
    "stage",
    "plot_type",
    "status",                  # pass | warn | fail
    "n_points_expected",
    "n_points_drawn",
    "n_labels",
    "label_overlap_count",
    "label_outside_axes_count",
    "legend_overlap",
    "stat_box_overlap",
    "clipped_text",
    "min_font_size",
    "svg_text_editable",
    "png_dpi",
    "has_source_data",
    "warnings",
)


# ---------------------------------------------------------------------------
# Plot metadata sidecar
# ---------------------------------------------------------------------------


def metadata_sidecar_path(plot_path: Path | str) -> Path:
    """Return the sibling ``.metadata.json`` path for a plot file.

    Both ``foo.png`` and ``foo.svg`` produced by the same rendering call
    share the same sidecar at ``foo.metadata.json`` so the validator
    only has to look once per plot stem.
    """
    plot_path = Path(plot_path)
    return plot_path.with_suffix(".metadata.json")


def write_plot_metadata(
    plot_path: Path | str,
    *,
    plot_type: str,
    stage: str,
    source_data: str | None = None,
    n_points_expected: int | None = None,
    n_points_drawn: int | None = None,
    n_labels: int | None = None,
    labels_drawn: int | None = None,
    label_policy: str | None = None,
    label_overlap_count: int | None = None,
    label_outside_axes_count: int | None = None,
    palette: str = "Okabe-Ito",
    formats: Sequence[str] = ("png", "svg"),
    dpi: int = 300,
    min_font_size: int | None = None,
    extra: dict | None = None,
) -> Path:
    """Write a ``.metadata.json`` sidecar next to *plot_path*.

    Returns the sidecar path. Existing sidecars are overwritten.

    The sidecar is the validator's contract: if ``n_points_expected``
    and ``n_points_drawn`` are both populated they must match, and if
    ``label_overlap_count`` is populated it must be zero. Sidecars
    without those fields are still valid but cannot be mechanically
    validated for that property.
    """
    plot_path = Path(plot_path)
    payload: dict = {
        "plot_type": plot_type,
        "stage": stage,
    }
    if source_data is not None:
        payload["source_data"] = source_data
    if n_points_expected is not None:
        payload["n_points_expected"] = int(n_points_expected)
    if n_points_drawn is not None:
        payload["n_points_drawn"] = int(n_points_drawn)
    if n_labels is not None:
        payload["n_labels"] = int(n_labels)
    if labels_drawn is not None:
        payload["labels_drawn"] = int(labels_drawn)
    if label_policy is not None:
        payload["label_policy"] = label_policy
    if label_overlap_count is not None:
        payload["label_overlap_count"] = int(label_overlap_count)
    if label_outside_axes_count is not None:
        payload["label_outside_axes_count"] = int(label_outside_axes_count)
    if min_font_size is not None:
        payload["min_font_size"] = int(min_font_size)
    payload["palette"] = palette
    payload["formats"] = list(formats)
    payload["dpi"] = int(dpi)
    if extra:
        payload.update(extra)

    sidecar = metadata_sidecar_path(plot_path)
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    sidecar.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return sidecar


def load_plot_metadata(plot_path: Path | str) -> dict | None:
    """Return the parsed sidecar JSON for *plot_path*, or ``None`` if missing."""
    sidecar = metadata_sidecar_path(plot_path)
    if not sidecar.is_file():
        return None
    try:
        return json.loads(sidecar.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


# ---------------------------------------------------------------------------
# Plot discovery
# ---------------------------------------------------------------------------


_PLOT_SUFFIXES: tuple[str, ...] = (".png", ".svg")


def discover_plots(run_root: Path | str) -> list[Path]:
    """Return every PNG / SVG plot under *run_root*, deduplicated by stem.

    When both ``foo.png`` and ``foo.svg`` exist for the same plot
    (the common case for publication-style outputs), only the PNG is
    returned; the SVG is checked alongside it during validation via
    :func:`metadata_sidecar_path` and the explicit SVG-text check.
    """
    run_root = Path(run_root)
    by_stem: dict[Path, Path] = {}
    for path in sorted(run_root.rglob("*")):
        if path.suffix.lower() not in _PLOT_SUFFIXES:
            continue
        if not path.is_file():
            continue
        stem_key = path.with_suffix("")
        # Prefer .png as the "primary" entry — its presence is the
        # publication-grade requirement (svg sidecar is checked
        # separately for editability).
        existing = by_stem.get(stem_key)
        if existing is None or path.suffix.lower() == ".png":
            by_stem[stem_key] = path
    return sorted(by_stem.values())


# ---------------------------------------------------------------------------
# Per-plot mechanical checks
# ---------------------------------------------------------------------------


def _png_dpi(path: Path) -> int | None:
    """Best-effort PNG dpi extraction from the pHYs chunk.

    Returns the x-axis dpi rounded to the nearest int. ``None`` when the
    file is not a valid PNG, the pHYs chunk is missing, or the unit is
    not 'metres' (PNG only encodes dpi via metres-per-unit).
    """
    try:
        with path.open("rb") as handle:
            sig = handle.read(8)
            if sig != b"\x89PNG\r\n\x1a\n":
                return None
            while True:
                length_bytes = handle.read(4)
                if len(length_bytes) < 4:
                    return None
                length = struct.unpack(">I", length_bytes)[0]
                chunk_type = handle.read(4)
                if not chunk_type:
                    return None
                data = handle.read(length)
                # Skip CRC.
                handle.read(4)
                if chunk_type == b"pHYs" and len(data) == 9:
                    px_per_unit_x = struct.unpack(">I", data[0:4])[0]
                    unit = data[8]
                    if unit == 1:  # metres
                        return int(round(px_per_unit_x * 0.0254))
                    return None
                if chunk_type == b"IEND":
                    return None
    except (OSError, struct.error):
        return None


def _svg_text_editable(path: Path) -> bool | None:
    """Return True iff the SVG contains real ``<text>`` elements.

    matplotlib's ``svg.fonttype="none"`` setting writes editable text
    elements; ``"path"`` (the matplotlib default outside our publication
    rc) converts text to glyph paths, breaking Illustrator editability.
    Returns ``None`` for plain (uncompressed) failures or when the file
    is not an SVG.
    """
    try:
        head = path.read_bytes()
    except OSError:
        return None
    # Handle gzipped svg sidecars (.svgz) too.
    if path.suffix.lower() == ".svgz" or head[:2] == b"\x1f\x8b":
        try:
            head = zlib.decompress(head, 16 + zlib.MAX_WBITS)
        except zlib.error:
            return None
    text = head.decode("utf-8", errors="replace")
    if "<svg" not in text:
        return None
    return "<text" in text


@dataclass
class FigureRecord:
    plot_path: Path
    stage: str = ""
    plot_type: str = ""
    status: str = "pass"
    n_points_expected: int | None = None
    n_points_drawn: int | None = None
    n_labels: int | None = None
    label_overlap_count: int | None = None
    label_outside_axes_count: int | None = None
    legend_overlap: int | None = None
    stat_box_overlap: int | None = None
    clipped_text: int | None = None
    min_font_size: int | None = None
    svg_text_editable: bool | None = None
    png_dpi: int | None = None
    has_source_data: bool = False
    warnings: list[str] = field(default_factory=list)

    def downgrade(self, target: str) -> None:
        """Move ``status`` toward fail (pass < warn < fail)."""
        order = {"pass": 0, "warn": 1, "fail": 2}
        if order.get(target, 0) > order.get(self.status, 0):
            self.status = target

    def add_warning(self, msg: str) -> None:
        self.warnings.append(msg)


def _stage_from_path(plot_path: Path, run_root: Path) -> str:
    try:
        rel = plot_path.relative_to(run_root)
    except ValueError:
        return ""
    parts = rel.parts
    if not parts:
        return ""
    head = parts[0]
    return head if head in {"align", "rpf", "rnaseq"} else ""


def _check_one_plot(
    plot_path: Path,
    *,
    run_root: Path,
    require_png_dpi: int = 300,
    strict: bool = False,
) -> FigureRecord:
    """Score a single plot path against the mechanical checks."""
    record = FigureRecord(plot_path=plot_path)
    record.stage = _stage_from_path(plot_path, run_root)

    # Sidecar - optional. Older RPF plot families do not write per-plot
    # metadata yet; missing sidecars should not swamp warnings.tsv.
    metadata = load_plot_metadata(plot_path)
    if metadata is not None:
        record.plot_type = str(metadata.get("plot_type") or "")
        record.stage = record.stage or str(metadata.get("stage") or "")
        record.n_points_expected = metadata.get("n_points_expected")
        record.n_points_drawn = metadata.get("n_points_drawn")
        record.n_labels = (
            metadata.get("n_labels") or metadata.get("labels_drawn")
        )
        record.label_overlap_count = metadata.get("label_overlap_count")
        record.label_outside_axes_count = metadata.get(
            "label_outside_axes_count"
        )
        record.min_font_size = metadata.get("min_font_size")
        record.has_source_data = bool(metadata.get("source_data"))

        # Hard contracts.
        n_exp = record.n_points_expected
        n_drw = record.n_points_drawn
        if n_exp is not None and n_drw is not None and n_exp != n_drw:
            record.add_warning(
                f"point-count mismatch: expected={n_exp} drawn={n_drw}"
            )
            record.downgrade("fail")
        if record.label_overlap_count and record.label_overlap_count > 0:
            record.add_warning(
                f"label_overlap_count={record.label_overlap_count}"
            )
            record.downgrade("fail")
        if (
            record.label_outside_axes_count
            and record.label_outside_axes_count > 0
        ):
            record.add_warning(
                f"label_outside_axes_count={record.label_outside_axes_count}"
            )
            record.downgrade("fail")
        if record.min_font_size is not None and record.min_font_size < 7:
            record.add_warning(
                f"min_font_size={record.min_font_size} below publication minimum (7pt)"
            )
            record.downgrade("warn")

    # PNG dpi check.
    if plot_path.suffix.lower() == ".png":
        dpi = _png_dpi(plot_path)
        record.png_dpi = dpi
        if dpi is None:
            record.add_warning("png dpi unreadable")
            record.downgrade("warn")
        elif dpi < require_png_dpi:
            record.add_warning(
                f"png_dpi={dpi} < required {require_png_dpi}"
            )
            record.downgrade("warn" if not strict else "fail")

    # SVG editability check (the SVG sibling of a PNG plot, or the SVG
    # itself when there is no PNG).
    svg_path = plot_path.with_suffix(".svg")
    if svg_path.is_file():
        editable = _svg_text_editable(svg_path)
        record.svg_text_editable = editable
        if editable is False:
            record.add_warning(
                f"svg text not editable: {svg_path.name}"
            )
            record.downgrade("warn" if not strict else "fail")

    return record


# ---------------------------------------------------------------------------
# Top-level validate + writers
# ---------------------------------------------------------------------------


def validate_figures(
    run_root: Path | str,
    *,
    strict: bool = False,
    require_png_dpi: int = 300,
) -> list[FigureRecord]:
    """Score every plot under *run_root* and return per-plot records.

    ``strict=True`` upgrades several warn-level conditions (low PNG dpi,
    SVG text non-editable) to fail so the CLI exits non-zero. Hard
    contract failures (point-count mismatch, label overlap > 0) are
    always fail regardless of ``strict``.
    """
    run_root = Path(run_root)
    records: list[FigureRecord] = []
    for plot_path in discover_plots(run_root):
        records.append(
            _check_one_plot(
                plot_path,
                run_root=run_root,
                require_png_dpi=require_png_dpi,
                strict=strict,
            )
        )
    return records


def build_figure_qc_rows(
    records: Iterable[FigureRecord],
    *,
    run_root: Path | str | None = None,
) -> list[dict[str, str]]:
    """Convert :class:`FigureRecord` instances to ``figure_qc.tsv`` rows.

    Paths are made run-root-relative when *run_root* is supplied.
    """
    rows: list[dict[str, str]] = []
    root = Path(run_root) if run_root is not None else None

    def _fmt_path(p: Path) -> str:
        if root is None:
            return str(p)
        try:
            return str(p.relative_to(root))
        except ValueError:
            return str(p)

    for r in records:
        rows.append({
            "plot_path": _fmt_path(r.plot_path),
            "stage": r.stage or "",
            "plot_type": r.plot_type or "",
            "status": r.status,
            "n_points_expected": _opt(r.n_points_expected),
            "n_points_drawn": _opt(r.n_points_drawn),
            "n_labels": _opt(r.n_labels),
            "label_overlap_count": _opt(r.label_overlap_count),
            "label_outside_axes_count": _opt(r.label_outside_axes_count),
            "legend_overlap": _opt(r.legend_overlap),
            "stat_box_overlap": _opt(r.stat_box_overlap),
            "clipped_text": _opt(r.clipped_text),
            "min_font_size": _opt(r.min_font_size),
            "svg_text_editable": _opt_bool(r.svg_text_editable),
            "png_dpi": _opt(r.png_dpi),
            "has_source_data": "1" if r.has_source_data else "0",
            "warnings": "; ".join(r.warnings),
        })
    return rows


def _opt(value) -> str:
    return "" if value is None else str(value)


def _opt_bool(value: bool | None) -> str:
    if value is None:
        return ""
    return "true" if value else "false"


def _scrub_cell(value: str) -> str:
    return value.replace("\t", " ").replace("\n", " ")


def write_figure_qc(
    rows: Iterable[dict[str, str]], path: Path | str
) -> Path:
    """Write *rows* to a ``figure_qc.tsv`` at *path*.

    The file is overwritten on each call. An empty input still produces
    a header-only file so ``mitoribopy validate-figures`` always leaves
    the path on disk.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(schema_header_line("figure_qc.tsv"))
        handle.write("\t".join(_FIGURE_QC_COLUMNS) + "\n")
        for row in rows:
            handle.write(
                "\t".join(_scrub_cell(row.get(c, "")) for c in _FIGURE_QC_COLUMNS)
                + "\n"
            )
    return path
