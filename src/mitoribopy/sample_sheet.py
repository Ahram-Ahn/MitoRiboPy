"""Unified per-project sample sheet.

A single TSV that replaces the historical pair of stage-specific tables
(`--sample-overrides` for align, `--condition-map` for rnaseq). The
sheet declares every sample's identity, assay, condition, FASTQ paths,
and (optionally) per-sample kit / UMI / strandedness overrides in one
place; downstream stages derive their per-stage inputs from it.

Schema
------

Required columns (header row, tab-separated)::

    sample_id    unique sample identifier; matches FASTQ basename for
                 align's per-sample resolver
    assay        'ribo' or 'rna'
    condition    free-form condition label (drives DE contrast)
    fastq_1      path to R1 FASTQ (or the SE FASTQ)

Optional columns (any subset; missing columns default to ``None``)::

    replicate              free-form replicate label ('1', '2', 'A', ...);
                           informational only — pairing is by sample_id
    biological_sample_id   logical biological sample id (groups Ribo+RNA
                           replicates that share a biological replicate)
    library_type           'single_end' | 'paired_end' | 'auto'
    fastq_2                path to R2 FASTQ for paired-end reads
    adapter                per-sample 3' adapter sequence override
    pretrimmed             'true' / 'false'; mark a sample's FASTQ as
                           already adapter-trimmed (cutadapt skips -a)
    umi_length             integer; per-sample UMI length override
    umi_position           '5p', '3p', or 'both' (dual-end UMI)
    umi_length_5p          integer; per-end 5' UMI length when
                           ``umi_position='both'``
    umi_length_3p          integer; per-end 3' UMI length when
                           ``umi_position='both'``
    strandedness           'forward' | 'reverse' | 'unstranded'
    dedup_strategy         'auto' | 'umi_coordinate' | 'skip'
                           ('umi-tools' / 'umi_tools' are accepted as
                           legacy aliases and normalised to
                           'umi_coordinate' at parse time)
    read_length_min        integer; per-sample minimum RPF length
    read_length_max        integer; per-sample maximum RPF length
    reference_id           logical reference identifier (e.g.
                           'human_mt_rCRS_v1') — informational, used
                           to flag drift across samples
    include                'true' / 'false' / blank — when present,
                           overrides ``exclude`` (true = keep the row).
                           Recognised so authors can use whichever
                           polarity they prefer.
    exclude                'true' / 'false' / blank — when true, the row
                           is skipped (use to drop a bad library without
                           deleting the row)
    notes                  free-form

Empty cells (``""``, ``"NA"``, ``"None"``, ``"-"``) are read as
``None``. Lookups are case-sensitive. ``assay`` is canonicalised to
lowercase.

Validation rules
----------------

* The header row must contain every required column.
* ``sample_id`` must be unique across the whole sheet (after exclusion
  rows are NOT yet dropped — exclusions run on a *materialised* row).
* ``assay`` must be ``'ribo'`` or ``'rna'``.
* ``umi_length`` must parse as a non-negative integer when set.
* ``strandedness`` and ``dedup_strategy``, when set, must use the
  literal values declared by :mod:`mitoribopy.align._types`.

The parser is deliberately strict: invalid sheets fail at load time
with a single error listing every problem, so users do not see partial
runs that silently dropped a sample.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Literal


__all__ = [
    "Assay",
    "SampleRow",
    "SampleSheet",
    "SampleSheetError",
    "check_sheet_conflicts",
    "format_sheet_conflict_error",
    "load_sample_sheet",
]


def check_sheet_conflicts(
    stage_cfg: dict | None,
    *,
    conflict_keys: Iterable[str],
) -> list[str]:
    """Return the list of conflict_keys present (with truthy values) in stage_cfg.

    Used by the orchestrator and per-stage CLIs to detect when a user
    declared the unified sample sheet AND a per-stage input that the
    sheet supersedes (e.g. ``samples:`` + ``align.fastq_dir:``). The
    caller decides how to format the error.

    Tolerates ``stage_cfg=None`` (an empty YAML stanza like ``align:``
    parses to ``None``) by reporting no conflicts.
    """
    if not stage_cfg:
        return []
    return [k for k in conflict_keys if stage_cfg.get(k)]


def format_sheet_conflict_error(
    stage_label: str,
    conflicts: list[str],
    *,
    sheet_label: str = "samples:",
) -> str:
    """Render a uniform error string for sheet-vs-per-stage conflicts.

    Used by the orchestrator (``mitoribopy all``) and standalone
    ``mitoribopy rnaseq`` so users see the same wording regardless of
    entry point.
    """
    return (
        f"[{stage_label}] ERROR: '{sheet_label}' is set, but the following "
        f"per-stage input(s) are also declared and would be shadowed: "
        + ", ".join(repr(c) for c in conflicts)
        + ". The unified sample sheet supersedes per-stage inputs; drop "
        "either the sheet or the conflicting flag(s)."
    )


Assay = Literal["ribo", "rna"]

_REQUIRED_COLUMNS: tuple[str, ...] = (
    "sample_id",
    "assay",
    "condition",
    "fastq_1",
)

_OPTIONAL_COLUMNS: tuple[str, ...] = (
    "replicate",
    "biological_sample_id",
    "library_type",
    "fastq_2",
    "adapter",
    "pretrimmed",
    "umi_length",
    "umi_position",
    "umi_length_5p",
    "umi_length_3p",
    "strandedness",
    "dedup_strategy",
    "read_length_min",
    "read_length_max",
    "reference_id",
    "include",
    "exclude",
    "notes",
)

_KNOWN_COLUMNS: frozenset[str] = frozenset(_REQUIRED_COLUMNS + _OPTIONAL_COLUMNS)

_VALID_ASSAYS: frozenset[str] = frozenset({"ribo", "rna"})
_VALID_STRANDEDNESS: frozenset[str] = frozenset(
    {"forward", "reverse", "unstranded"}
)
_VALID_DEDUP: frozenset[str] = frozenset(
    {
        "auto",
        # Canonical token (matches the align CLI's preferred form, see
        # ``mitoribopy.cli.align.build_parser`` --dedup-strategy choices).
        "umi_coordinate",
        # Legacy aliases — accepted for back-compat and normalised to
        # ``umi_coordinate`` by ``_normalise_dedup`` before the value
        # reaches the SampleRow / canonical_config / run_manifest.
        "umi-tools",
        "umi_tools",
        "skip",
    }
)


def _normalise_dedup(value: str | None) -> str | None:
    """Collapse legacy ``umi-tools`` / ``umi_tools`` to canonical ``umi_coordinate``.

    Other recognised tokens (``auto``, ``umi_coordinate``, ``skip``) are
    returned unchanged. Returning ``None`` for ``None`` keeps the
    column-not-set semantic intact.
    """
    if value is None:
        return None
    lowered = value.strip().lower()
    if lowered in {"umi-tools", "umi_tools"}:
        return "umi_coordinate"
    return lowered
_VALID_UMI_POS: frozenset[str] = frozenset({"5p", "3p", "both"})
_VALID_LIBRARY_TYPE: frozenset[str] = frozenset(
    {"single_end", "paired_end", "auto"}
)

_NULL_TOKENS: frozenset[str] = frozenset({"", "NA", "None", "none", "-", "null"})
_BOOL_TRUE: frozenset[str] = frozenset({"true", "yes", "1", "y", "t"})
_BOOL_FALSE: frozenset[str] = frozenset({"false", "no", "0", "n", "f"})


class SampleSheetError(ValueError):
    """Raised when a sample sheet fails validation at load time."""


@dataclass(frozen=True)
class SampleRow:
    """One materialised, validated row of the sample sheet.

    All fields are normalised: empty cells are ``None``, ``assay`` is
    lowercased, ``umi_length`` is an int, ``exclude`` is a bool.
    """

    sample_id: str
    assay: Assay
    condition: str
    fastq_1: Path
    replicate: str | None = None
    biological_sample_id: str | None = None
    library_type: str | None = None
    fastq_2: Path | None = None
    adapter: str | None = None
    pretrimmed: bool | None = None
    umi_length: int | None = None
    umi_position: str | None = None
    umi_length_5p: int | None = None
    umi_length_3p: int | None = None
    strandedness: str | None = None
    dedup_strategy: str | None = None
    read_length_min: int | None = None
    read_length_max: int | None = None
    reference_id: str | None = None
    exclude: bool = False
    notes: str | None = None


@dataclass(frozen=True)
class SampleSheet:
    """The full sheet, post-validation.

    Use :meth:`active` to iterate rows that should actually feed the
    pipeline (i.e. ``exclude=False``); use :attr:`rows` for the raw
    table including excluded rows.
    """

    path: Path
    rows: tuple[SampleRow, ...]
    columns: tuple[str, ...] = field(default_factory=tuple)

    def active(self) -> tuple[SampleRow, ...]:
        return tuple(row for row in self.rows if not row.exclude)

    def by_assay(self, assay: Assay) -> tuple[SampleRow, ...]:
        return tuple(row for row in self.active() if row.assay == assay)

    def condition_map(self) -> dict[str, str]:
        """``sample_id -> condition`` for active rows only."""
        return {row.sample_id: row.condition for row in self.active()}

    def fastq_paths(self, assay: Assay) -> list[Path]:
        """Flat list of R1 (and R2 when set) paths for the given assay,
        in sheet order, excluding rows marked ``exclude=true``."""
        out: list[Path] = []
        for row in self.by_assay(assay):
            out.append(row.fastq_1)
            if row.fastq_2 is not None:
                out.append(row.fastq_2)
        return out


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------


def load_sample_sheet(path: str | Path) -> SampleSheet:
    """Read and validate a unified sample sheet TSV.

    Raises
    ------
    SampleSheetError
        On any structural problem (missing required column, duplicate
        ``sample_id``, malformed ``umi_length``, unknown ``assay``,
        etc.). The error message lists every problem found in one
        pass so users can fix the sheet without iterating.
    """
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise SampleSheetError(
            f"sample sheet: cannot read {path}: {exc}"
        ) from exc

    lines = [line.rstrip("\r") for line in text.splitlines()]
    # Skip blank lines and full-line comments (`# ...`) for friendliness.
    data_lines = [
        line for line in lines
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not data_lines:
        raise SampleSheetError(f"sample sheet {path}: file is empty.")

    header = [col.strip() for col in data_lines[0].split("\t")]
    errors: list[str] = []

    missing = [c for c in _REQUIRED_COLUMNS if c not in header]
    if missing:
        errors.append(
            f"missing required column(s): {', '.join(missing)}; "
            f"required columns are {list(_REQUIRED_COLUMNS)}"
        )
    if "kit_preset" in header:
        errors.append(
            "the 'kit_preset' column was removed in v0.7.1; use 'adapter' "
            "(3' sequence) or 'pretrimmed' (true/false) instead"
        )
    unknown = [c for c in header if c not in _KNOWN_COLUMNS]
    if unknown:
        errors.append(
            f"unknown column(s): {', '.join(unknown)}; "
            f"recognised columns are {sorted(_KNOWN_COLUMNS)}"
        )
    if errors:
        raise SampleSheetError(
            f"sample sheet {path}: " + "; ".join(errors) + "."
        )

    column_index = {col: idx for idx, col in enumerate(header)}
    rows: list[SampleRow] = []
    seen_ids: dict[str, int] = {}

    for line_no, raw in enumerate(data_lines[1:], start=2):
        cells = [cell.strip() for cell in raw.split("\t")]
        if len(cells) < len(header):
            cells.extend([""] * (len(header) - len(cells)))

        def get(col: str) -> str | None:
            idx = column_index.get(col)
            if idx is None or idx >= len(cells):
                return None
            v = cells[idx]
            return None if v in _NULL_TOKENS else v

        sample_id = get("sample_id")
        assay_raw = (get("assay") or "").lower() or None
        condition = get("condition")
        fastq_1 = get("fastq_1")

        for required, value in (
            ("sample_id", sample_id),
            ("assay", assay_raw),
            ("condition", condition),
            ("fastq_1", fastq_1),
        ):
            if value is None:
                errors.append(
                    f"line {line_no}: required field {required!r} is empty"
                )

        if sample_id is not None:
            if sample_id in seen_ids:
                errors.append(
                    f"line {line_no}: duplicate sample_id {sample_id!r} "
                    f"(first seen on line {seen_ids[sample_id]})"
                )
            else:
                seen_ids[sample_id] = line_no

        if assay_raw is not None and assay_raw not in _VALID_ASSAYS:
            errors.append(
                f"line {line_no}: assay {assay_raw!r} must be one of "
                f"{sorted(_VALID_ASSAYS)}"
            )

        umi_raw = get("umi_length")
        umi_length: int | None = None
        if umi_raw is not None:
            try:
                umi_length = int(umi_raw)
            except ValueError:
                errors.append(
                    f"line {line_no}: umi_length {umi_raw!r} is not an integer"
                )
            else:
                if umi_length < 0:
                    errors.append(
                        f"line {line_no}: umi_length must be >= 0 "
                        f"(got {umi_length})"
                    )

        umi_pos = get("umi_position")
        if umi_pos is not None and umi_pos not in _VALID_UMI_POS:
            errors.append(
                f"line {line_no}: umi_position {umi_pos!r} must be one of "
                f"{sorted(_VALID_UMI_POS)}"
            )

        def _parse_optional_umi_end(colname: str) -> int | None:
            raw = get(colname)
            if raw is None:
                return None
            try:
                value = int(raw)
            except ValueError:
                errors.append(
                    f"line {line_no}: {colname} {raw!r} is not an integer"
                )
                return None
            if value < 0:
                errors.append(
                    f"line {line_no}: {colname} must be >= 0 (got {value})"
                )
                return None
            return value

        umi_length_5p = _parse_optional_umi_end("umi_length_5p")
        umi_length_3p = _parse_optional_umi_end("umi_length_3p")
        if umi_pos == "both":
            if not (umi_length_5p and umi_length_3p):
                errors.append(
                    f"line {line_no}: umi_position='both' requires both "
                    "umi_length_5p > 0 and umi_length_3p > 0"
                )
            elif umi_length is not None and umi_length != umi_length_5p + umi_length_3p:
                errors.append(
                    f"line {line_no}: umi_position='both' requires "
                    f"umi_length ({umi_length}) to equal umi_length_5p + "
                    f"umi_length_3p ({umi_length_5p + umi_length_3p})"
                )

        strandedness = get("strandedness")
        if strandedness is not None and strandedness not in _VALID_STRANDEDNESS:
            errors.append(
                f"line {line_no}: strandedness {strandedness!r} must be one "
                f"of {sorted(_VALID_STRANDEDNESS)}"
            )

        dedup = get("dedup_strategy")
        if dedup is not None and dedup not in _VALID_DEDUP:
            errors.append(
                f"line {line_no}: dedup_strategy {dedup!r} must be one of "
                f"{sorted(_VALID_DEDUP)}"
            )
        # Collapse legacy ``umi-tools`` / ``umi_tools`` to the canonical
        # ``umi_coordinate`` so the SampleRow + downstream
        # canonical_config / run_manifest record one stable token.
        dedup = _normalise_dedup(dedup)

        exclude_raw = get("exclude")
        include_raw = get("include")
        exclude_bool = False
        # Honour `include` as the priority polarity when both are set;
        # this lets users adopt either convention without surprise.
        if include_raw is not None:
            iv = include_raw.lower()
            if iv in _BOOL_TRUE:
                exclude_bool = False
            elif iv in _BOOL_FALSE:
                exclude_bool = True
            else:
                errors.append(
                    f"line {line_no}: include {include_raw!r} is not a "
                    "boolean (use 'true' / 'false' / blank)"
                )
        elif exclude_raw is not None:
            ev = exclude_raw.lower()
            if ev in _BOOL_TRUE:
                exclude_bool = True
            elif ev in _BOOL_FALSE:
                exclude_bool = False
            else:
                errors.append(
                    f"line {line_no}: exclude {exclude_raw!r} is not a "
                    "boolean (use 'true' / 'false' / blank)"
                )

        library_type = get("library_type")
        if library_type is not None and library_type not in _VALID_LIBRARY_TYPE:
            errors.append(
                f"line {line_no}: library_type {library_type!r} must be "
                f"one of {sorted(_VALID_LIBRARY_TYPE)}"
            )

        def _parse_optional_int(
            colname: str, allow_zero: bool = True,
        ) -> int | None:
            raw = get(colname)
            if raw is None:
                return None
            try:
                value = int(raw)
            except ValueError:
                errors.append(
                    f"line {line_no}: {colname} {raw!r} is not an integer"
                )
                return None
            if value < 0 or (not allow_zero and value == 0):
                errors.append(
                    f"line {line_no}: {colname} must be > 0 (got {value})"
                )
                return None
            return value

        read_length_min = _parse_optional_int("read_length_min", allow_zero=False)
        read_length_max = _parse_optional_int("read_length_max", allow_zero=False)
        if (
            read_length_min is not None
            and read_length_max is not None
            and read_length_min > read_length_max
        ):
            errors.append(
                f"line {line_no}: read_length_min ({read_length_min}) > "
                f"read_length_max ({read_length_max})"
            )

        # Bail before constructing the row when required fields are
        # missing — stop the per-row validation but keep collecting
        # other-row errors.
        if (
            sample_id is None
            or assay_raw is None
            or assay_raw not in _VALID_ASSAYS
            or condition is None
            or fastq_1 is None
        ):
            continue

        pretrimmed_raw = get("pretrimmed")
        pretrimmed_value: bool | None = None
        if pretrimmed_raw is not None:
            lowered = pretrimmed_raw.strip().lower()
            if lowered in _BOOL_TRUE:
                pretrimmed_value = True
            elif lowered in _BOOL_FALSE:
                pretrimmed_value = False
            else:
                errors.append(
                    f"line {line_no}: pretrimmed {pretrimmed_raw!r} is not "
                    "a boolean (use 'true' / 'false' / blank)"
                )

        rows.append(
            SampleRow(
                sample_id=sample_id,
                assay=assay_raw,  # type: ignore[arg-type]
                condition=condition,
                fastq_1=Path(fastq_1),
                replicate=get("replicate"),
                biological_sample_id=get("biological_sample_id"),
                library_type=library_type,
                fastq_2=Path(get("fastq_2")) if get("fastq_2") else None,
                adapter=get("adapter"),
                pretrimmed=pretrimmed_value,
                umi_length=umi_length,
                umi_position=umi_pos,
                umi_length_5p=umi_length_5p,
                umi_length_3p=umi_length_3p,
                strandedness=strandedness,
                dedup_strategy=dedup,
                read_length_min=read_length_min,
                read_length_max=read_length_max,
                reference_id=get("reference_id"),
                exclude=exclude_bool,
                notes=get("notes"),
            )
        )

    if errors:
        joined = "\n  - ".join(errors)
        raise SampleSheetError(
            f"sample sheet {path}: {len(errors)} validation error(s):\n  - "
            + joined
        )

    if not rows:
        raise SampleSheetError(
            f"sample sheet {path}: no data rows found."
        )

    return SampleSheet(
        path=path,
        rows=tuple(rows),
        columns=tuple(header),
    )
