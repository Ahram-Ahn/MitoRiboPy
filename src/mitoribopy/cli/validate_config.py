"""``mitoribopy validate-config <config>`` — pre-flight check (P1.9).

Parses a YAML / JSON / TOML config the way ``mitoribopy all`` would,
canonicalises legacy keys (via :mod:`mitoribopy.config.migrate`), and
runs every cheap structural check before any stage actually fires:

* file exists / parses
* unknown top-level keys are flagged
* mutually-exclusive sections (``de_table`` vs ``rna_fastq``) caught
* ``samples:`` block conflicts with explicit per-stage inputs caught
* every path-shaped value (FASTQs, references, indexes, sample sheet)
  is checked for existence on disk
* ``rnaseq.mode`` is resolved and validated against supplied inputs

Exit codes:
* 0 — config is valid (or has only warnings)
* 2 — at least one structural error
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

from . import common


VALIDATE_CONFIG_HELP = (
    "Pre-flight a MitoRiboPy YAML / JSON / TOML config: parse, "
    "canonicalise legacy keys, check file paths and mutually-exclusive "
    "sections, and resolve rnaseq.mode against supplied inputs. Exit "
    "code is 0 on success, 2 when at least one error was found."
)


# Stage section keys whose VALUE is a path (or list of paths) we should
# probe for existence. Missing files never fail parsing — they only
# fail this validator.
_ALIGN_PATH_KEYS = (
    "fastq", "fastq_dir", "contam_index", "mt_index", "sample_overrides",
)
_RPF_PATH_KEYS = (
    "fasta", "directory", "annotation_file", "codon_tables_file",
    "read_counts_file",
)
_RNASEQ_PATH_KEYS = (
    "rna_fastq", "ribo_fastq", "reference_fasta", "bowtie2_index",
    "de_table", "ribo_dir", "ribo_counts", "reference_gtf",
    "condition_map", "sample_sheet", "upstream_rpf_counts",
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mitoribopy validate-config",
        description=VALIDATE_CONFIG_HELP,
        formatter_class=common.MitoRiboPyHelpFormatter,
    )
    parser.add_argument(
        "config",
        metavar="PATH",
        help="Path to the YAML / JSON / TOML config to validate.",
    )
    parser.add_argument(
        "--no-path-checks",
        action="store_true",
        default=False,
        help=(
            "Skip the on-disk existence checks for path-shaped values. "
            "Useful when validating a config on a different host than "
            "where the run will execute (e.g. CI-side validation)."
        ),
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        default=False,
        help=(
            "Treat any legacy-key rewrite as a warning AND make the "
            "validator exit 2 when at least one warning fired. Use for "
            "publication-grade configs that should already be canonical."
        ),
    )
    return parser


def _check_paths(
    section: dict,
    keys: tuple[str, ...],
    *,
    section_label: str,
    errors: list[str],
    warnings: list[str],
) -> None:
    """Probe every path-shaped value in *section* under *keys* for existence."""
    for key in keys:
        value = section.get(key)
        if value is None:
            continue
        candidates: list[str] = []
        if isinstance(value, str):
            candidates = [value]
        elif isinstance(value, (list, tuple)):
            candidates = [str(v) for v in value]
        else:
            continue
        for candidate in candidates:
            p = Path(candidate)
            # bowtie2 indexes are PREFIXES, not files; existence is
            # signalled by sidecar files at "<prefix>.1.bt2" or
            # "<prefix>.1.bt2l". Fall back to the parent dir check.
            if key in {"contam_index", "mt_index", "bowtie2_index"}:
                if not (
                    Path(candidate + ".1.bt2").exists()
                    or Path(candidate + ".1.bt2l").exists()
                ):
                    errors.append(
                        f"{section_label}.{key}: bowtie2 index prefix "
                        f"{candidate!r} has no .1.bt2 / .1.bt2l sidecar."
                    )
                continue
            if not p.exists():
                errors.append(
                    f"{section_label}.{key}: path {candidate!r} does not exist."
                )


def _check_unknown_top_level_keys(
    cfg: dict, errors: list[str]
) -> None:
    # Keep in sync with
    # :data:`mitoribopy.config.canonical._KNOWN_TOP_LEVEL_SECTIONS`.
    known = {"samples", "align", "rpf", "rnaseq", "execution", "all", "shared", "resume"}
    unknown = sorted(set(cfg) - known)
    if unknown:
        errors.append(
            "unknown top-level key(s): "
            + ", ".join(repr(k) for k in unknown)
            + f". Recognised: {sorted(known)}."
        )


def _check_mutually_exclusive(
    cfg: dict, errors: list[str]
) -> None:
    """Catch the cross-section conflicts the orchestrator would also catch."""
    rnaseq = cfg.get("rnaseq")
    if not isinstance(rnaseq, dict):
        return
    has_de_table = bool(rnaseq.get("de_table") or rnaseq.get("de-table"))
    has_fastq = bool(
        rnaseq.get("rna_fastq")
        or rnaseq.get("rna-fastq")
        or rnaseq.get("sample_sheet")
        or rnaseq.get("sample-sheet")
    )
    if has_de_table and has_fastq:
        errors.append(
            "rnaseq: 'de_table' and 'rna_fastq' / 'sample_sheet' are "
            "mutually exclusive. Pick one, or set rnaseq.mode "
            "explicitly to disambiguate."
        )


def _check_sheet_conflicts(cfg: dict, errors: list[str]) -> None:
    """Sample-sheet vs per-stage input conflicts."""
    if not cfg.get("samples"):
        return
    align = cfg.get("align") or {}
    rnaseq = cfg.get("rnaseq") or {}
    align_conflicts = [
        k for k in ("fastq", "fastq_dir", "samples", "sample_overrides")
        if align.get(k)
    ]
    if align_conflicts:
        errors.append(
            "samples: top-level sample sheet conflicts with align section "
            "key(s): " + ", ".join(align_conflicts)
            + ". The unified sheet supersedes per-stage inputs; drop one."
        )
    rnaseq_conflicts = [
        k for k in ("rna_fastq", "ribo_fastq", "condition_map")
        if rnaseq.get(k)
    ]
    if rnaseq_conflicts:
        errors.append(
            "samples: top-level sample sheet conflicts with rnaseq section "
            "key(s): " + ", ".join(rnaseq_conflicts) + "."
        )


def _check_rnaseq_mode(cfg: dict, errors: list[str]) -> None:
    """Run the orchestrator's mode resolver and surface its error."""
    from .all_ import _resolve_rnaseq_mode_from_config

    rnaseq = cfg.get("rnaseq")
    if not isinstance(rnaseq, dict):
        return
    # The orchestrator's `_apply_top_level_samples` (cli/all_.py) wires a
    # top-level `samples:` block into `rnaseq.sample_sheet` before mode
    # resolution runs at execute time. Mirror that wiring here so the
    # validator does not falsely reject a `from_fastq` config whose RNA
    # inputs come exclusively from the unified sample sheet.
    if cfg.get("samples") and not (
        rnaseq.get("sample_sheet")
        or rnaseq.get("sample-sheet")
        or rnaseq.get("rna_fastq")
        or rnaseq.get("rna-fastq")
    ):
        rnaseq = {**rnaseq, "sample_sheet": "<from-top-level-samples>"}
    _, err = _resolve_rnaseq_mode_from_config(rnaseq)
    if err is not None:
        errors.append(f"rnaseq: {err}")


def _check_sample_sheet_loads(cfg: dict, errors: list[str]) -> None:
    raw = cfg.get("samples")
    sheet_path: str | None = None
    if isinstance(raw, str):
        sheet_path = raw
    elif isinstance(raw, dict):
        sheet_path = raw.get("table") or raw.get("path")
    if not sheet_path:
        return
    if not Path(sheet_path).exists():
        # Path check above will have caught this; no double-report.
        return
    from ..sample_sheet import SampleSheetError, load_sample_sheet

    try:
        load_sample_sheet(sheet_path)
    except SampleSheetError as exc:
        errors.append(f"samples: {exc}")


def run(argv: Iterable[str]) -> int:
    args = build_parser().parse_args(list(argv))

    try:
        raw = common.load_config_file(args.config)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"[mitoribopy validate-config] ERROR: {exc}", file=sys.stderr)
        return 2

    # Always canonicalise FIRST so downstream checks see canonical keys.
    from ..config.migrate import migrate

    canonical, change_log = migrate(raw)

    errors: list[str] = []
    warnings: list[str] = []

    if change_log and args.strict:
        for line in change_log:
            warnings.append(f"legacy key: {line}")
    elif change_log:
        for line in change_log:
            warnings.append(f"legacy key (auto-canonicalised): {line}")

    _check_unknown_top_level_keys(canonical, errors)
    _check_sheet_conflicts(canonical, errors)
    _check_mutually_exclusive(canonical, errors)
    _check_rnaseq_mode(canonical, errors)
    _check_sample_sheet_loads(canonical, errors)

    if not args.no_path_checks:
        align = canonical.get("align") or {}
        rpf = canonical.get("rpf") or {}
        rnaseq = canonical.get("rnaseq") or {}
        # Top-level samples sheet path.
        samples = canonical.get("samples")
        if isinstance(samples, str):
            if not Path(samples).exists():
                errors.append(
                    f"samples: path {samples!r} does not exist."
                )
        elif isinstance(samples, dict):
            sheet = samples.get("table") or samples.get("path")
            if sheet and not Path(sheet).exists():
                errors.append(f"samples.table: path {sheet!r} does not exist.")
        if isinstance(align, dict):
            _check_paths(
                align, _ALIGN_PATH_KEYS,
                section_label="align", errors=errors, warnings=warnings,
            )
        if isinstance(rpf, dict):
            _check_paths(
                rpf, _RPF_PATH_KEYS,
                section_label="rpf", errors=errors, warnings=warnings,
            )
        if isinstance(rnaseq, dict):
            _check_paths(
                rnaseq, _RNASEQ_PATH_KEYS,
                section_label="rnaseq", errors=errors, warnings=warnings,
            )

    for w in warnings:
        sys.stderr.write(f"[mitoribopy validate-config] WARNING: {w}\n")
    for e in errors:
        sys.stderr.write(f"[mitoribopy validate-config] ERROR: {e}\n")

    if errors:
        return 2
    if warnings and args.strict:
        return 2
    sys.stderr.write(
        f"[mitoribopy validate-config] OK ({args.config}).\n"
    )
    return 0
