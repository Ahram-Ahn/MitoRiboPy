"""``mitoribopy all`` subcommand - align -> rpf -> (optional) rnaseq.

End-to-end orchestrator. Runs the three per-stage subcommands in
sequence with a single shared config file and writes a composed
``run_manifest.json`` at the run root that records every parameter,
tool version, and input/output hash so a reviewer can reproduce the
full pipeline from the manifest alone.

Invocation pattern::

    mitoribopy all --config pipeline_config.yaml --output results/

The YAML (or JSON / TOML) config has three optional sections::

    align:   { kit_preset: truseq_smallrna, fastq_dir: fastqs/, ... }
    rpf:     { strain: h.sapiens, rpf: [29, 34], ... }
    rnaseq:  # one of:
             { rna_fastq: rna/, reference_fasta: tx.fa, condition_map: ..., ... }
             { de_table: de.tsv, ribo_dir: rpf/, reference_gtf: ..., ... }

Each section's keys correspond to the CLI flag names of the matching
subcommand with dashes replaced by underscores. ``mitoribopy all``
reconstructs a flag list from the section and calls the subcommand's
``run()`` entry point.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path
from typing import Iterable

from .. import __version__
from ..cli.common import load_config_file
from . import common
from ._resume_guard import (
    force_resume_requested,
    load_prior_manifest,
    validate_resume,
)


# Bumped whenever the run_manifest.json layout changes in a way that
# breaks downstream consumers (added fields are minor; renamed or
# removed fields are major). Read it from your own scripts to gate on a
# compatible manifest shape.
MANIFEST_SCHEMA_VERSION = "1.3.0"  # 1.3: + resource_plan (top-level execution audit, v0.6.2)


ALL_SUBCOMMAND_HELP = (
    "End-to-end orchestrator: align + rpf, plus rnaseq when the config "
    "carries an 'rnaseq' section configured for either flow "
    "(from-FASTQ via 'rna_fastq' + 'reference_fasta', or external-DE "
    "via 'de_table'). Writes a composed run_manifest.json with tool "
    "versions, parameters, and input/output hashes across all three stages."
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mitoribopy all",
        description=ALL_SUBCOMMAND_HELP,
        epilog=(
            "This subcommand owns only the orchestrator flags. Stage-specific options live\n"
            "inside the config file sections whose keys match the subcommand flags\n"
            "(with dashes replaced by underscores, e.g. '--kit-preset' -> 'kit_preset').\n"
            "  align:  keys for 'mitoribopy align --help'\n"
            "  rpf:    keys for 'mitoribopy rpf --help'\n"
            "  rnaseq: keys for 'mitoribopy rnaseq --help'\n"
            "\n"
            "Start a new project:\n"
            "  mitoribopy all --print-config-template > pipeline_config.yaml\n"
            "  # edit the file to point at your FASTQs / indexes, then:\n"
            "  mitoribopy all --config pipeline_config.yaml --output results/\n"
            "\n"
            "Inspect a stage's full flag list:\n"
            "  mitoribopy all --show-stage-help align\n"
            "  mitoribopy all --show-stage-help rpf\n"
            "  mitoribopy all --show-stage-help rnaseq"
        ),
        formatter_class=common.MitoRiboPyHelpFormatter,
    )
    common.add_common_arguments(parser)
    parser.add_argument(
        "--output",
        required=False,
        metavar="DIR",
        help=(
            "Run root directory. Each stage writes under "
            "<output>/align/, <output>/rpf/, and <output>/rnaseq/ when "
            "the stage is active."
        ),
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        default=False,
        help=(
            "Skip stages whose expected output already exists: align "
            "is skipped when <output>/align/read_counts.tsv is present; "
            "rpf when <output>/rpf/rpf_counts.tsv is present; rnaseq "
            "when <output>/rnaseq/delta_te.tsv is present. The skip "
            "decision is gated by a hash check against the prior "
            "run_manifest.json (config_source_sha256, sample_sheet_sha256, "
            "reference_checksum, mitoribopy_version, schema_version); "
            "edits to any of those fields force the affected stage(s) "
            "to re-run unless --force-resume is also set."
        ),
    )
    parser.add_argument(
        "--force-resume",
        action="store_true",
        default=False,
        help=(
            "Like --resume, but bypass the hash guard. Use only when "
            "you know the stage outputs are still valid for the new "
            "config (e.g. you edited a comment-only line). Also "
            "honoured via the MITORIBOPY_FORCE_RESUME=1 environment "
            "variable for CI scripts that cannot easily change argv."
        ),
    )
    parser.add_argument(
        "--skip-align",
        action="store_true",
        default=False,
        help="Skip the align stage even when an [align] section exists.",
    )
    parser.add_argument(
        "--skip-rpf",
        action="store_true",
        default=False,
        help="Skip the rpf stage.",
    )
    parser.add_argument(
        "--skip-rnaseq",
        action="store_true",
        default=False,
        help="Skip the rnaseq stage even when an [rnaseq] section exists.",
    )
    parser.add_argument(
        "--manifest",
        default="run_manifest.json",
        metavar="PATH",
        help="Manifest filename (relative to --output).",
    )
    parser.add_argument(
        "--show-stage-help",
        choices=["align", "rpf", "rnaseq"],
        default=None,
        metavar="STAGE",
        help=(
            "Print the full help for one stage and exit. Useful because "
            "'mitoribopy all --help' only shows orchestrator-level flags."
        ),
    )
    parser.add_argument(
        "--print-config-template",
        action="store_true",
        default=False,
        help=(
            "Print a commented YAML config template covering every stage "
            "(align / rpf / rnaseq) with sensible defaults, then exit. "
            "Pipe this into a file to start a new project: "
            "'mitoribopy all --print-config-template > pipeline_config.yaml'."
        ),
    )
    parser.add_argument(
        "--print-canonical-config",
        action="store_true",
        default=False,
        help=(
            "Load --config, apply every auto-wiring + sample-sheet "
            "expansion that 'mitoribopy all' would normally apply, then "
            "print the resulting canonical config to stdout (YAML if "
            "PyYAML is available, JSON otherwise) and exit. Useful for "
            "diffing your input config against what was actually "
            "executed: 'mitoribopy all --print-canonical-config "
            "--config pipeline_config.yaml --output results/'. The same "
            "blob is embedded in run_manifest.json under "
            "'config_canonical' on real runs."
        ),
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        default=False,
        help=(
            "Publication-safe mode. A single switch that forwards "
            "strictness to every stage and post-run validation:\n"
            "  * align: --strict-publication-mode (fail on non-default "
            "policies that would invalidate a publication run);\n"
            "  * config: --strict on the up-front validate-config pass "
            "(treat any deprecated-key rewrite as a hard error);\n"
            "  * figures: --strict on the post-run validate-figures pass "
            "(promote warn-only QC findings to fail);\n"
            "  * summary: warning rows in warnings.tsv are mirrored as "
            "WARN bullets in SUMMARY.md.\n"
            "Recommended for any run that backs a manuscript figure."
        ),
    )
    parser.add_argument(
        "--progress",
        choices=["auto", "plain", "bar", "rich", "jsonl", "off"],
        default="auto",
        metavar="MODE",
        help=(
            "Progress display mode. 'auto' (default) picks a tqdm bar on "
            "an interactive TTY and one-line plain logs on HPC / non-TTY "
            "streams. 'plain' forces stable [PROGRESS] log lines. 'bar' "
            "forces tqdm bars (falls back to plain if tqdm is missing). "
            "'rich' degrades to bar in this version. 'jsonl' streams "
            "machine-readable JSON events to stderr (or to "
            "--progress-file). 'off' silences all progress output, but "
            "--progress-file still attaches if set."
        ),
    )
    parser.add_argument(
        "--progress-file",
        default=None,
        metavar="PATH",
        help=(
            "Append every typed progress event to this file as JSONL "
            "(one JSON object per line). Independent of --progress: an "
            "explicit --progress-file always attaches, even when "
            "--progress=off. Default location is "
            "<output>/progress.jsonl when --output is set; pass "
            "explicitly to override."
        ),
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        default=False,
        help="Shortcut for --progress=off.",
    )
    return parser


# ---------------------------------------------------------------------------
# argv reconstruction from a config section
# ---------------------------------------------------------------------------


def _dict_to_argv(
    section: dict,
    *,
    flag_style: str = "hyphen",
    repeat_flags: set[str] | None = None,
    flag_overrides: dict[str, str] | None = None,
) -> list[str]:
    """Serialize a section dict into a CLI-style argv list.

    Rules:

    * ``flag_style="hyphen"`` (default, used for every stage as of
      v0.6.0): keys with underscores become ``--dashed-form`` flags.
      ``align``, ``rpf``, and ``rnaseq`` all accept the canonical
      hyphenated form; ``rpf`` also keeps the underscore-style aliases
      around for one transition cycle but they are no longer emitted
      from this serializer.
    * ``flag_style="underscore"`` (legacy): keys are emitted verbatim
      with a ``--`` prefix. Retained only for tests that pin the
      underscore-era serialisation; new call sites should always use
      ``"hyphen"``.
    * Boolean ``True`` emits just the flag; ``False`` emits nothing.
    * Lists / tuples are emitted as ``--flag v1 v2 v3`` (nargs="+" style).
      Keys listed in ``repeat_flags`` are emitted as repeated ``--flag v``
      pairs for argparse options that use ``action="append"``.
    * ``None`` values are skipped.
    * ``flag_overrides`` maps config keys to exact legacy flag spellings.
    """
    if flag_style not in {"hyphen", "underscore"}:
        raise ValueError(
            f"flag_style must be 'hyphen' or 'underscore', got {flag_style!r}."
        )

    repeat_flags = repeat_flags or set()
    flag_overrides = flag_overrides or {}
    argv: list[str] = []
    for key, value in section.items():
        if key in flag_overrides:
            flag = flag_overrides[key]
        elif flag_style == "underscore":
            flag = f"--{key.replace('-', '_')}"
        else:
            flag = f"--{key.replace('_', '-')}"
        if value is None:
            continue
        if isinstance(value, bool):
            if value:
                argv.append(flag)
            continue
        if isinstance(value, (list, tuple)):
            if key in repeat_flags:
                for item in value:
                    argv.extend([flag, str(item)])
            else:
                argv.append(flag)
                argv.extend(str(v) for v in value)
            continue
        argv.extend([flag, str(value)])
    return argv


def _normalize_align_inputs(
    align_cfg: dict, *, run_root: Path | None = None
) -> dict:
    """Make ``align.fastq`` polymorphic + materialise per-sample overrides.

    YAML users routinely want to point at a directory instead of listing
    every FASTQ. ``--fastq-dir`` already exists at the CLI; this helper
    lets the YAML key ``fastq:`` accept either a list (treated as
    explicit paths) or a single string (treated as ``fastq_dir``).

    A YAML ``samples:`` block (per-sample UMI / kit overrides) is
    materialised as a TSV under ``<run_root>/align/sample_overrides.tsv``
    and replaced in the returned dict with ``sample_overrides:`` so
    :func:`_dict_to_argv` serialises it as ``--sample-overrides PATH``.
    Pass ``run_root=None`` (e.g. in dry-run) to keep the ``samples:``
    block in the dict for plan rendering only.
    """
    if not align_cfg:
        return align_cfg
    cfg = dict(align_cfg)
    fastq_value = cfg.get("fastq")
    if isinstance(fastq_value, str):
        # String → directory shortcut. Promote to fastq_dir without
        # clobbering an existing fastq_dir if the user set both.
        cfg.setdefault("fastq_dir", fastq_value)
        cfg["fastq"] = None

    samples_block = cfg.pop("samples", None)
    if samples_block:
        if run_root is None:
            # Dry-run / non-materialising path: drop the block from the
            # argv (the align CLI does not accept it) and rely on the
            # caller to render it separately when planning.
            cfg["sample_overrides"] = None
        else:
            tsv_path = Path(run_root) / "align" / "sample_overrides.tsv"
            _write_samples_block_as_tsv(samples_block, tsv_path)
            cfg["sample_overrides"] = str(tsv_path)
    return cfg


def _write_samples_block_as_tsv(
    samples_block: list, tsv_path: Path
) -> Path:
    """Materialise an ``align.samples:`` YAML list as a sample_overrides TSV.

    Accepted item shape::

        - name: sampleA            # required
          kit_preset: …            # all five fields are optional
          adapter: …
          umi_length: 8
          umi_position: 5p
          dedup_strategy: skip

    The TSV columns match :data:`mitoribopy.align.sample_resolve._SAMPLE_OVERRIDE_COLUMNS`
    so :func:`mitoribopy.align.sample_resolve.read_sample_overrides_tsv`
    round-trips it.
    """
    tsv_path = Path(tsv_path)
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    columns = (
        "sample",
        "kit_preset",
        "adapter",
        "umi_length",
        "umi_position",
        "umi_length_5p",
        "umi_length_3p",
        "dedup_strategy",
    )
    rows: list[dict[str, str]] = []
    for index, entry in enumerate(samples_block):
        if not isinstance(entry, dict):
            raise ValueError(
                f"align.samples[{index}] must be a mapping with at least "
                f"a 'name' field, got {type(entry).__name__}."
            )
        sample = entry.get("name") or entry.get("sample")
        if not sample:
            raise ValueError(
                f"align.samples[{index}] is missing the required "
                "'name' field (the FASTQ basename without extension)."
            )
        rows.append(
            {
                "sample": str(sample),
                "kit_preset": _stringify(entry.get("kit_preset")),
                "adapter": _stringify(entry.get("adapter")),
                "umi_length": _stringify(entry.get("umi_length")),
                "umi_position": _stringify(entry.get("umi_position")),
                "umi_length_5p": _stringify(entry.get("umi_length_5p")),
                "umi_length_3p": _stringify(entry.get("umi_length_3p")),
                "dedup_strategy": _stringify(entry.get("dedup_strategy")),
            }
        )
    with tsv_path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(row[col] for col in columns) + "\n")
    return tsv_path


def _stringify(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def _apply_top_level_samples(config: dict, *, run_root: Path) -> int:
    """Materialise a top-level ``samples:`` block into stage configs.

    Accepted shapes::

        samples:
          table: samples.tsv     # path to a unified sample sheet TSV

        # OR shorthand:
        samples: samples.tsv

    When set, the sheet drives both stages:

    * ``align.fastq``           ← list of Ribo-seq FASTQ paths from the
                                  sheet (rows where ``assay='ribo'`` and
                                  ``exclude`` is false).
    * ``align.sample_overrides``← TSV materialised from the Ribo rows'
                                  per-sample kit / UMI columns under
                                  ``<run_root>/align/sample_overrides.tsv``.
    * ``rnaseq.sample_sheet``   ← path to the original sheet so the
                                  rnaseq stage CLI loads it directly.

    Reject (exit 2) when the sheet is set together with conflicting
    explicit per-stage inputs: the user has to pick one input style per
    stage. Returns 0 on success.
    """
    raw = config.get("samples")
    if raw in (None, "", {}):
        return 0

    if isinstance(raw, str):
        sheet_path = raw
    elif isinstance(raw, dict):
        sheet_path = raw.get("table") or raw.get("path")
        if not sheet_path:
            print(
                "[mitoribopy all] ERROR: top-level 'samples:' block must "
                "carry a 'table:' (or 'path:') key pointing at the "
                "sample-sheet TSV.",
                file=sys.stderr,
            )
            return 2
    else:
        print(
            "[mitoribopy all] ERROR: top-level 'samples:' must be either "
            "a path string or a mapping with a 'table:' key; "
            f"got {type(raw).__name__}.",
            file=sys.stderr,
        )
        return 2

    from ..sample_sheet import (
        SampleSheetError,
        check_sheet_conflicts,
        format_sheet_conflict_error,
        load_sample_sheet,
    )

    try:
        sheet = load_sample_sheet(sheet_path)
    except SampleSheetError as exc:
        print(f"[mitoribopy all] ERROR: {exc}", file=sys.stderr)
        return 2

    align_cfg = config.setdefault("align", {})
    # Only thread the sheet into rnaseq when the user explicitly opted
    # into the rnaseq stage. A bare `samples:` block must not silently
    # turn on rnaseq for align+rpf-only configs.
    rnaseq_explicit = config.get("rnaseq") is not None
    rnaseq_cfg = config.setdefault("rnaseq", {}) if rnaseq_explicit else {}

    # --- align side -----------------------------------------------------
    # Ribo rows feed `mitoribopy align`. Reject explicit conflicts so a
    # user does not get a silent shadow when they declared FASTQs in two
    # places.
    ribo_rows = sheet.by_assay("ribo")
    if ribo_rows:
        align_conflicts = check_sheet_conflicts(
            align_cfg,
            conflict_keys=(
                "fastq", "fastq_dir", "samples", "sample_overrides",
            ),
        )
        if align_conflicts:
            print(
                format_sheet_conflict_error(
                    "mitoribopy all/align",
                    [f"align.{k}" for k in align_conflicts],
                ),
                file=sys.stderr,
            )
            return 2
        align_cfg["fastq"] = [str(r.fastq_1) for r in ribo_rows]
        # Materialise per-sample overrides only when at least one Ribo
        # row carries per-sample kit / UMI fields — otherwise leave the
        # global align defaults in charge.
        if any(
            r.kit_preset or r.adapter or r.umi_length is not None
            or r.umi_position or r.dedup_strategy
            or getattr(r, "umi_length_5p", None) is not None
            or getattr(r, "umi_length_3p", None) is not None
            for r in ribo_rows
        ):
            tsv_path = Path(run_root) / "align" / "sample_overrides.tsv"
            tsv_path.parent.mkdir(parents=True, exist_ok=True)
            cols = (
                "sample", "kit_preset", "adapter",
                "umi_length", "umi_position",
                "umi_length_5p", "umi_length_3p",
                "dedup_strategy",
            )
            with tsv_path.open("w", encoding="utf-8") as handle:
                handle.write("\t".join(cols) + "\n")
                for r in ribo_rows:
                    row = {
                        "sample": r.sample_id,
                        "kit_preset": _stringify(r.kit_preset),
                        "adapter": _stringify(r.adapter),
                        "umi_length": _stringify(r.umi_length),
                        "umi_position": _stringify(r.umi_position),
                        "umi_length_5p": _stringify(
                            getattr(r, "umi_length_5p", None)
                        ),
                        "umi_length_3p": _stringify(
                            getattr(r, "umi_length_3p", None)
                        ),
                        "dedup_strategy": _stringify(r.dedup_strategy),
                    }
                    handle.write("\t".join(row[c] for c in cols) + "\n")
            align_cfg["sample_overrides"] = str(tsv_path)

    # --- rnaseq side ----------------------------------------------------
    # The rnaseq subcommand has its own --sample-sheet handler that
    # rejects mixing the sheet with --rna-fastq / --condition-map / etc.
    # We just hand the path through.
    if rnaseq_explicit:
        rnaseq_conflicts = check_sheet_conflicts(
            rnaseq_cfg,
            conflict_keys=(
                "rna_fastq", "rna-fastq", "ribo_fastq", "ribo-fastq",
                "condition_map", "condition-map",
            ),
        )
        if rnaseq_conflicts and rnaseq_cfg.get("sample_sheet") is None:
            print(
                format_sheet_conflict_error(
                    "mitoribopy all/rnaseq",
                    [f"rnaseq.{k}" for k in rnaseq_conflicts],
                ),
                file=sys.stderr,
            )
            return 2
        rnaseq_cfg.setdefault("sample_sheet", str(sheet.path))

    return 0


# ---------------------------------------------------------------------------
# Top-level `execution:` block (v0.6.2)
# ---------------------------------------------------------------------------


# Recognised keys inside the optional top-level `execution:` section.
# Unknown keys are reported by validate-config; nothing here is silently
# ignored. Matches the cascade rules documented in the README.
_EXECUTION_BLOCK_KEYS: tuple[str, ...] = (
    "threads",
    "memory_gb",
    "parallel_samples",
    "single_sample_mode",
    "min_threads_per_sample",
    "estimated_memory_per_sample_gb",
    "scheduler",
)


def _apply_execution_block(
    config: dict,
    *,
    cli_threads: int | None,
) -> dict:
    """Cascade the top-level ``execution:`` block into stage configs.

    The orchestrator-level CLI flag ``--threads N`` populates
    ``execution.threads`` when the YAML did not already pin it. Stage-
    specific keys (``align.threads``, ``align.max_parallel_samples``,
    ``align.single_sample_mode``, ``align.memory_gb``,
    ``rnaseq.align_threads``) win when set; otherwise the cascade fills
    them from ``execution.*``.

    Returns the (resolved) ``execution`` dict so the caller can persist
    it to ``resource_plan.json`` / the manifest. Mutates ``config`` in
    place to keep call-sites simple.
    """
    raw = config.get("execution") or {}
    if not isinstance(raw, dict):
        # Caller surfaces typed errors; here we just normalise.
        raw = {}
    execution: dict = {k: raw.get(k) for k in _EXECUTION_BLOCK_KEYS}

    # CLI --threads at the orchestrator level overrides nothing the user
    # explicitly pinned in YAML. Otherwise it becomes the global budget.
    if cli_threads is not None and execution.get("threads") in (None, "auto"):
        execution["threads"] = int(cli_threads)

    align_cfg = config.get("align")
    rnaseq_cfg = config.get("rnaseq")

    def _cascade(stage_cfg, stage_key, exec_key):
        if not isinstance(stage_cfg, dict):
            return
        if stage_cfg.get(stage_key) in (None, "auto"):
            value = execution.get(exec_key)
            if value is not None:
                stage_cfg[stage_key] = value

    if isinstance(align_cfg, dict):
        _cascade(align_cfg, "threads", "threads")
        _cascade(align_cfg, "max_parallel_samples", "parallel_samples")
        _cascade(align_cfg, "single_sample_mode", "single_sample_mode")
        _cascade(align_cfg, "memory_gb", "memory_gb")
    if isinstance(rnaseq_cfg, dict):
        # rnaseq's per-stage 'threads' is named align_threads and
        # only applies to its bowtie2 subprocess. Cascade if unset.
        _cascade(rnaseq_cfg, "align_threads", "threads")

    config["execution"] = {
        k: v for k, v in execution.items() if v is not None
    }
    return config["execution"]


# ---------------------------------------------------------------------------
# Periodicity config block
# ---------------------------------------------------------------------------

# Map of `periodicity.<key>` -> `rpf.periodicity_<key>` (the rpf parser
# flag). Spec defines a periodicity section; this mapping translates it
# into the rpf-stage CLI surface so the rpf parser stays the single
# source of truth for the actual flags.
_PERIODICITY_BLOCK_KEYS: dict[str, str] = {
    "enabled": "periodicity_enabled",
    "exclude_start_codons": "periodicity_exclude_start_codons",
    "exclude_stop_codons": "periodicity_exclude_stop_codons",
    "phase_score": "periodicity_phase_score",
    "fourier_spectrum": "periodicity_fourier_spectrum",
    "fourier_window_nt": "periodicity_fourier_window_nt",
    "metagene_nt": "periodicity_metagene_nt",
    "min_reads_per_length": "periodicity_min_reads_per_length",
}


def _apply_periodicity_block(config: dict) -> dict:
    """Cascade the top-level ``periodicity:`` block into ``rpf.periodicity_*``.

    Spec defines a top-level ``periodicity:`` section. The rpf parser
    owns the actual flags (``--periodicity-*``); this helper translates
    the section into rpf-stage keys so users can write the spec-shaped
    YAML and it Just Works through ``mitoribopy all``.

    Explicit ``rpf.periodicity_*`` keys win over the top-level section.
    Returns the resolved periodicity dict so callers can record it.
    Mutates ``config`` in place.
    """
    raw = config.get("periodicity") or {}
    if not isinstance(raw, dict):
        raw = {}

    # Spec also nests metric toggles + thresholds; flatten the ones we
    # surface as rpf flags. Unknown nested keys are silently ignored
    # so the YAML can carry future-spec keys without breaking older runs.
    metrics = raw.get("metrics") or {}
    if isinstance(metrics, dict):
        if "phase_score" in metrics and "phase_score" not in raw:
            raw["phase_score"] = metrics["phase_score"]
        if (
            "fourier_spectrum" in metrics
            and "fourier_spectrum" not in raw
        ):
            raw["fourier_spectrum"] = metrics["fourier_spectrum"]
    thresholds = raw.get("thresholds") or {}
    if isinstance(thresholds, dict):
        if "min_reads_per_length" in thresholds and "min_reads_per_length" not in raw:
            raw["min_reads_per_length"] = thresholds["min_reads_per_length"]

    rpf_cfg = config.get("rpf")
    if not isinstance(rpf_cfg, dict):
        rpf_cfg = {}
        config["rpf"] = rpf_cfg

    resolved: dict = {}
    for spec_key, rpf_key in _PERIODICITY_BLOCK_KEYS.items():
        if spec_key not in raw:
            continue
        # Explicit rpf-side override wins.
        if rpf_key in rpf_cfg:
            resolved[spec_key] = rpf_cfg[rpf_key]
            continue
        rpf_cfg[rpf_key] = raw[spec_key]
        resolved[spec_key] = raw[spec_key]

    if resolved:
        config["periodicity"] = resolved
    return resolved


def _resolve_and_write_run_root_resource_plan(
    config: dict,
    *,
    run_root: Path,
    n_samples: int,
) -> Path | None:
    """Plan resources at the orchestrator level and persist the result.

    Writes ``<run_root>/resource_plan.json`` derived from the resolved
    ``execution`` block (after cascade). Returns the path written, or
    ``None`` when the plan could not be computed (e.g. zero samples
    detected and nothing useful to record).
    """
    from ..pipeline.resource_plan import plan_parallelism, write_resource_plan

    execution = config.get("execution") or {}
    if n_samples < 1:
        # Still record the orchestrator's view so reviewers have an
        # artifact even on align-skip / dry-runs.
        n_samples = 1
    plan = plan_parallelism(
        n_samples=n_samples,
        requested_threads=execution.get("threads", "auto"),
        requested_parallel=execution.get("parallel_samples", "auto"),
        memory_gb=execution.get("memory_gb", "auto"),
        min_threads_per_sample=int(
            execution.get("min_threads_per_sample") or 2
        ),
        estimated_memory_per_sample_gb=float(
            execution.get("estimated_memory_per_sample_gb") or 4.0
        ),
    )
    return write_resource_plan(plan, run_root)


def _count_samples_from_config(config: dict) -> int:
    """Best-effort sample count for the run-root resource_plan.

    Looks at the unified `samples:` sheet first, falls back to
    `align.fastq:` (list or directory), and otherwise returns 1 so the
    plan is still written.
    """
    samples_block = config.get("samples")
    sheet_path: str | None = None
    if isinstance(samples_block, str):
        sheet_path = samples_block
    elif isinstance(samples_block, dict):
        sheet_path = samples_block.get("table") or samples_block.get("path")
    if sheet_path and Path(sheet_path).exists():
        try:
            from ..sample_sheet import load_sample_sheet

            sheet = load_sample_sheet(sheet_path)
            ribo_count = len(sheet.by_assay("ribo"))
            if ribo_count:
                return ribo_count
            return max(1, len(sheet.rows))
        except Exception:  # noqa: BLE001 — best-effort
            pass
    align_cfg = config.get("align") or {}
    fastq = align_cfg.get("fastq")
    if isinstance(fastq, list):
        return max(1, len(fastq))
    if isinstance(fastq, str):
        directory = Path(fastq)
        if directory.is_dir():
            patterns = ("*.fq.gz", "*.fastq.gz", "*.fq", "*.fastq")
            n = 0
            for pattern in patterns:
                n += sum(1 for _ in directory.glob(pattern))
            if n:
                return n
    return 1


# Per-stage config-template comment, also used by --print-config-template.
_CONFIG_TEMPLATE = """\
# mitoribopy all --config pipeline_config.yaml --output results/
#
# Every key below corresponds one-to-one to the matching subcommand's
# CLI flag (dashes replaced by underscores). Unset keys fall back to
# that subcommand's documented default. Omit a whole section to skip
# that stage; 'mitoribopy all' will record it as skipped in the manifest.
#
# Use 'mitoribopy all --show-stage-help {align,rpf,rnaseq}' for the
# full flag list with defaults.

# ---- samples (RECOMMENDED — single source of truth) -----------------------
# A unified per-project TSV declaring every sample's identity, assay,
# condition, FASTQ paths, and (optionally) per-sample kit / UMI /
# strandedness overrides. When this is set, you should NOT also set
# 'align.fastq' / 'align.fastq_dir' / 'align.samples' /
# 'rnaseq.rna_fastq' / 'rnaseq.condition_map' — the sheet supersedes
# those per-stage inputs and a conflict is a hard error.
#
# Required columns: sample_id, assay (ribo|rna), condition, fastq_1
# Optional columns: replicate, fastq_2, kit_preset, adapter, umi_length,
#                   umi_position, strandedness, dedup_strategy, exclude, notes
#
# samples:
#   table: samples.tsv

# ---- align -----------------------------------------------------------------
align:
  # Library chemistry. `auto` (default) detects the adapter family per
  # sample by scanning the head of each FASTQ; mixed-kit and mixed-UMI
  # batches are first-class. Set an explicit preset only when you need
  # a per-sample fallback for samples whose adapter cannot be
  # auto-identified.
  kit_preset: auto                # auto | illumina_smallrna | illumina_truseq |
                                  # illumina_truseq_umi | qiaseq_mirna |
                                  # pretrimmed | custom
                                  # (legacy aliases truseq_smallrna,
                                  #  nebnext_smallrna, nebnext_ultra_umi
                                  #  still accepted for back-compat)
  adapter: null                   # explicit fallback adapter; required when
                                  # kit_preset=custom and detection fails
  adapter_detection: auto         # auto | off | strict
  umi_length: null                # overrides kit preset's UMI length
  umi_position: null              # 5p | 3p

  # Adapter-detection tuning. Defaults are calibrated for typical
  # 36-150 nt sequencing libraries; touch only when detection is
  # producing surprises.
  adapter_detect_reads: 5000              # head-of-FASTQ scan size
  adapter_detect_min_rate: 0.30           # min fraction with adapter signal
  adapter_detect_min_len: 12              # adapter prefix length used
  adapter_detect_pretrimmed_threshold: 0.05  # all-kits below this -> pretrimmed
  # allow_pretrimmed_inference: true      # set false to hard-fail when
                                          #   adapter detection finds no kit

  # Inputs / reference indexes (bowtie2 prefixes built by bowtie2-build).
  # NOTE: when the top-level `samples:` block is set, `fastq:` and
  # `fastq_dir:` MUST be unset — the sheet declares per-sample FASTQs.
  # The keys below are for legacy / standalone runs without a sheet.
  # `fastq` accepts either a directory string OR an explicit list of paths.
  # Picked up patterns: *.fq, *.fq.gz, *.fastq, *.fastq.gz.
  fastq: input_data/              # legacy: only when not using top-level samples:
  # fastq:                        # OR list of explicit paths:
  #   - /path/to/sample_A.fq.gz
  #   - /path/to/sample_B.fq.gz
  # fastq_dir: null               # legacy alias for the directory form
  contam_index: null              # REQUIRED: rRNA/tRNA bowtie2 index prefix
  mt_index: null                  # REQUIRED: mt-transcriptome bowtie2 index prefix

  # Strandedness / length window / MAPQ / dedup
  library_strandedness: forward   # forward | reverse | unstranded
  min_length: 15
  max_length: 45
  quality: 20
  mapq: 10
  seed: 42
  dedup_strategy: auto            # auto | umi-tools | skip

  # Concurrency. >1 runs samples in parallel; the global --threads value
  # is divided across workers (each tool gets max(1, threads // N)
  # threads), so total CPU use stays around `threads`. The joint
  # `mitoribopy rpf` stage is unaffected -- offsets are selected across
  # all samples and remain serial there.
  # max_parallel_samples: 1

  # Per-sample overrides (mixed-UMI batches). Use this ONLY when not
  # using the top-level `samples:` block; the unified sheet's
  # umi_length / umi_position / kit_preset / dedup_strategy columns
  # are the recommended way to declare per-sample overrides.
  # The 'name' field must match the FASTQ basename with .fq[.gz] /
  # .fastq[.gz] stripped. Any unset override field falls through to
  # the globals above.
  # samples:
  #   - name: sampleA
  #     kit_preset: illumina_truseq_umi
  #     umi_length: 8
  #     umi_position: 5p
  #   - name: sampleB
  #     kit_preset: qiaseq_mirna
  #     umi_length: 12
  #     umi_position: 3p
  #   - name: sampleC               # SRA-deposited, already adapter-clipped
  #     kit_preset: pretrimmed
  #     umi_length: 0

# ---- rpf -------------------------------------------------------------------
rpf:
  # Organism + RPF length window.
  #   h.sapiens     human mt (default; ships annotation + codon table)
  #   s.cerevisiae  budding yeast mt (ships annotation + codon table)
  #   custom        any other organism; supply --annotation_file +
  #                 --codon_table_name (or --codon_tables_file).
  # `h` and `y` are still accepted as deprecated short aliases.
  strain: h.sapiens
  # `rpf: null` lets --footprint_class drive the default range.
  # short:    h.sapiens / s.cerevisiae 16-24 nt   (RNase truncation products)
  # monosome: h.sapiens 28-34 nt, s.cerevisiae 37-41 nt
  # disome:   h.sapiens 50-70 nt, s.cerevisiae 60-90 nt
  footprint_class: monosome
  rpf: [29, 34]
  fasta: null                     # reference FASTA (mt-transcriptome)
  directory: null                 # BED input dir (auto-wired from align/bed/)

  # Offset anchor + selection bounds
  align: stop                     # start | stop
  offset_type: "5"                # "5" | "3"  (report offsets from the 5' or 3' end)
  offset_site: p                  # p | a   coordinate space for the SELECTED
                                  # OFFSETS table only. Use analysis_sites to
                                  # control which downstream outputs are
                                  # produced (codon usage, coverage plots).
  offset_pick_reference: p_site   # p_site (default) | reported_site
                                  # p_site:        pick offset in canonical
                                  #                P-site space, then convert
                                  #                to --offset_site for output.
                                  # reported_site: pick offset directly in
                                  #                --offset_site space.
                                  # (Legacy alias 'selected_site' = 'reported_site'.)
  offset_mode: per_sample         # per_sample (default) | combined.
                                  # per_sample: each sample uses its own
                                  # offsets; combined: pool all samples and
                                  # apply one offset table (v0.3.x behaviour).
  analysis_sites: both            # both (default) | p | a
                                  # both writes parallel P-site and A-site
                                  # coverage plots and codon usage tables.
  min_5_offset: 10
  max_5_offset: 22
  min_3_offset: 10
  max_3_offset: 22
  offset_mask_nt: 5

  # Output / plotting
  plot_format: svg                # png | pdf | svg
  codon_density_window: true      # smooth codon-density with +/-1 nt window
                                  # (legacy key: 'merge_density')
  structure_density: false

  # Optional downstream modules
  cor_plot: false
  base_sample: null

# ---- rnaseq (optional) -----------------------------------------------------
# Uncomment ONE of the two mutually exclusive flows below to enable the
# translation-efficiency stage. Both produce te.tsv, delta_te.tsv, and plots.
#
# `rnaseq_mode` (RECOMMENDED): set explicitly to one of:
#   de_table   PUBLICATION: external DE table + prior rpf run.
#   from_fastq EXPLORATORY: in-tree pyDESeq2 on the mt-mRNA subset only.
#   none       stage section present but inert (no inputs).
# When omitted, the mode is inferred from supplied inputs.
#
# Publication-grade flow (mode=de_table): bring your own DE table.
# rnaseq:
#   rnaseq_mode: de_table
#   de_table: /path/to/de_table.tsv
#   gene_id_convention: hgnc
#   reference_gtf: /path/to/reference.fa       # must match align's --mt-index source
#   condition_map: /path/to/conditions.tsv
#   condition_a: control
#   condition_b: knockdown
#
# Exploratory flow (mode=from_fastq): from raw FASTQ (rna_fastq + reference_fasta).
# Recommended only for tutorials / smoke tests; for publication-grade DE,
# run DESeq2 / Xtail / Anota2Seq externally on the full transcriptome and
# come back via mode=de_table above.
# rnaseq:
#   rnaseq_mode: from_fastq
#   rna_fastq: /path/to/rnaseq/                # dir or list of FASTQs
#   ribo_fastq: /path/to/riboseq/              # optional; reuses align/ when omitted
#   reference_fasta: /path/to/transcriptome.fa # auto-wired from rpf.fasta if unset
#   condition_map: /path/to/conditions.tsv
#   condition_a: control
#   condition_b: knockdown
#   gene_id_convention: hgnc                   # hgnc | ensembl | refseq | bare
"""


def _print_config_template() -> None:
    """Write the commented config template to stdout."""
    sys.stdout.write(_CONFIG_TEMPLATE)


def _sha256_of(path: Path) -> str | None:
    try:
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for chunk in iter(lambda: handle.read(65536), b""):
                digest.update(chunk)
        return digest.hexdigest()
    except OSError:
        return None


def _read_stage_settings(stage_dir: Path) -> dict | None:
    path = stage_dir / "run_settings.json"
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def _git_commit() -> str | None:
    """Best-effort current git commit — never raises.

    Returns ``None`` when not in a repo, when git is missing, or when
    the call fails for any reason. Intentionally cheap; we never want
    a manifest write to fail because git is misconfigured.
    """
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=2,
            check=False,
        )
    except (FileNotFoundError, subprocess.SubprocessError, OSError):
        return None
    if out.returncode != 0:
        return None
    sha = out.stdout.strip()
    return sha or None


def _yaml_dump(payload: dict) -> str:
    """YAML dump if PyYAML is available, JSON fallback otherwise.

    The fallback is intentional: `--print-canonical-config` needs to
    work in lean install environments where PyYAML is not installed.
    JSON is a strict subset of YAML so the output remains valid YAML.
    """
    try:
        import yaml  # type: ignore[import-not-found]
    except ImportError:
        return json.dumps(payload, indent=2, sort_keys=True) + "\n"
    return yaml.safe_dump(
        payload, sort_keys=True, default_flow_style=False, allow_unicode=True
    )


def _build_stages_block(
    stages_run: list[str],
    stages_skipped: list[str],
    runtimes: dict[str, float],
    skip_reasons: dict[str, str],
) -> dict[str, dict]:
    """Reshape parallel run/skipped lists into a {stage: {status, ...}} map.

    Status values:
    * ``completed`` — the stage executed successfully.
    * ``skipped``   — the stage was not configured or was skipped via
                      --skip-* / --resume; ``reason`` carries the trigger.
    """
    out: dict[str, dict] = {}
    for stage in ("align", "rpf", "rnaseq"):
        if stage in stages_run:
            out[stage] = {"status": "completed"}
            if stage in runtimes:
                out[stage]["runtime_seconds"] = round(runtimes[stage], 3)
        elif stage in stages_skipped:
            entry: dict = {"status": "skipped"}
            reason = skip_reasons.get(stage)
            if reason:
                entry["reason"] = reason
            out[stage] = entry
        else:
            out[stage] = {"status": "not_configured"}
    return out


def _lift_tool_versions(
    align: dict | None, rpf: dict | None, rnaseq: dict | None
) -> dict[str, str]:
    """Pull tool versions out of per-stage run_settings.json files.

    Stages already record (varying subsets of) their tool versions in
    their own run_settings; this just lifts whatever's there into a
    flat top-level map. Missing keys stay missing — we don't probe.
    """
    out: dict[str, str] = {}
    out["python"] = platform.python_version()
    out["mitoribopy"] = __version__
    for settings in (align, rpf, rnaseq):
        if not isinstance(settings, dict):
            continue
        tools = settings.get("tools")
        if isinstance(tools, dict):
            for k, v in tools.items():
                if v is not None and k not in out:
                    out[k] = str(v)
        # Some stages stash individual tool versions at the top level.
        for legacy_key in (
            "cutadapt_version",
            "bowtie2_version",
            "umi_tools_version",
            "pysam_version",
            "pydeseq2_version",
        ):
            v = settings.get(legacy_key)
            if v is not None:
                tool_name = legacy_key.removesuffix("_version")
                out.setdefault(tool_name, str(v))
    return out


def _write_manifest(
    output_dir: Path,
    manifest_name: str,
    *,
    stages_run: list[str],
    stages_skipped: list[str],
    align_settings: dict | None,
    rpf_settings: dict | None,
    rnaseq_settings: dict | None,
    config_canonical: dict,
    config_source_path: str | None,
    sample_sheet_path: str | None,
    command_argv: list[str],
    runtimes: dict[str, float],
    skip_reasons: dict[str, str],
    total_runtime_seconds: float | None = None,
) -> Path:
    """Write the v1.0.0 ``run_manifest.json`` for an ``all`` run.

    Layout (top-level keys):
    * ``schema_version``       — version of THIS layout (see
                                 :data:`MANIFEST_SCHEMA_VERSION`).
    * ``mitoribopy_version``
    * ``git_commit``           — current commit when run inside a repo,
                                 ``null`` otherwise.
    * ``command``              — the original argv joined with spaces.
    * ``config_source``        — path to the user-supplied YAML.
    * ``config_source_sha256`` — SHA256 of that file as the user wrote
                                 it (useful for "is this the same input
                                 I used last time?" diffs).
    * ``config_canonical``     — the merged + auto-wired config that
                                 actually drove the run.
    * ``sample_sheet``         — path to the unified sheet, when set.
    * ``sample_sheet_sha256``  — its SHA256.
    * ``reference_checksum``   — promoted from rpf settings (legacy
                                 field; the rnaseq subcommand reads
                                 this directly).
    * ``stages``               — ``{align: {...}, rpf: {...}, rnaseq: {...}}``
                                 with ``status`` + optional
                                 ``runtime_seconds`` / ``reason``.
    * ``align`` / ``rpf`` / ``rnaseq``
                              — full per-stage run_settings.json
                                copies (kept as-is for back-compat).
    * ``tools``                — flat ``{tool: version}`` map lifted
                                 from per-stage settings.
    * ``warnings``             — structured warnings collected by
                                 ``mitoribopy.io.warnings_log`` during
                                 the run, mirrored to ``warnings.tsv``.
    * ``outputs``              — outputs_index rows (output_type,
                                 stage, path, description,
                                 recommended_for, schema_version).
                                 Same data is also written to
                                 ``outputs_index.tsv`` (P5.8).
    * ``runtime_seconds``      — total wall time for the orchestrator;
                                 per-stage values live under
                                 ``stages.<stage>.runtime_seconds`` (P5.8).
    * ``platform``             — :func:`platform.platform` string (P5.8).
    * ``python_version``       — running Python version (P5.8).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    from ..io.outputs_index import (
        build_outputs_index_rows,
        write_outputs_index,
    )
    from ..io.schema_versions import OUTPUT_SCHEMA_VERSIONS
    from ..io.warnings_log import collected as _collected_warnings
    from ..io.warnings_log import flush_tsv as _flush_warnings_tsv

    structured_warnings = [w.as_dict() for w in _collected_warnings()]
    _flush_warnings_tsv(output_dir / "warnings.tsv")
    # P5.8: write outputs_index.tsv before the manifest so the manifest
    # can reference the same row set under "outputs".
    write_outputs_index(output_dir)
    outputs_rows = build_outputs_index_rows(output_dir)

    # v0.6.2: lift resource_plan.json into the manifest so a downstream
    # script does not have to re-read the sidecar file. The path is
    # written by the orchestrator before any stage runs; we tolerate a
    # missing file (e.g. permissions issues) and just record null.
    resource_plan_blob: dict | None = None
    plan_path = output_dir / "resource_plan.json"
    if plan_path.is_file():
        try:
            resource_plan_blob = json.loads(plan_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            resource_plan_blob = None

    manifest: dict = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "subcommand": "all",
        "mitoribopy_version": __version__,
        "git_commit": _git_commit(),
        "command": "mitoribopy all " + " ".join(command_argv),
        "config_source": config_source_path,
        "config_source_sha256": (
            _sha256_of(Path(config_source_path)) if config_source_path else None
        ),
        "config_canonical": config_canonical,
        "sample_sheet": sample_sheet_path,
        "sample_sheet_sha256": (
            _sha256_of(Path(sample_sheet_path)) if sample_sheet_path else None
        ),
        "stages": _build_stages_block(
            stages_run, stages_skipped, runtimes, skip_reasons
        ),
        "align": align_settings,
        "rpf": rpf_settings,
        "rnaseq": rnaseq_settings,
        "tools": _lift_tool_versions(align_settings, rpf_settings, rnaseq_settings),
        # P1.12: schema versions for every advertised output TSV. A
        # downstream script can `jq .output_schemas.te_tsv` and gate on
        # a compatible major version.
        "output_schemas": dict(OUTPUT_SCHEMA_VERSIONS),
        # P1.11: structured warnings collected by mitoribopy.io.warnings_log
        # during the run; mirrored to <output>/warnings.tsv for diff-friendly
        # consumption.
        "warnings": structured_warnings,
        # P5.8: outputs_index — same rows as outputs_index.tsv. Lets a
        # script depending on the manifest alone enumerate every TSV /
        # plot / sidecar the run advertised.
        "outputs": outputs_rows,
        # P5.8: top-level runtime + platform metadata required by §8.
        "runtime_seconds": (
            float(total_runtime_seconds)
            if total_runtime_seconds is not None
            else sum(runtimes.values())
        ),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        # v0.6.2: orchestrator-level resource plan (mirrors
        # <run_root>/resource_plan.json). null when the orchestrator
        # could not write the sidecar.
        "resource_plan": resource_plan_blob,
    }

    # Promote rpf's reference_checksum so future rnaseq invocations
    # reading this manifest directly can find it without drilling into
    # the rpf section.
    if rpf_settings and rpf_settings.get("reference_checksum"):
        manifest["reference_checksum"] = rpf_settings["reference_checksum"]

    path = output_dir / manifest_name
    path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    return path


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


def _resolve_rnaseq_mode_from_config(rnaseq_section: dict) -> tuple[str, str | None]:
    """Resolve the rnaseq stage mode from a YAML rnaseq section.

    Mirrors :func:`mitoribopy.cli.rnaseq._resolve_rnaseq_mode` but
    operates on a config dict instead of an argparse Namespace, so the
    orchestrator can resolve the mode BEFORE dispatching to the rnaseq
    subcommand.

    Returns ``(mode, error_message)``; the caller exits 2 when
    ``error_message`` is set.
    """
    from .rnaseq import RNASEQ_MODES, _normalize_mode

    requested = _normalize_mode(
        rnaseq_section.get("rnaseq_mode")
        or rnaseq_section.get("rnaseq-mode")
        or rnaseq_section.get("mode")
    )
    has_de_table = bool(
        rnaseq_section.get("de_table") or rnaseq_section.get("de-table")
    )
    has_fastq_inputs = bool(
        rnaseq_section.get("rna_fastq")
        or rnaseq_section.get("rna-fastq")
        or rnaseq_section.get("sample_sheet")
        or rnaseq_section.get("sample-sheet")
    )
    inferred = (
        "de_table" if has_de_table
        else "from_fastq" if has_fastq_inputs
        else "none"
    )

    if requested is None:
        # Ambiguous: user supplied inputs for BOTH flows without picking
        # one. Fail loudly with the historical "mutually exclusive"
        # phrasing so existing scripts that grep for that string keep
        # working.
        if has_de_table and has_fastq_inputs:
            return inferred, (
                "rnaseq section has both 'de_table' and 'rna_fastq' / "
                "'sample_sheet'; the two flows are mutually exclusive. "
                "Drop one, or set rnaseq.mode={de_table|from_fastq} to "
                "disambiguate."
            )
        return inferred, None
    if requested not in RNASEQ_MODES:
        return inferred, (
            f"rnaseq.mode={requested!r} is not a valid mode; "
            f"choose one of {list(RNASEQ_MODES)}."
        )
    if requested == "de_table" and has_fastq_inputs:
        return requested, (
            "rnaseq.mode=de_table conflicts with rna_fastq / sample_sheet "
            "inputs; drop the FASTQ inputs or switch to mode=from_fastq."
        )
    if requested == "from_fastq" and has_de_table:
        return requested, (
            "rnaseq.mode=from_fastq conflicts with de_table input; drop "
            "de_table or switch to mode=de_table."
        )
    if requested == "de_table" and not has_de_table:
        return requested, "rnaseq.mode=de_table requires de_table: <path>."
    if requested == "from_fastq" and not has_fastq_inputs:
        return requested, (
            "rnaseq.mode=from_fastq requires either rna_fastq or sample_sheet."
        )
    return requested, None


def _auto_wire_paths(
    config: dict,
    *,
    run_root: Path,
    has_align: bool,
    has_rpf: bool,
    has_rnaseq: bool,
    has_de_table: bool,
    has_fastq_mode: bool,
) -> None:
    """Set stage-specific --output defaults and cross-stage wiring.

    After this call:
    * ``config['align']['output']``  defaults to ``<run_root>/align``.
    * ``config['rpf']['output']``    defaults to ``<run_root>/rpf``.
    * ``config['rpf']['directory']`` defaults to the align BED dir.
    * ``config['rnaseq']['output']`` defaults to ``<run_root>/rnaseq`` whenever
      either rnaseq mode is active (de_table or from-FASTQ).
    * ``config['rnaseq']['ribo_dir']`` defaults to the rpf output (de_table flow only).
    * ``config['rnaseq']['reference_fasta']`` defaults to ``rpf.fasta``
      when from-FASTQ mode is active and the user did not set one
      explicitly. Override by setting ``rnaseq.reference_fasta`` when
      RNA-seq uses a different transcriptome reference.
    """
    align_cfg = config.setdefault("align", {}) if has_align else {}
    rpf_cfg = config.setdefault("rpf", {}) if has_rpf else {}
    rnaseq_cfg = config.setdefault("rnaseq", {}) if has_rnaseq else {}

    if has_align:
        align_cfg.setdefault("output", str(run_root / "align"))

    if has_rpf:
        rpf_cfg.setdefault("output", str(run_root / "rpf"))
        # When align ran, its BED6 outputs live at <align>/bed/.
        if has_align:
            rpf_cfg.setdefault("directory", str(run_root / "align" / "bed"))
            rpf_cfg.setdefault(
                "read_counts_file", str(run_root / "align" / "read_counts.tsv")
            )

    if has_rnaseq and (has_de_table or has_fastq_mode):
        rnaseq_cfg.setdefault("output", str(run_root / "rnaseq"))
        # de_table flow consumes the rpf run dir; from-FASTQ flow uses
        # --ribo-fastq instead and ignores --ribo-dir.
        if has_de_table and has_rpf:
            rnaseq_cfg.setdefault("ribo-dir", str(run_root / "rpf"))
        if has_fastq_mode:
            user_reference = (
                rnaseq_cfg.get("reference_fasta")
                or rnaseq_cfg.get("reference-fasta")
            )
            rpf_fasta = config.get("rpf", {}).get("fasta") if has_rpf else None
            if rpf_fasta and not user_reference:
                rnaseq_cfg["reference_fasta"] = rpf_fasta
            # P0.2: when rpf will run as part of this `all` invocation,
            # default the rnaseq from-FASTQ stage to REUSE the rpf
            # stage's rpf_counts.tsv instead of independently re-aligning
            # the Ribo FASTQs. Users override by setting
            # `rnaseq.recount_ribo_fastq: true` in the config or by
            # passing --recount-ribo-fastq at the rnaseq subcommand
            # level. The wiring is a no-op when the user explicitly
            # set --recount-ribo-fastq to True.
            recount = bool(
                rnaseq_cfg.get("recount_ribo_fastq")
                or rnaseq_cfg.get("recount-ribo-fastq")
            )
            if has_rpf and not recount:
                rnaseq_cfg.setdefault(
                    "upstream_rpf_counts",
                    str(run_root / "rpf" / "rpf_counts.tsv"),
                )
            # Parallel align-only path: the orchestrator will run
            # `rnaseq --align-only` concurrently with `align`, producing
            # rna_counts.tsv before rpf even starts. The DE/TE
            # invocation at the end picks up the prebuilt matrix via
            # --upstream-rna-counts and skips re-alignment. Setting the
            # default here keeps the wiring symmetric with the Ribo
            # side and makes the canonical_config.yaml emitted at the
            # run root reflect the actual stage inputs.
            rnaseq_cfg.setdefault(
                "upstream_rna_counts",
                str(run_root / "rnaseq" / "rna_counts.tsv"),
            )


def _should_skip_align(output: Path) -> bool:
    return (output / "align" / "read_counts.tsv").is_file()


def _should_skip_rpf(output: Path) -> bool:
    return (output / "rpf" / "rpf_counts.tsv").is_file()


def _should_skip_rnaseq(output: Path) -> bool:
    return (output / "rnaseq" / "delta_te.tsv").is_file()


def run(argv: Iterable[str]) -> int:
    """Entry point for ``mitoribopy all <args>``."""
    parser = build_parser()
    args = parser.parse_args(list(argv))
    common.apply_common_arguments(args)

    if args.print_config_template:
        _print_config_template()
        return 0

    if args.print_canonical_config:
        if not args.config:
            print(
                "[mitoribopy all] ERROR: --print-canonical-config "
                "requires --config <path>.",
                file=sys.stderr,
            )
            return 2
        try:
            cfg = load_config_file(args.config)
        except (FileNotFoundError, RuntimeError, ValueError) as exc:
            print(f"[mitoribopy all] ERROR: {exc}", file=sys.stderr)
            return 2
        run_root = Path(args.output) if args.output else Path("results")
        rc = _apply_top_level_samples(cfg, run_root=run_root)
        if rc != 0:
            return rc
        rnaseq_section = cfg.get("rnaseq", {}) or {}
        rnaseq_mode, mode_err = _resolve_rnaseq_mode_from_config(rnaseq_section)
        if mode_err is not None:
            print(f"[mitoribopy all] ERROR: {mode_err}", file=sys.stderr)
            return 2
        if cfg.get("rnaseq") is not None:
            cfg["rnaseq"]["rnaseq_mode"] = rnaseq_mode
        _auto_wire_paths(
            cfg,
            run_root=run_root,
            has_align=bool(cfg.get("align")) and not args.skip_align,
            has_rpf=bool(cfg.get("rpf")) and not args.skip_rpf,
            has_rnaseq=bool(cfg.get("rnaseq")) and not args.skip_rnaseq,
            has_de_table=rnaseq_mode == "de_table",
            has_fastq_mode=rnaseq_mode == "from_fastq",
        )
        sys.stdout.write(_yaml_dump(cfg))
        return 0

    if args.show_stage_help:
        from . import align as align_cli
        from . import rnaseq as rnaseq_cli
        from ..config import DEFAULT_CONFIG
        from ..pipeline import runner as pipeline_runner

        stage_parsers = {
            "align": align_cli.build_parser,
            "rpf": lambda: pipeline_runner.build_parser(dict(DEFAULT_CONFIG)),
            "rnaseq": rnaseq_cli.build_parser,
        }
        stage_parsers[args.show_stage_help]().print_help()
        return 0

    if not args.config:
        if args.dry_run:
            return common.emit_dry_run(
                "all",
                [
                    "load --config YAML/JSON/TOML with align: / rpf: / rnaseq: sections",
                    "auto-wire stage --output and cross-stage --directory / --ribo-dir",
                    "run align -> rpf -> (optional) rnaseq, honoring --resume / --skip-*",
                    "write run_manifest.json with tool versions and reference_checksum",
                ],
            )
        print(
            "[mitoribopy all] ERROR: --config <path> is required; it holds "
            "the align / rpf / rnaseq sections.",
            file=sys.stderr,
        )
        return 2

    if not args.output and not args.dry_run:
        print(
            "[mitoribopy all] ERROR: --output <dir> is required (the run "
            "root for align/, rpf/, rnaseq/).",
            file=sys.stderr,
        )
        return 2

    try:
        config = load_config_file(args.config)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"[mitoribopy all] ERROR: {exc}", file=sys.stderr)
        return 2

    # Top-level --strict: run a validate-config preflight up front so a
    # deprecated-key rewrite or missing input never makes it as far as
    # the align stage. Mirrors `mitoribopy validate-config <path>
    # --strict` exactly. This is a hard gate -- if validate-config
    # rejects the config, the run aborts with the validator's exit
    # code so the user sees the same message they would in CI.
    if getattr(args, "strict", False):
        from . import validate_config as _validate_config_cli

        _vc_argv = [args.config, "--strict"]
        _vc_rc = _validate_config_cli.run(_vc_argv)
        if _vc_rc != 0:
            print(
                "[mitoribopy all] ERROR: --strict preflight "
                "(validate-config) rejected the config; aborting before "
                "any stage runs.",
                file=sys.stderr,
            )
            return _vc_rc

    run_root = Path(args.output) if args.output else Path(".")

    # Top-level `samples:` block (unified sample sheet) is the canonical
    # input description for the run. Apply it FIRST so it can populate
    # align.fastq / align.sample_overrides + rnaseq.sample_sheet before
    # the rest of the gating runs.
    rc = _apply_top_level_samples(config, run_root=run_root)
    if rc != 0:
        return rc

    has_align = bool(config.get("align")) and not args.skip_align
    has_rpf = bool(config.get("rpf")) and not args.skip_rpf
    has_rnaseq = bool(config.get("rnaseq")) and not args.skip_rnaseq
    rnaseq_section = config.get("rnaseq", {}) if has_rnaseq else {}
    rnaseq_mode, mode_err = _resolve_rnaseq_mode_from_config(rnaseq_section)
    if mode_err is not None:
        print(f"[mitoribopy all] ERROR: {mode_err}", file=sys.stderr)
        return 2
    has_de_table = has_rnaseq and rnaseq_mode == "de_table"
    has_fastq_mode = has_rnaseq and rnaseq_mode == "from_fastq"
    # Pin the resolved mode in the config so the downstream rnaseq
    # subcommand sees an explicit --rnaseq-mode flag and the manifest's
    # config_canonical blob records the resolved value.
    if has_rnaseq:
        config["rnaseq"]["rnaseq_mode"] = rnaseq_mode

    _auto_wire_paths(
        config,
        run_root=run_root,
        has_align=has_align,
        has_rpf=has_rpf,
        has_rnaseq=has_rnaseq,
        has_de_table=has_de_table,
        has_fastq_mode=has_fastq_mode,
    )

    # v0.6.2: cascade the top-level `execution:` block into stage
    # configs so a single thread / memory budget is honoured by every
    # stage. Stage-specific overrides still win.
    _apply_execution_block(
        config, cli_threads=getattr(args, "threads", None),
    )
    # Spec: top-level `periodicity:` block flows into rpf-stage flags
    # so users can write the spec-shaped YAML and have it Just Work.
    _apply_periodicity_block(config)

    # v0.6.2: write the run-root resource_plan.json EARLY so reviewers
    # have an audit artifact even if a downstream stage crashes. The
    # plan is recomputed by the per-stage align CLI under
    # <run_root>/align/resource_plan.json — those values may differ
    # because the align CLI sees the actual sample list resolution.
    if not args.dry_run and args.output:
        try:
            _resolve_and_write_run_root_resource_plan(
                config,
                run_root=run_root,
                n_samples=_count_samples_from_config(config),
            )
        except Exception as exc:  # pragma: no cover - belt-and-braces
            sys.stderr.write(
                "[mitoribopy all] WARNING: could not write "
                f"resource_plan.json: {exc}\n"
            )

    # Refactor-4 (report §3.5.C): canonical_config.yaml is now a
    # first-class output of every real run. Write it as soon as the
    # config is fully resolved (auto-wiring + samples-sheet expansion +
    # rnaseq_mode resolution) so even a crash mid-run leaves a record
    # of exactly what was about to be executed. The same blob is
    # embedded under run_manifest.json's `config_canonical` field.
    if not args.dry_run and args.output:
        try:
            run_root.mkdir(parents=True, exist_ok=True)
            (run_root / "canonical_config.yaml").write_text(
                _yaml_dump(config), encoding="utf-8"
            )
        except OSError as exc:  # pragma: no cover - filesystem corner cases
            sys.stderr.write(
                "[mitoribopy all] WARNING: could not write "
                f"canonical_config.yaml: {exc}\n"
            )

    if args.dry_run:
        plan: list[str] = []
        # When rnaseq is configured in from-FASTQ mode, the orchestrator
        # kicks off `rnaseq --align-only` in parallel with the align
        # stage so RNA cutadapt + bowtie2 happens concurrently with
        # Ribo align. The dry-run plan reflects this so users see the
        # new ordering at a glance.
        if has_rnaseq and has_fastq_mode:
            rnaseq_align_only_cfg = dict(config.get("rnaseq", {}))
            rnaseq_align_only_cfg["align_only"] = True
            for k in (
                "upstream_rpf_counts", "upstream-rpf-counts",
                "upstream_rna_counts", "upstream-rna-counts",
                "ribo_fastq", "ribo-fastq",
            ):
                rnaseq_align_only_cfg.pop(k, None)
            plan.append(
                "rnaseq (align-only, parallel with align): "
                + " ".join(
                    _dict_to_argv(rnaseq_align_only_cfg, flag_style="hyphen")
                )
            )
        if has_align:
            # Honour `mitoribopy all --strict` in the dry-run plan too,
            # so the printed argv matches what the real run would emit.
            if getattr(args, "strict", False):
                config["align"] = {
                    **config.get("align", {}),
                    "strict_publication_mode": True,
                }
            plan.append(
                "align: "
                + " ".join(
                    _dict_to_argv(
                        _normalize_align_inputs(
                            config.get("align", {}), run_root=None
                        ),
                        flag_style="hyphen",
                        repeat_flags={"fastq"},
                    )
                )
            )
        if has_rpf:
            plan.append(
                "rpf: "
                + " ".join(
                    _dict_to_argv(
                        config.get("rpf", {}),
                        flag_style="hyphen",
                        flag_overrides={"rpf": "-rpf"},
                    )
                )
            )
        if has_rnaseq and (has_de_table or has_fastq_mode):
            plan.append(
                "rnaseq: "
                + " ".join(_dict_to_argv(config.get("rnaseq", {}), flag_style="hyphen"))
            )
        if not plan:
            plan.append("(no stages selected after --skip-*/config evaluation)")
        plan.append(f"write manifest to {run_root / args.manifest}")
        return common.emit_dry_run("all", plan)

    stages_run: list[str] = []
    stages_skipped: list[str] = []
    runtimes: dict[str, float] = {}
    skip_reasons: dict[str, str] = {}
    # P5.8: top-level wall clock for the manifest's runtime_seconds.
    t_total_start = time.monotonic()

    # Resume hash guard.
    #
    # Only consulted when --resume is set. We compare the prior run's
    # recorded hashes (config_source_sha256, sample_sheet_sha256,
    # reference_checksum, mitoribopy_version, schema_version) against
    # the current values; on mismatch we refuse to skip the affected
    # stages so users do not re-use stale outputs that encode different
    # decisions. --force-resume bypasses the check, as does the
    # MITORIBOPY_FORCE_RESUME=1 env var.
    sample_sheet_for_guard: str | None = None
    samples_block = config.get("samples")
    if isinstance(samples_block, str):
        sample_sheet_for_guard = samples_block
    elif isinstance(samples_block, dict):
        sample_sheet_for_guard = samples_block.get("table") or samples_block.get("path")
    rpf_fasta_for_guard = config.get("rpf", {}).get("fasta") if has_rpf else None

    force_resume = force_resume_requested(cli_flag=args.force_resume)
    resume_active = args.resume or args.force_resume
    resume_report = None
    if resume_active:
        prior_manifest = load_prior_manifest(run_root, args.manifest)
        resume_report = validate_resume(
            prior_manifest=prior_manifest,
            config_path=args.config,
            sample_sheet_path=sample_sheet_for_guard,
            reference_fasta=rpf_fasta_for_guard,
            mitoribopy_version=__version__,
            manifest_schema_version=MANIFEST_SCHEMA_VERSION,
        )
        if not resume_report.ok and not force_resume:
            # Hard fail at the orchestrator level: we will not silently
            # produce a half-fresh, half-stale run. The error names
            # every field that drifted so the user can fix or override.
            print(
                "[mitoribopy all] ERROR: --resume cannot proceed: "
                + resume_report.render(),
                file=sys.stderr,
            )
            print(
                "[mitoribopy all] HINT: re-run without --resume to "
                "produce a fresh result, or pass --force-resume to "
                "bypass the hash guard (use only when you know the "
                "stale outputs are still valid).",
                file=sys.stderr,
            )
            return 2
        if not resume_report.ok and force_resume:
            sys.stderr.write(
                "[mitoribopy all] WARNING: --force-resume bypassing "
                "hash guard despite drift:\n"
                + resume_report.render()
                + "\n"
            )

    # Import subcommand entry points here to avoid a circular import at
    # module load time.
    from . import align as align_cli
    from . import rnaseq as rnaseq_cli
    from . import rpf as rpf_cli
    from ..io.outputs_index import write_outputs_index as _write_outputs_index
    from ..io.warnings_log import flush_tsv as _flush_warnings_tsv
    from ..progress import ProgressManager

    # Refactor-4: warnings.tsv and outputs_index.tsv must exist on disk
    # at the run root no matter what — even if a stage crashes before
    # the manifest write, downstream tooling (and the assessment §8
    # output contract) relies on the paths existing. Touch them now
    # with header-only content; they are re-written by _write_manifest
    # at the end of the run.
    run_root.mkdir(parents=True, exist_ok=True)
    _flush_warnings_tsv(run_root / "warnings.tsv")
    _write_outputs_index(run_root)

    # Event-driven progress system. The default --progress-file location
    # is <run-root>/progress.jsonl unless the user overrode it. The
    # manager is a context manager so file handles flush even if a
    # stage raises.
    progress_mode = "off" if args.no_progress else args.progress
    progress_file = args.progress_file
    if progress_file is None and progress_mode != "off":
        progress_file = run_root / "progress.jsonl"
    progress_mgr = ProgressManager.from_cli(
        mode=progress_mode,
        progress_file=progress_file,
    )
    progress_mgr.run_start(
        subcommand="all",
        config_source=args.config,
    )

    def _emit_resume_skip(stage: str, reason: str) -> None:
        progress_mgr.resume_skip(stage=stage, reason=reason)

    # --- align (Ribo) ‖ rnaseq align-only (RNA) -------------------------
    #
    # Historically these ran serially: align (Ribo only), then rpf, then
    # rnaseq (which did its own RNA cutadapt + bowtie2 internally before
    # DE/TE). For from-FASTQ rnaseq runs, the RNA alignment is the
    # second-largest CPU sink in the pipeline and there is no
    # data-dependency between Ribo align and RNA align — they read
    # disjoint FASTQs, write disjoint outputs, and only converge in the
    # DE/TE phase. Running them concurrently turns the wall time from
    # roughly (Ribo + RNA) into roughly max(Ribo, RNA), which matters
    # most when the pipeline is invoked on a multi-core box.
    #
    # When rnaseq is configured AND in from_fastq mode, we kick off a
    # `rnaseq --align-only` worker concurrently with the align stage.
    # The worker writes `rna_counts.tsv` (and optionally
    # `rpf_counts.tsv` / `rpf_counts_matrix.tsv` if Ribo FASTQs were
    # supplied directly to rnaseq, though the orchestrator's default
    # path defers Ribo counting to the rpf stage). The DE/TE phase then
    # runs at the end with --upstream-rna-counts pointing at the
    # prebuilt matrix.
    rna_align_only_active = (
        has_rnaseq
        and has_fastq_mode
        and not (resume_active and _should_skip_rnaseq(run_root))
    )
    rna_align_only_future: "Future | None" = None
    rna_align_only_executor: "ThreadPoolExecutor | None" = None
    rna_align_only_t0: float | None = None

    def _build_align_only_argv() -> list[str]:
        """Argv for the parallel `rnaseq --align-only` invocation."""
        rnaseq_align_only_cfg = dict(config["rnaseq"])
        rnaseq_align_only_cfg["output"] = str(run_root / "rnaseq")
        rnaseq_align_only_cfg["align_only"] = True
        # Do NOT propagate --upstream-rpf-counts to the align-only
        # invocation. We don't want the worker to try to load the rpf
        # counts file (it does not exist yet — rpf hasn't run) and there
        # is no Ribo work for the align-only worker to do anyway: the
        # orchestrator's default flow keeps Ribo counting in the rpf
        # stage. The worker only handles RNA FASTQs.
        rnaseq_align_only_cfg.pop("upstream_rpf_counts", None)
        rnaseq_align_only_cfg.pop("upstream-rpf-counts", None)
        # The auto-wire step set upstream_rna_counts to the same file
        # this worker is about to write. Pop it so the worker actually
        # runs alignment (otherwise reuse_upstream_rna would short-
        # circuit the worker and produce no counts).
        rnaseq_align_only_cfg.pop("upstream_rna_counts", None)
        rnaseq_align_only_cfg.pop("upstream-rna-counts", None)
        # Drop any Ribo FASTQ inputs the user happened to set on the
        # rnaseq stage; the rpf stage owns those.
        rnaseq_align_only_cfg.pop("ribo_fastq", None)
        rnaseq_align_only_cfg.pop("ribo-fastq", None)
        return _dict_to_argv(rnaseq_align_only_cfg, flag_style="hyphen")

    if rna_align_only_active:
        rna_align_only_executor = ThreadPoolExecutor(max_workers=1)
        rna_align_only_t0 = time.monotonic()
        rna_align_only_argv = _build_align_only_argv()
        progress_mgr.stage_start("rnaseq_prealign")
        rna_align_only_future = rna_align_only_executor.submit(
            rnaseq_cli.run, rna_align_only_argv,
        )

    if has_align:
        if resume_active and _should_skip_align(run_root):
            stages_skipped.append("align")
            skip_reasons["align"] = "resume: read_counts.tsv already exists"
            _emit_resume_skip("align", skip_reasons["align"])
            progress_mgr.stage_end(
                "align", status="skipped", reason=skip_reasons["align"]
            )
        else:
            align_cfg = _normalize_align_inputs(config["align"], run_root=run_root)
            if resume_active:
                # Propagate the orchestrator's --resume into the align
                # CLI so it skips per-sample work that already completed
                # (read_counts.tsv is missing -- we are running the
                # stage -- but individual .sample_done/<sample>.json
                # markers may exist from a previous crash).
                align_cfg = {**align_cfg, "resume": True}
            if getattr(args, "strict", False):
                # Top-level --strict propagates into the align stage as
                # --strict-publication-mode. The align CLI honours an
                # explicit value over the YAML's, so this still respects
                # an opt-OUT in the config.
                align_cfg = {**align_cfg, "strict_publication_mode": True}
            align_argv = _dict_to_argv(
                align_cfg,
                flag_style="hyphen",
                repeat_flags={"fastq"},
            )
            progress_mgr.stage_start("align")
            t0 = time.monotonic()
            rc = align_cli.run(align_argv)
            runtimes["align"] = time.monotonic() - t0
            if rc != 0:
                progress_mgr.stage_end(
                    "align", status="error",
                    elapsed_seconds=runtimes["align"],
                    reason=f"exit_code={rc}",
                )
                # Cancel / drain any in-flight RNA align-only worker so
                # we don't keep cutadapt + bowtie2 running for an aborted
                # pipeline.
                if rna_align_only_future is not None:
                    rna_align_only_future.cancel()
                if rna_align_only_executor is not None:
                    rna_align_only_executor.shutdown(wait=False, cancel_futures=True)
                # Re-flush warnings.tsv + outputs_index so they always
                # reflect the partial state the run produced before
                # the failure.
                _flush_warnings_tsv(run_root / "warnings.tsv")
                _write_outputs_index(run_root)
                progress_mgr.run_end(
                    status="error",
                    elapsed_seconds=time.monotonic() - t_total_start,
                )
                progress_mgr.close()
                print(
                    f"[mitoribopy all] align stage failed with exit code {rc}.",
                    file=sys.stderr,
                )
                return rc
            stages_run.append("align")
            progress_mgr.stage_end(
                "align", status="done", elapsed_seconds=runtimes["align"]
            )
    else:
        if args.skip_align:
            stages_skipped.append("align")
            skip_reasons["align"] = "--skip-align flag set"
        # Otherwise: no align section -> stage is "not_configured"
        # in the manifest (we leave it out of stages_skipped).

    # Wait for the parallel RNA align-only worker to finish before
    # entering rpf. Its outputs are not consumed by rpf, but joining
    # here gives the user a clean stage boundary in the progress log
    # and surfaces RNA-stage errors before the rpf stage runs and
    # blames itself for resource pressure.
    if rna_align_only_future is not None:
        try:
            rc_rna_align = rna_align_only_future.result()
        except Exception as exc:  # noqa: BLE001 — surface anything cleanly
            rc_rna_align = -1
            sys.stderr.write(
                f"[mitoribopy all] RNA align-only worker raised: {exc}\n"
            )
        finally:
            if rna_align_only_executor is not None:
                rna_align_only_executor.shutdown(wait=False)
        if rna_align_only_t0 is not None:
            runtimes["rnaseq_prealign"] = time.monotonic() - rna_align_only_t0
        if rc_rna_align != 0:
            progress_mgr.stage_end(
                "rnaseq_prealign", status="error",
                elapsed_seconds=runtimes.get("rnaseq_prealign"),
                reason=f"exit_code={rc_rna_align}",
            )
            _flush_warnings_tsv(run_root / "warnings.tsv")
            _write_outputs_index(run_root)
            progress_mgr.run_end(
                status="error",
                elapsed_seconds=time.monotonic() - t_total_start,
            )
            progress_mgr.close()
            print(
                "[mitoribopy all] RNA align-only worker failed with "
                f"exit code {rc_rna_align}.",
                file=sys.stderr,
            )
            return int(rc_rna_align) if rc_rna_align >= 0 else 2
        stages_run.append("rnaseq_prealign")
        progress_mgr.stage_end(
            "rnaseq_prealign", status="done",
            elapsed_seconds=runtimes.get("rnaseq_prealign"),
        )

    # --- rpf ------------------------------------------------------------
    if has_rpf:
        if resume_active and _should_skip_rpf(run_root):
            stages_skipped.append("rpf")
            skip_reasons["rpf"] = "resume: rpf_counts.tsv already exists"
            _emit_resume_skip("rpf", skip_reasons["rpf"])
            progress_mgr.stage_end(
                "rpf", status="skipped", reason=skip_reasons["rpf"]
            )
        else:
            rpf_argv = _dict_to_argv(
                config["rpf"],
                flag_style="hyphen",
                flag_overrides={"rpf": "-rpf"},
            )
            progress_mgr.stage_start("rpf")
            t0 = time.monotonic()
            rc = rpf_cli.run(rpf_argv)
            runtimes["rpf"] = time.monotonic() - t0
            if rc != 0:
                progress_mgr.stage_end(
                    "rpf", status="error",
                    elapsed_seconds=runtimes["rpf"],
                    reason=f"exit_code={rc}",
                )
                # Re-flush warnings.tsv + outputs_index so they always
                # reflect the partial state the run produced before
                # the failure.
                _flush_warnings_tsv(run_root / "warnings.tsv")
                _write_outputs_index(run_root)
                progress_mgr.run_end(
                    status="error",
                    elapsed_seconds=time.monotonic() - t_total_start,
                )
                progress_mgr.close()
                print(
                    f"[mitoribopy all] rpf stage failed with exit code {rc}.",
                    file=sys.stderr,
                )
                return rc
            stages_run.append("rpf")
            progress_mgr.stage_end(
                "rpf", status="done", elapsed_seconds=runtimes["rpf"]
            )
    else:
        if args.skip_rpf:
            stages_skipped.append("rpf")
            skip_reasons["rpf"] = "--skip-rpf flag set"
        # Otherwise: no rpf section -> not_configured.

    # --- rnaseq (optional) ---------------------------------------------
    # Run when either flow is configured: from-FASTQ (rna_fastq +
    # reference_fasta) or external DE table (de_table). The two are
    # mutually exclusive and validated above.
    if has_rnaseq and (has_de_table or has_fastq_mode):
        if resume_active and _should_skip_rnaseq(run_root):
            stages_skipped.append("rnaseq")
            skip_reasons["rnaseq"] = "resume: delta_te.tsv already exists"
            _emit_resume_skip("rnaseq", skip_reasons["rnaseq"])
            progress_mgr.stage_end(
                "rnaseq", status="skipped", reason=skip_reasons["rnaseq"]
            )
        else:
            rnaseq_argv = _dict_to_argv(config["rnaseq"], flag_style="hyphen")
            progress_mgr.stage_start("rnaseq")
            t0 = time.monotonic()
            rc = rnaseq_cli.run(rnaseq_argv)
            runtimes["rnaseq"] = time.monotonic() - t0
            if rc != 0:
                progress_mgr.stage_end(
                    "rnaseq", status="error",
                    elapsed_seconds=runtimes["rnaseq"],
                    reason=f"exit_code={rc}",
                )
                # Re-flush warnings.tsv + outputs_index so they always
                # reflect the partial state the run produced before
                # the failure.
                _flush_warnings_tsv(run_root / "warnings.tsv")
                _write_outputs_index(run_root)
                progress_mgr.run_end(
                    status="error",
                    elapsed_seconds=time.monotonic() - t_total_start,
                )
                progress_mgr.close()
                print(
                    f"[mitoribopy all] rnaseq stage failed with exit code {rc}.",
                    file=sys.stderr,
                )
                return rc
            stages_run.append("rnaseq")
            progress_mgr.stage_end(
                "rnaseq", status="done", elapsed_seconds=runtimes["rnaseq"]
            )
    else:
        if args.skip_rnaseq:
            stages_skipped.append("rnaseq")
            skip_reasons["rnaseq"] = "--skip-rnaseq flag set"
        elif has_rnaseq:
            # The rnaseq section was declared but neither flow's inputs
            # were supplied — treat as actively skipped (with a reason)
            # so users can see why their rnaseq config did not fire.
            stages_skipped.append("rnaseq")
            skip_reasons["rnaseq"] = (
                "rnaseq section present but neither --de-table nor "
                "--rna-fastq / sample_sheet input was configured"
            )
        # Otherwise: no rnaseq section at all -> not_configured.

    # --- manifest -------------------------------------------------------
    sample_sheet_path = None
    samples_block = config.get("samples")
    if isinstance(samples_block, str):
        sample_sheet_path = samples_block
    elif isinstance(samples_block, dict):
        sample_sheet_path = samples_block.get("table") or samples_block.get("path")

    _write_manifest(
        output_dir=run_root,
        manifest_name=args.manifest,
        stages_run=stages_run,
        stages_skipped=stages_skipped,
        align_settings=_read_stage_settings(run_root / "align"),
        rpf_settings=_read_stage_settings(run_root / "rpf"),
        rnaseq_settings=_read_stage_settings(run_root / "rnaseq"),
        config_canonical=config,
        config_source_path=args.config,
        sample_sheet_path=sample_sheet_path,
        command_argv=list(argv),
        runtimes=runtimes,
        skip_reasons=skip_reasons,
        total_runtime_seconds=time.monotonic() - t_total_start,
    )

    # P1.6 + P1.8: emit SUMMARY.md and summary_qc.tsv automatically so
    # users always get a one-glance view of the run. Failures here are
    # non-fatal (the manifest is the source of truth); we only log.
    try:
        from . import summarize as _summarize_cli

        _summarize_cli.run([str(run_root), "--manifest", args.manifest])
    except Exception as exc:  # pragma: no cover - belt-and-braces
        sys.stderr.write(
            f"[mitoribopy all] WARNING: SUMMARY.md generation failed: {exc}\n"
        )

    # Auto-run validate-figures so figure_qc.tsv is always produced.
    # The run's exit code is never elevated based on plot QC — the
    # contract for `all` is "did the pipeline finish", not "are the
    # figures publication-ready". Users who need that gate run
    # `mitoribopy validate-figures --strict` themselves.
    try:
        from . import validate_figures as _validate_figures_cli

        # Top-level --strict promotes warn-only QC findings to fail in
        # the post-run figure-QC pass. The pipeline still returns 0
        # overall — the gate is `mitoribopy validate-figures --strict`
        # itself, run from CI on the manifest.
        _vf_argv = [str(run_root)]
        if getattr(args, "strict", False):
            _vf_argv.append("--strict")
        _validate_figures_cli.run(_vf_argv)
    except Exception as exc:  # pragma: no cover - belt-and-braces
        sys.stderr.write(
            f"[mitoribopy all] WARNING: validate-figures run failed: {exc}\n"
        )

    # Re-write outputs_index.tsv so figure_qc.tsv and progress.jsonl
    # get advertised even when their writers ran AFTER the manifest
    # call above.
    try:
        from ..io.outputs_index import write_outputs_index as _wri

        _wri(run_root)
    except Exception:  # pragma: no cover
        pass

    # Emit a final RunEnd event with the warning / error totals lifted
    # from warnings_log so a downstream JSONL consumer can compute
    # success without re-parsing the manifest.
    try:
        from ..io.warnings_log import collected as _wl_collected

        _wl_records = _wl_collected()
        n_warn = sum(1 for r in _wl_records if r.severity == "warn")
        n_err = sum(1 for r in _wl_records if r.severity == "error")
    except Exception:  # pragma: no cover
        n_warn = 0
        n_err = 0
    progress_mgr.run_end(
        status="done",
        elapsed_seconds=time.monotonic() - t_total_start,
        n_warnings=n_warn,
        n_errors=n_err,
    )
    progress_mgr.close()
    return 0
