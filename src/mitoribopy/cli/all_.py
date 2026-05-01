"""``mitoribopy all`` subcommand - align -> rpf -> (optional) rnaseq.

Phase 6 of the v0.3.0 refactor. Runs the three per-stage subcommands in
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
from pathlib import Path
from typing import Iterable

from .. import __version__
from ..cli.common import load_config_file
from . import common


# Bumped whenever the run_manifest.json layout changes in a way that
# breaks downstream consumers (added fields are minor; renamed or
# removed fields are major). Read it from your own scripts to gate on a
# compatible manifest shape.
MANIFEST_SCHEMA_VERSION = "1.0.0"


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
            "when <output>/rnaseq/delta_te.tsv is present."
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

    * ``flag_style="hyphen"`` (default, for ``align`` and ``rnaseq``):
      keys with underscores become ``--dashed-form`` flags.
    * ``flag_style="underscore"`` (for ``rpf``, whose argparse declares
      its flags with underscores like ``--offset_type``): keys are
      emitted verbatim with a ``--`` prefix.
    * Boolean ``True`` emits just the flag; ``False`` emits nothing.
    * Lists / tuples are emitted as ``--flag v1 v2 v3`` (nargs="+" style).
      Keys listed in ``repeat_flags`` are emitted as repeated ``--flag v``
      pairs for argparse options that use ``action="append"``.
    * ``None`` values are skipped.
    * ``flag_overrides`` maps config keys to exact legacy flag spellings.

    The per-stage ``flag_style`` is necessary because the ``rpf`` stage
    parser (``pipeline.runner``) declares its flags in the legacy
    underscore form (``--offset_type``, ``--min_5_offset``, ...), so
    emitting hyphenated flags would be rejected with
    ``unrecognized arguments``.
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

    from ..sample_sheet import SampleSheetError, load_sample_sheet

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
        if align_cfg.get("fastq") or align_cfg.get("fastq_dir"):
            print(
                "[mitoribopy all] ERROR: top-level 'samples:' is set, "
                "but 'align.fastq' / 'align.fastq_dir' is also declared. "
                "Drop one — either let the sample sheet drive align "
                "inputs, or remove the 'samples:' block.",
                file=sys.stderr,
            )
            return 2
        if align_cfg.get("samples") or align_cfg.get("sample_overrides"):
            print(
                "[mitoribopy all] ERROR: top-level 'samples:' is set, "
                "but 'align.samples' / 'align.sample_overrides' is also "
                "declared. The unified sheet supersedes per-stage "
                "overrides; drop one.",
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
            for r in ribo_rows
        ):
            tsv_path = Path(run_root) / "align" / "sample_overrides.tsv"
            tsv_path.parent.mkdir(parents=True, exist_ok=True)
            cols = (
                "sample", "kit_preset", "adapter",
                "umi_length", "umi_position", "dedup_strategy",
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
                        "dedup_strategy": _stringify(r.dedup_strategy),
                    }
                    handle.write("\t".join(row[c] for c in cols) + "\n")
            align_cfg["sample_overrides"] = str(tsv_path)

    # --- rnaseq side ----------------------------------------------------
    # The rnaseq subcommand has its own --sample-sheet handler that
    # rejects mixing the sheet with --rna-fastq / --condition-map / etc.
    # We just hand the path through.
    if rnaseq_explicit:
        rnaseq_conflicts = [
            key for key in (
                "rna_fastq", "rna-fastq", "ribo_fastq", "ribo-fastq",
                "condition_map", "condition-map",
            )
            if rnaseq_cfg.get(key)
        ]
        if rnaseq_conflicts and rnaseq_cfg.get("sample_sheet") is None:
            print(
                "[mitoribopy all] ERROR: top-level 'samples:' is set, but "
                "'rnaseq' section also declares "
                + ", ".join(repr(c) for c in rnaseq_conflicts)
                + ". The unified sheet supersedes per-stage inputs; drop one.",
                file=sys.stderr,
            )
            return 2
        rnaseq_cfg.setdefault("sample_sheet", str(sheet.path))

    return 0


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
  # allow_pretrimmed_inference: true      # set false to keep v0.4.0 hard-fail
                                          #   behaviour when detection fails

  # Inputs / reference indexes (bowtie2 prefixes built by bowtie2-build).
  # `fastq` accepts either a directory string OR an explicit list of paths.
  # Picked up patterns: *.fq, *.fq.gz, *.fastq, *.fastq.gz.
  fastq: input_data/              # directory shortcut (recommended)
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

  # Per-sample overrides (mixed-UMI batches). Use this when each FASTQ
  # has a different UMI length / position / kit. The 'name' field must
  # match the FASTQ basename with .fq[.gz] / .fastq[.gz] stripped.
  # Any unset override field falls through to the globals above.
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
# Default flow: from raw FASTQ (rna_fastq + reference_fasta).
# rnaseq:
#   rna_fastq: /path/to/rnaseq/                # dir or list of FASTQs
#   ribo_fastq: /path/to/riboseq/              # optional; reuses align/ when omitted
#   reference_fasta: /path/to/transcriptome.fa # auto-wired from rpf.fasta if unset
#   condition_map: /path/to/conditions.tsv
#   condition_a: control
#   condition_b: knockdown
#   gene_id_convention: hgnc                   # hgnc | ensembl | refseq | bare
#
# Alternative flow: bring your own DE table.
# rnaseq:
#   de_table: /path/to/de_table.tsv
#   gene_id_convention: hgnc
#   reference_gtf: /path/to/reference.fa       # must match align's --mt-index source
#   condition_map: /path/to/conditions.tsv
#   condition_a: control
#   condition_b: knockdown
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
    * ``warnings``             — placeholder for future structured
                                 warnings; always an empty list in v1.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

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
        "warnings": [],
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
        _auto_wire_paths(
            cfg,
            run_root=run_root,
            has_align=bool(cfg.get("align")) and not args.skip_align,
            has_rpf=bool(cfg.get("rpf")) and not args.skip_rpf,
            has_rnaseq=bool(cfg.get("rnaseq")) and not args.skip_rnaseq,
            has_de_table=bool(
                rnaseq_section.get("de_table") or rnaseq_section.get("de-table")
            ),
            has_fastq_mode=bool(
                rnaseq_section.get("rna_fastq")
                or rnaseq_section.get("rna-fastq")
                or rnaseq_section.get("sample_sheet")
                or rnaseq_section.get("sample-sheet")
            ),
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
    has_de_table = has_rnaseq and bool(
        rnaseq_section.get("de_table") or rnaseq_section.get("de-table")
    )
    has_fastq_mode = has_rnaseq and bool(
        rnaseq_section.get("rna_fastq")
        or rnaseq_section.get("rna-fastq")
        or rnaseq_section.get("sample_sheet")
        or rnaseq_section.get("sample-sheet")
    )
    if has_de_table and has_fastq_mode:
        print(
            "[mitoribopy all] ERROR: rnaseq section has both 'de_table' and "
            "'rna_fastq'; the two flows are mutually exclusive. Drop one.",
            file=sys.stderr,
        )
        return 2

    _auto_wire_paths(
        config,
        run_root=run_root,
        has_align=has_align,
        has_rpf=has_rpf,
        has_rnaseq=has_rnaseq,
        has_de_table=has_de_table,
        has_fastq_mode=has_fastq_mode,
    )

    if args.dry_run:
        plan: list[str] = []
        if has_align:
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
                        flag_style="underscore",
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

    # Import subcommand entry points here to avoid a circular import at
    # module load time.
    from . import align as align_cli
    from . import rnaseq as rnaseq_cli
    from . import rpf as rpf_cli

    # --- align ----------------------------------------------------------
    if has_align:
        if args.resume and _should_skip_align(run_root):
            stages_skipped.append("align")
            skip_reasons["align"] = "resume: read_counts.tsv already exists"
        else:
            align_cfg = _normalize_align_inputs(config["align"], run_root=run_root)
            if args.resume:
                # Propagate the orchestrator's --resume into the align
                # CLI so it skips per-sample work that already completed
                # (read_counts.tsv is missing -- we are running the
                # stage -- but individual .sample_done/<sample>.json
                # markers may exist from a previous crash).
                align_cfg = {**align_cfg, "resume": True}
            align_argv = _dict_to_argv(
                align_cfg,
                flag_style="hyphen",
                repeat_flags={"fastq"},
            )
            t0 = time.monotonic()
            rc = align_cli.run(align_argv)
            runtimes["align"] = time.monotonic() - t0
            if rc != 0:
                print(
                    f"[mitoribopy all] align stage failed with exit code {rc}.",
                    file=sys.stderr,
                )
                return rc
            stages_run.append("align")
    else:
        if args.skip_align:
            stages_skipped.append("align")
            skip_reasons["align"] = "--skip-align flag set"
        # Otherwise: no align section -> stage is "not_configured"
        # in the manifest (we leave it out of stages_skipped).

    # --- rpf ------------------------------------------------------------
    if has_rpf:
        if args.resume and _should_skip_rpf(run_root):
            stages_skipped.append("rpf")
            skip_reasons["rpf"] = "resume: rpf_counts.tsv already exists"
        else:
            rpf_argv = _dict_to_argv(
                config["rpf"],
                flag_style="underscore",
                flag_overrides={"rpf": "-rpf"},
            )
            t0 = time.monotonic()
            rc = rpf_cli.run(rpf_argv)
            runtimes["rpf"] = time.monotonic() - t0
            if rc != 0:
                print(
                    f"[mitoribopy all] rpf stage failed with exit code {rc}.",
                    file=sys.stderr,
                )
                return rc
            stages_run.append("rpf")
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
        if args.resume and _should_skip_rnaseq(run_root):
            stages_skipped.append("rnaseq")
            skip_reasons["rnaseq"] = "resume: delta_te.tsv already exists"
        else:
            rnaseq_argv = _dict_to_argv(config["rnaseq"], flag_style="hyphen")
            t0 = time.monotonic()
            rc = rnaseq_cli.run(rnaseq_argv)
            runtimes["rnaseq"] = time.monotonic() - t0
            if rc != 0:
                print(
                    f"[mitoribopy all] rnaseq stage failed with exit code {rc}.",
                    file=sys.stderr,
                )
                return rc
            stages_run.append("rnaseq")
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
    )
    return 0
