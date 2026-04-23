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
    rpf:     { strain: h, rpf: [29, 34], ... }
    rnaseq:  { de_table: de.tsv, gene_id_convention: hgnc, ... }

Each section's keys correspond to the CLI flag names of the matching
subcommand with dashes replaced by underscores. ``mitoribopy all``
reconstructs a flag list from the section and calls the subcommand's
``run()`` entry point.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Iterable

from .. import __version__
from ..cli.common import load_config_file
from . import common


ALL_SUBCOMMAND_HELP = (
    "End-to-end orchestrator: align + rpf, plus rnaseq when the config "
    "carries an 'rnaseq' section with --de-table. Writes a composed "
    "run_manifest.json with tool versions, parameters, and input/output "
    "hashes across all three stages."
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mitoribopy all",
        description=ALL_SUBCOMMAND_HELP,
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
    return parser


# ---------------------------------------------------------------------------
# argv reconstruction from a config section
# ---------------------------------------------------------------------------


def _dict_to_argv(section: dict) -> list[str]:
    """Serialize a section dict into a CLI-style argv list.

    Rules:

    * Keys with underscores become ``--dashed-form`` flags.
    * Boolean ``True`` emits just the flag; ``False`` emits nothing.
    * Lists are emitted as ``--flag v1 v2 v3`` (nargs="+" style).
    * ``None`` values are skipped.
    """
    argv: list[str] = []
    for key, value in section.items():
        flag = f"--{key.replace('_', '-')}"
        if value is None:
            continue
        if isinstance(value, bool):
            if value:
                argv.append(flag)
            continue
        if isinstance(value, (list, tuple)):
            argv.append(flag)
            argv.extend(str(v) for v in value)
            continue
        argv.extend([flag, str(value)])
    return argv


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


def _write_manifest(
    output_dir: Path,
    manifest_name: str,
    *,
    stages_run: list[str],
    stages_skipped: list[str],
    align_settings: dict | None,
    rpf_settings: dict | None,
    rnaseq_settings: dict | None,
) -> Path:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "subcommand": "all",
        "mitoribopy_version": __version__,
        "stages_run": stages_run,
        "stages_skipped": stages_skipped,
        "align": align_settings,
        "rpf": rpf_settings,
        "rnaseq": rnaseq_settings,
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
) -> None:
    """Set stage-specific --output defaults and cross-stage wiring.

    After this call:
    * ``config['align']['output']``  defaults to ``<run_root>/align``.
    * ``config['rpf']['output']``    defaults to ``<run_root>/rpf``.
    * ``config['rpf']['directory']`` defaults to the align BED dir.
    * ``config['rnaseq']['output']`` defaults to ``<run_root>/rnaseq``.
    * ``config['rnaseq']['ribo_dir']`` defaults to the rpf output.
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
            rpf_cfg.setdefault(
                "directory", str(run_root / "align" / "bed")
            )

    if has_rnaseq and has_de_table:
        rnaseq_cfg.setdefault("output", str(run_root / "rnaseq"))
        if has_rpf:
            rnaseq_cfg.setdefault("ribo-dir", str(run_root / "rpf"))


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
    has_align = bool(config.get("align")) and not args.skip_align
    has_rpf = bool(config.get("rpf")) and not args.skip_rpf
    has_rnaseq = bool(config.get("rnaseq")) and not args.skip_rnaseq
    has_de_table = has_rnaseq and bool(config.get("rnaseq", {}).get("de_table")
                                       or config.get("rnaseq", {}).get("de-table"))

    _auto_wire_paths(
        config,
        run_root=run_root,
        has_align=has_align,
        has_rpf=has_rpf,
        has_rnaseq=has_rnaseq,
        has_de_table=has_de_table,
    )

    if args.dry_run:
        plan: list[str] = []
        if has_align:
            plan.append(
                "align: " + " ".join(_dict_to_argv(config.get("align", {})))
            )
        if has_rpf:
            plan.append(
                "rpf: " + " ".join(_dict_to_argv(config.get("rpf", {})))
            )
        if has_rnaseq and has_de_table:
            plan.append(
                "rnaseq: " + " ".join(_dict_to_argv(config.get("rnaseq", {})))
            )
        if not plan:
            plan.append("(no stages selected after --skip-*/config evaluation)")
        plan.append(f"write manifest to {run_root / args.manifest}")
        return common.emit_dry_run("all", plan)

    stages_run: list[str] = []
    stages_skipped: list[str] = []

    # Import subcommand entry points here to avoid a circular import at
    # module load time.
    from . import align as align_cli
    from . import rnaseq as rnaseq_cli
    from . import rpf as rpf_cli

    # --- align ----------------------------------------------------------
    if has_align:
        if args.resume and _should_skip_align(run_root):
            stages_skipped.append("align")
        else:
            align_argv = _dict_to_argv(config["align"])
            rc = align_cli.run(align_argv)
            if rc != 0:
                print(
                    f"[mitoribopy all] align stage failed with exit code {rc}.",
                    file=sys.stderr,
                )
                return rc
            stages_run.append("align")
    else:
        stages_skipped.append("align")

    # --- rpf ------------------------------------------------------------
    if has_rpf:
        if args.resume and _should_skip_rpf(run_root):
            stages_skipped.append("rpf")
        else:
            rpf_argv = _dict_to_argv(config["rpf"])
            rc = rpf_cli.run(rpf_argv)
            if rc != 0:
                print(
                    f"[mitoribopy all] rpf stage failed with exit code {rc}.",
                    file=sys.stderr,
                )
                return rc
            stages_run.append("rpf")
    else:
        stages_skipped.append("rpf")

    # --- rnaseq (optional) ---------------------------------------------
    if has_rnaseq and has_de_table:
        if args.resume and _should_skip_rnaseq(run_root):
            stages_skipped.append("rnaseq")
        else:
            rnaseq_argv = _dict_to_argv(config["rnaseq"])
            rc = rnaseq_cli.run(rnaseq_argv)
            if rc != 0:
                print(
                    f"[mitoribopy all] rnaseq stage failed with exit code {rc}.",
                    file=sys.stderr,
                )
                return rc
            stages_run.append("rnaseq")
    else:
        stages_skipped.append("rnaseq")

    # --- manifest -------------------------------------------------------
    _write_manifest(
        output_dir=run_root,
        manifest_name=args.manifest,
        stages_run=stages_run,
        stages_skipped=stages_skipped,
        align_settings=_read_stage_settings(run_root / "align"),
        rpf_settings=_read_stage_settings(run_root / "rpf"),
        rnaseq_settings=_read_stage_settings(run_root / "rnaseq"),
    )
    return 0
