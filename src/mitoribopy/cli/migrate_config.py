"""``mitoribopy migrate-config <old.yaml>`` — rewrite legacy keys to canonical.

P1.10. Reads a YAML / JSON / TOML config, applies the rewrites in
:mod:`mitoribopy.config.migrate`, and writes the canonical YAML to
stdout. Per-rewrite log lines go to stderr so users can audit what
changed.

Exit codes:
* 0 — config parsed and migrated (even when no rewrites were needed)
* 2 — config could not be loaded (missing file, parse error, etc.)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

from ..config.migrate import migrate
from . import common


MIGRATE_CONFIG_HELP = (
    "Rewrite legacy MitoRiboPy YAML keys to their canonical names. "
    "Input is read from a path; output is written to stdout (the change "
    "log goes to stderr). Use to upgrade old pipeline configs without "
    "manually hunting down every renamed key."
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mitoribopy migrate-config",
        description=MIGRATE_CONFIG_HELP,
        formatter_class=common.MitoRiboPyHelpFormatter,
    )
    parser.add_argument(
        "config",
        metavar="PATH",
        help="Path to the legacy YAML / JSON / TOML config file.",
    )
    parser.add_argument(
        "--in-place",
        action="store_true",
        default=False,
        help=(
            "Overwrite the input file in place (with a `.bak` backup) "
            "instead of writing the canonical config to stdout."
        ),
    )
    return parser


def _yaml_dump(payload: dict) -> str:
    """YAML dump if PyYAML is available; JSON fallback otherwise."""
    try:
        import yaml  # type: ignore[import-not-found]
    except ImportError:
        import json
        return json.dumps(payload, indent=2, sort_keys=True) + "\n"
    return yaml.safe_dump(
        payload, sort_keys=True, default_flow_style=False, allow_unicode=True
    )


def run(argv: Iterable[str]) -> int:
    args = build_parser().parse_args(list(argv))

    try:
        cfg = common.load_config_file(args.config)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"[mitoribopy migrate-config] ERROR: {exc}", file=sys.stderr)
        return 2

    canonical, log = migrate(cfg)

    for line in log:
        sys.stderr.write(f"[mitoribopy migrate-config] {line}\n")
    if not log:
        sys.stderr.write(
            "[mitoribopy migrate-config] no legacy keys found; "
            "config already canonical.\n"
        )

    body = _yaml_dump(canonical)
    if args.in_place:
        path = Path(args.config)
        backup = path.with_suffix(path.suffix + ".bak")
        backup.write_bytes(path.read_bytes())
        path.write_text(body, encoding="utf-8")
        sys.stderr.write(
            f"[mitoribopy migrate-config] wrote canonical config to "
            f"{path} (original backed up at {backup}).\n"
        )
    else:
        sys.stdout.write(body)
    return 0
