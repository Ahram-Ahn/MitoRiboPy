"""P0.4 — documentation-execution tests.

Three guarantees this test module enforces:

1. **Templates parse.** Every YAML under ``examples/templates/`` loads
   without error via the same ``load_config_file`` the CLIs use.
2. **Template keys map to real CLI flags.** Every key in each template
   resolves to an argparse ``dest`` on the matching subcommand's
   parser, so a template can never silently drift away from the
   CLI's actual surface.
3. **README dry-run examples are runnable.** Every ``mitoribopy ...``
   command shown in a fenced code block in ``README.md`` either
   carries ``--dry-run`` (or can have it appended safely) is
   executed via the in-process CLI; exit 0 is required and stderr
   must not contain ``ERROR:``.
4. **Output contract.** A toy ``mitoribopy all --dry-run`` lists the
   key documented output paths in the planned actions.
"""

from __future__ import annotations

import re
import shlex
from pathlib import Path

import pytest

from mitoribopy.cli.common import load_config_file


REPO_ROOT = Path(__file__).resolve().parents[1]
TEMPLATES_DIR = REPO_ROOT / "examples" / "templates"
README_PATH = REPO_ROOT / "README.md"


# Each template is bound to the subcommand whose parser declares its
# keys. The mapping drives both the parse and the key-to-flag tests.
# `pipeline_config.example.yaml` is the multi-section orchestrator
# config and is validated by section against the right per-stage
# parser instead.
_PER_STAGE_TEMPLATES = {
    "align_config.example.yaml": "align",
    "rpf_config.example.yaml": "rpf",
    "rnaseq_config.example.yaml": "rnaseq",
}


# Top-level keys the orchestrator owns directly (not delegated to a
# stage parser); these never need to map to a per-stage flag.
_ORCHESTRATOR_OWNED_KEYS: set[str] = {
    "samples",
    "align",
    "rpf",
    "rnaseq",
    "execution",   # v0.6.2: top-level resource-plan declaration
}


def _stage_parser(stage: str):
    """Return the argparse parser for one stage subcommand."""
    if stage == "align":
        from mitoribopy.cli import align as align_cli
        return align_cli.build_parser()
    if stage == "rpf":
        from mitoribopy.config import DEFAULT_CONFIG
        from mitoribopy.pipeline import runner as pipeline_runner
        return pipeline_runner.build_parser(dict(DEFAULT_CONFIG))
    if stage == "rnaseq":
        from mitoribopy.cli import rnaseq as rnaseq_cli
        return rnaseq_cli.build_parser()
    raise ValueError(f"Unknown stage {stage!r}")


def _parser_dests(stage: str) -> set[str]:
    """Set of argparse `dest` names for the stage's parser.

    Includes both underscore- and hyphen-form keys to match the YAML
    key normalization rules in `_dict_to_argv`. Also includes the
    long option strings (with the leading ``--`` stripped) so YAML
    keys that match the public flag name (instead of the shorter
    argparse `dest`) are accepted — e.g.
    ``allow_pseudo_replicates_for_demo_not_publication`` whose dest
    is the shorter ``allow_pseudo_replicates``.
    """
    parser = _stage_parser(stage)
    dests: set[str] = set()
    for action in parser._actions:
        if action.dest and action.dest != "help":
            dests.add(action.dest)
            dests.add(action.dest.replace("_", "-"))
        for opt in action.option_strings:
            if opt.startswith("--"):
                bare = opt[2:]
                dests.add(bare)
                dests.add(bare.replace("-", "_"))
    return dests


# ---------- 1. templates parse ---------------------------------------------


@pytest.mark.parametrize(
    "template_name",
    sorted(t for t in _PER_STAGE_TEMPLATES) + ["pipeline_config.example.yaml"],
)
def test_template_parses(template_name: str) -> None:
    """Every YAML template under examples/templates/ loads cleanly."""
    path = TEMPLATES_DIR / template_name
    assert path.is_file(), f"missing template: {path}"
    cfg = load_config_file(str(path))
    assert isinstance(cfg, dict), f"{template_name} did not parse to a dict"


# ---------- 2. template keys map to real CLI flags --------------------------


@pytest.mark.parametrize(
    "template_name,stage", sorted(_PER_STAGE_TEMPLATES.items())
)
def test_per_stage_template_keys_map_to_real_flags(
    template_name: str, stage: str
) -> None:
    """Every top-level key in a per-stage template must correspond to a
    real argparse `dest` on that stage's parser. This is the test that
    catches "I added a key to the YAML but never wired it up" drift."""
    cfg = load_config_file(str(TEMPLATES_DIR / template_name))
    dests = _parser_dests(stage)
    unknown = [k for k in cfg if k not in dests]
    assert not unknown, (
        f"{template_name} has key(s) not declared by mitoribopy {stage} "
        f"argparse: {unknown}. Either add the flag or remove the key."
    )


def test_pipeline_template_section_keys_map_to_per_stage_flags() -> None:
    """For the multi-section orchestrator template, validate each stage
    section against the correct stage parser's dests."""
    cfg = load_config_file(str(TEMPLATES_DIR / "pipeline_config.example.yaml"))
    # Top-level keys: orchestrator-owned only.
    unknown_top = [k for k in cfg if k not in _ORCHESTRATOR_OWNED_KEYS]
    assert not unknown_top, (
        "pipeline_config.example.yaml has unknown top-level key(s): "
        f"{unknown_top}. Recognised: {sorted(_ORCHESTRATOR_OWNED_KEYS)}"
    )
    for stage in ("align", "rpf", "rnaseq"):
        section = cfg.get(stage)
        if not section:
            continue
        dests = _parser_dests(stage)
        # `samples` inside the align section is the per-stage overrides
        # block (a list of mappings) that `mitoribopy all` materialises
        # into a sample_overrides TSV; it has no direct argparse dest.
        # `rnaseq_mode` / `mode` is a key the orchestrator pins onto the
        # rnaseq section before serialising it; the rnaseq parser knows
        # `--rnaseq-mode` so the mapped dest is `rnaseq_mode`. The bare
        # `mode` alias is accepted in YAML and treated as a synonym.
        # `recount_ribo_fastq` and `upstream_rpf_counts` come from
        # P0.2 wiring; both have argparse dests on the rnaseq parser.
        allowed_extra: set[str] = set()
        if stage == "align":
            allowed_extra.update({"samples"})
        if stage == "rnaseq":
            allowed_extra.update({"mode"})
        unknown = [
            k for k in section
            if k not in dests and k not in allowed_extra
        ]
        assert not unknown, (
            f"pipeline_config.example.yaml [{stage}] has key(s) not "
            f"declared by mitoribopy {stage} argparse: {unknown}."
        )


# ---------- 3. README dry-run examples are runnable -------------------------


_README_FENCE_RE = re.compile(
    r"```(?:bash|shell|sh)?\n(.*?)\n```",
    re.DOTALL,
)


def _extract_mitoribopy_commands(readme_text: str) -> list[str]:
    """Extract `mitoribopy ...` invocations from fenced code blocks.

    Filters to commands that:
    * start with `mitoribopy `
    * are single-line (multi-line backslash-joined commands are
      collapsed first)
    * mention `--dry-run` (we only execute documented dry-runs).
    """
    commands: list[str] = []
    for block in _README_FENCE_RE.findall(readme_text):
        # Collapse continuation lines (`\` at end of line).
        joined = re.sub(r"\\\n\s*", " ", block)
        for line in joined.splitlines():
            stripped = line.strip()
            if not stripped.startswith("mitoribopy "):
                continue
            if "--dry-run" not in stripped:
                continue
            # Strip leading `$ ` shell prompt if present.
            if stripped.startswith("$ "):
                stripped = stripped[2:]
            commands.append(stripped)
    return commands


def test_readme_has_dry_run_examples() -> None:
    """If this fails, README has lost its `--dry-run` smoke examples."""
    if not README_PATH.is_file():
        pytest.skip("README.md not present")
    cmds = _extract_mitoribopy_commands(README_PATH.read_text(encoding="utf-8"))
    # We do not assert a specific count — just that at least one
    # documented dry-run example exists for the harness to exercise.
    assert cmds, (
        "README.md contains no `mitoribopy ... --dry-run` examples; the "
        "doc-execution test cannot exercise documented commands."
    )


# README examples use placeholder paths (e.g. `pipeline_config.yaml`)
# that the user is expected to create. The harness substitutes those
# placeholders with the shipped template files so the shapes are
# exercised end-to-end without requiring the user to materialise the
# placeholder paths first.
_README_PATH_SUBSTITUTIONS = {
    "pipeline_config.yaml": str(TEMPLATES_DIR / "pipeline_config.example.yaml"),
    "align_config.yaml": str(TEMPLATES_DIR / "align_config.example.yaml"),
    "rpf_config.yaml": str(TEMPLATES_DIR / "rpf_config.example.yaml"),
    "rnaseq_config.yaml": str(TEMPLATES_DIR / "rnaseq_config.example.yaml"),
}


def _substitute_placeholders(tokens: list[str]) -> list[str]:
    return [_README_PATH_SUBSTITUTIONS.get(tok, tok) for tok in tokens]


def test_readme_dry_run_examples_execute_cleanly(capsys) -> None:
    """Every `mitoribopy ... --dry-run` command in README.md must:
    1. parse via the in-process CLI dispatcher (after placeholder
       paths like ``pipeline_config.yaml`` are substituted with shipped
       templates so the user-facing command shape is what is exercised),
    2. exit 0,
    3. not emit `ERROR:` on stderr."""
    if not README_PATH.is_file():
        pytest.skip("README.md not present")

    from mitoribopy import cli

    cmds = _extract_mitoribopy_commands(README_PATH.read_text(encoding="utf-8"))
    if not cmds:
        pytest.skip("no mitoribopy --dry-run examples in README.md")

    failures: list[str] = []
    for raw in cmds:
        # Drop the `mitoribopy` program token before handing to dispatch.
        try:
            tokens = shlex.split(raw)
        except ValueError as exc:
            failures.append(f"{raw!r}: shlex error {exc}")
            continue
        if tokens[:1] != ["mitoribopy"]:
            failures.append(f"{raw!r}: not a mitoribopy invocation")
            continue
        argv = _substitute_placeholders(tokens[1:])
        try:
            rc = cli.main(argv)
        except SystemExit as exc:
            rc = exc.code if isinstance(exc.code, int) else 1
        captured = capsys.readouterr()
        if rc != 0:
            failures.append(
                f"{raw!r}: exit {rc}; stderr={captured.err!r}"
            )
            continue
        if "ERROR:" in captured.err:
            failures.append(f"{raw!r}: stderr contains ERROR: {captured.err!r}")
    assert not failures, "\n".join(failures)


# ---------- 4. output contract: dry-run mentions key paths ------------------


# Documented top-level outputs that `mitoribopy all` advertises in
# README. The dry-run plan must mention each so any rename / move
# of an output path forces a doc update.
_REQUIRED_DRY_RUN_TOKENS = (
    "run_manifest.json",
    "align",
    "rpf",
)


def test_all_dry_run_lists_documented_output_paths(
    tmp_path: Path, capsys
) -> None:
    """A bare `mitoribopy all --dry-run` should mention the key files
    and directories the README advertises (run_manifest.json, align,
    rpf). Catches drift between README and orchestrator output layout."""
    from mitoribopy import cli

    rc = cli.main(["all", "--dry-run"])
    assert rc == 0
    out = capsys.readouterr().out
    missing = [tok for tok in _REQUIRED_DRY_RUN_TOKENS if tok not in out]
    assert not missing, (
        f"all --dry-run plan is missing documented token(s): {missing}; "
        f"plan output was:\n{out}"
    )


def test_print_canonical_config_round_trip_for_pipeline_template(
    tmp_path: Path, capsys
) -> None:
    """The shipped pipeline template must round-trip through
    --print-canonical-config without parser errors. Catches subtle
    schema drift (e.g., a renamed key in the template that the
    orchestrator's auto-wire path does not understand)."""
    from mitoribopy import cli

    template = TEMPLATES_DIR / "pipeline_config.example.yaml"
    # The template references paths that may not exist on disk; for
    # --print-canonical-config we only need it to PARSE + auto-wire,
    # not actually run any stage.
    rc = cli.main(
        [
            "all",
            "--print-canonical-config",
            "--config", str(template),
            "--output", str(tmp_path / "results"),
        ]
    )
    assert rc == 0, capsys.readouterr().err
    out = capsys.readouterr().out
    # Every documented section must show up in the canonical dump.
    assert "align" in out
    assert "rpf" in out
