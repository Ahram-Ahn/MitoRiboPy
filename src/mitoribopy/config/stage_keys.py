"""Stage-level config key allowlists.

`mitoribopy validate-config` already catches unknown keys at the top
level (`align`, `rpf`, `rnaseq`, ...). This module supplies the
allowlists for one level deeper: the keys *inside* each stage section.

A typo like ``rpf.offset_pick_refernce: p_site`` previously slipped
through canonicalisation and silently produced results under default
settings. Publication-grade software must fail early on that.

Sources of truth
----------------
* ``rpf`` keys are the keys of :data:`mitoribopy.config.runtime.DEFAULT_CONFIG`
  (the runtime config dict the rpf parser populates from).
* ``align`` and ``rnaseq`` keys are derived from the live argparse
  parsers in :mod:`mitoribopy.cli.align` and
  :mod:`mitoribopy.cli.rnaseq`. We introspect ``parser._actions``
  rather than maintaining a parallel list, so a new flag landing in
  the parser is automatically accepted by the validator.
* ``execution`` and ``periodicity`` keys are the small literal sets
  defined in :mod:`mitoribopy.cli.all_`
  (``_EXECUTION_BLOCK_KEYS`` and ``_PERIODICITY_BLOCK_KEYS``).

Allowed extras
--------------
A handful of YAML-only conventions are not argparse flags but are
documented in templates and accepted by the orchestrator (e.g.
``align.samples`` — a per-sample override list materialised as a TSV
before the align CLI sees it). Those live in
:data:`_STAGE_YAML_EXTRAS`.
"""

from __future__ import annotations

import argparse
import difflib
from typing import Any, Iterable


_STAGE_YAML_EXTRAS: dict[str, frozenset[str]] = {
    "align": frozenset(
        {
            # Materialised into align/sample_overrides.tsv by
            # mitoribopy.cli.all_._normalize_align_inputs.
            "samples",
            # Convenience YAML form: a string value is treated as fastq_dir.
            "fastq",
            "fastq_dir",
        }
    ),
    "rpf": frozenset(
        {
            # mitoribopy.cli.all_ injects this from the resolved
            # `samples:` block; it is not part of DEFAULT_CONFIG but is
            # documented in templates.
            "fasta",
            "directory",
            "output",
            "downstream_dir",
        }
    ),
    "rnaseq": frozenset(
        {
            # Top-level `samples:` cascades into rnaseq.sample_sheet.
            "sample_sheet",
            "rnaseq_mode",
            "de_table",
            # Explicit override for `mitoribopy all --strict` users who
            # genuinely want the exploratory mt-mRNA-only from-FASTQ
            # path. Not consumed by the rnaseq subcommand itself; only
            # the orchestrator reads it.
            "allow_exploratory_from_fastq_in_strict",
            # The CLI flag is
            # ``--allow-pseudo-replicates-for-demo-not-publication``
            # whose argparse ``dest`` is ``allow_pseudo_replicates``.
            # Both spellings are accepted as YAML keys by the
            # orchestrator (``cli.all_._run_from_fastq`` and
            # ``cli.rnaseq._run_from_fastq`` both check the long name
            # alongside the dest), and the flat-template ships the
            # long name as the publication-canonical YAML form so the
            # opt-in is visibly named in the config that produced the
            # run. Whitelisting the long form here keeps the validator
            # in lockstep with the orchestrator's accepted-keys set.
            "allow_pseudo_replicates_for_demo_not_publication",
        }
    ),
    "execution": frozenset(),
    "periodicity": frozenset(),
}


_EXECUTION_KEYS: frozenset[str] = frozenset(
    {
        "threads",
        "memory_gb",
        "parallel_samples",
        "single_sample_mode",
        "min_threads_per_sample",
        "estimated_memory_per_sample_gb",
        "scheduler",
    }
)

_PERIODICITY_KEYS: frozenset[str] = frozenset(
    {
        # Keep this set in lockstep with the orchestrator's
        # ``_PERIODICITY_BLOCK_KEYS`` mapping in
        # ``mitoribopy.cli.all_`` — every key the orchestrator
        # cascades into ``rpf.periodicity_*`` must validate here, or
        # the publication-profile template the package itself emits
        # fails ``validate-config --strict``.
        "enabled",
        "fourier_window_nt",
        "metagene_nt",
        "metagene_normalize",
        "fourier_bootstrap_n",
        "fourier_permutations_n",
        "fourier_ci_alpha",
        "fourier_random_seed",
        "no_fourier_stats",
    }
)


def _argparse_dests(parser: argparse.ArgumentParser) -> set[str]:
    """Collect every flag's ``dest`` attribute from a parser."""
    dests: set[str] = set()
    for action in parser._actions:
        if action.dest in (None, argparse.SUPPRESS, "help"):
            continue
        dests.add(action.dest)
    return dests


def _safe_align_dests() -> set[str]:
    from ..cli.align import build_parser

    return _argparse_dests(build_parser())


def _safe_rnaseq_dests() -> set[str]:
    from ..cli.rnaseq import build_parser

    return _argparse_dests(build_parser())


def _safe_rpf_keys() -> set[str]:
    from .runtime import DEFAULT_CONFIG, DEPRECATED_CONFIG_KEY_ALIASES

    keys = set(DEFAULT_CONFIG.keys())
    keys.update(DEPRECATED_CONFIG_KEY_ALIASES.keys())
    return keys


def allowed_keys_for(stage: str) -> set[str]:
    """Return the union of recognised keys for *stage*.

    Includes argparse dests, YAML-only conventions, and (for ``rpf``)
    deprecated aliases that ``migrate-config`` rewrites silently.
    """
    if stage == "align":
        return _safe_align_dests() | set(_STAGE_YAML_EXTRAS["align"])
    if stage == "rpf":
        return _safe_rpf_keys() | set(_STAGE_YAML_EXTRAS["rpf"])
    if stage == "rnaseq":
        return _safe_rnaseq_dests() | set(_STAGE_YAML_EXTRAS["rnaseq"])
    if stage == "execution":
        return set(_EXECUTION_KEYS)
    if stage == "periodicity":
        return set(_PERIODICITY_KEYS)
    raise ValueError(f"Unknown stage: {stage!r}")


def find_unknown_keys(
    section: dict[str, Any], stage: str
) -> list[tuple[str, list[str]]]:
    """Return ``[(unknown_key, suggestions), ...]`` for *section*.

    Suggestions come from :func:`difflib.get_close_matches` against the
    allowed-keys set for *stage*. An empty suggestions list means
    nothing close was found.
    """
    if not isinstance(section, dict):
        return []
    allowed = allowed_keys_for(stage)
    unknown: list[tuple[str, list[str]]] = []
    for key in section.keys():
        if key in allowed:
            continue
        suggestions = difflib.get_close_matches(key, allowed, n=2, cutoff=0.6)
        unknown.append((key, suggestions))
    return unknown


def format_unknown_key_message(
    stage: str, key: str, suggestions: Iterable[str]
) -> str:
    """Single-line error string with optional ``did you mean`` tail."""
    suggestions = list(suggestions)
    if suggestions:
        return (
            f"{stage}.{key} is not a valid key. "
            f"Did you mean {' or '.join(repr(s) for s in suggestions)}?"
        )
    return f"{stage}.{key} is not a valid key for the {stage!r} section."
