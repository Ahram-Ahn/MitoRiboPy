"""Legacy-key migration for MitoRiboPy YAML configs (P1.10).

The CLI accepts a number of legacy spellings for back-compat; over
time the documentation now leads with the canonical names. This
module rewrites a parsed YAML config dict to the canonical form,
returning the new dict alongside a human-readable change log.

Used by the standalone ``mitoribopy migrate-config`` subcommand and
also reusable by ``mitoribopy validate-config`` to print canonicalised
errors.
"""

from __future__ import annotations

from copy import deepcopy


__all__ = [
    "DEDUP_STRATEGY_VALUE_REWRITES",
    "LEGACY_ALIGN_TOPLEVEL",
    "LEGACY_RNASEQ_TOPLEVEL",
    "LEGACY_RPF_TOPLEVEL",
    "STRAIN_SHORTCUTS",
    "migrate",
]


# --- Legacy → canonical key maps -----------------------------------------


# rpf section (or top-level rpf-style configs).
LEGACY_RPF_TOPLEVEL: dict[str, str] = {
    "merge_density": "codon_density_window",
    "mrna_ref_patterns": "mt_mrna_substring_patterns",
}


# align section.
LEGACY_ALIGN_TOPLEVEL: dict[str, str] = {
    # `fastq_dir:` -> `fastq:` (string form). The string form of `fastq:`
    # is interpreted by `mitoribopy all` as a directory shortcut; see
    # `_normalize_align_inputs`.
    "fastq_dir": "fastq",
}


# Strain shortcuts.
STRAIN_SHORTCUTS: dict[str, str] = {
    "h": "h.sapiens",
    "y": "s.cerevisiae",
}


# `offset_pick_reference` value rewrite: `selected_site` -> `reported_site`.
# The KEY itself does not change; only the VALUE.
OFFSET_PICK_VALUE_REWRITES: dict[str, str] = {
    "selected_site": "reported_site",
}


# `dedup_strategy` value rewrite: implementation name -> statistical
# operation name. The KEY itself does not change; only the VALUE.
# Mirrors :data:`mitoribopy.align.dedup._DEDUP_STRATEGY_ALIASES`; keep
# the two in sync.
DEDUP_STRATEGY_VALUE_REWRITES: dict[str, str] = {
    "umi-tools": "umi_coordinate",
    "umi_tools": "umi_coordinate",
}


# Legacy rnaseq YAML keys -> canonical names. The
# `allow_pseudo_replicates` key is the only legacy spelling that
# requires special handling (it serialises to a never-existed CLI
# flag); migrate it to the long, intentionally-non-publication form.
LEGACY_RNASEQ_TOPLEVEL: dict[str, str] = {
    "allow_pseudo_replicates": "allow_pseudo_replicates_for_demo_not_publication",
}


# rnaseq mode value rewrites: hyphenated -> underscored.
RNASEQ_MODE_VALUE_REWRITES: dict[str, str] = {
    "from-fastq": "from_fastq",
    "de-table": "de_table",
}


# --- Migration ------------------------------------------------------------


def _rename_keys(
    section: dict, mapping: dict[str, str], *, log: list[str], path: str
) -> dict:
    """Rename keys in *section* per *mapping*, recording each change.

    Returns a NEW dict; the input is not mutated. When the canonical
    key already exists, the legacy entry is dropped (with a warning
    in the log) so the user does not get a silent shadow.
    """
    out: dict = {}
    for key, value in section.items():
        if key in mapping:
            new_key = mapping[key]
            if new_key in section:
                log.append(
                    f"{path}.{key}: legacy key dropped (canonical "
                    f"'{new_key}' already set)"
                )
                continue
            log.append(f"{path}.{key} -> {path}.{new_key}")
            out[new_key] = value
        else:
            out[key] = value
    return out


def _rewrite_strain(
    section: dict, *, log: list[str], path: str
) -> dict:
    out = dict(section)
    if "strain" in out and out["strain"] in STRAIN_SHORTCUTS:
        new = STRAIN_SHORTCUTS[out["strain"]]
        log.append(f"{path}.strain: '{out['strain']}' -> '{new}'")
        out["strain"] = new
    return out


def _drop_removed_kit_preset(
    section: dict, *, log: list[str], path: str
) -> dict:
    """Strip ``kit_preset`` keys (removed in v0.7.1) and record removals.

    Migration is intentionally lossy here: the user's previous kit name
    cannot be auto-translated into the new model (adapter sequence vs
    pretrimmed flag) without potentially picking the wrong default. We
    drop the key with a loud log entry so the user knows to add an
    explicit ``adapter:`` or ``pretrimmed:`` value if auto-detection
    would be insufficient.
    """
    out = dict(section)
    if "kit_preset" in out:
        log.append(
            f"{path}.kit_preset: REMOVED (no replacement); previous value "
            f"{out['kit_preset']!r}. Auto-detection now picks the kit; if "
            "you need to pin the adapter, add 'adapter: <SEQ>' (or "
            "'pretrimmed: true' for already-trimmed FASTQs)."
        )
        out.pop("kit_preset", None)
    samples = out.get("samples")
    if isinstance(samples, list):
        new_samples = []
        for index, entry in enumerate(samples):
            if isinstance(entry, dict) and "kit_preset" in entry:
                log.append(
                    f"{path}.samples[{index}].kit_preset: REMOVED "
                    f"(previous value {entry['kit_preset']!r})."
                )
                entry = {k: v for k, v in entry.items() if k != "kit_preset"}
            new_samples.append(entry)
        out["samples"] = new_samples
    return out


def _rewrite_offset_pick_reference(
    section: dict, *, log: list[str], path: str
) -> dict:
    out = dict(section)
    val = out.get("offset_pick_reference")
    if val in OFFSET_PICK_VALUE_REWRITES:
        new = OFFSET_PICK_VALUE_REWRITES[val]
        log.append(f"{path}.offset_pick_reference: '{val}' -> '{new}'")
        out["offset_pick_reference"] = new
    return out


def _rewrite_rnaseq_mode(
    section: dict, *, log: list[str], path: str
) -> dict:
    out = dict(section)
    for key in ("rnaseq_mode", "mode"):
        val = out.get(key)
        if val in RNASEQ_MODE_VALUE_REWRITES:
            new = RNASEQ_MODE_VALUE_REWRITES[val]
            log.append(f"{path}.{key}: '{val}' -> '{new}'")
            out[key] = new
    return out


def _rewrite_dedup_strategy(
    section: dict, *, log: list[str], path: str
) -> dict:
    out = dict(section)
    val = out.get("dedup_strategy")
    if isinstance(val, str) and val in DEDUP_STRATEGY_VALUE_REWRITES:
        new = DEDUP_STRATEGY_VALUE_REWRITES[val]
        log.append(f"{path}.dedup_strategy: '{val}' -> '{new}'")
        out["dedup_strategy"] = new
    # Per-sample overrides under align.samples[*].dedup_strategy too.
    samples = out.get("samples")
    if isinstance(samples, list):
        new_samples = []
        for index, entry in enumerate(samples):
            if (
                isinstance(entry, dict)
                and isinstance(entry.get("dedup_strategy"), str)
                and entry["dedup_strategy"] in DEDUP_STRATEGY_VALUE_REWRITES
            ):
                replacement = DEDUP_STRATEGY_VALUE_REWRITES[entry["dedup_strategy"]]
                log.append(
                    f"{path}.samples[{index}].dedup_strategy: "
                    f"'{entry['dedup_strategy']}' -> '{replacement}'"
                )
                entry = {**entry, "dedup_strategy": replacement}
            new_samples.append(entry)
        out["samples"] = new_samples
    return out


def migrate(raw_config: dict) -> tuple[dict, list[str]]:
    """Return ``(canonical_config, change_log)`` for a parsed YAML config.

    The input is deep-copied; the returned dict is independent.
    ``change_log`` is a list of human-readable lines, one per rewrite,
    suitable for printing to stderr by the CLI.
    """
    log: list[str] = []
    cfg = deepcopy(raw_config)

    # Sectioned (orchestrator) config.
    if "align" in cfg and isinstance(cfg["align"], dict):
        cfg["align"] = _rename_keys(
            cfg["align"], LEGACY_ALIGN_TOPLEVEL, log=log, path="align"
        )
        cfg["align"] = _drop_removed_kit_preset(
            cfg["align"], log=log, path="align"
        )
        cfg["align"] = _rewrite_dedup_strategy(
            cfg["align"], log=log, path="align"
        )
    if "rpf" in cfg and isinstance(cfg["rpf"], dict):
        cfg["rpf"] = _rename_keys(
            cfg["rpf"], LEGACY_RPF_TOPLEVEL, log=log, path="rpf"
        )
        cfg["rpf"] = _rewrite_strain(cfg["rpf"], log=log, path="rpf")
        cfg["rpf"] = _rewrite_offset_pick_reference(
            cfg["rpf"], log=log, path="rpf"
        )
    if "rnaseq" in cfg and isinstance(cfg["rnaseq"], dict):
        cfg["rnaseq"] = _rename_keys(
            cfg["rnaseq"], LEGACY_RNASEQ_TOPLEVEL, log=log, path="rnaseq"
        )
        cfg["rnaseq"] = _rewrite_rnaseq_mode(
            cfg["rnaseq"], log=log, path="rnaseq"
        )

    # Flat (per-stage) config: same rewrites, but applied at the top.
    # Heuristic: there is no align/rpf/rnaseq wrapper key.
    if not any(k in cfg for k in ("align", "rpf", "rnaseq")):
        cfg = _rename_keys(cfg, LEGACY_ALIGN_TOPLEVEL, log=log, path="")
        cfg = _rename_keys(cfg, LEGACY_RPF_TOPLEVEL, log=log, path="")
        cfg = _rename_keys(cfg, LEGACY_RNASEQ_TOPLEVEL, log=log, path="")
        cfg = _rewrite_strain(cfg, log=log, path="")
        cfg = _drop_removed_kit_preset(cfg, log=log, path="")
        cfg = _rewrite_offset_pick_reference(cfg, log=log, path="")
        cfg = _rewrite_rnaseq_mode(cfg, log=log, path="")
        cfg = _rewrite_dedup_strategy(cfg, log=log, path="")

    return cfg, log
