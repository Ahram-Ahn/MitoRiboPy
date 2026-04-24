"""cutadapt wrapper with kit-aware adapter and UMI handling.

Step B of the ``mitoribopy align`` pipeline. Runs cutadapt with:

* a kit-preset-resolved 3' adapter (never silently defaulted; ``custom``
  preset requires ``--adapter``),
* 5' or 3' UMI extraction into the read name when ``umi_length > 0``,
* a length filter whose defaults (15-45 nt) match published mt-RPF
  distributions and are configurable via CLI.

For 3' UMI kits (QIAseq style) we run cutadapt in two passes because
cutadapt's ``-u -N`` operates *before* adapter search in a single call;
pass 1 trims the adapter, pass 2 extracts the 3' UMI.
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import replace
from pathlib import Path

from ._types import (
    KIT_PRESETS,
    CutadaptResult,
    KitPreset,
    ResolvedKit,
    UmiPosition,
    resolve_kit_alias,
)


# ---------------------------------------------------------------------------
# Kit resolution
# ---------------------------------------------------------------------------


def resolve_kit_settings(
    kit: str,
    adapter: str | None = None,
    umi_length: int | None = None,
    umi_position: UmiPosition | None = None,
) -> ResolvedKit:
    """Apply the kit preset's defaults, then the explicit CLI overrides.

    Raises
    ------
    KeyError
        if *kit* is not a known preset name.
    ValueError
        if the effective configuration is internally inconsistent
        (notably: ``custom`` with no ``--adapter``).
    """
    # Translate any legacy vendor-specific alias to its canonical
    # adapter-family preset name before lookup.
    kit = resolve_kit_alias(kit)
    try:
        preset: KitPreset = KIT_PRESETS[kit]
    except KeyError as exc:
        known = ", ".join(sorted(KIT_PRESETS))
        raise KeyError(
            f"Unknown kit preset: {kit!r}. Known presets: {known}."
        ) from exc

    if kit == "auto":
        raise ValueError(
            "Cannot resolve kit settings against the 'auto' sentinel. The "
            "per-sample resolver in mitoribopy.align.sample_resolve must "
            "translate 'auto' to a real preset before calling this function."
        )

    effective_adapter = adapter if adapter is not None else preset.adapter
    # The 'pretrimmed' preset legitimately has no adapter — cutadapt will
    # skip the -a flag and only do length + quality filtering. Every
    # other adapter-less resolution is an error (most often: 'custom'
    # without --adapter).
    if effective_adapter is None and kit != "pretrimmed":
        named = ", ".join(
            name for name in KIT_PRESETS if name not in {"custom", "auto"}
        )
        raise ValueError(
            "The 'custom' kit preset requires an explicit --adapter <SEQ>. "
            "Either pass --adapter, switch to a named kit preset "
            f"({named}), or use --kit-preset pretrimmed for already-trimmed "
            "FASTQs."
        )

    effective_umi_length = umi_length if umi_length is not None else preset.umi_length
    if effective_umi_length < 0:
        raise ValueError("--umi-length must be zero or a positive integer.")

    effective_umi_position = umi_position if umi_position is not None else preset.umi_position
    if effective_umi_position not in ("5p", "3p"):
        raise ValueError(
            f"--umi-position must be '5p' or '3p'; got {effective_umi_position!r}."
        )

    return ResolvedKit(
        kit=preset.name,
        adapter=effective_adapter,
        umi_length=effective_umi_length,
        umi_position=effective_umi_position,
    )


# ---------------------------------------------------------------------------
# Command construction
# ---------------------------------------------------------------------------


def _build_pass1_command(
    *,
    fastq_in: Path,
    fastq_out: Path,
    resolved: ResolvedKit,
    min_length: int,
    max_length: int,
    quality: int,
    threads: int,
    log_json: Path,
) -> list[str]:
    """Assemble the primary cutadapt invocation.

    For 5' UMI libraries the UMI is extracted in this pass (``-u N`` is
    applied before adapter search by cutadapt's internal ordering).

    For 3' UMI libraries the UMI is NOT handled here; the caller runs a
    second cutadapt pass via :func:`_build_pass2_umi_command`.
    """
    cmd: list[str] = ["cutadapt"]

    if resolved.umi_length > 0 and resolved.umi_position == "5p":
        # Cut UMI from the 5' end; place the removed bases into the read name.
        cmd += [
            "-u",
            str(resolved.umi_length),
            "--rename={id}_{cut_prefix} {comment}",
        ]

    # Skip the -a flag entirely for the 'pretrimmed' kit so cutadapt
    # only enforces length + quality filtering. This is the supported
    # path for SRA-deposited or already-trimmed FASTQs.
    if resolved.adapter is not None:
        cmd += ["-a", resolved.adapter]

    cmd += [
        "--minimum-length",
        str(min_length),
        "--maximum-length",
        str(max_length),
        "-q",
        str(quality),
    ]

    if threads and threads > 1:
        cmd += ["--cores", str(threads)]

    cmd += [
        "--json",
        str(log_json),
        "-o",
        str(fastq_out),
        str(fastq_in),
    ]
    return cmd


def _build_pass2_umi_command(
    *,
    fastq_in: Path,
    fastq_out: Path,
    umi_length: int,
    threads: int,
    log_json: Path,
) -> list[str]:
    """Assemble the 3'-UMI-extraction cutadapt invocation.

    Used only when the kit has a 3' UMI (e.g. QIAseq). Runs AFTER adapter
    trimming so that the last ``umi_length`` bases of the trimmed read
    are the UMI (they were between the insert and the 3' adapter on the
    original read).
    """
    cmd: list[str] = [
        "cutadapt",
        "-u",
        f"-{umi_length}",
        "--rename={id}_{cut_suffix} {comment}",
    ]
    if threads and threads > 1:
        cmd += ["--cores", str(threads)]
    cmd += [
        "--json",
        str(log_json),
        "-o",
        str(fastq_out),
        str(fastq_in),
    ]
    return cmd


# ---------------------------------------------------------------------------
# JSON parsing
# ---------------------------------------------------------------------------


def parse_cutadapt_json(log_json_path: Path) -> dict[str, int]:
    """Extract stage counts from a cutadapt ``--json`` log.

    Returns a small dict with:

    * ``input_reads``        reads fed into cutadapt
    * ``reads_with_adapter`` reads in which the 3' adapter was detected
    * ``reads_passing_filters`` reads that survived all length/quality cuts

    The schema is the stable cutadapt >= 4.0 JSON layout. Missing keys
    default to 0 so partial logs still parse rather than crash.
    """
    with Path(log_json_path).open("r", encoding="utf-8") as handle:
        raw = json.load(handle)

    read_counts = raw.get("read_counts", {}) or {}
    input_reads = int(read_counts.get("input", 0) or 0)
    reads_passing_filters = int(read_counts.get("output", 0) or 0)

    adapters = raw.get("adapters_read1") or []
    reads_with_adapter = 0
    for entry in adapters:
        if not isinstance(entry, dict):
            continue
        total_matches = entry.get("total_matches")
        if isinstance(total_matches, (int, float)):
            reads_with_adapter += int(total_matches)

    return {
        "input_reads": input_reads,
        "reads_with_adapter": reads_with_adapter,
        "reads_passing_filters": reads_passing_filters,
    }


# ---------------------------------------------------------------------------
# Public runner
# ---------------------------------------------------------------------------


def run_cutadapt(
    *,
    fastq_in: Path,
    fastq_out: Path,
    resolved: ResolvedKit,
    min_length: int = 15,
    max_length: int = 45,
    quality: int = 20,
    threads: int = 1,
    log_json: Path,
    runner=subprocess.run,
) -> CutadaptResult:
    """Run cutadapt end-to-end for a single sample.

    Two-pass when ``resolved.umi_position == "3p"`` and ``umi_length > 0``;
    single-pass otherwise.

    Parameters
    ----------
    fastq_in, fastq_out:
        Input FASTQ (``.fq`` / ``.fq.gz``) and the trimmed output path.
    resolved:
        :class:`~mitoribopy.align._types.ResolvedKit` from
        :func:`resolve_kit_settings`.
    min_length, max_length:
        Inclusive read-length filter in nt. Defaults 15-45 for mt-RPF.
    quality:
        Phred+33 3' quality trim threshold.
    threads:
        Passed through as ``--cores`` to cutadapt.
    log_json:
        Path for the primary cutadapt JSON log. For two-pass runs the
        pass-2 log sits next to it with a ``.pass2.json`` suffix.
    runner:
        Injection point used by tests; defaults to :func:`subprocess.run`.

    Returns
    -------
    :class:`~mitoribopy.align._types.CutadaptResult`
        Counts aggregated across all passes; the final trimmed FASTQ path
        is always the caller-supplied ``fastq_out``.
    """
    fastq_in = Path(fastq_in)
    fastq_out = Path(fastq_out)
    log_json = Path(log_json)
    log_json.parent.mkdir(parents=True, exist_ok=True)

    two_pass = resolved.umi_length > 0 and resolved.umi_position == "3p"

    if two_pass:
        intermediate = fastq_out.with_suffix(fastq_out.suffix + ".prepass.fq.gz")
        pass1_log = log_json.with_name(log_json.stem + ".pass1.json")
        pass2_log = log_json.with_name(log_json.stem + ".pass2.json")

        pass1_cmd = _build_pass1_command(
            fastq_in=fastq_in,
            fastq_out=intermediate,
            resolved=replace(resolved, umi_length=0),  # skip UMI handling in pass 1
            min_length=min_length,
            max_length=max_length + resolved.umi_length,  # keep +UMI length for pass 2
            quality=quality,
            threads=threads,
            log_json=pass1_log,
        )
        _invoke(runner, pass1_cmd)

        pass2_cmd = _build_pass2_umi_command(
            fastq_in=intermediate,
            fastq_out=fastq_out,
            umi_length=resolved.umi_length,
            threads=threads,
            log_json=pass2_log,
        )
        _invoke(runner, pass2_cmd)

        pass1_counts = parse_cutadapt_json(pass1_log)
        pass2_counts = parse_cutadapt_json(pass2_log)
        # Use the union: input counted at the very first step, passing
        # counted at the final step so the Phase H read-counts table
        # reflects reads actually written to the trimmed FASTQ.
        return CutadaptResult(
            input_reads=pass1_counts["input_reads"],
            reads_with_adapter=pass1_counts["reads_with_adapter"],
            reads_passing_filters=pass2_counts["reads_passing_filters"],
            log_json_path=log_json,
        )

    pass1_cmd = _build_pass1_command(
        fastq_in=fastq_in,
        fastq_out=fastq_out,
        resolved=resolved,
        min_length=min_length,
        max_length=max_length,
        quality=quality,
        threads=threads,
        log_json=log_json,
    )
    _invoke(runner, pass1_cmd)

    counts = parse_cutadapt_json(log_json)
    return CutadaptResult(
        input_reads=counts["input_reads"],
        reads_with_adapter=counts["reads_with_adapter"],
        reads_passing_filters=counts["reads_passing_filters"],
        log_json_path=log_json,
    )


def _invoke(runner, cmd: list[str]) -> None:
    """Run a subprocess command and raise on failure with the full stderr."""
    completed = runner(cmd, check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        stderr = (getattr(completed, "stderr", "") or "").strip()
        raise RuntimeError(
            f"cutadapt failed with exit code {completed.returncode}: "
            f"{stderr or '<no stderr captured>'}"
        )
