"""cutadapt wrapper with adapter-detection-aware UMI handling.

Step B of the ``mitoribopy align`` pipeline. Runs cutadapt with:

* a 3' adapter resolved by auto-detection or supplied via ``--adapter``
  (never silently defaulted),
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
)


# ---------------------------------------------------------------------------
# Kit resolution
# ---------------------------------------------------------------------------


def resolve_kit_settings(
    kit: str,
    adapter: str | None = None,
    umi_length: int | None = None,
    umi_position: UmiPosition | None = None,
    umi_length_5p: int | None = None,
    umi_length_3p: int | None = None,
) -> ResolvedKit:
    """Apply the kit preset's defaults, then the explicit CLI overrides.

    *kit* is the internal preset name produced by adapter detection
    (``illumina_smallrna``, ``illumina_truseq``, …) or one of the
    sentinels ``pretrimmed`` / ``custom``. It is no longer a user input
    since v0.7.1; the per-sample resolver in :mod:`mitoribopy.align.
    sample_resolve` constructs it from the detector output or the
    user's ``--adapter`` / ``--pretrimmed`` flags.

    For ``umi_position == "both"`` libraries the caller must supply
    ``umi_length_5p`` and ``umi_length_3p`` (or have them in the kit
    preset). When *both* per-end lengths and ``umi_length`` are set,
    ``umi_length`` MUST equal ``umi_length_5p + umi_length_3p`` — that
    is the canonical QNAME UMI length the downstream umi_tools dedup
    parses after the last ``_`` separator. When ``umi_length`` is not
    set explicitly we derive it from the per-end sum.

    Raises
    ------
    KeyError
        if *kit* is not a known preset name.
    ValueError
        if the effective configuration is internally inconsistent
        (notably: ``custom`` with no adapter, dual-end UMI without
        per-end lengths, or a per-end / total length mismatch).
    """
    try:
        preset: KitPreset = KIT_PRESETS[kit]
    except KeyError as exc:
        known = ", ".join(sorted(KIT_PRESETS))
        raise KeyError(
            f"Unknown kit preset: {kit!r}. Known presets: {known}."
        ) from exc

    effective_adapter = adapter if adapter is not None else preset.adapter
    # The 'pretrimmed' preset legitimately has no adapter — cutadapt will
    # skip the -a flag and only do length + quality filtering. Every
    # other adapter-less resolution is an error (most often: 'custom'
    # was selected by the resolver without an accompanying adapter
    # sequence, which should never happen but is checked defensively).
    if effective_adapter is None and kit != "pretrimmed":
        raise ValueError(
            f"Resolved kit {kit!r} has no adapter sequence. Pass --adapter "
            "<SEQ>, or use --pretrimmed for already-trimmed FASTQs."
        )

    effective_umi_length = umi_length if umi_length is not None else preset.umi_length
    if effective_umi_length < 0:
        raise ValueError("--umi-length must be zero or a positive integer.")

    effective_umi_position = umi_position if umi_position is not None else preset.umi_position
    if effective_umi_position not in ("5p", "3p", "both"):
        raise ValueError(
            f"--umi-position must be '5p', '3p', or 'both'; got "
            f"{effective_umi_position!r}."
        )

    effective_umi_5p = (
        int(umi_length_5p) if umi_length_5p is not None else int(preset.umi_length_5p)
    )
    effective_umi_3p = (
        int(umi_length_3p) if umi_length_3p is not None else int(preset.umi_length_3p)
    )
    if effective_umi_5p < 0 or effective_umi_3p < 0:
        raise ValueError(
            "--umi-length-5p and --umi-length-3p must be zero or positive integers."
        )

    if effective_umi_position == "both":
        if effective_umi_5p == 0 or effective_umi_3p == 0:
            raise ValueError(
                "umi_position='both' requires both --umi-length-5p and "
                "--umi-length-3p to be > 0 (got umi_length_5p="
                f"{effective_umi_5p}, umi_length_3p={effective_umi_3p})."
            )
        per_end_total = effective_umi_5p + effective_umi_3p
        # When umi_length wasn't explicitly supplied, derive it from the
        # per-end sum so callers can pass only the 5p/3p values.
        if umi_length is None and preset.umi_length == 0:
            effective_umi_length = per_end_total
        elif effective_umi_length != per_end_total:
            raise ValueError(
                "umi_position='both' requires umi_length to equal "
                "umi_length_5p + umi_length_3p (got umi_length="
                f"{effective_umi_length}, umi_length_5p={effective_umi_5p}, "
                f"umi_length_3p={effective_umi_3p})."
            )

    return ResolvedKit(
        kit=preset.name,
        adapter=effective_adapter,
        umi_length=effective_umi_length,
        umi_position=effective_umi_position,
        umi_length_5p=effective_umi_5p,
        umi_length_3p=effective_umi_3p,
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

    For ``umi_position == "both"`` libraries pass 1 extracts ONLY the
    5' UMI (``umi_length_5p`` bases). The 3' UMI is taken in pass 2
    along with the standard adapter trim.

    For 3' UMI libraries the UMI is NOT handled here; the caller runs a
    second cutadapt pass via :func:`_build_pass2_umi_command`.
    """
    cmd: list[str] = ["cutadapt"]

    pass1_umi_length = 0
    if resolved.umi_position == "5p" and resolved.umi_length > 0:
        pass1_umi_length = resolved.umi_length
    elif resolved.umi_position == "both" and resolved.umi_length_5p > 0:
        pass1_umi_length = resolved.umi_length_5p

    if pass1_umi_length > 0:
        # Cut UMI from the 5' end; place the removed bases into the read name.
        # In dual-end ('both') mode we must omit the literal `` {comment}``
        # tail: cutadapt would otherwise append a trailing space to every
        # pass-1 header, and pass 2's ``{id}{cut_suffix}`` concatenation
        # would re-emit that space between the two UMI tokens (breaking
        # the umi_tools dedup contract that the substring after the last
        # ``_`` is the UMI). 5p-only and 3p-only modes preserve any
        # incoming FASTQ comment as before.
        if resolved.umi_position == "both":
            rename_template = "--rename={id}_{cut_prefix}"
        else:
            rename_template = "--rename={id}_{cut_prefix} {comment}"
        cmd += [
            "-u",
            str(pass1_umi_length),
            rename_template,
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
    append_to_existing_umi: bool = False,
) -> list[str]:
    """Assemble the 3'-UMI-extraction cutadapt invocation.

    Used when the kit has a 3' UMI (e.g. QIAseq) or when
    ``umi_position == "both"`` (dual-end UMIs).

    For ``append_to_existing_umi=True`` (the dual-end case) the rename
    template is ``{id}{cut_suffix}`` — note the lack of an underscore
    separator. Pass 1 already wrote ``<originalid>_<5pUMI>`` into the
    QNAME, so appending ``{cut_suffix}`` here yields
    ``<originalid>_<5pUMI><3pUMI>``. umi_tools dedup parses the substring
    after the LAST ``_`` as the UMI, so the concatenated string becomes
    the combined dedup token without us needing to encode the boundary.

    For ``append_to_existing_umi=False`` (single-end 3' UMI) the rename
    is ``{id}_<3pUMI>``, the original behaviour.
    """
    if append_to_existing_umi:
        # Dual-end mode: omit ``{comment}`` so the QNAME ends cleanly at
        # the concatenated UMI string (no trailing whitespace before the
        # umi_tools dedup token).
        rename_template = "--rename={id}{cut_suffix}"
    else:
        rename_template = "--rename={id}_{cut_suffix} {comment}"
    cmd: list[str] = [
        "cutadapt",
        "-u",
        f"-{umi_length}",
        rename_template,
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

    is_dual_end_umi = (
        resolved.umi_position == "both"
        and resolved.umi_length_5p > 0
        and resolved.umi_length_3p > 0
    )
    is_single_3p_umi = resolved.umi_length > 0 and resolved.umi_position == "3p"
    two_pass = is_single_3p_umi or is_dual_end_umi

    if two_pass:
        intermediate = fastq_out.with_suffix(fastq_out.suffix + ".prepass.fq.gz")
        pass1_log = log_json.with_name(log_json.stem + ".pass1.json")
        pass2_log = log_json.with_name(log_json.stem + ".pass2.json")

        if is_dual_end_umi:
            # Pass 1: extract the 5' UMI, trim adapter; pass 2 handles
            # the 3' UMI. The pass-1 length filter must leave room for
            # the 3' UMI bases that pass 2 will subsequently strip.
            pass2_umi_length = resolved.umi_length_3p
            pass1_resolved = resolved
        else:
            # Single-end 3' UMI: pass 1 skips UMI handling, pass 2 trims it.
            pass2_umi_length = resolved.umi_length
            pass1_resolved = replace(
                resolved, umi_length=0, umi_position="5p"
            )
        pass1_cmd = _build_pass1_command(
            fastq_in=fastq_in,
            fastq_out=intermediate,
            resolved=pass1_resolved,
            min_length=min_length + pass2_umi_length,
            max_length=max_length + pass2_umi_length,
            quality=quality,
            threads=threads,
            log_json=pass1_log,
        )
        _invoke(runner, pass1_cmd)

        pass2_cmd = _build_pass2_umi_command(
            fastq_in=intermediate,
            fastq_out=fastq_out,
            umi_length=pass2_umi_length,
            threads=threads,
            log_json=pass2_log,
            append_to_existing_umi=is_dual_end_umi,
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
