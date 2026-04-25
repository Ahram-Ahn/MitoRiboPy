"""Per-sample kit + dedup resolution for ``mitoribopy align``.

Pre-v0.4.0 the pipeline scanned only the first FASTQ and applied the
detected (or user-supplied) kit to every other sample. In real datasets
samples can come from different libraries (mixed UMI / no-UMI, mixed
kits) and need independent resolution. This module owns that work:

* one :class:`SampleResolution` per input FASTQ,
* a single pre-flight pass that fails fast before any cutadapt or
  bowtie2 subprocess starts,
* a TSV writer (:func:`write_kit_resolution_tsv`) that turns the
  resolution table into the per-run provenance file users see in the
  output directory.

Resolution rules per sample:

1. Detection mode ``off`` → trust the user's ``--kit-preset`` /
   ``--adapter`` verbatim (raises if ``custom`` without ``--adapter``).
2. Detection mode ``auto`` (default) →
   - scan the head of the FASTQ;
   - on success: use the detected preset (warn if user explicitly chose
     a different preset);
   - on failure: fall back to the user's preset / adapter when one is
     supplied; otherwise raise :class:`SampleResolutionError`.
3. Detection mode ``strict`` →
   - scan;
   - on success: use detected; HARD-FAIL if user supplied a conflicting
     preset;
   - on failure: HARD-FAIL.

The resolved dedup strategy follows the resolved kit's ``umi_length``
via :func:`mitoribopy.align.dedup.resolve_dedup_strategy`, so a UMI
sample naturally picks ``umi-tools`` and a no-UMI sample naturally
picks ``skip`` even when both are processed in the same run.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

from . import adapter_detect
from ._types import (
    KIT_PRESETS,
    DedupStrategy,
    ResolvedKit,
    SampleOverride,
    resolve_kit_alias,
)
from .dedup import resolve_dedup_strategy
from .trim import resolve_kit_settings


class SampleResolutionError(ValueError):
    """A sample could not be resolved to a usable kit / dedup combination."""


@dataclass(frozen=True)
class SampleResolution:
    """Frozen per-sample configuration produced by the pre-flight pass."""

    sample: str
    fastq: Path
    kit: ResolvedKit
    dedup_strategy: DedupStrategy
    detected_kit: str | None
    detection_match_rate: float
    detection_ambiguous: bool
    source: str  # "detected", "user_fallback", "explicit_off", "explicit_strict"


# The column order is what gets written to ``kit_resolution.tsv``. Keep
# it stable so downstream tooling can rely on positional access.
_KIT_RESOLUTION_COLUMNS: tuple[str, ...] = (
    "sample",
    "fastq",
    "applied_kit",
    "adapter",
    "umi_length",
    "umi_position",
    "dedup_strategy",
    "detected_kit",
    "detection_match_rate",
    "detection_ambiguous",
    "source",
)


def _label_source(base_source: str, override_applied: bool) -> str:
    """Compose the ``source`` column value for ``kit_resolution.tsv``.

    Keeps the original detection-path label (``detected``,
    ``user_fallback``, ``inferred_pretrimmed``, ``explicit_off``) so the
    user can still see *how* the resolver landed on the resolved kit, but
    prefixes ``per_sample_override:`` when an ``align.samples:`` /
    ``--sample-overrides`` entry was applied so the provenance file
    tells the truth about the user-supplied per-sample input.
    """
    if override_applied:
        return f"per_sample_override:{base_source}"
    return base_source


def _user_supplied_explicit_preset(kit_preset: str, adapter: str | None) -> bool:
    """True when the user did NOT leave the kit at the auto sentinel.

    Aliases for canonical presets count as explicit choices once
    resolved; ``resolve_kit_alias`` is run before this helper.
    """
    return kit_preset != "auto" or adapter is not None


def _resolve_one(
    *,
    sample: str,
    fastq: Path,
    kit_preset: str,
    adapter: str | None,
    umi_length: int | None,
    umi_position: str | None,
    dedup_strategy: DedupStrategy,
    confirm_mark_duplicates: bool,
    adapter_detection_mode: str,
    allow_pretrimmed_inference: bool = True,
    detector=adapter_detect.detect_adapter,
    override: SampleOverride | None = None,
) -> SampleResolution:
    """Resolve one sample. May raise :class:`SampleResolutionError`.

    When *override* is supplied, its non-``None`` fields swap in for the
    matching incoming kwargs BEFORE the existing ``auto`` / ``off`` /
    ``strict`` precedence runs. This makes a per-sample override look
    identical to a user-supplied explicit preset on the CLI, so all the
    existing detection / fallback / validation rules apply unchanged —
    the only visible difference is the ``source`` label on the resulting
    :class:`SampleResolution` (set to ``per_sample_override``).
    """
    override_applied = False
    if override is not None:
        if override.kit_preset is not None:
            kit_preset = resolve_kit_alias(override.kit_preset)
            override_applied = True
        if override.adapter is not None:
            adapter = override.adapter
            override_applied = True
        if override.umi_length is not None:
            umi_length = override.umi_length
            override_applied = True
        if override.umi_position is not None:
            umi_position = override.umi_position
            override_applied = True
        if override.dedup_strategy is not None:
            dedup_strategy = override.dedup_strategy
            override_applied = True

    user_explicit = _user_supplied_explicit_preset(kit_preset, adapter)

    if adapter_detection_mode == "off":
        if kit_preset == "auto":
            raise SampleResolutionError(
                f"{sample}: --adapter-detection off requires an explicit "
                "--kit-preset (or --adapter) because no auto-scan will run."
            )
        kit = resolve_kit_settings(
            kit_preset,
            adapter=adapter,
            umi_length=umi_length,
            umi_position=umi_position,
        )
        dedup = resolve_dedup_strategy(
            dedup_strategy,
            umi_length=kit.umi_length,
            confirm_mark_duplicates=confirm_mark_duplicates,
        )
        return SampleResolution(
            sample=sample,
            fastq=fastq,
            kit=kit,
            dedup_strategy=dedup,
            detected_kit=None,
            detection_match_rate=0.0,
            detection_ambiguous=False,
            source=_label_source("explicit_off", override_applied),
        )

    # auto / strict modes: scan the FASTQ.
    try:
        detection = detector(fastq)
    except OSError as exc:
        if adapter_detection_mode == "strict":
            raise SampleResolutionError(
                f"{sample}: cannot scan {fastq} for adapter detection "
                f"(strict mode): {exc}"
            ) from exc
        # Auto mode: treat the failed scan as "no detection" so the user
        # fallback path below can still produce a valid resolution.
        if not user_explicit:
            raise SampleResolutionError(
                f"{sample}: cannot scan {fastq} and no --kit-preset / "
                f"--adapter fallback was supplied: {exc}"
            ) from exc
        kit = resolve_kit_settings(
            kit_preset,
            adapter=adapter,
            umi_length=umi_length,
            umi_position=umi_position,
        )
        dedup = resolve_dedup_strategy(
            dedup_strategy,
            umi_length=kit.umi_length,
            confirm_mark_duplicates=confirm_mark_duplicates,
        )
        return SampleResolution(
            sample=sample,
            fastq=fastq,
            kit=kit,
            dedup_strategy=dedup,
            detected_kit=None,
            detection_match_rate=0.0,
            detection_ambiguous=False,
            source=_label_source("user_fallback", override_applied),
        )

    detected = detection.preset_name
    rates_str = adapter_detect.format_per_kit_rates(detection.per_kit_rates)

    if adapter_detection_mode == "strict":
        if detected is None:
            raise SampleResolutionError(
                f"{sample}: --adapter-detection strict could not identify a "
                f"known adapter ({detection.n_reads_scanned} reads scanned). "
                f"Per-kit match rates: {rates_str}."
            )
        if user_explicit and kit_preset != "auto" and detected != kit_preset:
            raise SampleResolutionError(
                f"{sample}: --adapter-detection strict: data looks like "
                f"'{detected}' ({detection.match_rate * 100:.1f}% match) but "
                f"--kit-preset is {kit_preset}. Per-kit match rates: {rates_str}."
            )
        applied_preset = detected
        source = "detected"

    else:
        # auto mode (default).
        if detected is not None:
            if user_explicit and kit_preset != "auto" and detected != kit_preset:
                # Honour the explicit user choice but record what we saw.
                applied_preset = kit_preset
                source = "user_fallback"
            else:
                applied_preset = detected
                source = "detected"
        elif user_explicit:
            applied_preset = kit_preset
            source = "user_fallback"
        elif detection.pretrimmed and allow_pretrimmed_inference:
            # Universal absence of every known adapter — best inference
            # is that the FASTQ is already trimmed (e.g. SRA-deposited).
            # cutadapt will skip the -a flag and only do length/quality.
            applied_preset = "pretrimmed"
            source = "inferred_pretrimmed"
        else:
            raise SampleResolutionError(
                f"{sample}: adapter auto-detection found no known kit "
                f"({detection.n_reads_scanned} reads scanned, "
                f"per-kit rates: {rates_str}) and no --kit-preset / "
                "--adapter fallback was supplied. Either pass --kit-preset "
                "explicitly, use --kit-preset pretrimmed for already-trimmed "
                "data, or use --adapter-detection off."
            )

    kit = resolve_kit_settings(
        applied_preset,
        adapter=adapter,
        umi_length=umi_length,
        umi_position=umi_position,
    )
    dedup = resolve_dedup_strategy(
        dedup_strategy,
        umi_length=kit.umi_length,
        confirm_mark_duplicates=confirm_mark_duplicates,
    )

    return SampleResolution(
        sample=sample,
        fastq=fastq,
        kit=kit,
        dedup_strategy=dedup,
        detected_kit=detected,
        detection_match_rate=detection.match_rate,
        detection_ambiguous=detection.ambiguous,
        source=_label_source(source, override_applied),
    )


def resolve_sample_resolutions(
    samples: Iterable[tuple[str, Path]],
    *,
    kit_preset: str,
    adapter: str | None,
    umi_length: int | None,
    umi_position: str | None,
    dedup_strategy: DedupStrategy,
    confirm_mark_duplicates: bool,
    adapter_detection_mode: str,
    allow_pretrimmed_inference: bool = True,
    detector=adapter_detect.detect_adapter,
    sample_overrides: dict[str, SampleOverride] | None = None,
) -> list[SampleResolution]:
    """Pre-flight resolution for every input sample.

    Failures raise :class:`SampleResolutionError` carrying a multi-line
    message that names every offending sample so the user can fix the
    config in one pass.

    *sample_overrides* maps the sample basename to a :class:`SampleOverride`
    whose non-``None`` fields take precedence over the corresponding
    global kwargs for that sample only. Samples without an override
    inherit every global. An override naming a sample that is not
    present in *samples* triggers a :class:`SampleResolutionError` so a
    typo'd sample name in the YAML / TSV is caught up front instead of
    silently producing the wrong resolution for every sample.
    """
    # Translate any legacy alias up front so the rest of the resolver
    # only deals with canonical preset names.
    kit_preset = resolve_kit_alias(kit_preset)
    if kit_preset not in KIT_PRESETS:
        known = ", ".join(sorted(KIT_PRESETS))
        raise SampleResolutionError(
            f"Unknown --kit-preset {kit_preset!r}. Known: {known}."
        )

    samples_list = list(samples)
    sample_overrides = dict(sample_overrides or {})
    sample_names_seen = {name for name, _ in samples_list}
    unknown_overrides = sorted(
        name for name in sample_overrides if name not in sample_names_seen
    )
    if unknown_overrides:
        raise SampleResolutionError(
            "Sample override(s) reference unknown sample name(s): "
            + ", ".join(unknown_overrides)
            + ". Known sample names (FASTQ basenames): "
            + ", ".join(sorted(sample_names_seen))
            + "."
        )

    resolutions: list[SampleResolution] = []
    errors: list[str] = []
    for sample, fastq in samples_list:
        try:
            resolutions.append(
                _resolve_one(
                    sample=sample,
                    fastq=fastq,
                    kit_preset=kit_preset,
                    adapter=adapter,
                    umi_length=umi_length,
                    umi_position=umi_position,
                    dedup_strategy=dedup_strategy,
                    confirm_mark_duplicates=confirm_mark_duplicates,
                    adapter_detection_mode=adapter_detection_mode,
                    allow_pretrimmed_inference=allow_pretrimmed_inference,
                    detector=detector,
                    override=sample_overrides.get(sample),
                )
            )
        except (SampleResolutionError, ValueError, KeyError) as exc:
            errors.append(str(exc))

    if errors:
        bullet = "\n  - ".join(errors)
        raise SampleResolutionError(
            "Could not resolve every sample's library configuration:\n  - "
            + bullet
        )
    return resolutions


# ---------------------------------------------------------------------------
# Sample-override TSV I/O
# ---------------------------------------------------------------------------


_SAMPLE_OVERRIDE_COLUMNS: tuple[str, ...] = (
    "sample",
    "kit_preset",
    "adapter",
    "umi_length",
    "umi_position",
    "dedup_strategy",
)


def read_sample_overrides_tsv(path: Path) -> dict[str, SampleOverride]:
    """Parse a per-sample override TSV into a dict keyed by sample name.

    The TSV must have a header row containing at least ``sample`` and at
    least one of the override columns. Empty cells (``""``, ``"NA"``,
    ``"None"``) are read as ``None`` so the resolver falls back to the
    global default for that field.

    Raises
    ------
    SampleResolutionError
        if the file is missing required columns, has a duplicate sample
        name, or contains a malformed integer ``umi_length``.
    """
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise SampleResolutionError(
            f"--sample-overrides: cannot read {path}: {exc}"
        ) from exc

    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        return {}

    header = [col.strip() for col in lines[0].split("\t")]
    if "sample" not in header:
        raise SampleResolutionError(
            f"--sample-overrides {path}: TSV is missing a 'sample' column. "
            f"Found columns: {header}."
        )
    overridable = set(_SAMPLE_OVERRIDE_COLUMNS) - {"sample"}
    if not overridable.intersection(header):
        raise SampleResolutionError(
            f"--sample-overrides {path}: TSV must contain at least one of "
            f"{sorted(overridable)} alongside 'sample'."
        )

    overrides: dict[str, SampleOverride] = {}
    for line_no, row_text in enumerate(lines[1:], start=2):
        cells = [cell.strip() for cell in row_text.split("\t")]
        if len(cells) < len(header):
            cells.extend([""] * (len(header) - len(cells)))
        row = dict(zip(header, cells))
        sample = row.get("sample", "").strip()
        if not sample:
            raise SampleResolutionError(
                f"--sample-overrides {path}: row {line_no} has an empty "
                "'sample' field."
            )
        if sample in overrides:
            raise SampleResolutionError(
                f"--sample-overrides {path}: duplicate sample name "
                f"{sample!r} on row {line_no}."
            )

        umi_length_raw = row.get("umi_length", "")
        umi_length: int | None
        if _is_blank(umi_length_raw):
            umi_length = None
        else:
            try:
                umi_length = int(umi_length_raw)
            except ValueError as exc:
                raise SampleResolutionError(
                    f"--sample-overrides {path}: row {line_no} has a "
                    f"non-integer umi_length {umi_length_raw!r}."
                ) from exc

        umi_position_raw = row.get("umi_position", "")
        umi_position = None if _is_blank(umi_position_raw) else umi_position_raw
        if umi_position is not None and umi_position not in ("5p", "3p"):
            raise SampleResolutionError(
                f"--sample-overrides {path}: row {line_no} has an invalid "
                f"umi_position {umi_position!r} (must be '5p' or '3p')."
            )

        kit_preset_raw = row.get("kit_preset", "")
        adapter_raw = row.get("adapter", "")
        dedup_raw = row.get("dedup_strategy", "")

        overrides[sample] = SampleOverride(
            sample=sample,
            kit_preset=None if _is_blank(kit_preset_raw) else kit_preset_raw,
            adapter=None if _is_blank(adapter_raw) else adapter_raw,
            umi_length=umi_length,
            umi_position=umi_position,  # type: ignore[arg-type]
            dedup_strategy=None if _is_blank(dedup_raw) else dedup_raw,  # type: ignore[arg-type]
        )

    return overrides


def write_sample_overrides_tsv(
    overrides: Iterable[SampleOverride], path: Path
) -> Path:
    """Write a per-sample override TSV (round-trips with
    :func:`read_sample_overrides_tsv`).

    Used by the ``mitoribopy all`` orchestrator to materialise an
    ``align.samples:`` YAML block as a file the align CLI can consume
    via ``--sample-overrides``.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(_SAMPLE_OVERRIDE_COLUMNS) + "\n")
        for override in overrides:
            row = {
                "sample": override.sample,
                "kit_preset": override.kit_preset or "",
                "adapter": override.adapter or "",
                "umi_length": (
                    "" if override.umi_length is None else str(override.umi_length)
                ),
                "umi_position": override.umi_position or "",
                "dedup_strategy": override.dedup_strategy or "",
            }
            handle.write(
                "\t".join(row[col] for col in _SAMPLE_OVERRIDE_COLUMNS) + "\n"
            )
    return path


def _is_blank(value: str) -> bool:
    return value.strip().lower() in {"", "na", "none", "null"}


def required_dedup_tools(resolutions: list[SampleResolution]) -> set[str]:
    """Union of external tools required by the resolved dedup strategies."""
    tools: set[str] = set()
    for resolution in resolutions:
        if resolution.dedup_strategy == "umi-tools":
            tools.add("umi_tools")
        elif resolution.dedup_strategy == "mark-duplicates":
            tools.add("picard")
    return tools


def write_kit_resolution_tsv(
    resolutions: list[SampleResolution],
    output_path: Path,
) -> Path:
    """Write the per-sample resolution table.

    Columns: see :data:`_KIT_RESOLUTION_COLUMNS`. One row per sample.
    The file is the spine the README troubleshooting section points to
    when explaining why a particular sample picked the kit it did.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(_KIT_RESOLUTION_COLUMNS) + "\n")
        for resolution in resolutions:
            row = {
                "sample": resolution.sample,
                "fastq": str(resolution.fastq),
                "applied_kit": resolution.kit.kit,
                "adapter": resolution.kit.adapter,
                "umi_length": resolution.kit.umi_length,
                "umi_position": resolution.kit.umi_position,
                "dedup_strategy": resolution.dedup_strategy,
                "detected_kit": resolution.detected_kit or "",
                "detection_match_rate": f"{resolution.detection_match_rate:.4f}",
                "detection_ambiguous": "true" if resolution.detection_ambiguous else "false",
                "source": resolution.source,
            }
            handle.write(
                "\t".join(str(row[col]) for col in _KIT_RESOLUTION_COLUMNS) + "\n"
            )
    return output_path


def resolution_summary_lines(resolutions: list[SampleResolution]) -> list[str]:
    """Human-readable per-sample summary used in dry-run plans."""
    lines: list[str] = []
    for resolution in resolutions:
        line = (
            f"{resolution.sample}: kit={resolution.kit.kit}, "
            f"adapter={resolution.kit.adapter}, "
            f"umi_length={resolution.kit.umi_length} ({resolution.kit.umi_position}), "
            f"dedup={resolution.dedup_strategy}, source={resolution.source}"
        )
        if resolution.detected_kit and resolution.detected_kit != resolution.kit.kit:
            line += (
                f", detected={resolution.detected_kit}@"
                f"{resolution.detection_match_rate * 100:.1f}%"
            )
        lines.append(line)
    return lines
