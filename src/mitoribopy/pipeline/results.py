"""Typed stage-result dataclasses for cross-stage handoff.

The orchestrator (``mitoribopy all``) wires stages together by passing
output paths from the previous stage into the next. Historically each
stage was wired through the raw ``argparse.Namespace``, which made the
contract implicit and easy to break: dropping a field on the producing
side silently sent ``None`` into the next stage.

These frozen dataclasses are the explicit contract:

* :class:`AlignResult`  — what ``align`` produces and ``rpf`` consumes.
* :class:`RpfResult`    — what ``rpf`` produces and ``rnaseq`` consumes.
* :class:`RnaseqResult` — what ``rnaseq`` produces (terminal).

Each class is a thin path manifest. None of the dataclasses own
DataFrames or open file handles; they are intentionally JSON-friendly
so they can be embedded in ``run_manifest.json`` if needed.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path


__all__ = [
    "AlignResult",
    "RpfResult",
    "RnaseqResult",
]


@dataclass(frozen=True)
class AlignResult:
    """Artifacts emitted by the ``align`` stage.

    Attributes
    ----------
    bed_dir
        Directory containing the per-sample ``<sample>.mt.bed`` BED6
        files that ``rpf`` ingests.
    read_counts
        Path to ``read_counts.tsv`` (per-sample read funnel summary).
    kit_resolution
        Path to ``kit_resolution.tsv`` (per-sample kit / adapter / UMI
        resolution).
    run_settings
        Path to ``run_settings.json`` for the align stage.
    bam_dir
        Optional directory of post-MAPQ BAMs kept under
        ``--keep-intermediates``. ``None`` when intermediates were
        discarded.
    """

    bed_dir: Path
    read_counts: Path
    kit_resolution: Path
    run_settings: Path
    bam_dir: Path | None = None

    def as_dict(self) -> dict:
        out = asdict(self)
        for key, value in list(out.items()):
            if isinstance(value, Path):
                out[key] = str(value)
        return out


@dataclass(frozen=True)
class RpfResult:
    """Artifacts emitted by the ``rpf`` stage.

    Attributes
    ----------
    rpf_counts
        Path to ``rpf_counts.tsv`` — the mt-RPF count matrix.
    offset_diagnostics_dir
        Directory holding per-sample offset selection diagnostics.
    run_settings
        Path to ``run_settings.json`` for the rpf stage.
    reference_checksum
        SHA256 of the mt-transcriptome FASTA that the rpf stage used.
        Re-checked by downstream stages (``rnaseq``) so a swapped
        reference fails fast rather than producing silent misalignment.
    rpf_counts_metadata
        Optional sidecar JSON next to ``rpf_counts.tsv`` recording the
        filter window and offset provenance.
    """

    rpf_counts: Path
    offset_diagnostics_dir: Path
    run_settings: Path
    reference_checksum: str
    rpf_counts_metadata: Path | None = None

    def as_dict(self) -> dict:
        out = asdict(self)
        for key, value in list(out.items()):
            if isinstance(value, Path):
                out[key] = str(value)
        return out


@dataclass(frozen=True)
class RnaseqResult:
    """Artifacts emitted by the ``rnaseq`` stage.

    The full table inventory varies with the stage mode (de_table /
    from_fastq / rna_only). Fields default to ``None`` so callers can
    populate only what their mode produced.

    Attributes
    ----------
    mode
        ``"de_table" | "from_fastq" | "rna_only" | "none"``.
    te_table
        Path to ``te.tsv`` when TE was computed.
    delta_te_table
        Path to ``delta_te.tsv`` when delta-TE was computed.
    rna_counts
        Path to ``rna_counts.tsv`` (always for rna_only / from_fastq).
    de_table_canonical
        Canonicalised DE table path used for de_table mode.
    run_settings
        Path to ``run_settings.json`` for the rnaseq stage.
    exploratory_sidecar
        Path to ``EXPLORATORY.md`` when the stage ran in pseudo-replicate
        or exploratory mode.
    """

    mode: str
    run_settings: Path
    te_table: Path | None = None
    delta_te_table: Path | None = None
    rna_counts: Path | None = None
    de_table_canonical: Path | None = None
    exploratory_sidecar: Path | None = None

    def as_dict(self) -> dict:
        out = asdict(self)
        for key, value in list(out.items()):
            if isinstance(value, Path):
                out[key] = str(value)
        return out
