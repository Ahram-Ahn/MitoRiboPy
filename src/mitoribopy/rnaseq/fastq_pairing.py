"""FASTQ enumeration and SE / PE pairing for the from-FASTQ rnaseq mode.

Users can pass a mix of files and directories via ``--rna-fastq`` /
``--ribo-fastq``. We expand directories to FASTQ files and group files
into samples by stripping a mate-pair token from the filename. The
supported tokens cover every encoding we have seen in real submissions:

* bcl2fastq:        ``_R1_001`` / ``_R2_001`` (e.g. ``A_S1_L001_R1_001.fastq.gz``)
* read1/read2:      ``_read1`` / ``_read2`` (case-insensitive)
* R1/R2:            ``_R1`` / ``_R2``
* dot-numbered:     ``.1.`` / ``.2.`` (interior, e.g. ``A.1.fastq.gz``)
* underscore-num:   ``_1`` / ``_2`` (suffix immediately before extension)

Tokens are matched longest-first so ``_R1_001`` does not get partially
stripped to ``_R1`` and lose the ``_001`` lane suffix.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path


_FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


@dataclass(frozen=True)
class FastqSample:
    """A single sample, single-end (``r2 is None``) or paired-end."""

    sample: str
    r1: Path
    r2: Path | None = None

    @property
    def paired(self) -> bool:
        return self.r2 is not None


# Patterns are evaluated in order. Each entry is
# (regex, mate_group, stem_replacement) where:
#   - regex matches a substring of the basename-without-extension
#   - mate_group captures the mate marker ('1' or '2')
#   - stem_replacement is what we put back in to get the canonical stem
# For tokens that are pure suffixes the replacement is the empty string;
# for the dot-numbered interior token we need to keep the surrounding
# context, hence the explicit ``"."``.
_PATTERNS: tuple[tuple[re.Pattern[str], str], ...] = (
    # bcl2fastq: ..._R1_001 or ..._R2_001 at end of stem
    (re.compile(r"_R([12])_001$", re.IGNORECASE), ""),
    # read1 / read2 at end of stem
    (re.compile(r"_read([12])$", re.IGNORECASE), ""),
    # R1 / R2 at end of stem
    (re.compile(r"_R([12])$", re.IGNORECASE), ""),
    # dot-numbered .1 / .2 at end of stem (originally interior in the
    # full filename — the trailing dot belonged to the .fq / .fastq
    # extension and is gone by the time we classify).
    (re.compile(r"\.([12])$"), ""),
    # _1 / _2 at end of stem
    (re.compile(r"_([12])$"), ""),
)


def _strip_fastq_suffix(name: str) -> str:
    """Return the basename with the longest known FASTQ extension stripped."""
    lower = name.lower()
    for suf in _FASTQ_SUFFIXES:
        if lower.endswith(suf):
            return name[: -len(suf)]
    return name


def _is_fastq(name: str) -> bool:
    return name.lower().endswith(_FASTQ_SUFFIXES)


def enumerate_fastqs(inputs: Iterable[Path | str]) -> list[Path]:
    """Expand a mixed list of files and directories into FASTQ paths.

    For each input that is a directory, every file directly inside it
    whose name ends in one of the supported FASTQ suffixes is collected
    (non-recursive — sequencing centres typically deliver a flat layout
    and recursing would surprise users with a deep ``work/`` tree).

    Unknown extensions (e.g. ``.bam``, ``.txt``) raise ``ValueError`` to
    surface mistakes early; silently ignoring unknown files is exactly
    how real runs end up missing a sample.

    Output is de-duplicated by absolute resolved path and stable-sorted.
    """
    out: list[Path] = []
    seen: set[Path] = set()

    def _add(path: Path) -> None:
        resolved = path.resolve()
        if resolved in seen:
            return
        seen.add(resolved)
        out.append(path)

    for raw in inputs:
        path = Path(raw)
        if path.is_dir():
            for child in sorted(path.iterdir()):
                if not child.is_file():
                    continue
                if not _is_fastq(child.name):
                    # A directory may legitimately contain unrelated files
                    # (READMEs, MD5 sums, etc.). Ignore them rather than
                    # rejecting — only top-level user inputs are validated
                    # for extension.
                    continue
                _add(child)
        elif path.is_file():
            if not _is_fastq(path.name):
                raise ValueError(
                    f"Unsupported FASTQ extension: {path}. "
                    f"Expected one of {', '.join(_FASTQ_SUFFIXES)}."
                )
            _add(path)
        else:
            raise FileNotFoundError(f"Input path does not exist: {path}")

    out.sort(key=lambda p: str(p))
    return out


def _classify(stem: str) -> tuple[str, str | None]:
    """Return ``(canonical_stem, mate)`` where mate is ``'1'``, ``'2'`` or ``None``."""
    for pattern, replacement in _PATTERNS:
        m = pattern.search(stem)
        if m:
            mate = m.group(1)
            canonical = stem[: m.start()] + replacement + stem[m.end():]
            # Trim any trailing underscore the strip may have left behind.
            canonical = canonical.rstrip("_.")
            return canonical, mate
    return stem, None


def detect_samples(paths: Iterable[Path]) -> list[FastqSample]:
    """Group ``paths`` into :class:`FastqSample` records.

    Pairing rules:

    * Both ``_R1`` and ``_R2`` for a stem present → paired sample.
    * Lone ``_R1`` → SE sample with ``sample = stem``.
    * Lone ``_R2`` → SE sample with ``sample = f"{stem}_R2"`` (so the
      user can see at a glance that R1 was missing).
    * Filename with no recognised mate token → SE sample, name = stem.
    * Two files producing the same ``(stem, mate)`` pair raise
      :class:`ValueError` listing every duplicated path.
    """
    by_stem_mate: dict[tuple[str, str | None], list[Path]] = {}
    no_mate: list[tuple[str, Path]] = []

    for raw in paths:
        path = Path(raw)
        stem = _strip_fastq_suffix(path.name)
        canonical, mate = _classify(stem)
        if mate is None:
            no_mate.append((canonical, path))
        else:
            by_stem_mate.setdefault((canonical, mate), []).append(path)

    duplicates: list[str] = []
    for (stem, mate), files in by_stem_mate.items():
        if len(files) > 1:
            duplicates.append(
                f"{stem!r} mate={mate}: " + ", ".join(str(p) for p in files)
            )
    if duplicates:
        raise ValueError(
            "Duplicate (sample, mate) entries detected:\n  "
            + "\n  ".join(duplicates)
        )

    samples: list[FastqSample] = []
    seen_stems: set[str] = set()
    # Iterate stems in stable sorted order so the output is reproducible.
    stems = sorted({stem for stem, _ in by_stem_mate.keys()})
    for stem in stems:
        r1_paths = by_stem_mate.get((stem, "1"), [])
        r2_paths = by_stem_mate.get((stem, "2"), [])
        r1 = r1_paths[0] if r1_paths else None
        r2 = r2_paths[0] if r2_paths else None
        if r1 and r2:
            samples.append(FastqSample(sample=stem, r1=r1, r2=r2))
        elif r1:
            samples.append(FastqSample(sample=stem, r1=r1, r2=None))
        elif r2:
            samples.append(FastqSample(sample=f"{stem}_R2", r1=r2, r2=None))
        seen_stems.add(stem)

    for canonical, path in sorted(no_mate, key=lambda x: str(x[1])):
        samples.append(FastqSample(sample=canonical, r1=path, r2=None))

    return samples
