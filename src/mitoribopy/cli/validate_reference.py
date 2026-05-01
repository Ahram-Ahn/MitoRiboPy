"""``mitoribopy validate-reference`` — pre-flight check for custom references.

Usage::

    mitoribopy validate-reference \\
        --fasta references/custom_mt.fa \\
        --annotation references/custom_mt.csv \\
        [--codon-table vertebrate_mitochondrial]

Verifies that a custom mitochondrial reference is internally
consistent BEFORE running ``mitoribopy rpf`` against it. The checks
mirror the failure modes that have actually bitten users:

* every transcript in the annotation has a FASTA record;
* every FASTA record is referenced by exactly one annotation row;
* annotation transcript lengths (``l_tr``) match the FASTA sequence
  length;
* the implied CDS length (``l_tr - l_utr5 - l_utr3``) is positive and
  divisible by 3;
* the CDS-encoded start codon is a valid start codon under the
  selected codon table (typically ``ATG``, plus the alternative
  starts allowed by NCBI translation tables — e.g. vertebrate mt
  uses ``ATG``, ``ATA``, ``ATT``, ``GTG``);
* the CDS-encoded stop codon is a valid stop codon under the same
  table.

Exit codes:
* ``0`` — clean (all checks passed)
* ``2`` — at least one rejection
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from . import common


VALIDATE_REFERENCE_HELP = (
    "Pre-flight a custom mitochondrial reference: check that the FASTA "
    "and annotation CSV are consistent (matching transcript IDs, "
    "matching lengths, CDS divisible by 3, valid start / stop codons "
    "under the selected codon table)."
)


# Conservative defaults. The vertebrate-mitochondrial code (NCBI table
# 2) is the most common case; we hard-code its start/stop codon set
# inline so the validator does not depend on a richer codon-table
# loader. ``--codon-table standard`` falls back to NCBI table 1.
_CODON_TABLE_STARTS: dict[str, set[str]] = {
    "vertebrate_mitochondrial": {"ATG", "ATA", "ATT", "GTG"},
    "yeast_mitochondrial": {"ATG", "ATA"},
    "standard": {"ATG"},
}
_CODON_TABLE_STOPS: dict[str, set[str]] = {
    "vertebrate_mitochondrial": {"TAA", "TAG", "AGA", "AGG"},
    "yeast_mitochondrial": {"TAA", "TAG"},
    "standard": {"TAA", "TAG", "TGA"},
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mitoribopy validate-reference",
        description=VALIDATE_REFERENCE_HELP,
        formatter_class=common.MitoRiboPyHelpFormatter,
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Path to the custom mt-transcriptome FASTA.",
    )
    parser.add_argument(
        "--annotation",
        required=True,
        help=(
            "Path to the per-transcript annotation CSV (the same "
            "schema accepted by `mitoribopy rpf --annotation_file`)."
        ),
    )
    parser.add_argument(
        "--codon-table",
        default="vertebrate_mitochondrial",
        choices=sorted(_CODON_TABLE_STARTS.keys()),
        help=(
            "Codon table for start / stop codon validation. Default "
            "is vertebrate_mitochondrial (NCBI translation table 2)."
        ),
    )
    return parser


def _parse_fasta(path: Path) -> dict[str, str]:
    """Return ``{header: sequence}`` from a plain FASTA.

    Header is taken from the first whitespace-separated token after
    ``>``. Sequence is uppercased and stripped of newlines. The reader
    is intentionally minimal — full validation of the FASTA itself is
    out of scope (the rpf stage will exercise more of the file).
    """
    records: dict[str, list[str]] = {}
    current: str | None = None
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = []
            elif current is not None:
                records[current].append(line.upper())
    return {name: "".join(parts) for name, parts in records.items()}


def _parse_annotation(path: Path) -> list[dict[str, str]]:
    """Return rows of the annotation CSV as ``{column: cell}`` dicts."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return [dict(row) for row in reader]


def _coerce_int(value: str | None, *, field: str, transcript: str) -> int:
    if value is None or value == "":
        raise ValueError(
            f"transcript {transcript!r}: missing required {field!r} column"
        )
    try:
        return int(value)
    except ValueError as exc:
        raise ValueError(
            f"transcript {transcript!r}: {field!r}={value!r} is not an integer"
        ) from exc


def validate_reference(
    *,
    fasta_path: Path,
    annotation_path: Path,
    codon_table: str,
) -> tuple[list[str], list[str]]:
    """Return ``(errors, warnings)`` for the given reference pair.

    ``errors`` block the run; ``warnings`` are surfaced but do not
    change the exit code.
    """
    errors: list[str] = []
    warnings: list[str] = []

    if not fasta_path.exists():
        errors.append(f"FASTA not found: {fasta_path}")
        return errors, warnings
    if not annotation_path.exists():
        errors.append(f"annotation CSV not found: {annotation_path}")
        return errors, warnings

    fasta = _parse_fasta(fasta_path)
    rows = _parse_annotation(annotation_path)
    if not rows:
        errors.append(f"annotation CSV is empty: {annotation_path}")
        return errors, warnings

    starts = _CODON_TABLE_STARTS[codon_table]
    stops = _CODON_TABLE_STOPS[codon_table]

    seen_seq_names: set[str] = set()

    for row in rows:
        transcript = row.get("transcript") or row.get("sequence_name") or ""
        if not transcript:
            errors.append("annotation row is missing 'transcript' / 'sequence_name'")
            continue

        seq_name = row.get("sequence_name") or transcript
        seen_seq_names.add(seq_name)

        # 1. FASTA header coverage.
        if seq_name not in fasta:
            errors.append(
                f"transcript {transcript!r}: sequence_name={seq_name!r} not "
                f"present in FASTA {fasta_path.name}"
            )
            continue

        seq = fasta[seq_name]

        # 2. Length match.
        try:
            l_tr = _coerce_int(row.get("l_tr"), field="l_tr", transcript=transcript)
            l_utr5 = _coerce_int(
                row.get("l_utr5"), field="l_utr5", transcript=transcript
            )
            l_utr3 = _coerce_int(
                row.get("l_utr3"), field="l_utr3", transcript=transcript
            )
        except ValueError as exc:
            errors.append(str(exc))
            continue

        if len(seq) != l_tr:
            errors.append(
                f"transcript {transcript!r}: l_tr={l_tr} but FASTA "
                f"sequence is {len(seq)} nt long"
            )
            continue

        # 3. CDS length and divisibility.
        l_cds = l_tr - l_utr5 - l_utr3
        if l_cds <= 0:
            errors.append(
                f"transcript {transcript!r}: implied CDS length "
                f"({l_cds}) is non-positive (l_tr={l_tr}, l_utr5="
                f"{l_utr5}, l_utr3={l_utr3})"
            )
            continue
        if l_cds % 3 != 0:
            errors.append(
                f"transcript {transcript!r}: CDS length {l_cds} is "
                "not divisible by 3"
            )
            # Continue so we can still flag start/stop codon issues
            # if they are also broken.

        # 4. Start codon.
        cds_start = seq[l_utr5 : l_utr5 + 3]
        if len(cds_start) == 3 and cds_start not in starts:
            errors.append(
                f"transcript {transcript!r}: CDS start codon "
                f"{cds_start!r} is not a valid start under "
                f"codon_table={codon_table!r} (allowed: "
                f"{sorted(starts)})"
            )

        # 5. Stop codon.
        if l_cds % 3 == 0:
            cds_stop = seq[l_utr5 + l_cds - 3 : l_utr5 + l_cds]
            if len(cds_stop) == 3 and cds_stop not in stops:
                errors.append(
                    f"transcript {transcript!r}: CDS stop codon "
                    f"{cds_stop!r} is not a valid stop under "
                    f"codon_table={codon_table!r} (allowed: "
                    f"{sorted(stops)})"
                )

    # 6. Orphan FASTA records (warnings, not errors — extra records
    # don't break a run, but the user typically wants to know).
    orphan = sorted(set(fasta.keys()) - seen_seq_names)
    if orphan:
        warnings.append(
            f"FASTA records present but not referenced by annotation: "
            + ", ".join(orphan)
        )

    return errors, warnings


def run(argv) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv))

    errors, warnings = validate_reference(
        fasta_path=Path(args.fasta),
        annotation_path=Path(args.annotation),
        codon_table=args.codon_table,
    )

    for w in warnings:
        sys.stderr.write(f"[mitoribopy validate-reference] WARNING: {w}\n")
    if errors:
        for e in errors:
            sys.stderr.write(f"[mitoribopy validate-reference] ERROR: {e}\n")
        sys.stderr.write(
            f"[mitoribopy validate-reference] {len(errors)} error(s) found.\n"
        )
        return 2

    sys.stderr.write(
        "[mitoribopy validate-reference] reference looks consistent "
        f"(codon_table={args.codon_table}).\n"
    )
    return 0
