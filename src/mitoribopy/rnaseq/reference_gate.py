"""Reference-consistency gate for ``mitoribopy rnaseq``.

Ribo-seq and RNA-seq TE comparisons are meaningful ONLY when both
datasets were aligned to the identical transcript reference. A silent
reference mismatch is the failure mode this gate exists to prevent:
it hashes the user-supplied reference (FASTA or GTF) and refuses to
proceed unless the hash matches the hash recorded by the prior
``mitoribopy rpf`` / ``mitoribopy align`` run.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


class ReferenceMismatchError(RuntimeError):
    """Raised when the rnaseq reference does not match the rpf reference."""


def compute_reference_checksum(path: Path, *, algorithm: str = "sha256") -> str:
    """Return a hex digest of the file at *path* using *algorithm*.

    ``sha256`` is the default; the hex form is written verbatim into
    the run manifest so text diffs of the manifest remain human-readable.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Reference file not found: {path}")

    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_recorded_checksum(ribo_dir: Path) -> str | None:
    """Return the reference checksum the prior rpf/align run recorded.

    Reads ``ribo_dir / 'run_settings.json'`` or
    ``ribo_dir / 'run_manifest.json'`` (the Phase-6 name); the first
    that exists wins. Missing file or missing key yields ``None``.
    """
    for name in ("run_manifest.json", "run_settings.json"):
        path = Path(ribo_dir) / name
        if path.is_file():
            try:
                data = json.loads(path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                continue
            checksum = data.get("reference_checksum") or data.get(
                "reference_sha256"
            )
            if checksum:
                return str(checksum)
    return None


def verify_reference_consistency(
    *,
    ribo_dir: Path,
    reference_path: Path | None = None,
    reference_checksum: str | None = None,
) -> str:
    """Verify the rnaseq reference matches the recorded rpf reference.

    Exactly one of *reference_path* / *reference_checksum* must be
    provided. When only a path is given we hash it; when only a
    checksum is given we use it verbatim (for users who cannot ship the
    reference file to the rnaseq host).

    Returns the verified hex digest. Raises
    :class:`ReferenceMismatchError` on mismatch, or when no recorded
    checksum is available to compare against (the rpf run must emit one;
    this is the loud failure that motivates the gate).
    """
    if not reference_path and not reference_checksum:
        raise ValueError(
            "verify_reference_consistency requires --reference-gtf or "
            "--reference-checksum."
        )
    if reference_path and reference_checksum:
        raise ValueError(
            "Pass EXACTLY ONE of --reference-gtf / --reference-checksum."
        )

    effective = (
        reference_checksum
        if reference_checksum is not None
        else compute_reference_checksum(Path(reference_path))
    )

    recorded = _read_recorded_checksum(Path(ribo_dir))
    if recorded is None:
        raise ReferenceMismatchError(
            "No 'reference_checksum' / 'reference_sha256' key found in "
            f"{ribo_dir}/run_manifest.json or run_settings.json. The prior "
            "mitoribopy rpf run must record the reference hash for rnaseq "
            "to verify consistency. Re-run rpf with a current MitoRiboPy "
            "(>=0.3.0) and ensure the manifest writer ran to completion."
        )

    if recorded.lower() != effective.lower():
        raise ReferenceMismatchError(
            "Reference MISMATCH between Ribo-seq and RNA-seq aligners.\n"
            f"  recorded (from rpf run at {ribo_dir}): {recorded}\n"
            f"  rnaseq-side reference digest:          {effective}\n"
            "TE comparisons across different references are invalid. "
            "Re-align one side so both use the identical transcript set, "
            "or verify you passed the correct --reference-gtf / "
            "--reference-checksum."
        )

    return effective
