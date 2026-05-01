"""Verify that external bioinformatics tools are available on PATH.

Every ``mitoribopy align`` invocation calls :func:`ensure_tools_available`
before touching user data. If any required tool is missing, we raise
:class:`ToolNotFoundError` with an actionable bioconda install command.
Optional tools (``fastqc``, ``umi_tools`` in no-UMI runs) are reported
but do not block a run.

``mitoribopy align`` is **fail-loud** by design: never silently skip,
never replace a missing tool with a default no-op, never continue with
partial data. A missing required tool aborts the run with a bioconda
install hint.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from typing import Iterable

from ._types import Strandedness, ToolInfo


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------


class ToolNotFoundError(RuntimeError):
    """Raised when a required external tool cannot be located on PATH."""


# ---------------------------------------------------------------------------
# Tool registry
# ---------------------------------------------------------------------------


# Command-line idiom used to read each tool's version. We pick a form that
# prints to stdout on success so subprocess capture is consistent.
_VERSION_CMDS: dict[str, list[str]] = {
    "cutadapt": ["cutadapt", "--version"],
    "bowtie2": ["bowtie2", "--version"],
    "bowtie2-build": ["bowtie2-build", "--version"],
    "samtools": ["samtools", "--version"],
    "umi_tools": ["umi_tools", "--version"],
    "bedtools": ["bedtools", "--version"],
    "fastqc": ["fastqc", "--version"],
}


# Bioconda install hint shown in the error message. Keeping the exact form
# in one place so the user sees the same line whether the failure is
# cutadapt or bowtie2.
_BIOCONDA_HINT = (
    "Install the required tools via bioconda:\n"
    "  conda install -c bioconda -c conda-forge cutadapt bowtie2 samtools "
    "umi_tools bedtools fastqc\n"
    "Or use the provided environment.yml under docs/environment/."
)


def _parse_version(tool_name: str, raw: str) -> str:
    """Extract a readable version string from the tool's stdout."""
    text = raw.strip().splitlines()[0] if raw.strip() else ""
    # Common shapes:
    #   "cutadapt 4.9"
    #   "/usr/local/bin/bowtie2-align-s version 2.5.4"
    #   "samtools 1.21"
    #   "UMI-tools version: 1.1.5"
    match = re.search(r"(\d+\.\d+(?:\.\d+)?)", text)
    if match:
        return match.group(1)
    # Fall back to the raw first line so we never claim "unknown" silently.
    return text or "unknown"


# ---------------------------------------------------------------------------
# Single-tool checks
# ---------------------------------------------------------------------------


def get_tool_version(tool_name: str) -> str:
    """Return a best-effort version string for *tool_name*.

    Returns ``"unknown"`` if the tool is on PATH but the ``--version``
    probe fails (rare; can happen on partially-installed packages).
    Raises :class:`ToolNotFoundError` if the tool is not on PATH.
    """
    path = shutil.which(tool_name)
    if path is None:
        raise ToolNotFoundError(
            f"Tool '{tool_name}' not found on PATH.\n{_BIOCONDA_HINT}"
        )

    cmd = _VERSION_CMDS.get(tool_name, [tool_name, "--version"])
    try:
        completed = subprocess.run(
            cmd,
            check=False,
            capture_output=True,
            text=True,
            timeout=15,
        )
    except (subprocess.SubprocessError, OSError):
        return "unknown"

    output = completed.stdout or completed.stderr or ""
    return _parse_version(tool_name, output)


def check_tool(tool_name: str) -> ToolInfo:
    """Return a :class:`ToolInfo` for *tool_name*; raise if not on PATH."""
    path = shutil.which(tool_name)
    if path is None:
        raise ToolNotFoundError(
            f"Required tool '{tool_name}' not found on PATH.\n{_BIOCONDA_HINT}"
        )
    return ToolInfo(name=tool_name, path=path, version=get_tool_version(tool_name))


# ---------------------------------------------------------------------------
# Batch check for an align run
# ---------------------------------------------------------------------------


def ensure_tools_available(
    required: Iterable[str],
    optional: Iterable[str] = (),
) -> dict[str, ToolInfo | None]:
    """Verify a batch of tools; raise if any required one is missing.

    Parameters
    ----------
    required:
        Tools that MUST be on PATH. A single missing tool raises
        :class:`ToolNotFoundError` with the complete missing-list so the
        user sees every problem at once rather than one per retry.
    optional:
        Tools that are nice-to-have (``fastqc``). Missing entries are
        returned as ``None`` in the result dict; no raise.

    Returns
    -------
    dict
        Mapping from tool name to :class:`ToolInfo` when present, or
        ``None`` when an *optional* tool is missing.
    """
    required_names = list(required)
    optional_names = list(optional)

    missing_required: list[str] = [
        name for name in required_names if shutil.which(name) is None
    ]
    if missing_required:
        listing = ", ".join(missing_required)
        raise ToolNotFoundError(
            "Required external tool(s) not found on PATH: "
            f"{listing}.\n{_BIOCONDA_HINT}"
        )

    resolved: dict[str, ToolInfo | None] = {}
    for name in required_names:
        resolved[name] = check_tool(name)
    for name in optional_names:
        if shutil.which(name) is None:
            resolved[name] = None
        else:
            resolved[name] = check_tool(name)
    return resolved


# ---------------------------------------------------------------------------
# Strandedness -> bowtie2 flag mapping
# ---------------------------------------------------------------------------


def bowtie2_strand_flag(strandedness: Strandedness) -> list[str]:
    """Return the bowtie2 flag list for a given library strandedness.

    This is surfaced here (rather than inside ``align.py``) so the
    strand-to-aligner mapping is auditable next to the other external-tool
    glue. Tests assert this mapping explicitly because a silently-wrong
    strand flag leaves the ND5 / ND6 antisense overlap ambiguous on
    genome references and produces transcript-contamination on
    transcriptome references.
    """
    if strandedness == "forward":
        return ["--norc"]
    if strandedness == "reverse":
        return ["--nofw"]
    if strandedness == "unstranded":
        return []
    raise ValueError(  # pragma: no cover - Strandedness literal guards this
        f"Unknown strandedness: {strandedness!r}"
    )
