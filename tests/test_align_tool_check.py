"""Unit tests for ``mitoribopy.align.tool_check``.

All tests use ``monkeypatch`` to stub out :func:`shutil.which` and
:func:`subprocess.run` so they run without any bioinformatics tool
installed. The separate integration suite in
``tests/test_align_integration.py`` exercises real tools under the
``requires_tools`` mark.
"""

from __future__ import annotations

import shutil
import subprocess

import pytest

from mitoribopy.align import tool_check
from mitoribopy.align._types import Strandedness, ToolInfo
from mitoribopy.align.tool_check import (
    ToolNotFoundError,
    bowtie2_strand_flag,
    check_tool,
    ensure_tools_available,
    get_tool_version,
)


# ---------- fixtures ---------------------------------------------------------


@pytest.fixture
def fake_path(monkeypatch):
    """Return a helper that registers fake binaries on the simulated PATH."""
    present: dict[str, str] = {}

    def fake_which(name):
        return present.get(name)

    monkeypatch.setattr(tool_check.shutil, "which", fake_which)
    return present


@pytest.fixture
def fake_subprocess_run(monkeypatch):
    """Return a helper that records subprocess.run calls and returns canned outputs."""
    version_map: dict[str, str] = {}
    recorded_calls: list[list[str]] = []

    class FakeCompleted:
        def __init__(self, stdout: str, stderr: str = "", returncode: int = 0) -> None:
            self.stdout = stdout
            self.stderr = stderr
            self.returncode = returncode

    def fake_run(cmd, **kwargs):
        recorded_calls.append(list(cmd))
        tool = cmd[0]
        out = version_map.get(tool, "unknown")
        return FakeCompleted(stdout=out)

    monkeypatch.setattr(tool_check.subprocess, "run", fake_run)
    return version_map, recorded_calls


# ---------- get_tool_version -------------------------------------------------


def test_get_tool_version_parses_semver_from_cutadapt(fake_path, fake_subprocess_run) -> None:
    fake_path["cutadapt"] = "/usr/local/bin/cutadapt"
    version_map, _ = fake_subprocess_run
    version_map["cutadapt"] = "cutadapt 4.9\n"

    assert get_tool_version("cutadapt") == "4.9"


def test_get_tool_version_parses_bowtie2_header_line(fake_path, fake_subprocess_run) -> None:
    fake_path["bowtie2"] = "/usr/local/bin/bowtie2"
    version_map, _ = fake_subprocess_run
    version_map["bowtie2"] = (
        "/usr/local/bin/bowtie2-align-s version 2.5.4\n"
        "64-bit\nBuilt on host ...\n"
    )

    assert get_tool_version("bowtie2") == "2.5.4"


def test_get_tool_version_falls_back_to_first_line_when_no_number(
    fake_path, fake_subprocess_run
) -> None:
    fake_path["umi_tools"] = "/usr/local/bin/umi_tools"
    version_map, _ = fake_subprocess_run
    version_map["umi_tools"] = "umi_tools\n"

    # No digits in the stdout -> we return the raw first line rather than
    # claim "unknown", so the user still sees something traceable.
    assert get_tool_version("umi_tools") == "umi_tools"


def test_get_tool_version_raises_when_tool_absent(fake_path) -> None:
    with pytest.raises(ToolNotFoundError) as exc:
        get_tool_version("cutadapt")
    assert "cutadapt" in str(exc.value)
    assert "bioconda" in str(exc.value)


def test_get_tool_version_returns_unknown_on_subprocess_error(
    monkeypatch, fake_path
) -> None:
    fake_path["samtools"] = "/usr/local/bin/samtools"

    def raiser(*_args, **_kwargs):
        raise subprocess.SubprocessError("exploded")

    monkeypatch.setattr(tool_check.subprocess, "run", raiser)

    assert get_tool_version("samtools") == "unknown"


# ---------- check_tool -------------------------------------------------------


def test_check_tool_returns_tool_info(fake_path, fake_subprocess_run) -> None:
    fake_path["samtools"] = "/usr/local/bin/samtools"
    version_map, _ = fake_subprocess_run
    version_map["samtools"] = "samtools 1.21\n"

    info = check_tool("samtools")

    assert isinstance(info, ToolInfo)
    assert info.name == "samtools"
    assert info.path == "/usr/local/bin/samtools"
    assert info.version == "1.21"


def test_check_tool_raises_when_tool_absent(fake_path) -> None:
    with pytest.raises(ToolNotFoundError):
        check_tool("bowtie2")


# ---------- ensure_tools_available -------------------------------------------


def test_ensure_tools_available_all_present(fake_path, fake_subprocess_run) -> None:
    version_map, _ = fake_subprocess_run
    for tool, ver in [
        ("cutadapt", "cutadapt 4.9"),
        ("bowtie2", "bowtie2 version 2.5.4"),
        ("samtools", "samtools 1.21"),
    ]:
        fake_path[tool] = f"/usr/local/bin/{tool}"
        version_map[tool] = ver

    resolved = ensure_tools_available(
        required=("cutadapt", "bowtie2", "samtools"),
        optional=("fastqc",),  # missing; should be returned as None
    )

    assert set(resolved) == {"cutadapt", "bowtie2", "samtools", "fastqc"}
    assert resolved["fastqc"] is None
    for name in ("cutadapt", "bowtie2", "samtools"):
        assert resolved[name] is not None
        assert resolved[name].version


def test_ensure_tools_available_reports_all_missing_required_at_once(
    fake_path, fake_subprocess_run
) -> None:
    fake_path["samtools"] = "/usr/local/bin/samtools"
    version_map, _ = fake_subprocess_run
    version_map["samtools"] = "samtools 1.21"

    with pytest.raises(ToolNotFoundError) as exc:
        ensure_tools_available(
            required=("cutadapt", "bowtie2", "samtools"),
        )

    message = str(exc.value)
    # The user sees EVERY missing tool in one error, not one-at-a-time.
    assert "cutadapt" in message
    assert "bowtie2" in message
    assert "samtools" not in message.split("\n")[0]  # samtools is present; not listed
    assert "bioconda" in message


def test_ensure_tools_available_passes_optional_through_when_present(
    fake_path, fake_subprocess_run
) -> None:
    version_map, _ = fake_subprocess_run
    for tool in ("cutadapt", "fastqc"):
        fake_path[tool] = f"/usr/local/bin/{tool}"
        version_map[tool] = f"{tool} 1.0"

    resolved = ensure_tools_available(
        required=("cutadapt",),
        optional=("fastqc",),
    )

    assert resolved["cutadapt"] is not None
    assert resolved["fastqc"] is not None
    assert resolved["fastqc"].name == "fastqc"


# ---------- bowtie2_strand_flag ---------------------------------------------


@pytest.mark.parametrize(
    "strandedness, expected",
    [
        ("forward", ["--norc"]),
        ("reverse", ["--nofw"]),
        ("unstranded", []),
    ],
)
def test_bowtie2_strand_flag_maps_every_strandedness(
    strandedness: Strandedness, expected: list[str]
) -> None:
    assert bowtie2_strand_flag(strandedness) == expected


def test_bowtie2_strand_flag_rejects_unknown_value() -> None:
    with pytest.raises(ValueError):
        bowtie2_strand_flag("sideways")  # type: ignore[arg-type]
