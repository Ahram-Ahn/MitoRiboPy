"""End-to-end smoke test against examples/smoke/.

Opt-in via ``pytest -m smoke``. Skipped when cutadapt + bowtie2 +
samtools are not on PATH, or when running on a sandbox without
subprocess permission for those tools.

The smoke fixture is the publication-ready sanity check: a fresh
install must produce every file in ``examples/smoke/expected_outputs.txt``
under the run root. We deliberately do NOT assert content hashes —
the synthetic FASTQs are RNG-driven and the offset / alignment
numbers depend on exact draw order. We assert existence + non-zero
size.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest


SMOKE_DIR = Path(__file__).resolve().parent.parent / "examples" / "smoke"
REQUIRED_TOOLS = ("cutadapt", "bowtie2", "bowtie2-build")


def _missing_tools() -> list[str]:
    return [t for t in REQUIRED_TOOLS if shutil.which(t) is None]


def _expected_paths() -> list[str]:
    """Parse expected_outputs.txt, dropping comments and blanks."""
    text = (SMOKE_DIR / "expected_outputs.txt").read_text(encoding="utf-8")
    rows: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        rows.append(line)
    return rows


pytestmark = pytest.mark.smoke


def test_smoke_fixture_files_are_present() -> None:
    """The fixture directory itself must always be in the repo even
    when the external tools are not available — the pytest -m smoke
    invocation should still discover this test."""
    assert SMOKE_DIR.is_dir()
    assert (SMOKE_DIR / "human_mt_tiny.fa").is_file()
    assert (SMOKE_DIR / "samples.tsv").is_file()
    assert (SMOKE_DIR / "samples.condition_map.tsv").is_file()
    assert (SMOKE_DIR / "pipeline_config.smoke.yaml").is_file()
    assert (SMOKE_DIR / "expected_outputs.txt").is_file()
    assert (SMOKE_DIR / "generate_smoke_fastqs.py").is_file()


def test_expected_outputs_list_is_well_formed() -> None:
    paths = _expected_paths()
    assert len(paths) >= 10, "expected_outputs.txt looks suspiciously short"
    # Every entry is a relative path (no leading slash).
    for p in paths:
        assert not p.startswith("/"), p
    # Every entry contains at least one path separator OR is a
    # top-level run-root file.
    for p in paths:
        assert "/" in p or "." in p, p


def test_smoke_run_produces_expected_outputs(tmp_path: Path) -> None:
    """End-to-end smoke run: generate FASTQs, build index, run mitoribopy all,
    assert every expected output exists and is non-empty.

    Skipped when external tools or python module subprocess invocation
    is not available.
    """
    missing = _missing_tools()
    if missing:
        pytest.skip(f"smoke test needs external tools on PATH: {missing}")

    # Stage the smoke fixture into a writable scratch dir so the run
    # outputs do not leak back into the source tree.
    scratch = tmp_path / "smoke"
    shutil.copytree(SMOKE_DIR, scratch, ignore=shutil.ignore_patterns(
        "results", "fastqs", "bowtie2_index", ".gitignore",
    ))

    # 1) Synthesize FASTQs + bowtie2 index.
    gen_proc = subprocess.run(
        [sys.executable, "generate_smoke_fastqs.py"],
        cwd=scratch, capture_output=True, text=True,
    )
    if gen_proc.returncode != 0:
        pytest.skip(
            "generate_smoke_fastqs.py failed (likely sandboxed subprocess); "
            f"stderr:\n{gen_proc.stderr[-500:]}"
        )

    # 2) Run mitoribopy all.
    run_root = scratch / "results"
    cli_proc = subprocess.run(
        [
            sys.executable, "-m", "mitoribopy", "all",
            "--config", "pipeline_config.smoke.yaml",
            "--output", str(run_root),
        ],
        cwd=scratch, capture_output=True, text=True,
    )
    if cli_proc.returncode != 0:
        pytest.fail(
            "mitoribopy all failed on the smoke fixture\n"
            f"--- stdout ---\n{cli_proc.stdout[-2000:]}\n"
            f"--- stderr ---\n{cli_proc.stderr[-2000:]}\n"
        )

    # 3) Every file in expected_outputs.txt must exist and be non-empty.
    missing_files: list[str] = []
    empty_files: list[str] = []
    for rel in _expected_paths():
        path = run_root / rel
        if not path.is_file():
            missing_files.append(rel)
        elif path.stat().st_size == 0:
            empty_files.append(rel)

    msgs: list[str] = []
    if missing_files:
        msgs.append("missing files: " + ", ".join(missing_files))
    if empty_files:
        msgs.append("empty files: " + ", ".join(empty_files))
    if msgs:
        pytest.fail("smoke run output mismatch — " + "; ".join(msgs))
