"""Tests for dual-end UMI support (``umi_position == "both"``).

Covers:
* :func:`mitoribopy.align.trim.resolve_kit_settings` validation that
  ``umi_length`` equals ``umi_length_5p + umi_length_3p`` when position
  is ``both``.
* :func:`mitoribopy.align.trim._build_pass1_command` /
  :func:`_build_pass2_umi_command` shape: pass 1 carries
  ``-u <umi_length_5p>`` with ``{cut_prefix}`` rename; pass 2 carries
  ``-u -<umi_length_3p>`` with ``{cut_suffix}`` rename and the special
  append-to-existing-umi template.
* End-to-end synthetic FASTQ trimming where we know the inserted 5' /
  3' UMI bases — after the two-pass run the QNAME contains
  ``_<5pUMI><3pUMI>`` and the trimmed insert no longer carries those
  UMIs.
"""

from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest

from mitoribopy.align import trim
from mitoribopy.align._types import ResolvedKit
from mitoribopy.align.dedup import recommend_umi_method
from mitoribopy.align.trim import (
    _build_pass1_command,
    _build_pass2_umi_command,
    resolve_kit_settings,
    run_cutadapt,
)


# ---------------------------------------------------------------------------
# (a) KitPreset / ResolvedKit validation
# ---------------------------------------------------------------------------


def test_resolve_kit_both_requires_per_end_lengths() -> None:
    with pytest.raises(ValueError) as exc:
        resolve_kit_settings(
            "custom",
            adapter="AAAA",
            umi_position="both",
        )
    assert "umi_length_5p" in str(exc.value)


def test_resolve_kit_both_derives_umi_length_from_per_end_sum() -> None:
    resolved = resolve_kit_settings(
        "custom",
        adapter="AAAA",
        umi_position="both",
        umi_length_5p=6,
        umi_length_3p=6,
    )
    assert resolved.umi_position == "both"
    assert resolved.umi_length_5p == 6
    assert resolved.umi_length_3p == 6
    # Auto-derivation: 6 + 6 = 12.
    assert resolved.umi_length == 12


def test_resolve_kit_both_rejects_umi_length_mismatch() -> None:
    with pytest.raises(ValueError) as exc:
        resolve_kit_settings(
            "custom",
            adapter="AAAA",
            umi_position="both",
            umi_length=10,
            umi_length_5p=6,
            umi_length_3p=6,
        )
    msg = str(exc.value)
    assert "umi_length_5p + umi_length_3p" in msg


def test_resolve_kit_both_accepts_consistent_explicit_umi_length() -> None:
    resolved = resolve_kit_settings(
        "custom",
        adapter="AAAA",
        umi_position="both",
        umi_length=14,
        umi_length_5p=6,
        umi_length_3p=8,
    )
    assert resolved.umi_length == 14
    assert resolved.umi_length_5p == 6
    assert resolved.umi_length_3p == 8


def test_resolve_kit_single_end_unaffected_by_per_end_defaults() -> None:
    # Single-end UMI presets must keep working without per-end fields.
    resolved = resolve_kit_settings("nebnext_ultra_umi")
    assert resolved.umi_length == 8
    assert resolved.umi_position == "5p"
    assert resolved.umi_length_5p == 0
    assert resolved.umi_length_3p == 0


# ---------------------------------------------------------------------------
# (b) Build-command shape
# ---------------------------------------------------------------------------


def test_build_pass1_command_dual_end_extracts_only_5p_umi(tmp_path: Path) -> None:
    resolved = ResolvedKit(
        kit="custom",
        adapter="AAAA",
        umi_length=12,
        umi_position="both",
        umi_length_5p=6,
        umi_length_3p=6,
    )
    cmd = _build_pass1_command(
        fastq_in=tmp_path / "in.fq.gz",
        fastq_out=tmp_path / "intermediate.fq.gz",
        resolved=resolved,
        min_length=15,
        max_length=45,
        quality=20,
        threads=1,
        log_json=tmp_path / "log.json",
    )
    assert "-u" in cmd
    ui = cmd.index("-u")
    assert cmd[ui + 1] == "6"  # 5p length only, NOT the combined 12.
    rename = next(t for t in cmd if t.startswith("--rename="))
    assert "{cut_prefix}" in rename
    # Pass 1 must still trim the adapter so pass 2 sees a clean tail.
    assert "-a" in cmd and cmd[cmd.index("-a") + 1] == "AAAA"


def test_build_pass2_command_dual_end_appends_to_existing_umi(tmp_path: Path) -> None:
    cmd = _build_pass2_umi_command(
        fastq_in=tmp_path / "intermediate.fq.gz",
        fastq_out=tmp_path / "out.fq.gz",
        umi_length=6,
        threads=1,
        log_json=tmp_path / "pass2.json",
        append_to_existing_umi=True,
    )
    assert "-u" in cmd
    ui = cmd.index("-u")
    assert cmd[ui + 1] == "-6"
    rename = next(t for t in cmd if t.startswith("--rename="))
    # Concatenation: '{id}{cut_suffix}' (no underscore separator).
    assert "{cut_suffix}" in rename
    assert "{id}{cut_suffix}" in rename
    assert "{id}_{cut_suffix}" not in rename


def test_build_pass2_command_legacy_3p_uses_underscore_separator(
    tmp_path: Path,
) -> None:
    cmd = _build_pass2_umi_command(
        fastq_in=tmp_path / "intermediate.fq.gz",
        fastq_out=tmp_path / "out.fq.gz",
        umi_length=12,
        threads=1,
        log_json=tmp_path / "pass2.json",
        append_to_existing_umi=False,
    )
    rename = next(t for t in cmd if t.startswith("--rename="))
    assert "{id}_{cut_suffix}" in rename


def test_run_cutadapt_dual_end_issues_two_passes(tmp_path: Path) -> None:
    recorded: list[list[str]] = []

    def runner(cmd, **kwargs):
        recorded.append(list(cmd))
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    log_json = tmp_path / "sample.cutadapt.json"
    pass1_log = log_json.with_name(log_json.stem + ".pass1.json")
    pass2_log = log_json.with_name(log_json.stem + ".pass2.json")
    pass1_log.write_text(
        json.dumps(
            {
                "read_counts": {"input": 1000, "output": 950},
                "adapters_read1": [{"total_matches": 920}],
            }
        )
    )
    pass2_log.write_text(
        json.dumps(
            {
                "read_counts": {"input": 950, "output": 920},
                "adapters_read1": [],
            }
        )
    )

    resolved = ResolvedKit(
        kit="custom",
        adapter="AAAA",
        umi_length=12,
        umi_position="both",
        umi_length_5p=5,
        umi_length_3p=7,
    )
    result = run_cutadapt(
        fastq_in=tmp_path / "in.fq.gz",
        fastq_out=tmp_path / "out.fq.gz",
        resolved=resolved,
        log_json=log_json,
        runner=runner,
    )

    assert len(recorded) == 2
    pass1_cmd, pass2_cmd = recorded

    # Pass 1: 5p UMI extraction (5 nt) + adapter trim.
    assert pass1_cmd[pass1_cmd.index("-u") + 1] == "5"
    assert pass1_cmd[pass1_cmd.index("-a") + 1] == "AAAA"

    # Pass 2: 3p UMI extraction (-7) with the appending rename template.
    assert pass2_cmd[pass2_cmd.index("-u") + 1] == "-7"
    rename = next(t for t in pass2_cmd if t.startswith("--rename="))
    assert "{id}{cut_suffix}" in rename

    # Pass 1 length filter must reserve room for the 3p UMI bases pass 2 strips.
    assert pass1_cmd[pass1_cmd.index("--minimum-length") + 1] == "22"
    assert pass1_cmd[pass1_cmd.index("--maximum-length") + 1] == "52"

    assert result.input_reads == 1000
    assert result.reads_passing_filters == 920


# ---------------------------------------------------------------------------
# (c) End-to-end synthetic FASTQ — exercises a real cutadapt binary
# ---------------------------------------------------------------------------


_INSERT = "ACGTACGTACGTACGTACGT"  # 20 nt biological insert
_UMI_5P = "AAAAA"                 # 5 nt 5' UMI
_UMI_3P = "TTTTTTT"               # 7 nt 3' UMI
_ADAPTER = "GGGGGGGGGGGGGGGG"     # 16 nt 3' adapter

_HAS_CUTADAPT = shutil.which("cutadapt") is not None


def _write_fastq(path: Path) -> None:
    """Write a one-record FASTQ that carries 5'UMI + insert + 3'UMI + adapter."""
    seq = _UMI_5P + _INSERT + _UMI_3P + _ADAPTER
    qual = "I" * len(seq)
    record = f"@read1\n{seq}\n+\n{qual}\n"
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as handle:
            handle.write(record)
    else:
        path.write_text(record)


def _read_first_record(fastq_out: Path) -> tuple[str, str]:
    if fastq_out.suffix == ".gz":
        with gzip.open(fastq_out, "rt") as handle:
            lines = handle.read().splitlines()
    else:
        lines = fastq_out.read_text().splitlines()
    assert len(lines) >= 4, f"unexpected FASTQ output: {lines}"
    qname = lines[0][1:]  # strip the leading '@'
    seq = lines[1]
    return qname, seq


@pytest.mark.skipif(not _HAS_CUTADAPT, reason="cutadapt not installed")
def test_dual_end_umi_end_to_end_qname_and_insert(tmp_path: Path) -> None:
    fastq_in = tmp_path / "in.fq.gz"
    fastq_out = tmp_path / "out.fq.gz"
    _write_fastq(fastq_in)

    resolved = resolve_kit_settings(
        "custom",
        adapter=_ADAPTER,
        umi_position="both",
        umi_length_5p=len(_UMI_5P),
        umi_length_3p=len(_UMI_3P),
    )
    run_cutadapt(
        fastq_in=fastq_in,
        fastq_out=fastq_out,
        resolved=resolved,
        min_length=10,
        max_length=80,
        quality=0,
        log_json=tmp_path / "cutadapt.json",
    )

    qname, insert = _read_first_record(fastq_out)
    # QNAME must be ``read1_<5pUMI><3pUMI>`` (concatenated, no
    # separator between the two UMIs).
    assert qname.split()[0] == f"read1_{_UMI_5P}{_UMI_3P}"
    # The trimmed insert no longer carries the UMI bases.
    assert insert == _INSERT


# ---------------------------------------------------------------------------
# (d) recommend_umi_method picks based on combined length
# ---------------------------------------------------------------------------


def test_recommend_umi_method_uses_combined_length_for_dual_end() -> None:
    # Per-end 5 nt + 5 nt = 10 nt combined; high duplication -> directional.
    method, warn = recommend_umi_method(
        umi_length=10,
        duplicate_fraction=0.5,
        umi_length_5p=5,
        umi_length_3p=5,
    )
    assert method == "directional"
    assert warn == "none"


def test_recommend_umi_method_short_combined_warns() -> None:
    # 3 + 3 = 6 nt combined: below threshold (8), warn + force unique.
    method, warn = recommend_umi_method(
        umi_length=6,
        duplicate_fraction=0.5,
        umi_length_5p=3,
        umi_length_3p=3,
    )
    assert method == "unique"
    assert warn == "short_umi_collision_risk"


# ---------------------------------------------------------------------------
# (e) validate-config rejects a 'both' config missing per-end lengths
# ---------------------------------------------------------------------------


def test_validate_config_rejects_both_without_per_end_lengths(
    tmp_path: Path, capsys
) -> None:
    from mitoribopy import cli

    cfg = tmp_path / "dual.yaml"
    cfg.write_text(
        "align:\n"
        "  kit_preset: custom\n"
        "  adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\n"
        "  umi_position: both\n"
        # umi_length_5p / umi_length_3p deliberately omitted
        "rpf:\n"
        "  strain: h.sapiens\n"
        "  fasta: /tmp/tx.fa\n",
        encoding="utf-8",
    )
    rc = cli.main(["validate-config", str(cfg), "--no-path-checks"])
    assert rc == 2
    err = capsys.readouterr().err
    assert "umi_position='both'" in err
    assert "umi_length_5p" in err and "umi_length_3p" in err


def test_validate_config_accepts_both_with_consistent_lengths(
    tmp_path: Path, capsys
) -> None:
    from mitoribopy import cli

    cfg = tmp_path / "dual_ok.yaml"
    cfg.write_text(
        "align:\n"
        "  kit_preset: custom\n"
        "  adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\n"
        "  umi_position: both\n"
        "  umi_length: 12\n"
        "  umi_length_5p: 6\n"
        "  umi_length_3p: 6\n"
        "rpf:\n"
        "  strain: h.sapiens\n"
        "  fasta: /tmp/tx.fa\n",
        encoding="utf-8",
    )
    rc = cli.main(["validate-config", str(cfg), "--no-path-checks"])
    assert rc == 0
    err = capsys.readouterr().err
    assert "OK" in err
