"""MitoRiboPy command-line entrypoint with subcommand dispatch.

Top-level synopsis::

    mitoribopy [--version] [--help]
    mitoribopy align              [options]   # FASTQ -> BAM + BED + read counts
    mitoribopy rpf                [options]   # Ribo-seq analysis from BED / BAM
    mitoribopy rnaseq             [options]   # DE table + rpf -> TE / dTE
    mitoribopy all                [options]   # align + rpf + (optional) rnaseq
    mitoribopy periodicity        [options]   # 3-nt periodicity QC bundle
    mitoribopy migrate-config     <path>      # rewrite legacy YAML keys
    mitoribopy validate-config    <path>      # schema + cross-stage checks
    mitoribopy validate-reference [options]   # FASTA / annotation sanity checks
    mitoribopy validate-figures   <run_dir>   # figure_qc.tsv producer
    mitoribopy summarize          <run_dir>   # SUMMARY.md / outputs_index.tsv
    mitoribopy benchmark          [options]   # wall / RSS benchmark harness

A subcommand is required. The pre-v0.6.0 fallback that re-routed
``mitoribopy -s h -f ref.fa ...`` to ``mitoribopy rpf ...`` was
removed in v0.6.0 (publication freeze). Use the explicit form::

    mitoribopy rpf -s h.sapiens -f ref.fa ...
"""

from __future__ import annotations

import importlib
import sys
from typing import Callable, Iterable

from .. import __version__

__all__ = ["main", "run_pipeline_cli"]


def __getattr__(name: str):
    """Lazy attribute access for the back-compat ``run_pipeline_cli``
    re-export.

    Eagerly importing ``mitoribopy.pipeline.runner`` at module-init
    time creates a cli<->pipeline circular when callers import
    ``mitoribopy.pipeline.context`` (or any other pipeline submodule)
    before they ever reach the ``cli`` package. Deferring the import
    until the attribute is actually accessed keeps the convenience
    while letting pipeline-only callers (tests, scripts) avoid the
    cli load.
    """
    if name == "run_pipeline_cli":
        from ..pipeline.runner import run_pipeline_cli as _impl
        return _impl
    raise AttributeError(f"module 'mitoribopy.cli' has no attribute {name!r}")


_SubcommandRun = Callable[[Iterable[str]], int]


def _lazy_run(module_name: str) -> _SubcommandRun:
    """Load subcommand modules only when the user invokes that subcommand."""
    def _run(argv: Iterable[str]) -> int:
        module = importlib.import_module(f"{__name__}.{module_name}")
        return module.run(argv)

    return _run


_SUBCOMMAND_DEFS: list[tuple[str, str, str]] = [
    (
        "align",
        "align",
        "Preprocess FASTQ inputs: cutadapt trim + bowtie2 contaminant "
        "subtraction + bowtie2 mt-transcriptome alignment + MAPQ filter "
        "+ dedup + BAM->BED6. Produces drop-in inputs for 'mitoribopy rpf'.",
    ),
    (
        "rpf",
        "rpf",
        "Run the Ribo-seq analysis pipeline from BED or BAM inputs. "
        "See 'mitoribopy rpf --help' for the full flag list.",
    ),
    (
        "rnaseq",
        "rnaseq",
        "Translation efficiency (TE / delta-TE) from paired RNA-seq + "
        "Ribo-seq. Default flow: pass --rna-fastq + --ribo-fastq + "
        "--reference-fasta and the subcommand runs trimming, bowtie2 "
        "alignment, per-transcript counting, and pyDESeq2 itself before "
        "emitting te.tsv, delta_te.tsv, and plots. Alternative: pass "
        "--de-table from a prior external DESeq2 / Xtail / Anota2Seq run "
        "together with --ribo-dir; this path is mutually exclusive with "
        "--rna-fastq and enforces a SHA256 reference-consistency gate.",
    ),
    (
        "all",
        "all_",
        "End-to-end orchestrator: align + rpf, plus rnaseq when the config "
        "carries an 'rnaseq' section configured for either flow "
        "(from-FASTQ via 'rna_fastq' + 'reference_fasta', or external-DE "
        "via 'de_table'). Writes a composed run_manifest.json with tool "
        "versions, parameters, and input/output hashes across all three stages.",
    ),
    (
        "periodicity",
        "periodicity",
        "Quantify 3-nt periodicity by running the metagene Fourier "
        "analysis on a pre-assigned site table.",
    ),
    (
        "migrate-config",
        "migrate_config",
        "Rewrite legacy MitoRiboPy YAML keys to their canonical names. "
        "Input is read from a path; output is written to stdout (the change "
        "log goes to stderr). Use to upgrade old pipeline configs without "
        "manually hunting down every renamed key.",
    ),
    (
        "validate-config",
        "validate_config",
        "Pre-flight a MitoRiboPy YAML / JSON / TOML config: parse, "
        "canonicalise legacy keys, check file paths and mutually-exclusive "
        "sections, and resolve rnaseq.mode against supplied inputs. Exit "
        "code is 0 on success, 2 when at least one error was found.",
    ),
    (
        "validate-reference",
        "validate_reference",
        "Pre-flight a custom mitochondrial reference: check that the FASTA "
        "and annotation CSV are consistent (matching transcript IDs, "
        "matching lengths, CDS divisible by 3, valid start / stop codons "
        "under the selected codon table).",
    ),
    (
        "validate-figures",
        "validate_figures",
        "Mechanically validate every plot under a finished MitoRiboPy run "
        "root: check label / legend / stat-box overlap, label clipping, "
        "point counts vs source TSV, SVG text editability, and PNG dpi. "
        "Writes <RUN_DIR>/figure_qc.tsv. Exit "
        "0 / 1 / 2 (all pass / warn-only / fail; --strict upgrades warn to fail).",
    ),
    (
        "summarize",
        "summarize",
        "Regenerate SUMMARY.md and summary_qc.tsv from a finished MitoRiboPy "
        "run by reading the run_manifest.json and per-stage TSV outputs. "
        "Useful for re-rendering summaries on archival runs without "
        "re-executing any pipeline stage.",
    ),
    (
        "benchmark",
        "benchmark",
        "Time and disk-measure a full `mitoribopy all` run, optionally "
        "after pre-subsampling each FASTQ to N reads. Produces "
        "benchmark.tsv and benchmark_summary.md at the run root for "
        "tuning thread counts, disk budgets, and per-stage wall time.",
    ),
]


_SUBCOMMANDS: dict[str, _SubcommandRun] = {
    name: _lazy_run(module_name) for name, module_name, _ in _SUBCOMMAND_DEFS
}


_SUBCOMMAND_SUMMARIES: list[tuple[str, str]] = [
    (name, summary) for name, _, summary in _SUBCOMMAND_DEFS
]


_TOP_HELP = """\
usage: mitoribopy [--version] <subcommand> [options]

MitoRiboPy - mitochondrial ribosome profiling analysis pipeline.

Subcommands:
{rows}

Top-level options:
  --version, -V  Print MitoRiboPy version and exit.
  --help,    -h  Print this help message.

Run 'mitoribopy <subcommand> --help' for per-subcommand options.
"""


def _normalize_args(argv: Iterable[str]) -> list[str]:
    """Strip legacy 'run' / '--' separators at the front of *argv*.

    Preserved verbatim from the pre-subcommand CLI so that existing tests
    and user scripts invoking ``mitoribopy run -- <args>`` continue to work.
    """
    args = list(argv)
    if args and args[0] == "run":
        args = args[1:]
    if args and args[0] == "--":
        args = args[1:]
    return args


def _print_top_help() -> None:
    name_w = max(len(name) for name, _ in _SUBCOMMAND_SUMMARIES) + 2
    rows = "\n".join(
        f"  {name:<{name_w}}{summary}" for name, summary in _SUBCOMMAND_SUMMARIES
    )
    print(_TOP_HELP.format(rows=rows))


def main(argv: Iterable[str] | None = None) -> int:
    """Dispatch the top-level ``mitoribopy`` entry point."""
    raw_args = sys.argv[1:] if argv is None else list(argv)
    args = _normalize_args(raw_args)

    if args in (["--version"], ["-V"]):
        print(f"MitoRiboPy {__version__}")
        return 0

    if not args or args[0] in ("-h", "--help"):
        _print_top_help()
        return 0

    first, rest = args[0], args[1:]

    if first in _SUBCOMMANDS:
        return _SUBCOMMANDS[first](rest)

    if first.startswith("-"):
        sys.stderr.write(
            f"mitoribopy: error: missing subcommand. Got '{first}' as the "
            "first token, but a subcommand is required as of v0.6.0. "
            f"Did you mean 'mitoribopy rpf {first} ...'? "
            "Run 'mitoribopy --help' for the full list.\n"
        )
        return 2

    sys.stderr.write(
        f"mitoribopy: error: unknown subcommand '{first}'. "
        "Known subcommands: " + ", ".join(_SUBCOMMANDS) + ". "
        "Try 'mitoribopy --help'.\n"
    )
    return 2
