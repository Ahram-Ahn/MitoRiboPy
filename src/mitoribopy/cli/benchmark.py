"""``mitoribopy benchmark`` (P2.13) — measured `mitoribopy all` runs.

Wraps :mod:`mitoribopy.cli.all_` with wall-time + peak-RSS + disk-
footprint measurement and emits ``benchmark.tsv`` plus
``benchmark_summary.md`` at the run root. Optionally pre-subsamples
each Ribo / RNA FASTQ to ``--subsample N`` reads for fast tuning runs.

Output columns of ``benchmark.tsv``:

* ``stage``               — align | rpf | rnaseq | total
* ``status``              — completed | skipped | not_configured
* ``wall_time_sec``       — runtime captured by the orchestrator
                            (``manifest.stages.<stage>.runtime_seconds``)
* ``cumulative_wall_sec`` — wall-clock sum up to and including the row
* ``max_rss_mb_total``    — peak RSS (MB) across the WHOLE
                            ``mitoribopy benchmark`` invocation; only
                            populated on the ``total`` row
* ``disk_mb``             — output-directory size after the run, on
                            the ``total`` row
* ``threads``             — value of ``--threads`` for this run
* ``subsample_reads``     — value of ``--subsample`` (or empty)

The MD summary mirrors the TSV in human-readable form.
"""

from __future__ import annotations

import argparse
import json
import resource
import shutil
import sys
import time
from copy import deepcopy
from pathlib import Path
from typing import Iterable

from . import common


BENCHMARK_SUBCOMMAND_HELP = (
    "Time and disk-measure a full `mitoribopy all` run, optionally "
    "after pre-subsampling each FASTQ to N reads. Produces "
    "benchmark.tsv and benchmark_summary.md at the run root for "
    "tuning thread counts, disk budgets, and per-stage wall time."
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mitoribopy benchmark",
        description=BENCHMARK_SUBCOMMAND_HELP,
        formatter_class=common.MitoRiboPyHelpFormatter,
    )
    parser.add_argument(
        "--config",
        required=True,
        metavar="PATH",
        help="Pipeline YAML / JSON / TOML, exactly as accepted by "
             "`mitoribopy all --config`.",
    )
    parser.add_argument(
        "--output",
        required=True,
        metavar="DIR",
        help="Run root for the benchmarked execution. Will be created "
             "if missing; existing contents are not deleted.",
    )
    parser.add_argument(
        "--subsample",
        type=int,
        default=None,
        metavar="N",
        help="If set, reservoir-sample each Ribo / RNA FASTQ to N reads "
             "before running. Subsampled copies are written under "
             "<output>/.benchmark_subsamples/ and the canonical config "
             "is rewritten to point at them. Requires reservoir-friendly "
             "FASTQ inputs (no streaming sources).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        metavar="N",
        help="Forwarded to `mitoribopy all --threads`. Recorded in the "
             "benchmark.tsv `threads` column.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for the FASTQ subsampler (no effect when "
             "--subsample is omitted).",
    )
    return parser


# ---------- subsampler ----------------------------------------------------


def _subsample_fastq_inputs(
    config: dict, *, n_reads: int, seed: int, work_dir: Path
) -> dict:
    """Return a copy of *config* with every FASTQ path swapped for a
    reservoir-sampled copy under *work_dir*.

    Touches: ``align.fastq`` (string or list), ``align.fastq_dir``
    (when set, every *.fq[.gz] / *.fastq[.gz] inside is sampled),
    ``rnaseq.rna_fastq``, ``rnaseq.ribo_fastq``. Top-level ``samples:``
    is left as-is — sample-sheet driven runs are pre-subsampled by
    swapping the FASTQ paths inside the sheet, which is more involved
    and not yet supported.
    """
    from ..tools.subsample_fastq import reservoir_sample_fastq

    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    cfg = deepcopy(config)

    def _sample_one(src: Path) -> Path:
        dst = work_dir / src.name
        if dst.exists():
            return dst
        sampled = reservoir_sample_fastq(src, n_reads, seed)
        # Replicate _open_text from subsample_fastq inline to keep
        # this module independent of the script's main().
        if dst.suffix == ".gz":
            import gzip

            with gzip.open(dst, "wt", encoding="utf-8") as out:
                for header, seq, plus, qual in sampled:
                    out.write(header)
                    out.write(seq)
                    out.write(plus)
                    out.write(qual)
        else:
            with dst.open("w", encoding="utf-8") as out:
                for header, seq, plus, qual in sampled:
                    out.write(header)
                    out.write(seq)
                    out.write(plus)
                    out.write(qual)
        return dst

    def _expand_dir(d: Path) -> list[Path]:
        patterns = ("*.fq.gz", "*.fastq.gz", "*.fq", "*.fastq")
        out: list[Path] = []
        for pat in patterns:
            out.extend(sorted(d.glob(pat)))
        return out

    align = cfg.get("align") or {}
    if isinstance(align, dict):
        if isinstance(align.get("fastq"), str):
            d = Path(align["fastq"])
            if d.is_dir():
                paths = [_sample_one(p) for p in _expand_dir(d)]
                align["fastq"] = [str(p) for p in paths]
            elif d.is_file():
                align["fastq"] = str(_sample_one(d))
        elif isinstance(align.get("fastq"), list):
            align["fastq"] = [str(_sample_one(Path(p))) for p in align["fastq"]]
        if align.get("fastq_dir"):
            d = Path(align["fastq_dir"])
            if d.is_dir():
                paths = [_sample_one(p) for p in _expand_dir(d)]
                align["fastq"] = [str(p) for p in paths]
                align["fastq_dir"] = None

    rnaseq = cfg.get("rnaseq") or {}
    for key in ("rna_fastq", "ribo_fastq"):
        value = rnaseq.get(key)
        if isinstance(value, str):
            d = Path(value)
            if d.is_dir():
                rnaseq[key] = [str(_sample_one(p)) for p in _expand_dir(d)]
            elif d.is_file():
                rnaseq[key] = [str(_sample_one(d))]
        elif isinstance(value, list):
            rnaseq[key] = [str(_sample_one(Path(p))) for p in value]

    return cfg


# ---------- measurement ---------------------------------------------------


def _measure_max_rss_mb() -> float:
    """Peak RSS for THIS process (and its children on Linux)."""
    info = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    # ru_maxrss is kilobytes on Linux, bytes on macOS. Normalise to MB.
    factor = 1024.0 if sys.platform == "darwin" else 1.0  # macOS bytes→KB
    raw = max(info.ru_maxrss, children.ru_maxrss)
    return raw / factor / 1024.0


def _directory_size_mb(path: Path) -> float:
    total = 0
    for p in path.rglob("*"):
        try:
            if p.is_file():
                total += p.stat().st_size
        except OSError:
            continue
    return total / 1024.0 / 1024.0


def _yaml_dump(payload: dict) -> str:
    try:
        import yaml  # type: ignore[import-not-found]
    except ImportError:
        return json.dumps(payload, indent=2, sort_keys=True) + "\n"
    return yaml.safe_dump(
        payload, sort_keys=True, default_flow_style=False, allow_unicode=True
    )


# ---------- runner --------------------------------------------------------


def _write_benchmark_tsv(
    rows: list[dict], path: Path
) -> Path:
    from ..io.schema_versions import OUTPUT_SCHEMA_VERSIONS

    columns = (
        "stage", "status", "wall_time_sec", "cumulative_wall_sec",
        "max_rss_mb_total", "disk_mb", "threads", "subsample_reads",
    )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # benchmark.tsv schema version is implicit (not yet in the registry)
    # but we still prepend a marker line so consumers can rely on the
    # comment-line convention.
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# schema_version: 1.0\n")
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write(
                "\t".join(str(row.get(col, "") or "") for col in columns)
                + "\n"
            )
    return path


def _write_benchmark_summary_md(
    rows: list[dict],
    *,
    config_source: str,
    threads: int | None,
    subsample: int | None,
    path: Path,
) -> Path:
    lines = [
        "# MitoRiboPy benchmark summary",
        "",
        f"- config: `{config_source}`",
        f"- threads: `{threads if threads is not None else 'default'}`",
        f"- subsample: "
        + (f"`{subsample}` reads/FASTQ" if subsample else "(off)"),
        "",
        "## Per-stage timing",
        "",
        "| stage | status | wall (s) | cumulative (s) |",
        "|---|---|---|---|",
    ]
    for row in rows:
        if row["stage"] == "total":
            continue
        lines.append(
            f"| {row['stage']} | {row['status']} | "
            f"{row['wall_time_sec']} | {row['cumulative_wall_sec']} |"
        )
    total = next((r for r in rows if r["stage"] == "total"), None)
    if total is not None:
        lines.extend([
            "",
            "## Resource totals",
            "",
            f"- total wall: **{total['wall_time_sec']} s**",
            f"- peak RSS: **{total['max_rss_mb_total']} MB**",
            f"- disk footprint: **{total['disk_mb']} MB**",
        ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def run(argv: Iterable[str]) -> int:
    args = build_parser().parse_args(list(argv))
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        cfg = common.load_config_file(args.config)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"[mitoribopy benchmark] ERROR: {exc}", file=sys.stderr)
        return 2

    config_for_run = cfg
    if args.subsample is not None and args.subsample > 0:
        sub_dir = output_dir / ".benchmark_subsamples"
        try:
            config_for_run = _subsample_fastq_inputs(
                cfg,
                n_reads=args.subsample,
                seed=args.seed,
                work_dir=sub_dir,
            )
        except (FileNotFoundError, ValueError) as exc:
            print(
                f"[mitoribopy benchmark] ERROR: subsample failed: {exc}",
                file=sys.stderr,
            )
            return 2
        sys.stderr.write(
            f"[mitoribopy benchmark] subsampled FASTQ inputs to "
            f"{args.subsample} reads under {sub_dir}.\n"
        )

    # Persist the (possibly subsampled) config so the wrapped `all`
    # invocation has a real --config path to read.
    benchmark_cfg = output_dir / "benchmark_config.yaml"
    benchmark_cfg.write_text(_yaml_dump(config_for_run), encoding="utf-8")

    forwarded_argv = ["--config", str(benchmark_cfg), "--output", str(output_dir)]
    if args.threads is not None:
        forwarded_argv.extend(["--threads", str(args.threads)])

    from . import all_ as _all_cli

    t0 = time.monotonic()
    rc = _all_cli.run(forwarded_argv)
    wall_total = time.monotonic() - t0
    if rc != 0:
        sys.stderr.write(
            f"[mitoribopy benchmark] inner `mitoribopy all` exited {rc}; "
            "benchmark.tsv will record a partial run.\n"
        )

    # Read per-stage times back out of the manifest the orchestrator wrote.
    manifest_path = output_dir / "run_manifest.json"
    stage_times: dict[str, dict] = {}
    if manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            stage_times = manifest.get("stages") or {}
        except (OSError, json.JSONDecodeError):
            stage_times = {}

    rows: list[dict] = []
    cumulative = 0.0
    for stage in ("align", "rpf", "rnaseq"):
        info = stage_times.get(stage, {})
        wall = info.get("runtime_seconds")
        if wall is not None:
            cumulative += float(wall)
        rows.append(
            {
                "stage": stage,
                "status": info.get("status", "not_configured"),
                "wall_time_sec": (
                    f"{wall:.3f}" if isinstance(wall, (int, float)) else ""
                ),
                "cumulative_wall_sec": f"{cumulative:.3f}",
                "max_rss_mb_total": "",
                "disk_mb": "",
                "threads": str(args.threads) if args.threads is not None else "",
                "subsample_reads": str(args.subsample) if args.subsample else "",
            }
        )
    rows.append(
        {
            "stage": "total",
            "status": "ok" if rc == 0 else f"exit:{rc}",
            "wall_time_sec": f"{wall_total:.3f}",
            "cumulative_wall_sec": f"{wall_total:.3f}",
            "max_rss_mb_total": f"{_measure_max_rss_mb():.1f}",
            "disk_mb": f"{_directory_size_mb(output_dir):.1f}",
            "threads": str(args.threads) if args.threads is not None else "",
            "subsample_reads": str(args.subsample) if args.subsample else "",
        }
    )

    tsv_path = _write_benchmark_tsv(rows, output_dir / "benchmark.tsv")
    md_path = _write_benchmark_summary_md(
        rows,
        config_source=args.config,
        threads=args.threads,
        subsample=args.subsample,
        path=output_dir / "benchmark_summary.md",
    )
    sys.stderr.write(
        f"[mitoribopy benchmark] wrote {tsv_path} and {md_path}.\n"
    )
    return rc
