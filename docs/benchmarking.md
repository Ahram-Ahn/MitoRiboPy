# Benchmarking protocol

Reproducible runtime / memory / disk numbers for `mitoribopy all`.
The numbers in this page back the "Benchmark" table in the
manuscript Methods; CI does not maintain them, but the benchmarking
*command* below is how to regenerate them on a comparable machine.

> **Status (v0.7.1).** Numbers below are placeholders pending the
> publication-grade dataset run on the maintainer's reference HPC
> tier. The benchmarking *protocol* (commands, table schema,
> required artefacts) is finalised; only the row entries are TBD.
> See P1.4 in
> [`docs/developer/release_checklist.md`](developer/release_checklist.md)
> for the readiness gate this page closes.

---

## Command

`mitoribopy benchmark` is a thin wrapper around `mitoribopy all`
that records wall time, peak RSS, and disk footprint per stage:

```bash
mitoribopy benchmark \
  --config <pipeline_config.yaml> \
  --output <run_root> \
  --threads <N> \
  --subsample <K>     # optional: reservoir-sample every Ribo / RNA FASTQ to K reads
```

The wrapper writes two artefacts at the run root:

* `benchmark.tsv` — one row per stage with `wall_seconds`,
  `peak_rss_gb`, `disk_delta_gb`, `threads`, `parallel_samples`,
  `cpu_user_seconds`, `cpu_system_seconds`.
* `benchmark_summary.md` — human-readable rollup with the totals
  and a recommended thread / parallel-samples plan for the same
  hardware.

For tuning runs, the `--subsample N` flag pre-reservoir-samples
every Ribo / RNA FASTQ to N reads (default seed 42; configurable
via `--seed`). Subsampled FASTQs land under
`<output>/.benchmark_subsamples/` and the canonical config is
rewritten to point at them, so re-runs of `mitoribopy all` against
the same `--output` reproduce the subsampled run.

Subsampled rows are tagged in the benchmark table — never compare a
subsampled wall-time to a full-data wall-time except as an upper /
lower bound.

---

## Required table schema (manuscript-grade)

The published Methods table follows this exact column order so a
reviewer can compare runs across labs / hardware:

| Column | Units | Source | Notes |
|---|---|---|---|
| `dataset` | — | hand-filled | Free-form name (e.g. "TACO1-KO biological", "GSE############ subset") |
| `organism` | — | hand-filled | `h.sapiens` / `s.cerevisiae` / custom strain label |
| `assay` | — | hand-filled | `ribo` / `rna` / `ribo+rna` |
| `n_samples` | int | sample sheet | Active rows after `exclude=true` filter |
| `reads_per_sample` | reads | `align/read_counts.tsv` | Median `total_reads` across active samples |
| `input_size_gb` | GB | `du -sh <fastq_dir>` | Sum of compressed FASTQ sizes |
| `threads` | int | `--threads` arg | Total worker threads |
| `parallel_samples` | int | resolved plan | From `resource_plan.json` |
| `wall_time_min` | minutes | `benchmark.tsv` | Sum across stages |
| `peak_rss_gb` | GB | `benchmark.tsv` | Maximum observed across stages |
| `disk_peak_gb` | GB | `benchmark.tsv` | Peak `disk_delta_gb` (excludes subsamples) |
| `mitoribopy_version` | — | `mitoribopy --version` | Must be the PyPI / Zenodo-archived version |
| `commit_sha` | — | `git rev-parse HEAD` | Captured into `run_manifest.json` automatically |
| `cutadapt_version` | — | `cutadapt --version` | Mirrored in `run_manifest.json` `tool_versions` |
| `bowtie2_version` | — | `bowtie2 --version` | Same |
| `umi_tools_version` | — | `umi_tools --version` | Same; empty when no sample needs it |
| `pysam_version` | — | `pysam.__version__` | Same |
| `notes` | — | hand-filled | "subsampled to N reads", "rerun from cache", etc. |

The first 11 columns are objective; the last 6 lock in
reproducibility. The `notes` column is the only place to mention
caveats — anything that should change how a reader interprets the
runtime (subsampling, filesystem class, background load) goes there.

---

## Required benchmark cases

The Methods section lists at minimum these five cases:

1. **Synthetic smoke dataset** — `examples/smoke/` with the
   default 6 k reads × 2 samples. The "lower bound on wall time"
   reference.
2. **Small human public mt-Ribo-seq subset** — a single sample
   from a public GEO library, subsampled to 1 M reads. Calibrates
   the per-million-reads rate.
3. **Small yeast public mt-Ribo-seq subset** — same as above for
   `s.cerevisiae`. Confirms strain alignment / annotation paths.
4. **TACO1 WT/KO biological** — three biological replicates per
   condition, full reads, full TE / ΔTE pipeline. The
   "publication-realistic" reference.
5. **One full realistic non-subsampled dataset** — a current-lab
   library at typical sequencing depth, no subsample. The
   "what does a real run cost" anchor.

Subsampled and non-subsampled rows are clearly separated in the
table (the `notes` column).

---

## Reference numbers (TBD)

| dataset | organism | assay | n_samples | reads_per_sample | input_size_gb | threads | parallel_samples | wall_time_min | peak_rss_gb | disk_peak_gb | mitoribopy_version | commit_sha | cutadapt_version | bowtie2_version | umi_tools_version | pysam_version | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| smoke synthetic | h.sapiens (synthetic) | ribo | 2 | 3 000 | 0.001 | 2 | 1 | TBD | TBD | TBD | 0.7.1 | TBD | TBD | TBD | n/a | TBD | smoke fixture; opt-in via `pytest -m smoke` |
| public-human-1M | h.sapiens | ribo | 1 | 1 000 000 | TBD | 8 | 1 | TBD | TBD | TBD | 0.7.1 | TBD | TBD | TBD | TBD | TBD | subsampled (`--subsample 1000000`) |
| public-yeast-1M | s.cerevisiae | ribo | 1 | 1 000 000 | TBD | 8 | 1 | TBD | TBD | TBD | 0.7.1 | TBD | TBD | TBD | TBD | TBD | subsampled |
| TACO1 WT/KO | h.sapiens | ribo | 6 | TBD | TBD | 16 | 4 | TBD | TBD | TBD | 0.7.1 | TBD | TBD | TBD | TBD | TBD | full reads, biological replicates |
| current-lab-realistic | h.sapiens | ribo+rna | TBD | TBD | TBD | 16 | 4 | TBD | TBD | TBD | 0.7.1 | TBD | TBD | TBD | TBD | TBD | non-subsampled, end-to-end |

> Numbers will land in this table when the maintainer's HPC
> benchmark batch completes. The protocol above is what reviewers
> need to know **how** they were produced; the numbers themselves
> are the publication artefact.

---

## See also

* [`docs/reference/output_schema.md`](reference/output_schema.md)
  — what every benchmarked output file actually contains.
* [`docs/developer/release_checklist.md`](developer/release_checklist.md)
  — release-day gates including the link-check and manuscript
  version-pin sections.
* [`docs/tutorials/05_hpc_cluster_run.md`](tutorials/05_hpc_cluster_run.md)
  — HPC scheduler integration; the same recipe drives the
  benchmark cases.
