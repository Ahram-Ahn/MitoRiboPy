# Benchmarking Protocol

`mitoribopy benchmark` wraps a real `mitoribopy all` run, records
runtime and resource totals, and writes benchmark artifacts into the
same run root. Use it for local cluster sizing and for producing
manuscript runtime tables from measured runs.

This page documents the implemented command and output schema. It does
not contain placeholder benchmark numbers.

---

## Command

```bash
mitoribopy benchmark \
  --config <pipeline_config.yaml> \
  --output <run_root> \
  --threads <N> \
  --subsample <K>     # optional
```

The wrapper:

* optionally reservoir-samples explicit FASTQ inputs to `K` reads;
* writes the rewritten config to `<run_root>/benchmark_config.yaml`;
* runs `mitoribopy all --config <run_root>/benchmark_config.yaml`;
* reads stage runtimes back from `run_manifest.json`;
* writes `<run_root>/benchmark.tsv` and
  `<run_root>/benchmark_summary.md`.

`--subsample` currently rewrites explicit FASTQ paths in `align.fastq`,
`align.fastq_dir`, `rnaseq.rna_fastq`, and `rnaseq.ribo_fastq`.
Top-level `samples:` sheets are not rewritten yet; use an explicit
FASTQ config for subsampled tuning runs.

---

## `benchmark.tsv`

The file starts with `# schema_version: 1.0`, followed by these columns:

| Column | Source | Notes |
|---|---|---|
| `stage` | benchmark wrapper | `align`, `rpf`, `rnaseq`, then `total` |
| `status` | `run_manifest.json` for stage rows; wrapper exit for `total` | `not_configured` when a stage did not run |
| `wall_time_sec` | stage manifest or wrapper timer | stage rows use `manifest.stages.<stage>.runtime_seconds`; total row uses wall-clock wrapper time |
| `cumulative_wall_sec` | benchmark wrapper | running sum for stage rows; total wall for `total` |
| `max_rss_mb_total` | `resource.getrusage` | populated on the `total` row only |
| `disk_mb` | output directory walk | populated on the `total` row only |
| `threads` | CLI argument | empty when `--threads` was omitted |
| `subsample_reads` | CLI argument | empty when `--subsample` was omitted |

`benchmark_summary.md` mirrors the same information in a short
human-readable report.

---

## Manuscript Table

For a publication table, derive one row per measured run from
`benchmark.tsv`, `run_manifest.json`, `resource_plan.json`,
`align/read_counts.tsv`, and the original input-file inventory.

Suggested columns:

| Column | Units | Source |
|---|---|---|
| `dataset` | text | hand-filled dataset label |
| `organism` | text | config / methods |
| `assay` | text | `ribo`, `rna`, or `ribo+rna` |
| `n_samples` | count | sample sheet or active FASTQ list |
| `reads_per_sample` | reads | median `total_reads` from `align/read_counts.tsv` |
| `input_size_gb` | GB | compressed input FASTQ inventory |
| `threads` | count | `benchmark.tsv` |
| `parallel_samples` | count | `resource_plan.json` |
| `wall_time_min` | minutes | `benchmark.tsv` total `wall_time_sec` / 60 |
| `peak_rss_gb` | GB | `benchmark.tsv` total `max_rss_mb_total` / 1024 |
| `disk_gb` | GB | `benchmark.tsv` total `disk_mb` / 1024 |
| `mitoribopy_version` | version | `run_manifest.json` / `mitoribopy --version` |
| `commit_sha` | SHA | `run_manifest.json` |
| `cutadapt_version` | version | `run_manifest.json` `tool_versions` |
| `bowtie2_version` | version | `run_manifest.json` `tool_versions` |
| `umi_tools_version` | version | `run_manifest.json` `tool_versions`, empty when unused |
| `pysam_version` | version | `run_manifest.json` `tool_versions` |
| `notes` | text | subsampling, filesystem class, rerun/cache caveats |

Do not mix subsampled and non-subsampled rows without stating the
subsampling level in `notes`.

---

## Reference Cases To Measure

The Methods section should include measured rows for:

1. `examples/smoke/` synthetic fixture.
2. One public human mt-Ribo-seq sample subsampled to 1M reads.
3. One public yeast mt-Ribo-seq sample subsampled to 1M reads.
4. TACO1 WT/KO biological replicates at full depth.
5. One realistic non-subsampled lab dataset, with RNA-seq included if
   that workflow is being claimed.

Keep the committed docs number-free until those rows come from real
run artifacts.

---

## See Also

* [`reference/output_schema.md`](reference/output_schema.md) — output
  table schemas and provenance files.
* [`developer/roadmap.md`](developer/roadmap.md) — remaining benchmark
  documentation work.
* [`tutorials/05_hpc_cluster_run.md`](tutorials/05_hpc_cluster_run.md)
  — scheduler and scratch-space guidance.
