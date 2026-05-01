# Running MitoRiboPy on HPC clusters (SLURM / LSF)

This guide covers the practical knobs that matter when scaling
`mitoribopy all` from a laptop to a shared cluster: thread accounting,
parallelism, scratch storage, and disk budgeting.

---

## 1. Thread accounting

MitoRiboPy has **two** thread settings that control different things:

| Flag | Effect |
|---|---|
| `--threads N` | Caps BLAS / pyDESeq2 / OpenMP pool size. Forwarded to `OMP_NUM_THREADS` and `OPENBLAS_NUM_THREADS`. |
| `align.max_parallel_samples M` | Number of samples processed in parallel during the align stage. |

The total CPU budget is roughly `M × (N // M)` — every per-sample
worker gets `max(1, N // M)` threads. Set both, with the cluster
allocation in mind:

```yaml
align:
  max_parallel_samples: 4   # 4 workers (1 per sample)
                            # cluster job below requests 16 cores
```

```bash
# SLURM job request: 16 cores -> each worker gets 4 threads
srun --cpus-per-task=16 \
    mitoribopy all \
    --config pipeline_config.yaml \
    --output /scratch/$USER/mitoribopy/run42 \
    --threads 16
```

The `rpf` and `rnaseq` stages are single-process; `--threads N` only
caps their BLAS / NumPy pool.

---

## 2. Local scratch vs network storage

* **Inputs (FASTQs, references)** — staged once on shared storage.
* **Outputs (`<run>/align/`, `<run>/rpf/`, `<run>/rnaseq/`)** — write
  to **local scratch first** (e.g. `/scratch/$USER/`), then `rsync`
  the final results back to lab storage. Intermediate files
  (`trimmed/`, `contam_filtered/`, `aligned/`) are deleted by the
  orchestrator unless `--keep-intermediates` is set, but the
  per-sample BAM files in `aligned/<sample>.mapq.bam` and the
  per-sample BED files in `bed/` are kept and can be GBs.
* **Bowtie2 indexes** — keep in shared storage; bowtie2 mmap's the
  index, so accessing it from network storage costs only a one-time
  page-in per worker.

The `mitoribopy rnaseq` from-FASTQ flow caches its bowtie2 index at
`<workdir>/bt2_cache/transcriptome_<sha12>`; pass `--workdir` to point
this at local scratch when the same FASTA is reused across runs.

---

## 3. Disk budget

Approximate per-sample on-disk footprint:

| Stage | Output | Size for ~20M-read RPF library |
|---|---|---|
| align | `aligned/<sample>.mapq.bam` | ~150 MB |
| align | `bed/<sample>.bed` | ~80 MB |
| align | `trimmed/<sample>.cutadapt.json` | ~5 KB |
| rpf | `rpf_counts.tsv` (long) | ~10 KB / sample |
| rpf | per-sample BED windows | ~50 MB |

The `mitoribopy benchmark --subsample N` mode is the recommended way
to estimate the real budget for your data — it runs the full pipeline
on a reservoir-sampled FASTQ and writes `benchmark.tsv` with wall
time, peak RSS, and the actual on-disk footprint.

```bash
mitoribopy benchmark \
    --config pipeline_config.yaml \
    --output /scratch/$USER/bench \
    --subsample 200000 \
    --threads 8
cat /scratch/$USER/bench/benchmark_summary.md
```

---

## 4. SLURM example

See `examples/templates/run_pipeline.slurm.sh` for a complete script.

```bash
#!/usr/bin/env bash
#SBATCH --job-name=mitoribopy
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=mitoribopy.%j.out

set -euo pipefail

RUN_ROOT=/scratch/$USER/mitoribopy/$SLURM_JOB_ID
mkdir -p "$RUN_ROOT"

# Validate before launching the long-running stages.
mitoribopy validate-config pipeline_config.yaml --strict

mitoribopy all \
    --config pipeline_config.yaml \
    --output "$RUN_ROOT" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --resume      # safe: hash-validated against the prior manifest

# Resume after a crash:
#   sbatch run_pipeline.slurm.sh
# The orchestrator detects existing sentinel files + matches the
# config / sample-sheet / reference SHA256 from the prior run's
# run_manifest.json, then skips completed stages.

# Optional: ship the run back to shared storage.
rsync -a "$RUN_ROOT"/ "/lab_storage/$USER/mitoribopy/$SLURM_JOB_ID/"
```

---

## 5. LSF example

See `examples/templates/run_pipeline.lsf.sh`.

```bash
#!/usr/bin/env bash
#BSUB -J mitoribopy
#BSUB -n 16
#BSUB -M 32000
#BSUB -W 06:00
#BSUB -o mitoribopy.%J.out

set -euo pipefail

RUN_ROOT=/scratch/$USER/mitoribopy/$LSB_JOBID
mkdir -p "$RUN_ROOT"

mitoribopy validate-config pipeline_config.yaml --strict
mitoribopy all \
    --config pipeline_config.yaml \
    --output "$RUN_ROOT" \
    --threads "$LSB_DJOB_NUMPROC" \
    --resume

rsync -a "$RUN_ROOT"/ "/lab_storage/$USER/mitoribopy/$LSB_JOBID/"
```

---

## 6. Resume on the cluster

`mitoribopy all --resume` is a hash-validated resume: the orchestrator
re-hashes the current config / sample-sheet / reference and compares
them to the values recorded in the prior run's `run_manifest.json`.
On any mismatch, the resume FAILS rather than re-using stale outputs:

```
[mitoribopy all] ERROR: --resume cannot proceed: resume hash mismatch
(1 field(s) differ from the prior run_manifest.json):
  - config_source_sha256: prior='abc...', current='def...'
[mitoribopy all] HINT: re-run without --resume to produce a fresh
result, or pass --force-resume to bypass the hash guard ...
```

For CI scripts that intentionally re-run with edited configs, set the
env var `MITORIBOPY_FORCE_RESUME=1` (equivalent to `--force-resume`)
without changing the argv.

---

## 7. Common pitfalls

* **`--threads` ignored for cutadapt / bowtie2.** These are external
  tools; the orchestrator forwards `--threads` to them via
  `align.max_parallel_samples`'s per-worker share. Setting a global
  `--threads 16` with `max_parallel_samples: 1` means each tool gets
  16 threads sequentially.
* **OOM on the rpf stage.** `mitoribopy rpf` loads BED data per
  sample; if your samples are very deep (>50M reads each), bump the
  cluster memory request and consider `--max_parallel_samples 1` (no
  parallelism within rpf, but the next sample only starts after the
  previous frees memory).
* **Mixed-UMI batches**: declare `umi_length` / `umi_position` per
  sample in the sample sheet. Without an explicit declaration, the
  rnaseq from-FASTQ path *infers* UMI length by entropy and emits a
  `UMI_INFERRED_NO_DECLARATION` warning into `warnings.tsv`.

See also:
* `docs/reference/cli.md` — full CLI flag reference.
* `docs/tutorials/01_end_to_end_fastq.md` — single-machine quick start.
