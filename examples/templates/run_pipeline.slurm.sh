#!/usr/bin/env bash
#SBATCH --job-name=mitoribopy
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=mitoribopy.%j.out

# MitoRiboPy on SLURM. See docs/tutorials/05_hpc_cluster_run.md for the
# full guide. Edit pipeline_config.yaml + the SBATCH directives above
# for your cluster, then `sbatch run_pipeline.slurm.sh`.

set -euo pipefail

CONFIG="${CONFIG:-pipeline_config.yaml}"
RUN_ROOT="${RUN_ROOT:-/scratch/$USER/mitoribopy/$SLURM_JOB_ID}"
mkdir -p "$RUN_ROOT"

echo "[runner] config=$CONFIG output=$RUN_ROOT cpus=$SLURM_CPUS_PER_TASK"

# Pre-flight: parse + canonicalise + path-check the config so a typo
# fails in seconds, not after hours of bowtie2 work.
mitoribopy validate-config "$CONFIG"

# Run with hash-validated resume. Edits to the config / sample sheet /
# reference between runs cause this to fail loudly (rather than silently
# re-using stale outputs). Pass --force-resume only when you know the
# stale outputs are still valid.
mitoribopy all \
    --config "$CONFIG" \
    --output "$RUN_ROOT" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --resume

# Optional: ship results back to lab storage.
# rsync -a "$RUN_ROOT"/ "/lab_storage/$USER/mitoribopy/$SLURM_JOB_ID/"
