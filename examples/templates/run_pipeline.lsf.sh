#!/usr/bin/env bash
#BSUB -J mitoribopy
#BSUB -n 16
#BSUB -M 32000
#BSUB -W 06:00
#BSUB -o mitoribopy.%J.out

# MitoRiboPy on LSF. See docs/tutorials/05_hpc_cluster_run.md.
# Adjust pipeline_config.yaml + the BSUB directives, then
# `bsub < run_pipeline.lsf.sh`.

set -euo pipefail

CONFIG="${CONFIG:-pipeline_config.yaml}"
RUN_ROOT="${RUN_ROOT:-/scratch/$USER/mitoribopy/$LSB_JOBID}"
mkdir -p "$RUN_ROOT"

echo "[runner] config=$CONFIG output=$RUN_ROOT cpus=$LSB_DJOB_NUMPROC"

mitoribopy validate-config "$CONFIG"

mitoribopy all \
    --config "$CONFIG" \
    --output "$RUN_ROOT" \
    --threads "$LSB_DJOB_NUMPROC" \
    --resume

# Optional: ship results back to lab storage.
# rsync -a "$RUN_ROOT"/ "/lab_storage/$USER/mitoribopy/$LSB_JOBID/"
