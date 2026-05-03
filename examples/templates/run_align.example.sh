#!/usr/bin/env bash
# MitoRiboPy align -- exhaustive shell-script template (align stage only).
# Compatible with: MitoRiboPy 0.7.1+
#
# Runs only the align stage:
#
#   FASTQ -> cutadapt trim (+ UMI extraction into QNAME)
#         -> bowtie2 contam subtract
#         -> bowtie2 mt-transcriptome align
#         -> MAPQ filter
#         -> umi-tools dedup (or skip)
#         -> BAM -> BED6
#
# Output: <OUTPUT_DIR>/{bed/, kit_resolution.tsv, read_counts.tsv,
# run_settings.json, .sample_done/, ...}. The BED and read_counts.tsv
# are the inputs to `mitoribopy rpf` (see run_rpf.example.sh).
#
# Quick start:
#   1. Edit the variables in the "Edit me" block below.
#   2. chmod +x run_align.example.sh
#   3. ./run_align.example.sh
#
# Default values match `mitoribopy align --help` defaults; delete any
# line you do not need to override. Lines starting with `#` are
# commented-out optional flags; uncomment to enable.
set -euo pipefail

# ============================================================================
# Edit me
# ============================================================================
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Inputs (REQUIRED)
FASTQ_DIR="${PROJECT_ROOT}/input_data/seq"
CONTAM_INDEX="${PROJECT_ROOT}/input_data/indexes/rrna_contam"
MT_INDEX="${PROJECT_ROOT}/input_data/indexes/mt_tx"

# Output directory.
OUTPUT_DIR="${PROJECT_ROOT}/results/align"

# Compute knobs.
THREADS=8                 # global thread budget for cutadapt + bowtie2.
MAX_PARALLEL_SAMPLES=4    # >1 runs samples concurrently in a thread pool;
                          # THREADS is auto-divided across workers, so each
                          # tool sees max(1, THREADS // MAX_PARALLEL_SAMPLES)
                          # threads (here: 4 workers x 2 threads = ~8 cores).
                          # Set to 1 for the original serial behaviour.
LOG_LEVEL=INFO            # DEBUG | INFO | WARNING | ERROR


# ============================================================================
# align  --  FASTQ -> BED6 + read_counts.tsv
# ============================================================================
ALIGN_OPTS=(
  # --- Inputs / outputs --------------------------------------------------
  --fastq-dir "${FASTQ_DIR}"
  --contam-index "${CONTAM_INDEX}"
  --mt-index "${MT_INDEX}"
  --output "${OUTPUT_DIR}"
  # Add per-FASTQ paths as needed:
  # --fastq /path/to/extra_sample.fq.gz
  # --fastq /path/to/another.fq.gz

  # --- Library prep / adapter --------------------------------------------
  # Auto-detection is the default; pin --adapter only when detection
  # cannot identify your library. --pretrimmed declares already-trimmed
  # FASTQs and is mutually exclusive with --adapter.
  # --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
  # --pretrimmed                   # already-trimmed FASTQs
  # --umi-length 8                 # global UMI override
  # --umi-position 5p              # 5p | 3p | both
  # --sample-overrides "${PROJECT_ROOT}/per_sample_overrides.tsv"
                                   # for mixed-UMI / mixed-adapter batches;
                                   # columns: sample, adapter, pretrimmed,
                                   # umi_length, umi_position, dedup_strategy

  # --- Adapter detection -------------------------------------------------
  --adapter-detection auto         # auto | off | strict
  --adapter-detect-reads 5000      # head-of-FASTQ scan size
  --adapter-detect-min-rate 0.30   # min match rate to call a kit detected
  --adapter-detect-min-len 12      # adapter prefix length used as the scan needle
  --adapter-detect-pretrimmed-threshold 0.05
                                   # all-kits below this -> classify as pretrimmed
  # --no-pretrimmed-inference      # disable the pretrimmed auto-fallback

  # --- Trim / align ------------------------------------------------------
  --library-strandedness forward   # forward | reverse | unstranded
  --min-length 15                  # cutadapt --minimum-length
  --max-length 45                  # cutadapt --maximum-length
  --quality 20                     # cutadapt 3' quality trim threshold
  --mapq 10                        # MAPQ threshold (NUMT suppression)
  --seed 42                        # bowtie2 --seed (deterministic output)

  # --- Deduplication -----------------------------------------------------
  --dedup-strategy auto            # auto | umi-tools | skip
                                   # The legacy 'mark-duplicates' option was removed
                                   # in v0.4.5 because coordinate-only dedup destroys
                                   # codon-occupancy signal on mt-Ribo-seq.
  --umi-dedup-method unique        # unique | percentile | cluster | adjacency | directional

  # --- Intermediate files / resume --------------------------------------
  # --keep-intermediates           # default off; keep trimmed.fq.gz / nocontam.fq.gz / pre-MAPQ .bam
  # --resume                       # skip samples whose .sample_done/<sample>.json marker exists

  # --- Execution / concurrency -----------------------------------------
  --max-parallel-samples "${MAX_PARALLEL_SAMPLES}"
                                   # Run multiple samples concurrently via a
                                   # ThreadPoolExecutor. Per-sample work is
                                   # independent (cutadapt -> bowtie2 ->
                                   # MAPQ -> dedup -> BAM->BED). Resume-cached
                                   # samples short-circuit the pool. Failures
                                   # are fail-fast: pending futures are
                                   # cancelled and the first exception is
                                   # re-raised. Default 1 (serial).

  # --- Shared --------------------------------------------------------------
  --threads "${THREADS}"
  --log-level "${LOG_LEVEL}"
)

mitoribopy align "${ALIGN_OPTS[@]}"


# ============================================================================
# Done
# ============================================================================
echo
echo "MitoRiboPy align finished."
echo "Inspect outputs under: ${OUTPUT_DIR}/"
echo
echo "Key files to look at first:"
echo "  ${OUTPUT_DIR}/kit_resolution.tsv      # per-sample kit + dedup decisions"
echo "  ${OUTPUT_DIR}/read_counts.tsv         # per-stage drop-off invariants:"
echo "                                        #   rrna_aligned + post_rrna_filter == post_trim"
echo "                                        #   mt_aligned   + unaligned_to_mt == post_rrna_filter"
echo "  ${OUTPUT_DIR}/bed/<sample>.bed        # strand-aware BED6 -- input to mitoribopy rpf"
echo "  ${OUTPUT_DIR}/run_settings.json       # full per-sample resolution + tool versions"
echo "  ${OUTPUT_DIR}/mitoribopy.log          # per-stage timing per sample, plus the"
echo "                                        # end-of-run [ALIGN] Timing summary table"
echo "                                        # (total / mean / max per stage + wall)."
