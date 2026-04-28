#!/usr/bin/env bash
# MitoRiboPy -- exhaustive shell-script template (end-to-end).
# Compatible with: MitoRiboPy 0.5.0+
#
# Calls each subcommand (align -> rpf -> optional rnaseq) directly with
# every available flag spelled out, so you can see and edit exactly
# what is being passed. Mirrors the YAML template at
# ./pipeline_config.example.yaml -- you only need ONE of the two:
#
#   * Use the YAML if you prefer a single declarative file managed
#     under version control:
#         mitoribopy all --config pipeline_config.example.yaml \
#                        --output results/
#
#   * Use this shell script if you prefer per-stage commands you can
#     tweak interactively or split into separate cluster jobs.
#
# Quick start:
#   1. Edit the variables in the "Edit me" block below.
#   2. chmod +x run_pipeline.example.sh
#   3. ./run_pipeline.example.sh
#
# Default values match `mitoribopy <stage> --help` defaults, so you can
# delete any line you do not need to override -- the CLI will fall back
# to the documented default. Lines starting with `#` are commented-out
# optional flags; uncomment to enable.
#
# Strict bash so a single failed command stops the pipeline instead of
# silently moving on. Drop `pipefail` if you intentionally want partial
# runs.
set -euo pipefail

# ============================================================================
# Edit me
# ============================================================================
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Inputs (REQUIRED for align)
FASTQ_DIR="${PROJECT_ROOT}/input_data/seq"
CONTAM_INDEX="${PROJECT_ROOT}/input_data/indexes/rrna_contam"
MT_INDEX="${PROJECT_ROOT}/input_data/indexes/mt_tx"

# Reference (REQUIRED for rpf)
MT_FASTA="${PROJECT_ROOT}/input_data/human-mt-mRNA.fasta"

# Top-level output directory; per-stage subdirs are created under it.
OUTPUT_DIR="${PROJECT_ROOT}/results"

# Optional rnaseq inputs (only needed when the rnaseq stage runs).
# `mitoribopy rnaseq` has TWO mutually-exclusive flows:
#   * DEFAULT (from-FASTQ) -- pass raw RNA-seq + Ribo-seq FASTQs and a
#     transcriptome FASTA; the subcommand runs cutadapt + bowtie2 +
#     pyDESeq2 itself. Recommended starting point. Requires the
#     `[fastq]` extra: pip install 'mitoribopy[fastq]'.
#   * ALTERNATIVE (--de-table) -- bring an existing DE table from an
#     external DESeq2 / Xtail / Anota2Seq run. SHA256 reference-
#     consistency gate is enforced against the rpf run above. Use this
#     for publication-grade DE statistics over the full transcriptome.
RNASEQ_FLOW="default"      # "default" (from-FASTQ) | "de-table" (alternative)

# Default-flow inputs:
RNA_FASTQ_DIR="${PROJECT_ROOT}/input_data/rna_seq"
RIBO_FASTQ_DIR="${PROJECT_ROOT}/input_data/ribo_seq"
RNASEQ_ALIGN_THREADS=4

# Alternative-flow inputs:
DE_TABLE="${PROJECT_ROOT}/de.tsv"

# Both flows:
CONDITION_MAP="${PROJECT_ROOT}/samples.tsv"
CONDITION_A="control"
CONDITION_B="knockdown"

# Compute knobs.
THREADS=8                 # global thread budget for cutadapt + bowtie2.
MAX_PARALLEL_SAMPLES=4    # align-only: >1 runs samples concurrently; THREADS
                          # is auto-divided across workers (each tool gets
                          # max(1, THREADS // MAX_PARALLEL_SAMPLES) threads).
                          # Set to 1 for the pre-parallelism serial behaviour.
                          # The joint `mitoribopy rpf` stage is unaffected.
LOG_LEVEL=INFO            # DEBUG | INFO | WARNING | ERROR
RUN_RNASEQ=false          # set true to run the optional translation-efficiency stage


# ============================================================================
# Stage 1: align  --  FASTQ -> BED6 + read_counts.tsv
# ============================================================================
ALIGN_OPTS=(
  # --- Inputs / outputs --------------------------------------------------
  --fastq-dir "${FASTQ_DIR}"
  --contam-index "${CONTAM_INDEX}"
  --mt-index "${MT_INDEX}"
  --output "${OUTPUT_DIR}/align"
  # Add per-FASTQ paths as needed:
  # --fastq /path/to/extra_sample.fq.gz
  # --fastq /path/to/another.fq.gz

  # --- Library prep / kit ------------------------------------------------
  --kit-preset auto                # auto | illumina_smallrna | illumina_truseq |
                                   # illumina_truseq_umi | qiaseq_mirna |
                                   # pretrimmed | custom
  # --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA   # required when kit_preset=custom
  # --umi-length 8                 # global UMI override
  # --umi-position 5p              # 5p | 3p
  # --sample-overrides "${PROJECT_ROOT}/per_sample_overrides.tsv"
                                   # for mixed-UMI batches; columns: sample,
                                   # kit_preset, adapter, umi_length,
                                   # umi_position, dedup_strategy

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
                                   # ThreadPoolExecutor. THREADS is divided
                                   # across workers, so total CPU stays
                                   # around ${THREADS} regardless of N.
                                   # Resume-cached samples skip the pool.
                                   # Default 1 (serial). The joint rpf stage
                                   # is unaffected by this flag.

  # --- Shared --------------------------------------------------------------
  --threads "${THREADS}"
  --log-level "${LOG_LEVEL}"
)

mitoribopy align "${ALIGN_OPTS[@]}"


# ============================================================================
# Stage 2: rpf  --  BED -> offsets + translation profile + coverage plots
# ============================================================================
RPF_OPTS=(
  # --- Core inputs / outputs ---------------------------------------------
  --fasta "${MT_FASTA}"
  --directory "${OUTPUT_DIR}/align/bed"             # auto-wired when running via `mitoribopy all`
  --read_counts_file "${OUTPUT_DIR}/align/read_counts.tsv"
  --output "${OUTPUT_DIR}/rpf"
  --bam_mapq 10                    # MAPQ threshold for BAM inputs (BAM->BED6)

  # --- Organism / strain -------------------------------------------------
  --strain h.sapiens               # h.sapiens | s.cerevisiae | custom
                                   # Synonyms `h` and `y` are accepted.
  # --annotation_file "${PROJECT_ROOT}/my_annotation.csv"
                                   # REQUIRED for --strain custom
  # --codon_tables_file "${PROJECT_ROOT}/my_codons.json"
                                   # only when your code is not in the bundled NCBI list
  # --codon_table_name vertebrate_mitochondrial
                                   # NCBI Genetic Code name (e.g. mold_mitochondrial,
                                   # invertebrate_mitochondrial, ...)
  # --start_codons ATG ATA         # override strain default
  --atp8_atp6_baseline ATP6        # ATP6 | ATP8
  --nd4l_nd4_baseline ND4          # ND4  | ND4L

  # --- RPF length window -------------------------------------------------
  --footprint_class monosome       # short | monosome | disome | custom
                                   # short:    h.sapiens / s.cerevisiae 16-24 nt
                                   # monosome: h.sapiens 28-34, s.cerevisiae 37-41
                                   # disome:   h.sapiens 50-70, s.cerevisiae 60-90
  -rpf 29 34                       # explicit [min, max] (overrides footprint_class)
  --unfiltered_read_length_range 15 50

  # --- Offset enrichment / selection -------------------------------------
  --align stop                     # start | stop
  --range 20                       # plot offsets from -range to +range
  --min_offset 11                  # shared min absolute offset (legacy fallback)
  --max_offset 20                  # shared max absolute offset (legacy fallback)
  --rpf_min_count_frac 0.20        # auto-prune read-length bins below 20% of the
                                   # most-enriched length's count (set 0 to disable).
  --min_5_offset 10                # end-specific 5' min (preferred)
  --max_5_offset 22                # end-specific 5' max
  --min_3_offset 10                # end-specific 3' min
  --max_3_offset 22                # end-specific 3' max
  --offset_mask_nt 5               # mask -N..-1 and +1..+N near the anchor codon
  --offset_pick_reference p_site   # p_site | reported_site
  --offset_type 5                  # 5 | 3 -- which read end the offset is measured from
  --offset_site p                  # p | a -- coordinate space of the SELECTED OFFSETS table
  --analysis_sites both            # p | a | both
  --codon_overlap_mode full        # full | any
  # --psite_offset 12              # fixed offset for every read length (bypass enrichment)
  --offset_mode per_sample         # per_sample | combined

  # --- Outputs / plotting ------------------------------------------------
  --downstream_dir footprint_density
  --plot_dir offset_diagnostics    # CSVs go to <plot_dir>/csv/, plots to <plot_dir>/plots/.
                                   # Per-sample offsets at <plot_dir>/csv/per_sample_offset/.
  --plot_format svg                # png | pdf | svg
  --line_plot_style combined       # combined | separate
  --cap_percentile 0.999
  --codon_density_window           # smooth codon-density with +/-1 nt window
  # --x_breaks 10 20 30            # custom x-axis tick marks
  # --order_samples ctrl_1 ctrl_2 kd_1 kd_2

  # --- Read-count normalization -----------------------------------------
  --rpm_norm_mode total            # total | mt_mrna
  --mt_mrna_substring_patterns mt_genome mt-mrna mt_mrna
  # --read_counts_sample_col Sample
  # --read_counts_reads_col Reads
  # --read_counts_reference_col Reference

  # --- Optional modules --------------------------------------------------
  # --structure_density
  # --structure_density_norm_perc 0.99
  --cor_plot                       # publication-quality codon-correlation scatter,
                                   # one folder per requested site (P, A) when both
                                   # are generated.
  --base_sample WT_R1              # required by --cor_plot; change to YOUR reference.
  # --cor_mask_method percentile   # percentile | fixed | none
  # --cor_mask_percentile 0.99
  # --cor_mask_threshold 1.5

  # --- Read-coverage outputs (default: both raw and rpm) ------------------
  --no-read_coverage_raw           # skip read_coverage_raw[_codon]/.
  --read_coverage_rpm              # write read_coverage_rpm[_codon]/.
  --igv_export                     # write igv_tracks/<sample>/<sample>_{p,a}_site.bedgraph.

  # --- Shared --------------------------------------------------------------
  --threads "${THREADS}"
  --log-level "${LOG_LEVEL}"
)

mitoribopy rpf "${RPF_OPTS[@]}"


# ============================================================================
# Stage 3: rnaseq (optional)  --  DE table + rpf -> TE / dTE
# ============================================================================
if [[ "${RUN_RNASEQ}" == "true" ]]; then
  RNASEQ_COMMON=(
    # --- Gene identifier convention (REQUIRED, no default) -------------
    --gene-id-convention hgnc      # ensembl | refseq | hgnc | mt_prefixed | bare
    --organism h.sapiens           # h.sapiens | s.cerevisiae

    # --- Conditions (REQUIRED in Mode B; required for replicate ΔTE in A) ---
    --condition-map "${CONDITION_MAP}"
    --condition-a "${CONDITION_A}"
    --condition-b "${CONDITION_B}"

    # --- Output ---------------------------------------------------------
    --output "${OUTPUT_DIR}/rnaseq"

    # --- Shared ---------------------------------------------------------
    --threads "${THREADS}"
    --log-level "${LOG_LEVEL}"
  )

  if [[ "${RNASEQ_FLOW}" == "default" ]]; then
    RNASEQ_OPTS=(
      "${RNASEQ_COMMON[@]}"

      # --- From-FASTQ (default; mutually exclusive with --de-table) ----
      --rna-fastq "${RNA_FASTQ_DIR}"
      --ribo-fastq "${RIBO_FASTQ_DIR}"
      --reference-fasta "${MT_FASTA}"
      # --bowtie2-index "${PROJECT_ROOT}/cache/transcriptome_5ca397907373"
                                   # pre-built prefix; skips bowtie2-build
      --workdir "${OUTPUT_DIR}/rnaseq/work"
      --align-threads "${RNASEQ_ALIGN_THREADS}"
      # --no-auto-pseudo-replicate
                                   # default OFF; when set, n=1
                                   # conditions are NOT split into
                                   # rep1/rep2 and the run will fail at
                                   # pyDESeq2 dispersion.
    )
  elif [[ "${RNASEQ_FLOW}" == "de-table" ]]; then
    RNASEQ_OPTS=(
      "${RNASEQ_COMMON[@]}"

      # --- DE table (alternative flow) ---------------------------------
      --de-table "${DE_TABLE}"
      --de-format auto             # auto | deseq2 | xtail | anota2seq | custom
      # --de-gene-col gene_id      # only when --de-format custom
      # --de-log2fc-col log2FoldChange
      # --de-padj-col padj
      # --de-basemean-col baseMean

      # --- Ribo-seq inputs (alternative flow) --------------------------
      --ribo-dir "${OUTPUT_DIR}/rpf"
      # --ribo-counts "${OUTPUT_DIR}/rpf/rpf_counts.tsv"   # explicit override

      # --- Reference-consistency gate (alternative flow; exactly one) -
      --reference-gtf "${MT_FASTA}"
      # --reference-checksum <sha256>
    )
  else
    echo "RNASEQ_FLOW must be 'default' or 'de-table' (got: ${RNASEQ_FLOW})" >&2
    exit 2
  fi

  mitoribopy rnaseq "${RNASEQ_OPTS[@]}"
fi


# ============================================================================
# Done
# ============================================================================
echo
echo "MitoRiboPy pipeline finished."
echo "Inspect outputs under: ${OUTPUT_DIR}/"
echo
echo "Key files to look at first:"
echo "  ${OUTPUT_DIR}/align/kit_resolution.tsv      # per-sample kit + dedup decisions"
echo "  ${OUTPUT_DIR}/align/read_counts.tsv         # per-stage drop-off invariants"
echo "  ${OUTPUT_DIR}/align/mitoribopy.log          # per-stage timing per sample +"
echo "                                              # end-of-run [ALIGN] Timing summary"
echo "  ${OUTPUT_DIR}/rpf/mitoribopy.log            # per-step timing for rpf +"
echo "                                              # end-of-run [PIPELINE] Timing summary"
echo "  ${OUTPUT_DIR}/rpf/offset_diagnostics/plots/offset_drift_*.svg"
echo "  ${OUTPUT_DIR}/rpf/offset_diagnostics/csv/per_sample_offset/<sample>/offset_applied.csv"
echo "                                              # exact offsets row applied per sample"
echo "  ${OUTPUT_DIR}/rpf/coverage_profile_plots/p_site_density_rpm_frame/  # frame-coloured QC"
echo "  ${OUTPUT_DIR}/rpf/translation_profile/<sample>/codon_usage/p_site_codon_usage_total.csv"
echo "  ${OUTPUT_DIR}/rpf/codon_correlation/{p_site,a_site}/<base>_vs_<sample>_*.{svg,png}"
echo "  ${OUTPUT_DIR}/rpf/igv_tracks/<sample>/<sample>_{p_site,a_site}.bedgraph  # IGV tracks"
if [[ "${RUN_RNASEQ}" == "true" ]]; then
echo "  ${OUTPUT_DIR}/rnaseq/te.tsv  +  delta_te.tsv  +  plots/"
fi
