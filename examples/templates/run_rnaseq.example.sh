#!/usr/bin/env bash
# MitoRiboPy rnaseq -- exhaustive shell-script template (rnaseq stage only).
# Compatible with: MitoRiboPy 0.5.0+
#
# `mitoribopy rnaseq` has TWO mutually-exclusive modes, picked by which
# input flag you set:
#
#   * Mode A (PRE-COMPUTED DE) -- you already ran DESeq2 / Xtail /
#     Anota2Seq externally on the full transcriptome and have the
#     results table. Pass `--de-table`. SHA256 reference-consistency
#     gate is enforced against a prior `mitoribopy rpf` run.
#
#   * Mode B (FROM-FASTQ, added in v0.5.0) -- pass raw RNA-seq +
#     Ribo-seq FASTQs and a transcriptome FASTA. The subcommand auto-
#     detects SE vs PE from filename mate tokens, runs cutadapt +
#     bowtie2 per sample, counts reads per transcript, runs pyDESeq2,
#     then falls through into the same TE / ΔTE / plot path as Mode A.
#     Requires the `[fastq]` extra: pip install 'mitoribopy[fastq]'.
#
# Quick start:
#   1. Set MODE below (a or b) and edit the variables in its block.
#   2. chmod +x run_rnaseq.example.sh
#   3. ./run_rnaseq.example.sh
#
# Default values match `mitoribopy rnaseq --help` defaults; delete any
# line you do not need to override. Lines starting with `#` are
# commented-out optional flags; uncomment to enable.
set -euo pipefail

# ============================================================================
# Edit me
# ============================================================================
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Pick the mode: "a" = pre-computed DE table, "b" = from-FASTQ.
MODE=a

# Common knobs.
THREADS=8                 # global thread budget (BLAS / pyDESeq2 caps)
LOG_LEVEL=INFO            # DEBUG | INFO | WARNING | ERROR
OUTPUT_DIR="${PROJECT_ROOT}/results/rnaseq"

# Reference (required by both modes).
MT_FASTA="${PROJECT_ROOT}/input_data/human-mt-mRNA.fasta"

# Gene-ID convention used in the DE table / FASTA records (REQUIRED).
GENE_ID_CONVENTION="hgnc"      # ensembl | refseq | hgnc | mt_prefixed | bare
ORGANISM="h.sapiens"           # h.sapiens | s.cerevisiae

# Conditions (REQUIRED in Mode B; required for replicate-based ΔTE in Mode A).
CONDITION_MAP="${PROJECT_ROOT}/samples.tsv"
CONDITION_A="control"
CONDITION_B="knockdown"


# ============================================================================
# Mode A: pre-computed DE table
# ============================================================================
# Inputs (REQUIRED for Mode A):
DE_TABLE="${PROJECT_ROOT}/de.tsv"
RIBO_DIR="${PROJECT_ROOT}/results/rpf"


# ============================================================================
# Mode B: from-FASTQ
# ============================================================================
# Inputs (REQUIRED for Mode B):
RNA_FASTQ_DIR="${PROJECT_ROOT}/input_data/rna_seq"
RIBO_FASTQ_DIR="${PROJECT_ROOT}/input_data/ribo_seq"

# Per-stage compute knobs (Mode B only).
ALIGN_THREADS=4                 # threads passed to cutadapt + bowtie2
WORKDIR="${OUTPUT_DIR}/work"    # scratch (defaults to <output>/work if unset)


# ============================================================================
# Build per-mode argv
# ============================================================================
COMMON_OPTS=(
  # --- Gene identifier convention (REQUIRED, no default) -----------------
  --gene-id-convention "${GENE_ID_CONVENTION}"
  --organism "${ORGANISM}"

  # --- Output ------------------------------------------------------------
  --output "${OUTPUT_DIR}"

  # --- Conditions --------------------------------------------------------
  --condition-map "${CONDITION_MAP}"
  --condition-a "${CONDITION_A}"
  --condition-b "${CONDITION_B}"

  # --- Shared ------------------------------------------------------------
  --threads "${THREADS}"
  --log-level "${LOG_LEVEL}"
)

if [[ "${MODE}" == "a" ]]; then
  RNASEQ_OPTS=(
    "${COMMON_OPTS[@]}"

    # --- DE table (Mode A) ---------------------------------------------
    --de-table "${DE_TABLE}"
    --de-format auto              # auto | deseq2 | xtail | anota2seq | custom
    # --de-gene-col gene_id       # only when --de-format custom
    # --de-log2fc-col log2FoldChange
    # --de-padj-col padj
    # --de-basemean-col baseMean

    # --- Ribo-seq inputs (Mode A) --------------------------------------
    --ribo-dir "${RIBO_DIR}"
    # --ribo-counts "${RIBO_DIR}/rpf_counts.tsv"   # explicit override

    # --- Reference-consistency gate (Mode A; pass exactly one) --------
    --reference-gtf "${MT_FASTA}"
    # --reference-checksum 5ca397907373...        # OR pre-computed SHA256
  )
elif [[ "${MODE}" == "b" ]]; then
  RNASEQ_OPTS=(
    "${COMMON_OPTS[@]}"

    # --- From-FASTQ inputs (Mode B; mutually exclusive with --de-table) -
    --rna-fastq "${RNA_FASTQ_DIR}"
    --ribo-fastq "${RIBO_FASTQ_DIR}"
    --reference-fasta "${MT_FASTA}"
    # --bowtie2-index "${PROJECT_ROOT}/cache/transcriptome_5ca397907373"
                                  # pre-built bowtie2 prefix; skips bowtie2-build
    --workdir "${WORKDIR}"
    --align-threads "${ALIGN_THREADS}"

    # --- Pseudo-replicate auto-split (default ON) ----------------------
    # When a condition has only 1 sample, the FASTQ is auto-split by
    # record parity into rep1 / rep2 so pyDESeq2 has n>=2 for dispersion
    # estimation. The augmented condition map is written to
    # <output>/condition_map.augmented.tsv. Pass the line below to
    # disable -- the run will then fail at the DESeq2 dispersion step
    # on n=1-per-condition designs.
    # --no-auto-pseudo-replicate
  )
else
  echo "MODE must be 'a' or 'b' (got: ${MODE})" >&2
  exit 2
fi


# ============================================================================
# Run
# ============================================================================
mitoribopy rnaseq "${RNASEQ_OPTS[@]}"


# ============================================================================
# Done
# ============================================================================
echo
echo "MitoRiboPy rnaseq finished (Mode ${MODE})."
echo "Inspect outputs under: ${OUTPUT_DIR}/"
echo
echo "Key files to look at first:"
echo "  ${OUTPUT_DIR}/te.tsv                     # per-(sample,gene) TE"
echo "  ${OUTPUT_DIR}/delta_te.tsv               # per-gene ΔTE log2 + padj + note"
echo "  ${OUTPUT_DIR}/run_settings.json          # full provenance + per-sample stats"
echo "  ${OUTPUT_DIR}/plots/mrna_vs_rpf.png      # log2FC RPF vs mRNA scatter"
echo "  ${OUTPUT_DIR}/plots/delta_te_volcano.png # ΔTE volcano"
echo "  ${OUTPUT_DIR}/plots/ma.png               # DESeq2 MA plot"
echo "  ${OUTPUT_DIR}/plots/te_bar_by_condition.png  # primary biological readout"
echo "  ${OUTPUT_DIR}/plots/te_heatmap.png       # gene × sample log2(TE) heatmap"
if [[ "${MODE}" == "b" ]]; then
echo "  ${OUTPUT_DIR}/plots/sample_pca.png       # sample PCA (Mode B only)"
echo "  ${OUTPUT_DIR}/de_table.tsv               # pyDESeq2 result (Mode B only)"
echo "  ${OUTPUT_DIR}/rna_counts.tsv             # wide RNA-seq counts (Mode B only)"
echo "  ${OUTPUT_DIR}/rpf_counts.tsv             # long-format Ribo-seq counts (Mode B only)"
echo "  ${OUTPUT_DIR}/condition_map.augmented.tsv  # rep1/rep2 names (when auto-split fired)"
fi
