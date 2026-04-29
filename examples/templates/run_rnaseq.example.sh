#!/usr/bin/env bash
# MitoRiboPy rnaseq -- exhaustive shell-script template (rnaseq stage only).
# Compatible with: MitoRiboPy 0.5.1+
#
# `mitoribopy rnaseq` has TWO mutually-exclusive flows:
#
#   * DEFAULT (from-FASTQ) -- pass raw RNA-seq + Ribo-seq FASTQs and a
#     transcriptome FASTA. The subcommand auto-detects SE vs PE from
#     filename mate tokens, runs cutadapt + bowtie2 per sample, counts
#     reads per transcript, runs pyDESeq2, then emits TE / ΔTE / plots.
#     This is the default and the recommended starting point. Requires
#     the `[fastq]` extra: pip install 'mitoribopy[fastq]'.
#
#   * ALTERNATIVE (--de-table) -- you already ran DESeq2 / Xtail /
#     Anota2Seq externally on the full transcriptome and have the
#     results table. Pass `--de-table`. SHA256 reference-consistency
#     gate is enforced against a prior `mitoribopy rpf` run. Use this
#     for publication-grade DE statistics (the default flow runs DE on
#     just the 13 mt-mRNAs, which is fine for exploration but noisy
#     for inference).
#
# Quick start:
#   1. Set FLOW below ("default" or "de-table") and edit the
#      corresponding variables block.
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

# Pick the flow: "default" = from-FASTQ (recommended), "de-table" = bring your own DE.
FLOW="default"

# Common knobs.
THREADS=8                 # global thread budget (BLAS / pyDESeq2 caps)
LOG_LEVEL=INFO            # DEBUG | INFO | WARNING | ERROR
OUTPUT_DIR="${PROJECT_ROOT}/results/rnaseq"

# Reference (required by both modes).
MT_FASTA="${PROJECT_ROOT}/input_data/human-mt-mRNA.fasta"

# Gene-ID convention used in the DE table / FASTA records (REQUIRED).
GENE_ID_CONVENTION="hgnc"      # ensembl | refseq | hgnc | mt_prefixed | bare
ORGANISM="h.sapiens"           # h.sapiens | s.cerevisiae

# Conditions (REQUIRED in default flow; optional / replicate-based ΔTE in --de-table flow).
# BASE_SAMPLE = reference condition (denominator of the WT-vs-X contrast); seeds the
# baseline label on every comparison plot. COMPARE_SAMPLE = comparison condition.
# These map to the canonical --base-sample / --compare-sample CLI flags.
# `--condition-a` / `--condition-b` are still accepted as legacy aliases.
CONDITION_MAP="${PROJECT_ROOT}/samples.tsv"
BASE_SAMPLE="control"
COMPARE_SAMPLE="knockdown"


# ============================================================================
# Default flow: from-FASTQ
# ============================================================================
# Inputs (REQUIRED for the default flow):
RNA_FASTQ_DIR="${PROJECT_ROOT}/input_data/rna_seq"
RIBO_FASTQ_DIR="${PROJECT_ROOT}/input_data/ribo_seq"

# Per-stage compute knobs (default flow only).
ALIGN_THREADS=4                 # threads passed to cutadapt + bowtie2
WORKDIR="${OUTPUT_DIR}/work"    # scratch (defaults to <output>/work if unset)


# ============================================================================
# Alternative flow: pre-computed DE table
# ============================================================================
# Inputs (REQUIRED for the --de-table flow):
DE_TABLE="${PROJECT_ROOT}/de.tsv"
RIBO_DIR="${PROJECT_ROOT}/results/rpf"


# ============================================================================
# Build per-flow argv
# ============================================================================
COMMON_OPTS=(
  # --- Gene identifier convention (REQUIRED, no default) -----------------
  --gene-id-convention "${GENE_ID_CONVENTION}"
  --organism "${ORGANISM}"

  # --- Output ------------------------------------------------------------
  --output "${OUTPUT_DIR}"

  # --- Conditions --------------------------------------------------------
  --condition-map "${CONDITION_MAP}"
  --base-sample "${BASE_SAMPLE}"
  --compare-sample "${COMPARE_SAMPLE}"

  # --- Shared ------------------------------------------------------------
  --threads "${THREADS}"
  --log-level "${LOG_LEVEL}"
)

if [[ "${FLOW}" == "default" ]]; then
  RNASEQ_OPTS=(
    "${COMMON_OPTS[@]}"

    # --- From-FASTQ inputs (default; mutually exclusive with --de-table)
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
elif [[ "${FLOW}" == "de-table" ]]; then
  RNASEQ_OPTS=(
    "${COMMON_OPTS[@]}"

    # --- DE table (alternative flow) -----------------------------------
    --de-table "${DE_TABLE}"
    --de-format auto              # auto | deseq2 | xtail | anota2seq | custom
    # --de-gene-col gene_id       # only when --de-format custom
    # --de-log2fc-col log2FoldChange
    # --de-padj-col padj
    # --de-basemean-col baseMean

    # --- Ribo-seq inputs (alternative flow) ----------------------------
    --ribo-dir "${RIBO_DIR}"
    # --ribo-counts "${RIBO_DIR}/rpf_counts.tsv"   # explicit override

    # --- Reference-consistency gate (alternative flow; exactly one) ----
    --reference-gtf "${MT_FASTA}"
    # --reference-checksum 5ca397907373...        # OR pre-computed SHA256
  )
else
  echo "FLOW must be 'default' or 'de-table' (got: ${FLOW})" >&2
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
echo "MitoRiboPy rnaseq finished (FLOW=${FLOW})."
echo "Inspect outputs under: ${OUTPUT_DIR}/"
echo
echo "Key files to look at first:"
echo "  ${OUTPUT_DIR}/te.tsv                       # per-(sample,gene) TE"
echo "  ${OUTPUT_DIR}/delta_te.tsv                 # per-gene ΔTE log2 + padj + note"
echo "  ${OUTPUT_DIR}/run_settings.json            # full provenance + per-sample stats"
echo
echo "Plots — every PNG ships with a sibling .svg (300 dpi PNG + editable-text SVG):"
echo "  ${OUTPUT_DIR}/plots/mrna_vs_rpf.png        # 4-quadrant log2FC scatter (mRNA vs RPF)"
echo "  ${OUTPUT_DIR}/plots/delta_te_volcano.png   # ΔTE volcano with stat box"
echo "  ${OUTPUT_DIR}/plots/ma.png                 # DESeq2 MA plot"
echo "  ${OUTPUT_DIR}/plots/de_volcano_mrna.png    # mRNA DE volcano (WT-vs-X)"
echo "  ${OUTPUT_DIR}/plots/te_bar_by_condition.png  # primary readout — bars + replicate dots"
echo "  ${OUTPUT_DIR}/plots/te_heatmap.png         # gene × sample log2(TE) with condition strip"
echo "  ${OUTPUT_DIR}/plots/te_compare_scatter.png # base vs compare TE scatter (when --condition-map +"
echo "                                             # --base-sample + --compare-sample are all set)"
echo "  ${OUTPUT_DIR}/plots/te_log2fc_bar.png      # sorted log2(TE_compare/TE_base) per gene"
if [[ "${FLOW}" == "default" ]]; then
echo "  ${OUTPUT_DIR}/plots/sample_pca.png         # sample PCA (default flow only)"
echo "  ${OUTPUT_DIR}/plots/de_volcano_rpf.png     # Ribo-seq DE volcano (default flow only;"
echo "                                             # skipped when Ribo subset has < 2 condition levels)"
echo "  ${OUTPUT_DIR}/de_table.tsv                 # mRNA pyDESeq2 result"
echo "  ${OUTPUT_DIR}/rpf_de_table.tsv             # Ribo-seq pyDESeq2 result (default flow only)"
echo "  ${OUTPUT_DIR}/rna_counts.tsv               # wide RNA-seq counts"
echo "  ${OUTPUT_DIR}/rpf_counts.tsv               # long-format Ribo-seq counts"
echo "  ${OUTPUT_DIR}/rpf_counts_matrix.tsv        # wide Ribo-seq counts"
echo "  ${OUTPUT_DIR}/condition_map.augmented.tsv  # rep1/rep2 names (when auto-split fired)"
fi
