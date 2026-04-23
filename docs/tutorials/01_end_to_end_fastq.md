# Tutorial 01 &mdash; End-to-end mt-Ribo-seq from FASTQ to codon usage

This tutorial walks through a full mt-Ribo-seq analysis with MitoRiboPy
v0.3.0, starting from raw FASTQ and ending with per-sample
translation-profile outputs ready for codon-occupancy interpretation.

Every command below is a _shell_ command (lines starting with `$`). The
outputs are sketched so you know what to expect at each stage. The
tutorial uses a human NEBNext Small RNA library as the running example;
adjust `--kit-preset`, `--library-strandedness`, and `-s` / `-rpf` to
match your own library chemistry and organism.

> **Prerequisites.** Install MitoRiboPy (see the repo README) and the
> external tools it shells out to. The quickest way is the bioconda
> environment shipped under `docs/environment/environment.yml`:
>
> ```bash
> $ conda env create -f docs/environment/environment.yml
> $ conda activate mitoribopy
> $ mitoribopy --version   # expect "MitoRiboPy 0.3.0"
> ```
>
> `cutadapt`, `bowtie2`, `bowtie2-build`, `samtools`, and `umi_tools`
> must be on `$PATH`. `fastqc` and `picard` are optional.

## Step 0 &mdash; Data layout

Create a workspace with the following shape:

```text
<project_root>/
  fastqs/
    ctrl_1.fq.gz
    ctrl_2.fq.gz
    kd_1.fq.gz
    kd_2.fq.gz
  references/
    human_mt_transcriptome.fa      # one FASTA record per mt-mRNA
    human_rrna_contam.fa            # 12S, 16S, cytoplasmic rRNA spike
  pipeline_config.yaml
```

Your FASTA headers in `human_mt_transcriptome.fa` must match the
`sequence_name` column of the annotation CSV. MitoRiboPy ships
annotation CSVs for human and yeast under
`src/mitoribopy/data/`. For an in-house transcriptome reference
inspect those files and mirror the naming, or build your own
annotation CSV and pass it with `--annotation_file`.

## Step 1 &mdash; Build the bowtie2 indexes

```bash
$ mkdir -p indexes
$ bowtie2-build --quiet references/human_rrna_contam.fa       indexes/rrna_contam
$ bowtie2-build --quiet references/human_mt_transcriptome.fa  indexes/mt_tx
```

Each command produces six `*.bt2` files sharing a common prefix. The
prefix (e.g. `indexes/mt_tx`) is what you pass to MitoRiboPy, not an
individual `.bt2` file.

## Step 2 &mdash; Write a shared pipeline config

Config files can be YAML, JSON, or TOML. The following YAML is the
most common shape for `mitoribopy all`:

```yaml
align:
  kit_preset: nebnext_smallrna
  library_strandedness: forward
  fastq_dir: fastqs/
  contam_index: indexes/rrna_contam
  mt_index: indexes/mt_tx
  mapq: 10
  min_length: 15
  max_length: 45
  dedup_strategy: auto   # -> skip because nebnext_smallrna has umi_length == 0

rpf:
  strain: h
  fasta: references/human_mt_transcriptome.fa
  rpf: [29, 34]
  align: stop
  offset_type: 5
  offset_site: p
  offset_pick_reference: p_site
  min_5_offset: 10
  max_5_offset: 22
  min_3_offset: 10
  max_3_offset: 22
  offset_mask_nt: 5
  plot_format: svg
  merge_density: true

# Optional rnaseq section - see Tutorial 02.
# rnaseq:
#   de_table: de.tsv
#   gene_id_convention: hgnc
#   reference_gtf: references/human_mt_transcriptome.fa
```

> Every key under a section maps to the corresponding subcommand's CLI
> flag with dashes turned to underscores (so `kit_preset` -> `--kit-preset`,
> `library_strandedness` -> `--library-strandedness`). Booleans emit the
> bare flag (`true`) or are omitted entirely (`false`).

## Step 3 &mdash; Dry-run the full pipeline

```bash
$ mitoribopy all --config pipeline_config.yaml --output results/ --dry-run
[all] dry-run: planned actions
  1. align: --kit-preset nebnext_smallrna --library-strandedness forward --fastq-dir fastqs/ ...
  2. rpf: --strain h --fasta references/human_mt_transcriptome.fa --rpf 29 34 ...
  3. write manifest to results/run_manifest.json
```

The dry-run prints the exact argv each subcommand will receive, and
auto-wires `rpf`'s `--directory` to `results/align/bed/` so `rpf`
ingests the BED6 produced by `align`.

## Step 4 &mdash; Run end-to-end

```bash
$ mitoribopy all --config pipeline_config.yaml --output results/
```

On completion, inspect:

```text
results/
  align/
    trimmed/                *.trimmed.fq.gz, *.cutadapt.json
    contam_filtered/        *.nocontam.fq.gz
    aligned/                *.bam, *.mapq.bam
    deduped/                *.dedup.bam (or hardlink of *.mapq.bam when dedup=skip)
    bed/                    *.bed              <-- strand-aware BED6
    read_counts.tsv
    run_settings.json
  rpf/
    mitoribopy.log
    plots_and_csv/          offset_*.csv, offset_*.svg, p_site_offsets_*.csv
    <sample>/
      footprint_density/    codon-level P/A/E-site tables
      translating_frame/    frame usage
      codon_usage/          per-codon tables
      debug_csv/
    coverage_profile_plots/
    rpf_counts.tsv          per-sample per-gene RPF counts   <-- feeds rnaseq
    run_settings.json       includes reference_checksum
  run_manifest.json          composed provenance
```

## Step 5 &mdash; Interpret the outputs

- `align/read_counts.tsv` tells you how many reads survived each stage.
  The invariants `rrna_aligned + post_rrna_filter == post_trim` and
  `mt_aligned + unaligned_to_mt == post_rrna_filter` must hold; if they
  don't, something is wrong with the alignment step.
- `rpf/plots_and_csv/offset_stop.csv` (and the SVG next to it) is the
  offset enrichment around the stop codon. A sharp peak at the
  expected P-site offset (~12-15 nt from the 5' end for most mt-Ribo-seq
  libraries) confirms the library quality.
- `rpf/<sample>/footprint_density/` holds P/A/E-site codon-resolution
  tables; these are the primary input to pause-site analysis.

## Step 6 &mdash; Resume a failed run

If `align` finished but `rpf` bailed on a config typo, fix the config
and re-run with `--resume`:

```bash
$ mitoribopy all --config pipeline_config.yaml --output results/ --resume
```

`--resume` skips a stage when its sentinel output file already exists
(`align/read_counts.tsv`, `rpf/rpf_counts.tsv`, `rnaseq/delta_te.tsv`).
Individual stages can also be force-skipped via `--skip-align`,
`--skip-rpf`, `--skip-rnaseq`.

## Strandedness + ND5 / ND6 note

Human mt-ND5 and mt-ND6 are transcribed from opposite strands and their
3' ends overlap at the genome level. On Path A (transcriptome
reference, the MitoRiboPy default) every mt-mRNA is its own FASTA
record, so the overlap is not a concern at alignment time.
`--library-strandedness forward` (default) additionally passes
`--norc` to bowtie2 so reverse-complement hits are rejected. The
resulting `bed/<sample>.bed` should show only `+` strand rows; any
`-` strand rows are worth investigating because they suggest either a
library-prep mismatch or a reverse-stranded kit that should have been
run with `--library-strandedness reverse`.
