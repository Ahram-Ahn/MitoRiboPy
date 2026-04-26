# MitoRiboPy CLI reference

Concise per-subcommand flag list. For every flag's full help text and current default, run:

```bash
mitoribopy --help
mitoribopy <subcommand> --help
mitoribopy all --show-stage-help <align|rpf|rnaseq>
```

The [`README.md`](../../README.md) at the repo root has the long-form description and per-flag explanation; this file is the at-a-glance grouped index.

---

## Top level

```
usage: mitoribopy [--version] <subcommand> [options]
```

| Subcommand | Purpose |
|---|---|
| `align` | FASTQ → BAM → BED6 + per-sample read counts |
| `rpf` | BED/BAM → offsets, translation profile, codon usage, coverage plots |
| `rnaseq` | DE table + rpf → TE / ΔTE with SHA256 reference-consistency gate |
| `all` | End-to-end orchestrator (align + rpf + optional rnaseq) |

Top-level flags: `--version, -V`, `--help, -h`.

### Shared options (every subcommand)

- `--config PATH` — JSON / YAML / TOML configuration file; CLI flags win on conflict.
- `--dry-run` — print planned actions and exit 0 without executing.
- `--threads N` — exports `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `MITORIBOPY_THREADS`.
- `--log-level {DEBUG,INFO,WARNING,ERROR}` (default `INFO`).

---

## `mitoribopy align`

### Inputs

- `--fastq-dir DIR` — directory of `*.fq(.gz)` / `*.fastq(.gz)`.
- `--fastq PATH` — individual FASTQ; repeatable.
- `--contam-index BT2_PREFIX` — contaminant (rRNA + tRNA + spikes) bowtie2 index.
- `--mt-index BT2_PREFIX` — mt-transcriptome bowtie2 index (one record per mt-mRNA).
- `--output DIR`.

### Library prep / kit

- `--kit-preset {auto, illumina_smallrna, illumina_truseq, illumina_truseq_umi, qiaseq_mirna, pretrimmed, custom}` (default `auto`). Vendor synonyms (`truseq_smallrna`, `nebnext_smallrna`, `nebnext_ultra_umi`, `truseq_stranded_total`, `smarter_pico_v3`, `sequoia_express`) are also accepted.
- `--adapter SEQ` — explicit 3' adapter (required for `--kit-preset custom`).
- `--umi-length N`, `--umi-position {5p,3p}` — kit-preset overrides.
- `--sample-overrides TSV` — per-sample overrides for `kit_preset` / `adapter` / `umi_length` / `umi_position` / `dedup_strategy`. Required header column `sample` plus at least one of the override columns. The orchestrator (`mitoribopy all`) materialises an `align.samples:` YAML block as this TSV automatically.

### Adapter detection

- `--adapter-detection {auto,off,strict}` (default `auto`).
- `--adapter-detect-reads N` (default 5000).
- `--adapter-detect-min-rate FRAC` (default 0.30).
- `--adapter-detect-min-len N` (default 12).
- `--adapter-detect-pretrimmed-threshold FRAC` (default 0.05).
- `--no-pretrimmed-inference` — disable the auto-fallback to `pretrimmed`.

### Trim / align

- `--library-strandedness {forward,reverse,unstranded}` (default `forward` → bowtie2 `--norc`).
- `--min-length 15`, `--max-length 45`, `--quality 20`.
- `--mapq 10`, `--seed 42`.

### Deduplication

- `--dedup-strategy {auto,umi-tools,skip}` (default `auto`, per sample). When the resolved strategy is `skip`, the orchestrator wires `aligned/<sample>.mapq.bam` straight into BED conversion — no duplicate `deduped/<sample>.dedup.bam` is written. The legacy `mark-duplicates` (picard) option was removed in v0.4.5 because coordinate-only dedup destroys codon-occupancy signal on mt-Ribo-seq libraries.
- `--umi-dedup-method {unique,percentile,cluster,adjacency,directional}` (default `unique`).

### Intermediate files + resume

- `--keep-intermediates` (default off) — keep `trimmed/<sample>.trimmed.fq.gz`, `contam_filtered/<sample>.nocontam.fq.gz`, and the pre-MAPQ `aligned/<sample>.bam`. Without this flag they are deleted as soon as the next step consumes them.
- `--resume` — skip samples whose `<output>/.sample_done/<sample>.json` marker is already present (auto-set by `mitoribopy all --resume` when the stage's `read_counts.tsv` is missing).

---

## `mitoribopy rpf`

### Core inputs

- `-f, --fasta REF_FASTA` (**required**).
- `-s, --strain {h.sapiens, s.cerevisiae, custom}` (default `h.sapiens`). Synonyms: `h`, `y` (printed deprecation notice on first use).
- `-d, --directory <BED/BAM dir>`, `--bam_mapq Q` (default 10 for BAM inputs).
- `-rpf MIN_LEN MAX_LEN` (default depends on strain × footprint_class).
- `--footprint_class {short, monosome, disome, custom}` (default `monosome`). `short` covers 16–24 nt RNase truncation products; `disome` covers collided ribosomes (h.sapiens 50–70, s.cerevisiae 60–90).
- `--annotation_file ANNOTATION.csv` — required for `--strain custom`. See the README's [Custom organisms](../../README.md#custom-organisms) section for the schema.
- `--codon_tables_file CODON_TABLES.json` — only when your organism's genetic code is not in the built-in NCBI list.
- `--codon_table_name TABLE_NAME` — pick a built-in NCBI Genetic Code (e.g. `vertebrate_mitochondrial`, `mold_mitochondrial`, `invertebrate_mitochondrial`). Run `mitoribopy rpf --help` for the full list.
- `--start_codons CODON [CODON ...]` — override the strain default.
- `--atp8_atp6_baseline {ATP6,ATP8}`, `--nd4l_nd4_baseline {ND4,ND4L}`.

### Offset enrichment + selection

- `-a, --align {start,stop}`, `-r, --range NT`.
- `--min_offset`, `--max_offset` — shared bounds (legacy fallback).
- `--min_5_offset`, `--max_5_offset`, `--min_3_offset`, `--max_3_offset` — end-specific bounds (preferred).
- `--offset_mask_nt NT` (default 5).
- `--offset_pick_reference {p_site, reported_site}` (default `p_site`). `p_site` picks in canonical P-site space then converts to the reported space; `reported_site` picks directly in the `--offset_site` space.
- `--offset_type {5,3}` (default `5`).
- `--offset_site {p,a}` (default `p`) — controls the SELECTED OFFSETS table coordinate space only.
- `--analysis_sites {p,a,both}` (default `both`) — controls which downstream sites are generated.
- `--codon_overlap_mode {full,any}` (default `full`).
- `-p, --psite_offset NT` — fixed offset for every read length (bypasses enrichment).
- `--offset_mode {per_sample,combined}` (default `per_sample`).

### Outputs + plotting

- `-o, --output DIR`, `--downstream_dir NAME`, `--plot_dir NAME`.
- `-fmt, --plot_format {png,pdf,svg}`.
- `--x_breaks`, `--line_plot_style {combined,separate}`, `--cap_percentile`.
- `-m, --codon_density_window` — sum coverage at the codon centre and ±1 nt before writing codon-density tables / plots (smooths short-window noise; does **not** collapse reading frames).
- `--order_samples NAME [NAME ...]`.

### Read-count normalization

- `--read_counts_file PATH`, `--read_counts_sample_col`, `--read_counts_reads_col`, `--read_counts_reference_col`.
- `--unfiltered_read_length_range MIN MAX` (default depends on `--footprint_class`: `short` → 10 30, `monosome` → 15 50, `disome` → 40 100).
- `--rpm_norm_mode {total,mt_mrna}`, `--mt_mrna_substring_patterns PATTERN [PATTERN ...]`.

### Optional modules

- `--structure_density`, `--structure_density_norm_perc`.
- `--cor_plot`, `--base_sample`, `--cor_mask_method {percentile,fixed,none}`, `--cor_mask_percentile`, `--cor_mask_threshold`.

---

## `mitoribopy rnaseq`

### DE table

- `--de-table PATH` (**required**).
- `--de-format {auto,deseq2,xtail,anota2seq,custom}` (default `auto`).
- `--de-gene-col`, `--de-log2fc-col`, `--de-padj-col`, `--de-basemean-col` (used with `--de-format custom`).

### Gene identifiers

- `--gene-id-convention {ensembl,refseq,hgnc,mt_prefixed,bare}` (**required, no default**).
- `--organism {h.sapiens, s.cerevisiae}` (default `h.sapiens`). Synonyms: `h`, `y`, `human`, `yeast`.

### Ribo-seq inputs

- `--ribo-dir DIR` — output of a prior `mitoribopy rpf` run.
- `--ribo-counts PATH` — defaults to `<ribo-dir>/rpf_counts.tsv`.

### Reference-consistency gate (exactly one)

- `--reference-gtf PATH`
- `--reference-checksum SHA256`

### Conditions (optional, required for replicate-based ΔTE)

- `--condition-map PATH` (TSV with `sample` and `condition` columns).
- `--condition-a NAME`, `--condition-b NAME`.

### Output

- `--output DIR`.

---

## `mitoribopy all`

- `--config PATH` (**required** unless using `--print-config-template` or `--show-stage-help`).
- `--output DIR` (**required** for non-dry-run).
- `--resume` — skip stages whose sentinel output already exists; also propagates `--resume` into the align stage so per-sample `.sample_done/<sample>.json` markers are honoured.
- `--skip-align`, `--skip-rpf`, `--skip-rnaseq`.
- `--manifest PATH` (default `run_manifest.json`).
- `--show-stage-help {align,rpf,rnaseq}` — print the full per-stage help and exit.
- `--print-config-template` — print a curated commented YAML template covering every stage and exit. The exhaustive copy-and-edit starter lives at the repo root as [`pipeline_config.example.yaml`](../../pipeline_config.example.yaml).

### Auto-wiring

When `align` and `rpf` both run:
- `rpf.directory` → `<run_root>/align/bed`
- `rpf.read_counts_file` → `<run_root>/align/read_counts.tsv`

When `rpf` and `rnaseq` both run:
- `rnaseq.ribo_dir` → `<run_root>/rpf`

---

## Synonyms

A few short / legacy spellings are accepted as synonyms. They emit a single `[mitoribopy] DEPRECATED:` line on first use and continue. The full table lives in the README's [Synonyms](../../README.md#synonyms) subsection.

| Synonym | Canonical |
|---|---|
| `-s h` / `-s y` | `-s h.sapiens` / `-s s.cerevisiae` |
| `--merge_density` | `--codon_density_window` |
| `--mrna_ref_patterns` | `--mt_mrna_substring_patterns` |
| `--offset_pick_reference selected_site` | `--offset_pick_reference reported_site` |

Plain `mitoribopy <flags>` (no subcommand) routes to `mitoribopy rpf` with a one-line notice; prefer the explicit subcommand form.
