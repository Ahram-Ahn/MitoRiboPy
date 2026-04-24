# MitoRiboPy CLI reference (v0.4.1)

Concise per-subcommand flag list. For every flag's full help text and current default, run:

```bash
mitoribopy --help
mitoribopy <subcommand> --help
mitoribopy all --show-stage-help <align|rpf|rnaseq>
```

The `README.md` at the repo root has the long-form description and per-flag explanation; this file is the at-a-glance grouped index.

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

- `--kit-preset {auto, illumina_smallrna, illumina_truseq, illumina_truseq_umi, qiaseq_mirna, pretrimmed, custom}` (default `auto`). Legacy aliases (`truseq_smallrna`, `nebnext_smallrna`, `nebnext_ultra_umi`, `truseq_stranded_total`, `smarter_pico_v3`, `sequoia_express`) are still accepted.
- `--adapter SEQ` — explicit 3' adapter (required for `--kit-preset custom`).
- `--umi-length N`, `--umi-position {5p,3p}` — kit-preset overrides.

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

- `--dedup-strategy {auto,umi-tools,skip,mark-duplicates}` (default `auto`, per sample).
- `--umi-dedup-method {unique,percentile,cluster,adjacency,directional}` (default `unique`).
- `--i-understand-mark-duplicates-destroys-mt-ribo-seq-signal` — required confirmation for `--dedup-strategy mark-duplicates`.

---

## `mitoribopy rpf`

### Core inputs

- `-f, --fasta REF_FASTA` (**required**), `-s, --strain {h,y,vm,ym,custom}` (default `y`).
- `-d, --directory <BED/BAM dir>`, `--bam_mapq Q` (default 10 for BAM inputs).
- `-rpf MIN_LEN MAX_LEN` (default depends on strain × footprint_class).
- `--footprint_class {monosome,disome,custom}` (default `monosome`).
- `--annotation_file`, `--codon_tables_file`, `--codon_table_name`, `--start_codons`.
- `--atp8_atp6_baseline {ATP6,ATP8}`, `--nd4l_nd4_baseline {ND4,ND4L}`.

### Offset enrichment + selection

- `-a, --align {start,stop}`, `-r, --range NT`.
- `--min_offset`, `--max_offset` — shared bounds (legacy fallback).
- `--min_5_offset`, `--max_5_offset`, `--min_3_offset`, `--max_3_offset` — end-specific bounds (preferred).
- `--offset_mask_nt NT` (default 5).
- `--offset_pick_reference {p_site, selected_site}` (default `p_site`).
- `--offset_type {5,3}` (default `5`).
- `--offset_site {p,a}` (default `p`) — controls SELECTED OFFSETS table coordinate space only.
- `--analysis_sites {p,a,both}` (default `both`) — controls which downstream sites are generated.
- `--codon_overlap_mode {full,any}` (default `full`).
- `-p, --psite_offset NT` — fixed offset for every read length (bypasses enrichment).
- `--offset_mode {per_sample,combined}` (default `per_sample`).

### Outputs + plotting

- `-o, --output DIR`, `--downstream_dir NAME`, `--plot_dir NAME`.
- `-fmt, --plot_format {png,pdf,svg}`.
- `--x_breaks`, `--line_plot_style {combined,separate}`, `--cap_percentile`.
- `-m, --merge_density`, `--order_samples`.

### Read-count normalization

- `--read_counts_file PATH`, `--read_counts_sample_col`, `--read_counts_reads_col`, `--read_counts_reference_col`.
- `--unfiltered_read_length_range MIN MAX` (default 15 50).
- `--rpm_norm_mode {total,mt_mrna}`, `--mrna_ref_patterns`.

### Optional modules

- `--structure_density`, `--structure_density_norm_perc`.
- `--cor_plot`, `--base_sample`, `--cor_mask_method {percentile,fixed,none}`, `--cor_mask_percentile`, `--cor_mask_threshold`.
- `--use_rna_seq` — **DEPRECATED** in v0.3.0; **removed** in v0.4.0. Use `mitoribopy rnaseq`.

---

## `mitoribopy rnaseq`

### DE table

- `--de-table PATH` (**required**).
- `--de-format {auto,deseq2,xtail,anota2seq,custom}` (default `auto`).
- `--de-gene-col`, `--de-log2fc-col`, `--de-padj-col`, `--de-basemean-col` (used with `--de-format custom`).

### Gene identifiers

- `--gene-id-convention {ensembl,refseq,hgnc,mt_prefixed,bare}` (**required, no default**).
- `--organism {h,y}` (default `h`).

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
- `--resume` — skip stages whose sentinel output already exists.
- `--skip-align`, `--skip-rpf`, `--skip-rnaseq`.
- `--manifest PATH` (default `run_manifest.json`).
- `--show-stage-help {align,rpf,rnaseq}` — print the full per-stage help and exit.
- `--print-config-template` — print a commented YAML template covering every stage and exit.

### Auto-wiring

When `align` and `rpf` both run:
- `rpf.directory` → `<run_root>/align/bed`
- `rpf.read_counts_file` → `<run_root>/align/read_counts.tsv`

When `rpf` and `rnaseq` both run:
- `rnaseq.ribo_dir` → `<run_root>/rpf`

---

## Backward compatibility

`mitoribopy <flags>` without a subcommand (the v0.2.x invocation style) is still accepted in v0.3+ and routed to `mitoribopy rpf <flags>` with a stderr `DEPRECATION` warning. The fallback is scheduled for removal in a future major version; new code should use the explicit subcommand form.
