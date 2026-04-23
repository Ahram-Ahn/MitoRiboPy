# MitoRiboPy CLI reference

Auto-generated summary of every subcommand's flag surface at v0.3.0.
For the canonical, always-current help text, run:

```bash
mitoribopy --help
mitoribopy <subcommand> --help
```

## Top level

```
usage: mitoribopy [--version] <subcommand> [options]

Subcommands:
  align    Preprocess FASTQ inputs (cutadapt + bowtie2 + dedup + BAM->BED6)
  rpf      Run the Ribo-seq analysis pipeline from BED/BAM inputs
  rnaseq   DE + rpf -> TE / delta-TE with SHA256 reference-consistency gate
  all      End-to-end orchestrator (align + rpf + optional rnaseq)

Top-level options:
  --version, -V
  --help,    -h
```

### Shared options (every subcommand)

- `--config PATH` &mdash; JSON / YAML / TOML configuration file; CLI flags win on conflict.
- `--dry-run` &mdash; Print planned actions and exit 0 without executing.
- `--threads N` &mdash; Exports `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `MITORIBOPY_THREADS`.
- `--log-level {DEBUG,INFO,WARNING,ERROR}`

## `mitoribopy align`

### Inputs
- `--fastq-dir DIR` &mdash; directory of `*.fq(.gz)` / `*.fastq(.gz)`
- `--fastq PATH` (repeatable)
- `--contam-index BT2_PREFIX` &mdash; rRNA / tRNA / spike bowtie2 index
- `--mt-index BT2_PREFIX` &mdash; mt-transcriptome bowtie2 index
- `--output DIR`

### Library prep
- `--kit-preset {truseq_smallrna, nebnext_smallrna, nebnext_ultra_umi, qiaseq_mirna, custom}` (default `custom`)
- `--adapter SEQ` (required when `--kit-preset custom`)
- `--umi-length N`, `--umi-position {5p,3p}`
- `--library-strandedness {forward,reverse,unstranded}` (default `forward` &rarr; bowtie2 `--norc`)
- `--min-length 15`, `--max-length 45`, `--quality 20`

### Alignment
- `--mapq 10`, `--seed 42`

### Deduplication
- `--dedup-strategy {auto,umi-tools,skip,mark-duplicates}` (default `auto`)
- `--umi-dedup-method {unique,percentile,cluster,adjacency,directional}` (default `unique`)
- `--i-understand-mark-duplicates-destroys-mt-ribo-seq-signal` (required for `mark-duplicates`)

## `mitoribopy rpf`

### Core inputs
- `-f / --fasta REF_FASTA` (required), `-s / --strain {y,h,custom}`, `-d / --directory <BED/BAM dir>`
- `-rpf MIN_LEN MAX_LEN`, `--bam_mapq Q` (default 10 for BAM inputs)
- `--annotation_file`, `--codon_tables_file`, `--codon_table_name`, `--start_codons`
- `--atp8_atp6_baseline {ATP6,ATP8}`, `--nd4l_nd4_baseline {ND4,ND4L}`

### Offsets
- `-a / --align {start,stop}`, `-r / --range NT`
- `--min_offset`, `--max_offset`, `--min_5_offset`, `--max_5_offset`, `--min_3_offset`, `--max_3_offset`
- `--offset_mask_nt NT`
- `--offset_pick_reference {p_site, selected_site}`
- `--offset_type {5,3}`, `--offset_site {p,a}`
- `--codon_overlap_mode {full,any}`
- `-p / --psite_offset NT`

### Outputs + normalization
- `-o / --output DIR`, `--downstream_dir NAME`, `--plot_dir NAME`, `-fmt / --plot_format {png,pdf,svg}`
- `--x_breaks`, `--line_plot_style {combined,separate}`, `--cap_percentile`
- `-m / --merge_density`, `--order_samples`
- `--read_counts_file`, `--read_counts_sample_col`, `--read_counts_reads_col`, `--read_counts_reference_col`
- `--rpm_norm_mode {total,mt_mrna}`, `--mrna_ref_patterns`

### Optional modules
- `--structure_density`, `--structure_density_norm_perc`
- `--cor_plot`, `--base_sample`, `--cor_mask_method`, `--cor_mask_percentile`, `--cor_mask_threshold`
- `--use_rna_seq` (**DEPRECATED** in v0.3.0; removed in v0.4.0; use `mitoribopy rnaseq`)

## `mitoribopy rnaseq`

### DE table
- `--de-table PATH` (required)
- `--de-format {auto,deseq2,xtail,anota2seq,custom}` (default `auto`)
- `--de-gene-col`, `--de-log2fc-col`, `--de-padj-col`, `--de-basemean-col` (for `--de-format custom`)

### Gene identifiers
- `--gene-id-convention {ensembl,refseq,hgnc,mt_prefixed,bare}` (**required**; no default)
- `--organism {h,y,human,yeast}` (default `h`)

### Ribo-seq inputs
- `--ribo-dir DIR` (required) &mdash; output of a prior `mitoribopy rpf` run
- `--ribo-counts PATH` (defaults to `<ribo-dir>/rpf_counts.tsv`)

### Reference-consistency gate (exactly one)
- `--reference-gtf PATH` or `--reference-checksum SHA256`

### Conditions (optional, required for replicate-based &Delta;TE)
- `--condition-map PATH`, `--condition-a NAME`, `--condition-b NAME`

### Output
- `--output DIR`

## `mitoribopy all`

- `--config PATH` (required; YAML / JSON / TOML with `align:` / `rpf:` / `rnaseq:` sections)
- `--output DIR` (required; run root)
- `--resume` &mdash; skip stages whose sentinel output already exists
- `--skip-align / --skip-rpf / --skip-rnaseq`
- `--manifest PATH` (default `run_manifest.json`, relative to `--output`)

## Backward compatibility

`mitoribopy <flags>` without a subcommand (the v0.2.x invocation
style) is still accepted in v0.3.x and routed to `mitoribopy rpf
<flags>` with a stderr `DEPRECATION` warning. The fallback is
scheduled for removal in v0.4.0.
