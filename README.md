# MitoRiboPy

MitoRiboPy is a Python package for mitochondrial ribosome profiling (mt-Ribo-seq) analysis. Starting in v0.3.0 it spans the full pipeline from raw FASTQ through translation-efficiency integration with paired RNA-seq:

- `mitoribopy align` &mdash; FASTQ &rarr; BAM &rarr; BED6 + per-sample read counts (cutadapt + bowtie2 + umi_tools + pysam)
- `mitoribopy rpf` &mdash; BED/BAM &rarr; offsets, translation-profile, codon usage, coverage plots
- `mitoribopy rnaseq` &mdash; DE table (DESeq2 / Xtail / Anota2Seq) + rpf outputs &rarr; TE and &Delta;TE tables + plots, SHA256 reference-consistency gate
- `mitoribopy all` &mdash; end-to-end orchestrator with a shared config file and a composed `run_manifest.json`

## Highlights

- Subcommand CLI (`align` / `rpf` / `rnaseq` / `all`) with shared `--config`, `--dry-run`, `--threads`, `--log-level`
- Config files in JSON, YAML, or TOML (auto-detected by path suffix)
- Kit-aware FASTQ trimming: `truseq_smallrna`, `nebnext_smallrna`, `nebnext_ultra_umi`, `qiaseq_mirna`, or explicit `--adapter`
- Strand-aware mt-transcriptome alignment (`--library-strandedness forward` by default) so ND5 / ND6 antisense overlap is resolved by construction on Path A (transcriptome reference)
- Deduplication safe by default: `--dedup-strategy auto` picks UMI-aware when UMIs are present and skips otherwise; `mark-duplicates` is behind a long confirmation flag because coordinate-only dedup destroys codon-occupancy signal on low-complexity mt-Ribo-seq libraries
- BAM input to `rpf` via pysam (no samtools / bedtools PATH dependency)
- SHA256 reference-consistency gate on `rnaseq`: Ribo-seq and RNA-seq sides must be aligned to the identical transcript reference; mismatches are a hard fail
- Built-in human and yeast reference data (annotation CSVs + codon tables)
- End-specific 5'/3' offset selection, P-site vs A-site workflows, bicistronic ATP8/ATP6 and ND4L/ND4 handling
- Custom organism support via `--annotation_file`, `--codon_tables_file`, `--codon_table_name`, `--start_codons`
- Persistent per-run logging in `<output>/mitoribopy.log`
- Consistent terminal + file progress reporting for `align` and `rpf`
- Provenance: every stage writes a `run_settings.json`; `mitoribopy all` composes them into `run_manifest.json`

## Installation

From the repository root:

```bash
python -m pip install -e .
```

For development and tests:

```bash
python -m pip install -e ".[dev]"
```

Then confirm the CLI is available:

```bash
mitoribopy --help
```

If you prefer not to install the package yet:

```bash
PYTHONPATH=src python -m mitoribopy --help
```

## Quick Start

### `mitoribopy rpf` &mdash; BED/BAM through the analysis pipeline

```bash
mitoribopy rpf \
  -s h \
  -f <reference.fa> \
  --directory <ribo_bed_dir> \
  -rpf 29 34 \
  --output <results_dir>
```

Plain `mitoribopy <flags>` still works in v0.3.x but routes to `rpf` with a deprecation warning. Use the explicit subcommand form.

### `mitoribopy align` &mdash; FASTQ &rarr; BAM + BED

```bash
mitoribopy align \
  --kit-preset nebnext_smallrna \
  --library-strandedness forward \
  --fastq-dir <fastqs_dir> \
  --contam-index <bowtie2_rRNA_index_prefix> \
  --mt-index <bowtie2_mt_transcriptome_index_prefix> \
  --output <align_results_dir>
```

Use `--kit-preset custom --adapter <SEQ>` when your library isn't one of the built-in presets. External tools (`cutadapt`, `bowtie2`, `umi_tools`) must be on `$PATH`; see [docs/environment/environment.yml](docs/environment/environment.yml) for a ready-made bioconda env.

### `mitoribopy rnaseq` &mdash; DE table + rpf &rarr; TE / &Delta;TE

```bash
mitoribopy rnaseq \
  --de-table <deseq2_or_xtail_or_anota2seq_output.tsv> \
  --gene-id-convention hgnc \
  --ribo-dir <rpf_results_dir> \
  --reference-gtf <shared_reference.fa> \
  --condition-map <samples_to_conditions.tsv> \
  --condition-a control --condition-b knockdown \
  --output <rnaseq_results_dir>
```

`--gene-id-convention` is required (no default). The reference-consistency gate will hard-fail unless the hash of `--reference-gtf` matches the hash the prior `rpf` run recorded.

### `mitoribopy all` &mdash; end-to-end orchestrator

```bash
mitoribopy all --config pipeline_config.yaml --output <run_root>
```

Where `pipeline_config.yaml` has `align:`, `rpf:`, and optional `rnaseq:` sections; each section's keys correspond to the subcommand's CLI flag names. See [docs/tutorials/01_end_to_end_fastq.md](docs/tutorials/01_end_to_end_fastq.md) for a worked example.

Useful details for `mitoribopy all`:

- `mitoribopy all --help` shows only orchestrator-level flags. For full stage help, use:
  - `mitoribopy all --show-stage-help align`
  - `mitoribopy all --show-stage-help rpf`
  - `mitoribopy all --show-stage-help rnaseq`
- When `align` and `rpf` both run, `all` auto-wires:
  - `rpf.directory -> <run_root>/align/bed`
  - `rpf.read_counts_file -> <run_root>/align/read_counts.tsv`
- When `rpf` and `rnaseq` both run, `all` auto-wires:
  - `rnaseq.ribo_dir -> <run_root>/rpf`

### Logs and progress

- `mitoribopy align` writes `<output>/mitoribopy.log` and emits per-sample stage updates for trim, contaminant filtering, mt alignment, MAPQ filtering, deduplication, and BED export.
- `mitoribopy rpf` writes `<output>/mitoribopy.log` and emits numbered pipeline-step progress plus downstream plotting/profile progress.
- The same status lines are written to both the terminal and the log file.

## Built-In References

MitoRiboPy ships with packaged reference data for:

- Human mitochondrial translation using the `vertebrate_mitochondrial` codon table
- Yeast mitochondrial translation using the `yeast_mitochondrial` codon table

Built-in annotation tables are stored as CSV and built-in codon tables are stored as JSON under [src/mitoribopy/data](src/mitoribopy/data).

For bicistronic transcript regions:

- Titles stay consistent as `ATP8/ATP6` and `ND4L/ND4`
- The default sequence baselines are `ATP6` and `ND4`
- You can switch them with `--atp8_atp6_baseline ATP8|ATP6` and `--nd4l_nd4_baseline ND4L|ND4`

Legacy FASTA/BED identifiers such as `ATP86` and `ND4L4` are still recognized through built-in aliases.

## Custom Organisms

Custom organisms are supported through:

- `--annotation_file`
- `--codon_tables_file`
- `--codon_table_name`
- `--start_codons`

For `--strain custom`, provide an explicit RPF range as well:

```bash
mitoribopy \
  -s custom \
  -f <reference.fa> \
  --directory <ribo_bed_dir> \
  -rpf 28 34 \
  --annotation_file examples/custom_reference/annotation_template.csv \
  --codon_tables_file examples/custom_reference/codon_tables_template.json \
  --codon_table_name custom_example \
  --start_codons ATG GTG \
  --output <results_dir>
```

Example templates are included here:

- [examples/custom_reference/annotation_template.csv](examples/custom_reference/annotation_template.csv)
- [examples/custom_reference/codon_tables_template.json](examples/custom_reference/codon_tables_template.json)
- [examples/custom_reference/README.md](examples/custom_reference/README.md)

## CLI Parameters

### Required parameters

- `-f, --fasta`: reference FASTA

### Usually required for a normal run

These are not all technically mandatory in the parser, but they are the recommended minimum for a reproducible run:

- `-s, --strain`
- `--directory`
- `-rpf <min> <max>`
- `--output`

### Additional required parameters for `--strain custom`

- `--annotation_file`
- `--codon_tables_file` or `--codon_table_name`
- `-rpf <min> <max>`

### Common optional parameters

- `--align start|stop`
- `--offset_type 5|3`
- `--offset_site p|a`
- `--offset_pick_reference p_site|selected_site`
- `--min_5_offset`, `--max_5_offset`
- `--min_3_offset`, `--max_3_offset`
- `--offset_mask_nt`
- `--read_counts_file`
- `--read_counts_sample_col`
- `--read_counts_reference_col`
- `--read_counts_reads_col`
- `--unfiltered_read_length_range <min> <max>`
- `--rpm_norm_mode total|mt_mrna`
- `--plot_format png|pdf|svg`
- `-m, --merge_density`
- `--structure_density`
- `--cor_plot`
- `--use_rna_seq`

## Example Usage

### Human or yeast with default-style analysis

```bash
mitoribopy rpf \
  -s h \
  -f <reference.fa> \
  --directory <ribo_bed_dir> \
  -rpf 29 34 \
  --align stop \
  --offset_type 5 \
  --offset_site p \
  --offset_pick_reference p_site \
  --offset_mask_nt 5 \
  --min_5_offset 10 \
  --max_5_offset 22 \
  --min_3_offset 10 \
  --max_3_offset 22 \
  --plot_format svg \
  --output <results_dir> \
  -m
```

### Run with read-count normalization

```bash
mitoribopy \
  -s h \
  -f <reference.fa> \
  --directory <ribo_bed_dir> \
  -rpf 29 34 \
  --read_counts_file <read_counts.csv> \
  --read_counts_sample_col sample \
  --read_counts_reads_col reads \
  --read_counts_reference_col reference \
  --rpm_norm_mode mt_mrna \
  --mrna_ref_patterns mt_genome \
  --output <results_dir>
```

### Inspect broader read-length QC ranges

```bash
mitoribopy rpf \
  -s h \
  -f <reference.fa> \
  --directory <ribo_bed_dir> \
  -rpf 29 34 \
  --unfiltered_read_length_range 15 60 \
  --output <results_dir>
```

This keeps the filtered analysis range at `29-34 nt` while broadening the unfiltered QC tables and heatmaps so longer footprints remain visible.

### Run optional downstream modules

```bash
mitoribopy \
  -s h \
  -f <reference.fa> \
  --directory <ribo_bed_dir> \
  -rpf 29 34 \
  --structure_density \
  --cor_plot \
  --base_sample <sample_name> \
  --output <results_dir>
```

### Run a custom organism

```bash
mitoribopy \
  -s custom \
  -f <reference.fa> \
  --directory <ribo_bed_dir> \
  -rpf 28 34 \
  --annotation_file <annotation.csv> \
  --codon_tables_file <codon_tables.json> \
  --codon_table_name <table_name> \
  --start_codons ATG GTG \
  --output <results_dir>
```

## Input Files

### BED

Expected columns:

1. `chrom`
2. `start`
3. `end`

Additional BED columns are tolerated. Coordinates are treated as standard 0-based, end-exclusive intervals.

### FASTA

FASTA headers should match the annotation `sequence_name` or one of its `sequence_aliases`.

### Annotation CSV

Required columns:

- `transcript`
- `l_tr`
- `l_utr5`
- `l_utr3`

Optional columns:

- `l_cds`
- `sequence_name`
- `sequence_aliases`
- `display_name`

Meaning:

- `transcript` is the logical CDS name used in frame and codon outputs
- `sequence_name` is the FASTA/BED sequence ID that the row maps onto
- `sequence_aliases` contains alternate FASTA/BED names separated by semicolons
- `display_name` controls plot titles and grouped transcript labels

If `l_cds` is omitted, it is computed as `l_tr - l_utr5 - l_utr3`.

### Codon-Table JSON

Two formats are supported:

- One flat 64-codon mapping
- A dictionary of named 64-codon mappings

When multiple named tables are present, choose one with `--codon_table_name`.

### Read-Count Table

`.csv`, `.tsv`, and `.txt` are supported. Column matching is flexible and case-insensitive, with fallback to positional matching:

- first column: sample
- second column: reference
- third column: read count

## Output Overview

Typical output structure:

```text
<output>/
  mitoribopy.log
  plots_and_csv/
  <sample>/
    footprint_density/
    translating_frame/
    codon_usage/
    debug_csv/
  coverage_profile_plots/
  structure_density/      # if --structure_density
  codon_correlation/      # if --cor_plot
  rna_seq_results/        # if --use_rna_seq
```

Key outputs include:

- offset enrichment CSVs and plots
- selected offset tables by read length
- footprint-density CSVs for P-site, A-site, and E-site
- frame-usage summaries
- transcript-level and total codon-usage summaries
- RPM and raw coverage-profile plots
- CDS-aware codon-binned coverage plots (`*_codon/`, 3 nt combined per codon)
- optional structure-density exports from footprint-density tables

## Important Runtime Notes

- `--offset_type 5|3`: downstream site placement from the read 5' or 3' end
- `--offset_site p|a`: whether reported offsets represent P-site or A-site positions
- `--offset_pick_reference p_site|selected_site`: how the best offset is chosen
- `--min_5_offset`, `--max_5_offset`, `--min_3_offset`, `--max_3_offset`: recommended end-specific selection bounds
- `--offset_mask_nt`: mask near-anchor bins from enrichment summaries and plots
- `--rpm_norm_mode total|mt_mrna`: read-count normalization mode
- `--structure_density`: export log2 and scaled density values from footprint-density tables

For the full interface, run:

```bash
mitoribopy --help
```

## Development

Run the test suite with:

```bash
PYTHONPATH=src pytest
```

This repository also includes package migration notes and release materials under [docs/README.md](docs/README.md).

## License

MIT. See [LICENSE](LICENSE).
