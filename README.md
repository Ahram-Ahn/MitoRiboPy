# MitoRiboPy

MitoRiboPy is a package for mitochondrial ribosome profiling analysis. It runs a package-native pipeline from BED inputs through offset selection, translation-profile analysis, codon usage, coverage-profile plotting, and optional downstream modules such as structure-density export, codon correlation, and RNA-seq integration.

## Highlights

- Standalone package CLI with no runtime dependency on legacy pipeline scripts
- Built-in human and yeast reference data loaded from packaged CSV and JSON files
- End-specific offset selection with separate 5' and 3' bounds
- P-site and A-site workflows with explicit offset-picking behavior
- In-memory BED filtering with no duplicated filtered BED output files
- Persistent per-run logging in `<output>/mitoribopy.log`
- Custom organism support through user-supplied annotation CSV and codon-table JSON files
- Bicistronic transcript handling for ATP8/ATP6 and ND4L/ND4 with configurable baseline sequence IDs

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

Minimal example:

```bash
mitoribopy \
  -s h \
  -f <reference.fa> \
  --directory <ribo_bed_dir> \
  -rpf 29 34 \
  --output <results_dir>
```

Compatibility wrapper:

```bash
python main.py --help
```

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
- `--rpm_norm_mode total|mt_mrna`
- `--plot_format png|pdf|svg`
- `-m, --merge_density`
- `--structure_density`
- `--cor_plot`
- `--use_rna_seq`

## Example Usage

### Human or yeast with default-style analysis

```bash
mitoribopy \
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

This repository also includes package migration notes and release materials under [docs](docs).

## License

MIT. See [LICENSE](LICENSE).
