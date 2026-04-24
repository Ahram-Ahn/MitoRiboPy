# MitoRiboPy

MitoRiboPy is a Python package for mitochondrial ribosome profiling (mt-Ribo-seq) analysis. Starting in v0.3.0 it spans the full pipeline from raw FASTQ through translation-efficiency integration with paired RNA-seq; v0.4.0 makes every per-sample decision (kit, dedup, offsets) independent so mixed-library batches just work:

- `mitoribopy align` &mdash; FASTQ &rarr; BAM &rarr; BED6 + per-sample read counts (cutadapt + bowtie2 + umi_tools + pysam). **Per-sample kit detection and dedup** in v0.4.0.
- `mitoribopy rpf` &mdash; BED/BAM &rarr; offsets, translation-profile, codon usage, coverage plots. **Per-sample offsets by default** in v0.4.0; combined still available as a diagnostic.
- `mitoribopy rnaseq` &mdash; DE table (DESeq2 / Xtail / Anota2Seq) + rpf outputs &rarr; TE and &Delta;TE tables + plots, SHA256 reference-consistency gate
- `mitoribopy all` &mdash; end-to-end orchestrator with a shared YAML config and a composed `run_manifest.json`

## Highlights

- **Per-sample resolution everywhere**: each FASTQ gets its own auto-detected kit, its own UMI-aware-or-skip dedup decision, and its own offset table. Mixed-kit / mixed-UMI runs are first-class. Per-sample provenance lands in `kit_resolution.tsv` and `run_settings.json -> per_sample`.
- **Per-sample offsets by default** (`--offset_mode per_sample`): downstream codon usage and coverage plots use each sample's own offsets, so inter-sample drift no longer biases the combined output. The combined offset table is still emitted for QC, and `offset_drift_<align>.svg` makes drift visible at a glance. Pass `--offset_mode combined` to force the v0.3.x global behaviour.
- **Both P-site AND A-site outputs by default** (`--analysis_sites both`): codon usage and coverage plots are generated for both sites side by side under per-site subdirectories. `--analysis_sites p` or `--analysis_sites a` restricts to one site.
- **`--kit-preset auto`** (the new default) lets per-sample detection pick the kit; an explicit preset becomes a per-sample fallback only when detection fails.
- Subcommand CLI (`align` / `rpf` / `rnaseq` / `all`) with shared `--config`, `--dry-run`, `--threads`, `--log-level`.
- Config files in JSON, YAML, or TOML (auto-detected by path suffix); YAML is the recommended invocation surface.
- **Polymorphic `align.fastq` YAML key**: pass a directory string OR an explicit list of paths; the directory is glob-scanned for `*.fq`, `*.fq.gz`, `*.fastq`, `*.fastq.gz`.
- Kit-aware FASTQ trimming: `truseq_smallrna`, `nebnext_smallrna`, `nebnext_ultra_umi`, `qiaseq_mirna`, or explicit `--adapter`. Strict mode (`--adapter-detection strict`) hard-fails any sample whose data disagrees with an explicit preset.
- Strand-aware mt-transcriptome alignment (`--library-strandedness forward` by default) so ND5 / ND6 antisense overlap is resolved by construction on Path A (transcriptome reference).
- Deduplication safe by default: `--dedup-strategy auto` resolves to `umi-tools` per sample whenever a UMI is present and `skip` otherwise; `mark-duplicates` is gated behind a long confirmation flag because coordinate-only dedup destroys codon-occupancy signal on low-complexity mt-Ribo-seq libraries.
- BAM input to `rpf` via pysam (no samtools / bedtools PATH dependency).
- SHA256 reference-consistency gate on `rnaseq`: Ribo-seq and RNA-seq sides must be aligned to the identical transcript reference; mismatches are a hard fail.
- Strain presets (`-s h` / `-s y` / `-s vm` / `-s ym` / `-s custom`): human + yeast ship a built-in annotation; `vm` / `ym` / `custom` pick up the matching codon table but require user-supplied `--annotation_file` and an explicit `-rpf` range.
- Footprint-class defaults (`--footprint_class monosome|disome|custom`): monosome uses the canonical 28-34 nt (vertebrate) / 37-41 nt (yeast) RPF window; disome widens to 60-90 nt / 65-95 nt for collided-ribosome studies.
- End-specific 5'/3' offset selection, P-site vs A-site workflows, bicistronic ATP8/ATP6 and ND4L/ND4 handling.
- Custom organism support via `--annotation_file`, `--codon_tables_file`, `--codon_table_name`, `--start_codons`.
- Persistent per-run logging in `<output>/mitoribopy.log`.
- Consistent terminal + file progress reporting for `align` and `rpf`.
- Provenance: every stage writes a `run_settings.json` (with per-sample resolution embedded for `align`); `mitoribopy all` composes them into `run_manifest.json`.

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

### Recommended path: YAML config + `mitoribopy all`

The shortest route from raw FASTQ to translation-profile and coverage outputs is a single YAML file plus one command. This is the path the README demonstrates first; the bash-wrapper alternative further down is for batch / cluster jobs where every flag should be explicit.

```bash
# 1. Conda env with cutadapt / bowtie2 / umi_tools / samtools / pysam.
conda env create -f docs/environment/environment.yml
conda activate mitoribopy

# 2. Drop a working YAML template next to your data and fill in the paths.
mitoribopy all --print-config-template > pipeline_config.yaml
$EDITOR pipeline_config.yaml

# 3. Optional: dry-run prints the per-stage argv so you can review.
mitoribopy all --config pipeline_config.yaml --output results/ --dry-run

# 4. Run.
mitoribopy all --config pipeline_config.yaml --output results/
```

A working `pipeline_config.yaml` for a typical human mt-Ribo-seq run looks like this (annotated):

```yaml
align:
  # Per-sample auto detection; explicit kit_preset becomes the fallback.
  kit_preset: auto
  adapter_detection: auto         # auto | off | strict (per-sample)
  library_strandedness: forward
  # Pass a DIRECTORY (string) and every *.fq.gz / *.fastq.gz / *.fq /
  # *.fastq inside is picked up. To restrict to specific files, pass
  # `fastq:` as a list instead.
  fastq: input_data/
  contam_index: input_data/indexes/rrna_contam
  mt_index: input_data/indexes/mt_tx
  mapq: 10
  min_length: 15
  max_length: 45
  dedup_strategy: auto            # umi-tools per sample if UMI, else skip

rpf:
  strain: h                       # human mt-mRNA reference + codon table
  fasta: input_data/human-mt-mRNA.fasta
  rpf: [29, 34]                   # filtered RPF range
  align: stop                     # anchor offsets at the stop codon
  offset_type: "5"                # offsets reported from the read 5' end
  offset_site: p                  # selection coordinate space (P-site)
  offset_pick_reference: p_site
  offset_mode: per_sample         # per-sample offsets drive downstream
  analysis_sites: both            # write BOTH P-site and A-site outputs
  min_5_offset: 10
  max_5_offset: 22
  min_3_offset: 10
  max_3_offset: 22
  offset_mask_nt: 5
  plot_format: svg
  merge_density: true
```

### Alternative: bash wrapper for batch / cluster jobs

When every flag should be visible in the job script (e.g. for a cluster scheduler that captures stdout / stderr per task), wrap the same YAML invocation in a thin shell wrapper. The bash form composes well with `set -euo pipefail`, exit-code capture, and per-job log rotation.

```bash
#!/usr/bin/env bash
set -uo pipefail

ENV_BIN=/path/to/conda/envs/mitoribopy/bin
export PATH="$ENV_BIN:$PATH"
ROOT=/path/to/your/run
cd "$ROOT"

OUT=results/full_run
mkdir -p "$OUT"

mitoribopy all \
  --config pipeline_config.yaml \
  --output "$OUT" \
  --threads 8
RC=$?

echo
echo "================ kit_resolution.tsv ================"
cat "$OUT/align/kit_resolution.tsv"

echo
echo "================ read_counts.tsv ================"
cat "$OUT/align/read_counts.tsv"

echo "$RC" > "$OUT.exitcode"
exit $RC
```

The YAML form is canonical because it is self-documenting and loads exactly the same way `mitoribopy all` reads it programmatically. The bash wrapper only exists for environments that prefer one-file job scripts.

### Inspecting and overriding the config

```bash
# Print the full per-stage --help (orchestrator --help shows only its own flags).
mitoribopy all --show-stage-help align
mitoribopy all --show-stage-help rpf
mitoribopy all --show-stage-help rnaseq
```

### Strain presets

| `-s` | Organism / codon table | Ships annotation? | Ships `-rpf` default? |
|------|------------------------|:-:|:-:|
| `h`      | Human mt (`vertebrate_mitochondrial`) | yes | yes (28-34 nt monosome) |
| `y`      | Yeast mt (`yeast_mitochondrial`)      | yes | yes (37-41 nt monosome) |
| `vm`     | Any vertebrate mt (`vertebrate_mitochondrial`) | no | no &mdash; pass `--annotation_file` + `-rpf` |
| `ym`     | Any fungus with yeast-mito code (`yeast_mitochondrial`) | no | no &mdash; pass `--annotation_file` + `-rpf` |
| `custom` | Fully user-specified                  | no | no &mdash; also requires `--codon_tables_file` or `--codon_table_name` |

Pair `-s` with `--footprint_class`:

| `--footprint_class` | RPF window default | `--unfiltered_read_length_range` default | Use for |
|---------------------|--------------------|------------------------------------------|---------|
| `monosome` (default) | h/vm: 28-34, y/ym: 37-41 | 15-50 | Standard single-ribosome footprints |
| `disome`             | h/vm: 60-90, y/ym: 65-95 | 40-110 | Collided-ribosome studies (e.g. eIF5A depletion, stalling) |
| `custom`             | user must pass `-rpf` | unchanged | Any non-standard footprint class |

An explicit `-rpf MIN MAX` or `--unfiltered_read_length_range MIN MAX` always wins over the footprint-class default.

### `mitoribopy rpf` &mdash; BED/BAM through the analysis pipeline

```bash
mitoribopy rpf \
  -s h \
  -f <reference.fa> \
  --directory <ribo_bed_dir> \
  -rpf 29 34 \
  --offset_mode per_sample \
  --analysis_sites both \
  --output <results_dir>
```

`--offset_mode per_sample` (the default) selects offsets independently for each sample so inter-sample offset drift no longer biases pooled output. The combined-across-samples table is still written as a diagnostic and surfaced in `offset_drift_<align>.svg`. Pass `--offset_mode combined` to recover the v0.3.x global behaviour.

`--analysis_sites both` (the default) writes parallel P-site and A-site downstream outputs into per-site subdirectories (`translation_profile_p/`, `translation_profile_a/`, `coverage_profile_plots_p/`, `coverage_profile_plots_a/`). Pass `--analysis_sites p` or `--analysis_sites a` to restrict to one site (legacy single-site layout).

Plain `mitoribopy <flags>` still works in v0.3.x but routes to `rpf` with a deprecation warning. Use the explicit subcommand form.

### `mitoribopy align` &mdash; FASTQ &rarr; BAM + BED

```bash
mitoribopy align \
  --library-strandedness forward \
  --fastq-dir <fastqs_dir> \
  --contam-index <bowtie2_rRNA_index_prefix> \
  --mt-index <bowtie2_mt_transcriptome_index_prefix> \
  --output <align_results_dir>
```

`--kit-preset` defaults to `auto`: each input FASTQ is scanned at run start and its kit is resolved independently. To pin an explicit fallback for samples whose data cannot be auto-identified, pass `--kit-preset truseq_smallrna` (or one of the other presets); to refuse anything but a hard match, pass `--adapter-detection strict`. Use `--kit-preset custom --adapter <SEQ>` when your library isn't one of the built-in presets. The per-sample resolution table lands at `<align_results_dir>/kit_resolution.tsv` and is also embedded under `run_settings.json -> per_sample`.

External tools (`cutadapt`, `bowtie2`, `umi_tools` for UMI samples) must be on `$PATH`; see [docs/environment/environment.yml](docs/environment/environment.yml) for a ready-made bioconda env. `umi_tools` is only required when at least one resolved sample has a UMI.

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

#### `align` stage

- `--kit-preset {auto,truseq_smallrna,nebnext_smallrna,nebnext_ultra_umi,qiaseq_mirna,custom}` &mdash; `auto` (default) runs per-sample detection; explicit preset becomes a per-sample fallback.
- `--adapter <SEQ>` &mdash; explicit adapter sequence (required for `--kit-preset custom`); otherwise an optional fallback used only when detection fails.
- `--adapter-detection {auto,off,strict}` &mdash; per-sample policy. `strict` hard-fails any sample whose detection disagrees with `--kit-preset` or yields no match.
- `--dedup-strategy {auto,umi-tools,skip,mark-duplicates}` &mdash; `auto` resolves per sample (umi-tools when UMIs present, skip otherwise).
- `--library-strandedness {forward,reverse,unstranded}`
- `--min-length`, `--max-length`, `--quality`, `--mapq`, `--seed`

#### `rpf` stage

- `--align start|stop`
- `--offset_type 5|3` &mdash; report offsets from the read 5' or 3' end.
- `--offset_site p|a` &mdash; coordinate space for the SELECTED OFFSETS table only. Use `--analysis_sites` to control which downstream outputs are generated.
- `--offset_pick_reference p_site|selected_site`
- `--offset_mode per_sample|combined` &mdash; per-sample offsets drive downstream by default; `combined` matches v0.3.x.
- `--analysis_sites p|a|both` &mdash; which downstream sites to generate (default `both` writes parallel P and A outputs).
- `--min_5_offset`, `--max_5_offset`, `--min_3_offset`, `--max_3_offset`, `--offset_mask_nt`
- `--read_counts_file`, `--read_counts_sample_col`, `--read_counts_reference_col`, `--read_counts_reads_col`
- `--unfiltered_read_length_range <min> <max>`
- `--rpm_norm_mode total|mt_mrna`
- `--plot_format png|pdf|svg`
- `-m, --merge_density`
- `--structure_density`
- `--cor_plot`
- `--use_rna_seq` (deprecated; use the `rnaseq` subcommand)

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

### FASTQ (primary input to the `align` stage)

Accepted file extensions: `*.fq`, `*.fq.gz`, `*.fastq`, `*.fastq.gz`. Both gzipped and uncompressed are auto-detected.

Two ways to point the pipeline at your FASTQs:

1. **Directory** (recommended): pass a directory containing every input FASTQ.
   - CLI: `--fastq-dir input_data/`
   - YAML: `align.fastq: input_data/` (a single string is treated as a directory)
2. **Explicit list**: name each FASTQ.
   - CLI: `--fastq sampleA.fq.gz --fastq sampleB.fq.gz` (repeatable)
   - YAML: `align.fastq: [sampleA.fq.gz, sampleB.fq.gz]` (a list is treated as explicit paths)

Sample names are derived from the FASTQ filename with the extension stripped (`WT_R1.fq.gz` → `WT_R1`). The same name flows through every per-sample table (`read_counts.tsv`, `kit_resolution.tsv`, downstream profile and codon usage subdirs).

#### Per-sample kit detection

When `--kit-preset` is `auto` (the default), the pipeline scans the head of every input FASTQ at run start, picks the matching kit per sample, and records its decision in `<output>/kit_resolution.tsv`. Mixed-kit batches (TruSeq + NEBNext + qiaseq in one run) and mixed UMI / non-UMI batches just work: the dedup strategy is also resolved per sample so UMI-bearing samples take the umi-tools path while non-UMI samples skip dedup.

Pass an explicit preset (e.g. `--kit-preset truseq_smallrna`) when you want a deterministic fallback for samples whose data cannot be auto-identified, or `--adapter-detection strict` to refuse to continue on any disagreement.

### BED (input to the `rpf` stage if you already have aligned BEDs)

Expected columns:

1. `chrom`
2. `start`
3. `end`

Additional BED columns are tolerated. Coordinates are treated as standard 0-based, end-exclusive intervals.

When the `align` stage runs first (or you use `mitoribopy all`), `mitoribopy rpf` consumes its `<align>/bed/` output automatically &mdash; you never need to handle BED files by hand.

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

Typical output structure for an `mitoribopy all` run with the v0.4.0 defaults (`--offset_mode per_sample`, `--analysis_sites both`):

```text
<output>/
  run_manifest.json
  align/
    mitoribopy.log
    read_counts.tsv               # per-stage counts
    kit_resolution.tsv            # per-sample kit + dedup decisions
    run_settings.json             # includes per_sample[] block
    trimmed/, contam_filtered/, aligned/, deduped/, bed/
  rpf/
    mitoribopy.log
    rpf_counts.tsv                # per-sample per-gene RPF counts
    run_settings.json             # includes reference_checksum
    plots_and_csv/
      offset_<align>.csv                 # combined diagnostic
      p_site_offsets_<align>.csv         # combined diagnostic
      offset_drift_<align>.svg           # per-sample offset comparison
      per_sample/
        <sample>/
          offset_<align>.csv             # per-sample enrichment
          p_site_offsets_<align>.csv     # per-sample selected offsets
    translation_profile_p/        # P-site outputs (when analysis_sites=both)
      <sample>/
        footprint_density/
        translating_frame/
        codon_usage/
        debug_csv/
    translation_profile_a/        # A-site outputs (when analysis_sites=both)
      <sample>/...
    coverage_profile_plots_p/     # P-site coverage plots
    coverage_profile_plots_a/     # A-site coverage plots
    structure_density/            # if --structure_density (uses P-site by default)
    codon_correlation/            # if --cor_plot
  rnaseq/                         # if rnaseq config supplied
```

When `--analysis_sites=p` or `=a` (single site), the per-site directory names collapse back to the legacy v0.3 layout: `<output>/rpf/<sample>/...` and `<output>/rpf/coverage_profile_plots/`.

Key outputs:

- **Per-sample enrichment + selection tables** under `plots_and_csv/per_sample/<sample>/`.
- **`offset_drift_<align>.svg`** &mdash; per-sample 5'/3' offset by read length plus the combined diagnostic dashed line. Read this first to spot inter-sample drift.
- **`kit_resolution.tsv`** in `<align>/` &mdash; per-sample detected kit, applied kit, match rate, and dedup decision. The first thing to check when a particular sample looks off.
- **Footprint-density CSVs** for P-site, A-site, and E-site under each sample's `footprint_density/`.
- **Frame-usage summaries** and **transcript-level + total codon-usage summaries** per site.
- **RPM and raw coverage-profile plots** per site.
- **CDS-aware codon-binned coverage plots** (`*_codon/`, 3 nt combined per codon).
- Optional **structure-density exports** from footprint-density tables.

## Important Runtime Notes

- `--offset_type 5|3`: measure offsets from the read 5' or 3' end.
- `--offset_site p|a`: coordinate space for the SELECTED OFFSETS table only. Does not control which downstream outputs are generated.
- `--offset_pick_reference p_site|selected_site`: how the best offset is chosen.
- `--offset_mode per_sample|combined`: per-sample is the default; combined matches v0.3.x behaviour.
- `--analysis_sites p|a|both`: which downstream sites to generate. Default `both` writes parallel P and A outputs.
- `--min_5_offset`, `--max_5_offset`, `--min_3_offset`, `--max_3_offset`: recommended end-specific selection bounds.
- `--offset_mask_nt`: mask near-anchor bins from enrichment summaries and plots.
- `--rpm_norm_mode total|mt_mrna`: read-count normalization mode.
- `--structure_density`: export log2 and scaled density values from footprint-density tables.

For the full interface, run:

```bash
mitoribopy --help
```

## Troubleshooting

**~99 % of reads disappear at the trim step for one sample, `post_trim` is a tiny fraction of `total_reads` in `read_counts.tsv`.**
That sample's auto-detected kit was wrong, or auto-detection silently
fell back to the global `--kit-preset` for a sample whose adapter is
something else. Open `<output>/align/kit_resolution.tsv`: the
`detected_kit`, `match_rate`, and `applied_kit` columns reveal what the
detector saw and what was actually applied. Re-run with `--kit-preset
<right_kit>` (a per-sample fallback), or pass `--adapter-detection
strict` to refuse to continue on any disagreement.

**Two of my samples have UMIs and two don't &mdash; how do I configure dedup?**
You don't. With the default `--dedup-strategy auto` and `--kit-preset
auto`, the per-sample resolver picks `umi-tools` for the UMI samples
and `skip` for the non-UMI samples in the same run. The dedup decisions
are recorded per sample in `kit_resolution.tsv` and embedded under
`run_settings.json -> per_sample`.

**Per-sample offsets diverge by more than 1-2 nt &mdash; which sample is the outlier?**
Open `<rpf_output>/plots_and_csv/offset_drift_<align>.svg`. Each sample
has its own per-read-length 5' and 3' offset bar; the combined
diagnostic is overlaid as a dashed line. Outliers are visible by eye in
seconds. If a single sample drifts by more than 2 nt at most read
lengths, inspect its enrichment heatmap under
`plots_and_csv/per_sample/<sample>/` &mdash; usually the cause is low
coverage / a bad library at that sample, not a real biological signal.

**A-site vs P-site &mdash; which downstream files do I read?**
`--analysis_sites both` (the default) writes parallel outputs:
- P-site results live in `translation_profile_p/` and `coverage_profile_plots_p/`.
- A-site results live in `translation_profile_a/` and `coverage_profile_plots_a/`.

Restrict to one site with `--analysis_sites p` or `--analysis_sites a`
to suppress the other and use the legacy single-directory layout.
`--offset_site` only affects the SELECTED OFFSETS table's coordinate
space, not which downstream outputs exist; use `--analysis_sites` to
control the latter.

**Filtered BED is empty &rarr; "no data remained after BED filtering".**
Either your RPF range does not cover the actual read-length
distribution (open the per-sample `*_read_length_distribution.svg`;
the shaded band shows the currently selected window) or every
mapped read has been filtered earlier by MAPQ / contaminant
subtraction. Widen `-rpf MIN MAX`, try `--footprint_class disome` if
you are studying collided ribosomes, or lower `--mapq`.

**Offset selection produced no rows &rarr; `p_site_offsets_*.csv` is empty.**
The `--min_5_offset` / `--max_5_offset` window (default 10-22 nt) did
not overlap the enrichment peak. Re-open the
`offset_enrichment_heatmap_*.svg` and widen the window explicitly.

**`RPM` is 0 for every sample.**
Either `--read_counts_file` was not passed, or the file does not
contain entries for the sample name(s) the pipeline inferred from
the BED filenames. The `[QC] WARNING: No total read-count entry found
for sample(s): ...` log line lists the samples that missed the
lookup. Add a matching row to the counts file and re-run.

**`--show-stage-help` output is too dense.**
It is the full argparse `--help` for that stage. Pair it with
`mitoribopy all --print-config-template` to get a pre-populated YAML
and only override the keys you care about.

**Reference-consistency gate failure in `rnaseq`.**
The reference FASTA you just passed does not hash-match the one the
prior `rpf` stage recorded in its `run_settings.json`. You must
re-align both sides against the identical transcript set.

## Development

Run the test suite with:

```bash
PYTHONPATH=src pytest
```

This repository also includes package migration notes and release materials under [docs/README.md](docs/README.md).

## License

MIT. See [LICENSE](LICENSE).
