# MitoRiboPy

Mitochondrial ribosome profiling (mt-Ribo-seq) analysis, end to end.

MitoRiboPy is a Python package + CLI for analysing mt-Ribo-seq data from raw FASTQ all the way through translation-efficiency integration with paired RNA-seq. v0.4.1 makes every per-sample decision (kit, dedup, offsets) independent so mixed-library batches just work.

The package is built around four subcommands:

| Subcommand | What it does |
|---|---|
| `mitoribopy align` | FASTQ → BAM → BED6 + per-sample read counts (cutadapt + bowtie2 + umi_tools + pysam) |
| `mitoribopy rpf` | BED/BAM → offsets, translation profile, codon usage, coverage plots |
| `mitoribopy rnaseq` | DE table (DESeq2 / Xtail / Anota2Seq) + rpf outputs → TE and ΔTE tables + plots, with a SHA256 reference-consistency gate |
| `mitoribopy all` | End-to-end orchestrator that runs align + rpf + (optional) rnaseq from one YAML config and writes a composed `run_manifest.json` |

---

## Table of contents

1. [What MitoRiboPy is for](#what-mitoribopy-is-for)
2. [Pipeline overview](#pipeline-overview)
3. [Installation](#installation)
4. [Quick start](#quick-start)
5. [Input files](#input-files)
6. [How to run — YAML vs shell wrapper](#how-to-run--yaml-vs-shell-wrapper)
7. [Strain presets and footprint classes](#strain-presets-and-footprint-classes)
8. [Subcommand reference](#subcommand-reference)
   - [`mitoribopy align`](#mitoribopy-align)
   - [`mitoribopy rpf`](#mitoribopy-rpf)
   - [`mitoribopy rnaseq`](#mitoribopy-rnaseq)
   - [`mitoribopy all`](#mitoribopy-all)
9. [Output overview](#output-overview)
10. [Custom organisms](#custom-organisms)
11. [Built-in references](#built-in-references)
12. [Examples](#examples)
13. [Logs and provenance](#logs-and-provenance)
14. [Development](#development)
15. [License](#license)

---

## What MitoRiboPy is for

MitoRiboPy is a focused tool for the 13 mt-mRNAs of human mitochondria (or 8 mt-mRNAs in yeast, plus configurable codon tables for any other mitochondrion). It ships:

- **Per-sample adapter detection** with auto-fallback to the right kit. Mixed-kit and mixed-UMI batches resolve each sample independently. Pre-trimmed FASTQs (e.g. SRA-deposited data) are auto-detected and routed through cutadapt with no `-a` flag.
- **Per-sample offset selection** so inter-sample drift in the canonical 12–15 nt 5' P-site offset doesn't bias your downstream codon-usage tables. A combined-across-samples diagnostic is still emitted and an `offset_drift_<align>.svg` plot makes drift visible at a glance.
- **Both P-site and A-site downstream outputs** by default, side by side under per-site subdirectories. No more ambiguity about which output corresponds to which site.
- **Strict reference-consistency gate** for the optional translation-efficiency (TE / ΔTE) integration with paired RNA-seq: Ribo-seq and RNA-seq must hash to the identical transcript reference or the run aborts.
- **Strain-aware defaults**: built-in human (`-s h`) and yeast (`-s y`) annotations + codon tables, plus generic `vm` (vertebrate-mito) and `ym` (yeast-mito-codon-code) presets, and `custom` for any other organism.

What MitoRiboPy is **not**:

- Not a general-purpose nuclear Ribo-seq pipeline. The defaults, references, and dedup heuristics are calibrated for the low-complexity 13-mRNA mt universe.
- Not a DE engine. RNA-seq DE (DESeq2 / Xtail / Anota2Seq) is run externally on the full transcriptome; MitoRiboPy consumes the resulting table and does TE / ΔTE on the mt-mRNA subset.

---

## Pipeline overview

```
                        FASTQ (one or more samples; mixed kits OK)
                           │
                           ▼
      ┌────────────────────────────────────────────────┐
      │   align stage  (mitoribopy align)              │
      │                                                │
      │   per-sample adapter detection → cutadapt trim │
      │   → bowtie2 contam subtract                    │
      │   → bowtie2 mt-transcriptome align             │
      │   → MAPQ filter (NUMT suppression)             │
      │   → per-sample dedup (umi-tools | skip)        │
      │   → BAM → BED6                                 │
      └────────────────────────────────────────────────┘
                           │
              BED6 + read_counts.tsv + kit_resolution.tsv
                           │
                           ▼
      ┌────────────────────────────────────────────────┐
      │   rpf stage  (mitoribopy rpf)                  │
      │                                                │
      │   BED filter on RPF length window              │
      │   → offset enrichment (per sample + combined)  │
      │   → offset selection (per sample + combined)   │
      │   → translation profile (P + A site)           │
      │   → codon usage  (per sample, per site)        │
      │   → coverage plots (per sample, per site)      │
      │   → optional: structure density, codon         │
      │     correlation                                │
      └────────────────────────────────────────────────┘
                           │
        rpf_counts.tsv + run_settings.json (with reference_checksum)
                           │
                           ▼
      ┌────────────────────────────────────────────────┐
      │   rnaseq stage  (mitoribopy rnaseq) — optional │
      │                                                │
      │   DE table (DESeq2 / Xtail / Anota2Seq)        │
      │   + rpf_counts.tsv                             │
      │   → SHA256 reference-consistency gate          │
      │   → te.tsv (per sample × gene)                 │
      │   → delta_te.tsv (per gene)                    │
      │   → mrna_vs_rpf scatter, ΔTE volcano           │
      └────────────────────────────────────────────────┘
```

---

## Installation

From the repository root:

```bash
python -m pip install -e .
```

For development and tests:

```bash
python -m pip install -e ".[dev]"
```

Confirm the CLI is available:

```bash
mitoribopy --version    # should print "MitoRiboPy 0.4.1"
mitoribopy --help
```

If you prefer not to install yet:

```bash
PYTHONPATH=src python -m mitoribopy --help
```

### External tool dependencies

MitoRiboPy shells out to a small set of standard bioinformatics tools. All of them must be on `$PATH` for a real run:

| Tool | Used by | Required when |
|---|---|---|
| `cutadapt` | `align` | always (length + quality filter even for pre-trimmed data) |
| `bowtie2` + `bowtie2-build` | `align` | always |
| `umi_tools` | `align` | at least one sample's resolved kit has UMIs |
| `picard` | `align` | only when `--dedup-strategy mark-duplicates` is opted into |
| `pysam` (Python lib) | `align`, `rpf` | always (installed automatically via `pip`) |
| `samtools` | optional | recommended for inspecting outputs; not required |

The bioconda environment under [docs/environment/environment.yml](docs/environment/environment.yml) installs everything in one command:

```bash
conda env create -f docs/environment/environment.yml
conda activate mitoribopy
```

---

## Quick start

The shortest path from raw FASTQ to translation-profile + coverage outputs is one YAML file plus one command.

```bash
# 1. Drop a working YAML template next to your data and fill in the paths.
mitoribopy all --print-config-template > pipeline_config.yaml
$EDITOR pipeline_config.yaml

# 2. Optional: dry-run prints the per-stage argv so you can review.
mitoribopy all --config pipeline_config.yaml --output results/ --dry-run

# 3. Run.
mitoribopy all --config pipeline_config.yaml --output results/
```

A working `pipeline_config.yaml` for a typical human mt-Ribo-seq run looks like this (annotated):

```yaml
align:
  # Per-sample auto detection; explicit kit_preset becomes a fallback.
  kit_preset: auto                # auto | illumina_smallrna | illumina_truseq |
                                  # illumina_truseq_umi | qiaseq_mirna |
                                  # pretrimmed | custom
  adapter_detection: auto         # auto | off | strict
  library_strandedness: forward
  # Pass a directory string (auto-glob of *.fq, *.fq.gz, *.fastq, *.fastq.gz)
  # OR an explicit list of paths.
  fastq: input_data/seq
  contam_index: input_data/indexes/rrna_contam
  mt_index: input_data/indexes/mt_tx
  mapq: 10
  min_length: 15
  max_length: 45
  dedup_strategy: auto            # umi-tools per sample if UMI, else skip

rpf:
  strain: h                       # human mt-mRNA reference + codon table
  fasta: input_data/human-mt-mRNA.fasta
  rpf: [29, 34]                   # filtered RPF length range
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

After the run, you'll have:

```text
results/
  align/    bed/, deduped/, kit_resolution.tsv, read_counts.tsv, run_settings.json
  rpf/      plots_and_csv/, translation_profile_p/, translation_profile_a/,
            coverage_profile_plots_p/, coverage_profile_plots_a/, rpf_counts.tsv
  run_manifest.json
```

See [Output overview](#output-overview) for the full directory tree.

---

## Input files

### FASTQ (primary input to `align`)

Accepted file extensions: `*.fq`, `*.fq.gz`, `*.fastq`, `*.fastq.gz`. Both gzipped and uncompressed are auto-detected.

Two ways to point the pipeline at your FASTQs:

1. **Directory** (recommended): pass a directory containing every input FASTQ.
   - CLI: `--fastq-dir input_data/`
   - YAML: `align.fastq: input_data/` (a single string is treated as a directory)
2. **Explicit list**: name each FASTQ.
   - CLI: `--fastq sampleA.fq.gz --fastq sampleB.fq.gz` (repeatable)
   - YAML: `align.fastq: [sampleA.fq.gz, sampleB.fq.gz]` (a list is treated as explicit paths)

Sample names are derived from the FASTQ filename with the extension stripped (`WT_R1.fq.gz` → `WT_R1`). The same name flows through every per-sample table (`read_counts.tsv`, `kit_resolution.tsv`, downstream profile and codon usage subdirs).

### BED (input to `rpf` if you already have aligned BEDs)

Expected columns:

1. `chrom`
2. `start`
3. `end`

Additional BED columns are tolerated. Coordinates are 0-based, end-exclusive intervals (standard BED).

When the `align` stage runs first (or you use `mitoribopy all`), `mitoribopy rpf` consumes `<align>/bed/` automatically — you never need to handle BED files by hand.

### BAM (alternative input to `rpf`)

`mitoribopy rpf --directory <dir>` accepts BAM files mixed with BED. BAMs are auto-converted to BED6 under `<output>/bam_converted/` via pysam. The `--bam_mapq` flag (default 10) filters BAM reads on MAPQ before conversion to suppress NUMT cross-talk; set to 0 to disable.

### Reference FASTA

One FASTA record per mt-mRNA. Headers must match the `sequence_name` column of the annotation CSV (or any of its `sequence_aliases`). The built-in human and yeast annotations cover the canonical mt-mRNAs and ship under `src/mitoribopy/data/`.

For total-genome FASTAs (rare in mt-Ribo-seq), use `--annotation_file` to map FASTA records onto your own annotation rows.

### Annotation CSV (custom organisms only)

Required columns:

- `transcript` — logical CDS name used in frame and codon outputs
- `l_tr` — full transcript length
- `l_utr5` — 5' UTR length
- `l_utr3` — 3' UTR length

Optional columns:

- `l_cds` — if omitted, computed as `l_tr − l_utr5 − l_utr3`
- `sequence_name` — FASTA/BED sequence ID this row maps to
- `sequence_aliases` — alternate FASTA/BED names separated by semicolons
- `display_name` — controls plot titles and grouped transcript labels

Templates: [examples/custom_reference/annotation_template.csv](examples/custom_reference/annotation_template.csv).

### Codon-table JSON (custom organisms only)

Two formats are supported:

- One flat 64-codon mapping (`{"AAA": "K", "AAC": "N", ...}`)
- A dictionary of named 64-codon mappings (`{"my_table": {"AAA": "K", ...}, "another": {...}}`)

Choose with `--codon_table_name` when multiple are present. Templates: [examples/custom_reference/codon_tables_template.json](examples/custom_reference/codon_tables_template.json).

### Read-count table (optional, for RPM normalization)

`.csv`, `.tsv`, and `.txt` accepted; delimiter is auto-detected. Column matching is flexible and case-insensitive, with positional fallback:

- column 1: sample name
- column 2: reference (used when `--rpm_norm_mode mt_mrna`)
- column 3: read count

When `mitoribopy all` runs `align` first, the read-count table is auto-wired from `<run_root>/align/read_counts.tsv`.

---

## How to run — YAML vs shell wrapper

### Recommended: YAML config

```bash
mitoribopy all --config pipeline_config.yaml --output results/ --threads 8
```

This is the canonical invocation. The YAML is self-documenting, version-controllable, and loads exactly the same way the CLI reads it programmatically.

### Alternative: bash wrapper for batch / cluster jobs

When every flag should be visible in the job script (e.g. for a cluster scheduler that captures stdout/stderr per task), wrap the YAML invocation in a thin shell wrapper:

```bash
#!/usr/bin/env bash
set -uo pipefail

ENV_BIN=/path/to/conda/envs/mitoribopy/bin
export PATH="$ENV_BIN:$PATH"

ROOT=/path/to/project
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

### Direct subcommand invocation

You can run any single stage directly without going through `mitoribopy all`. This is useful when you only have BED inputs (skip `align`), or when you want to iterate on `rpf` parameters without re-running alignment.

```bash
# Just align
mitoribopy align --kit-preset auto --fastq-dir fastqs/ \
  --contam-index idx/rrna --mt-index idx/mt --output results/align/

# Just rpf, against an existing BED dir
mitoribopy rpf -s h -f ref.fa --directory bed/ -rpf 29 34 --output results/rpf/

# Just rnaseq, against existing rpf output + a DE table
mitoribopy rnaseq --de-table de.tsv --gene-id-convention hgnc \
  --ribo-dir results/rpf --reference-gtf ref.fa --output results/rnaseq/
```

---

## Strain presets and footprint classes

### Strain (`-s` / `--strain`)

| Value | Organism / codon table | Ships annotation? | Ships `-rpf` default? |
|---|---|:-:|:-:|
| `h` | Human mt (`vertebrate_mitochondrial`) | ✓ | ✓ (28–34 nt monosome) |
| `y` | Yeast mt (`yeast_mitochondrial`) | ✓ | ✓ (37–41 nt monosome) |
| `vm` | Any vertebrate mt (`vertebrate_mitochondrial`) | ✗ | ✗ — pass `--annotation_file` + `-rpf` |
| `ym` | Any fungus with yeast-mito code (`yeast_mitochondrial`) | ✗ | ✗ — pass `--annotation_file` + `-rpf` |
| `custom` | Fully user-specified | ✗ | ✗ — also requires `--codon_tables_file` or `--codon_table_name` |

### Footprint class (`--footprint_class`)

Pair `-s` with `--footprint_class` to pick sensible RPF and unfiltered-length defaults:

| Value | RPF window default | `--unfiltered_read_length_range` default | Use for |
|---|---|---|---|
| `monosome` (default) | h/vm: 28–34, y/ym: 37–41 | 15–50 | Standard single-ribosome footprints |
| `disome` | h/vm: 60–90, y/ym: 65–95 | 40–110 | Collided-ribosome studies (eIF5A depletion, stalling) |
| `custom` | user must pass `-rpf` | unchanged | Any non-standard footprint class |

An explicit `-rpf MIN MAX` or `--unfiltered_read_length_range MIN MAX` always wins over the footprint-class default.

---

## Subcommand reference

Every subcommand inherits these shared options:

| Flag | Default | Description |
|---|---|---|
| `--config PATH` | — | Configuration file (.json, .yaml, .yml, or .toml). CLI flags override values from the file. |
| `--dry-run` | off | Print planned actions and exit 0 without executing. |
| `--threads N` | 1 | Preferred thread count; exports `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `MITORIBOPY_THREADS`. |
| `--log-level {DEBUG,INFO,WARNING,ERROR}` | `INFO` | Python logging level for console output. |

---

### `mitoribopy align`

Preprocesses FASTQ inputs into BAM + BED6 + per-sample read counts. Pipeline (per sample):

1. Adapter detection (head-of-FASTQ scan)
2. cutadapt trim (kit-aware; optional UMI extraction)
3. bowtie2 contaminant subtraction
4. bowtie2 mt-transcriptome alignment (Path A: per-mRNA FASTA records)
5. MAPQ filter (NUMT suppression)
6. Deduplication (`umi-tools` for UMI samples, `skip` for no-UMI)
7. BAM → BED6 (strand-aware)

#### Inputs

| Flag | Default | Description |
|---|---|---|
| `--fastq-dir DIR` | — | Directory of `*.fq(.gz)` / `*.fastq(.gz)`. |
| `--fastq PATH` | — | Individual FASTQ; repeatable. Pass `--fastq-dir` OR `--fastq` (or both). |
| `--contam-index BT2_PREFIX` | — | bowtie2 index prefix for the contaminant panel (rRNA + tRNA + spike-ins). Build with `bowtie2-build contaminants.fa <prefix>`. **Required** for non-dry-run. |
| `--mt-index BT2_PREFIX` | — | bowtie2 index prefix for the mt-transcriptome (one FASTA record per mt-mRNA). **Required** for non-dry-run. |
| `--output DIR` | — | Output directory. **Required**. |

#### Library prep

| Flag | Default | Description |
|---|---|---|
| `--kit-preset PRESET` | `auto` | Library-prep adapter family. See [Kit presets](#kit-presets) below for the canonical list. |
| `--adapter SEQ` | — | Explicit 3' adapter sequence. **Required** when `--kit-preset custom`; otherwise an optional fallback used only when detection fails. |
| `--umi-length N` | from preset | Override the kit preset's UMI length. |
| `--umi-position {5p,3p}` | from preset | Override the kit preset's UMI position. |
| `--adapter-detection MODE` | `auto` | `auto` (default), `strict`, or `off`. See [Adapter detection](#adapter-detection) below. |
| `--adapter-detect-reads N` | 5000 | FASTQ reads scanned per sample during detection. |
| `--adapter-detect-min-rate FRAC` | 0.30 | Minimum fraction of scanned reads with adapter signal for the kit to be considered detected. |
| `--adapter-detect-min-len N` | 12 | Adapter prefix length used as the search needle (nt). |
| `--adapter-detect-pretrimmed-threshold FRAC` | 0.05 | When EVERY kit's match rate is at or below this, classify as already-trimmed. |
| `--no-pretrimmed-inference` | off | Disable the auto-fallback to `pretrimmed`; restore the v0.4.0 hard-fail behaviour. |
| `--library-strandedness {forward,reverse,unstranded}` | `forward` | `forward` enforces bowtie2 `--norc`; `reverse` enforces `--nofw`; `unstranded` leaves bowtie2 permissive. |
| `--min-length NT` | 15 | Minimum read length kept after trimming. |
| `--max-length NT` | 45 | Maximum read length kept after trimming. |
| `--quality Q` | 20 | cutadapt `-q` Phred+33 3' quality trim threshold. |

##### Kit presets

| Preset | 3' adapter | UMI | Covers (representative) |
|---|---|---|---|
| `auto` | per-sample detection | per-sample | default; scans every input FASTQ |
| `pretrimmed` | none | none | Already-trimmed FASTQs (SRA-deposited, prior trim step). cutadapt skips `-a` |
| `custom` | user-supplied via `--adapter` | configurable | anything not in the table |
| `illumina_smallrna` | `TGGAATTCTCGGGTGCCAAGG` | none | Illumina TruSeq Small RNA |
| `illumina_truseq` | `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` | none | NEBNext Multiplex Small RNA, TruSeq Stranded Total RNA Gold, Takara SMARTer Stranded Total v3 Pico, Bio-Rad SEQuoia Express Standard, … (any Illumina R1 adapter without a UMI) |
| `illumina_truseq_umi` | `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` | 8 nt 5' | NEBNext Ultra II UMI, Bio-Rad SEQuoia Complete UMI, … |
| `qiaseq_mirna` | `AACTGTAGGCACCATCAAT` | 12 nt 3' | QIAseq miRNA Library Kit |

Legacy vendor names from v0.4.0 (`truseq_smallrna`, `nebnext_smallrna`, `nebnext_ultra_umi`, plus the short-lived `truseq_stranded_total`, `smarter_pico_v3`, `sequoia_express`) are still accepted as aliases — existing YAML configs need no changes.

##### Adapter detection

| Mode | Behaviour |
|---|---|
| `auto` (default) | Scan each FASTQ; pick the matching preset per sample. Samples whose scan fails fall back to the user's `--kit-preset` / `--adapter` if supplied; otherwise fall through to `pretrimmed` when the data shows no adapter signal at all (the v0.4.1 auto-inference). |
| `strict` | Scan; HARD-FAIL any sample whose scan disagrees with an explicit `--kit-preset` or yields no match. Use for batch / CI runs where silent surprises are unacceptable. |
| `off` | Skip the scan entirely; trust `--kit-preset` / `--adapter` for every sample (requires an explicit non-`auto` preset). |

The per-sample resolution table (sample, detected_kit, applied_kit, match_rate, dedup_strategy, source) is written to `<output>/kit_resolution.tsv` and embedded under `run_settings.json -> per_sample`. The `source` column distinguishes `detected`, `user_fallback`, `inferred_pretrimmed`, `explicit_off`, and `dry_run_*`.

#### Alignment

| Flag | Default | Description |
|---|---|---|
| `--mapq Q` | 10 | MAPQ threshold for the post-alignment filter (NUMT suppression). |
| `--seed N` | 42 | bowtie2 `--seed` value (deterministic output). |

#### Deduplication

| Flag | Default | Description |
|---|---|---|
| `--dedup-strategy {auto,umi-tools,skip,mark-duplicates}` | `auto` | Per-sample resolved. `auto` → `umi-tools` for UMI samples, `skip` otherwise. `mark-duplicates` is coordinate-only and destroys codon-occupancy signal on mt-Ribo-seq; gated behind the long confirmation flag. |
| `--umi-dedup-method {unique,percentile,cluster,adjacency,directional}` | `unique` | umi_tools `--method`. `unique` collapses only on exact coord+UMI match; other methods may over-collapse in low-complexity mt regions. |
| `--i-understand-mark-duplicates-destroys-mt-ribo-seq-signal` | — | Required confirmation to opt into `--dedup-strategy mark-duplicates`. |

---

### `mitoribopy rpf`

Runs the Ribo-seq analysis pipeline against a directory of BED (or BAM) files. Pipeline:

1. Load + filter BED on the RPF length window
2. Compute offset enrichment (combined + per sample)
3. Select offsets per read length (combined + per sample, depending on `--offset_mode`)
4. Render diagnostic plots (heatmaps, line plots, drift plot)
5. Translation profile: A-/P-/E-site footprint density, frame usage, codon usage (per sample, per requested site)
6. Coverage profile plots (per sample, per requested site)
7. Optional modules: structure-density export, codon correlation

Plain `mitoribopy <flags>` (no subcommand) still routes to `mitoribopy rpf` with a deprecation warning for v0.2.x compatibility.

#### Core inputs

| Flag | Default | Description |
|---|---|---|
| `-f, --fasta REF_FASTA` | — | Reference FASTA for the transcript/annotation context. **Required**. |
| `-s, --strain {h,y,vm,ym,custom}` | `y` | Reference preset; see [Strain presets](#strain-presets-and-footprint-classes). |
| `-d, --directory BED_DIR` | cwd | Directory of `.bed` and/or `.bam` files. BAMs are auto-converted to BED6 via pysam. |
| `--bam_mapq Q` | 10 | MAPQ threshold for BAM inputs before BAM→BED6 conversion. Set 0 to disable. |
| `-rpf MIN_LEN MAX_LEN` | strain/class default | Inclusive read-length filter range, e.g. `-rpf 29 34`. |
| `--footprint_class {monosome,disome,custom}` | `monosome` | Picks `-rpf` and `--unfiltered_read_length_range` defaults. |
| `--annotation_file ANNOTATION.csv` | — | Annotation CSV override. Required for `vm`, `ym`, `custom`. |
| `--codon_tables_file CODON_TABLES.json` | — | Codon-table JSON override. |
| `--codon_table_name TABLE_NAME` | strain default | Codon table to load (built-in or from `--codon_tables_file`). 27 built-in tables ship; run `--help` for the full list. |
| `--start_codons CODON [CODON ...]` | strain default | Allowed start codons. Defaults: `y` → `ATG`, `h` → `ATG ATA`, `custom` → `ATG`. |
| `--atp8_atp6_baseline {ATP6,ATP8}` | `ATP6` | Baseline name for the bicistronic ATP8/ATP6 transcript. Titles remain ATP8/ATP6. |
| `--nd4l_nd4_baseline {ND4,ND4L}` | `ND4` | Baseline name for the bicistronic ND4L/ND4 transcript. Titles remain ND4L/ND4. |

#### Offset enrichment + selection

| Flag | Default | Description |
|---|---|---|
| `-a, --align {start,stop}` | `start` | Anchor offset enrichment around the start or stop codon. |
| `-r, --range NT` | 20 | Plot offsets from -range to +range around the anchor codon. |
| `--min_offset NT` | 11 | Shared minimum absolute offset (used only when end-specific bounds are not provided). |
| `--max_offset NT` | 20 | Shared maximum absolute offset (used only when end-specific bounds are not provided). |
| `--min_5_offset NT`, `--max_5_offset NT` | from `--min_offset` / `--max_offset` | End-specific 5' selection bounds. **Preferred** over the shared bounds. |
| `--min_3_offset NT`, `--max_3_offset NT` | from `--min_offset` / `--max_offset` | End-specific 3' selection bounds. **Preferred**. |
| `--offset_mask_nt NT` | 5 | Mask near-anchor bins from -N..-1 and +1..+N in summaries and plots. |
| `--offset_pick_reference {selected_site,p_site}` | `p_site` | `p_site`: pick in canonical P-site space, then convert to the reported space. `selected_site`: pick directly in the final `--offset_site` space (legacy). |
| `--offset_type {5,3}` | `5` | Which read end the offset is measured from. |
| `--offset_site {p,a}` | `p` | Coordinate space for the SELECTED OFFSETS table (the values in `p_site_offsets_<align>.csv`). Does NOT control which downstream outputs are generated; use `--analysis_sites` for that. |
| `--analysis_sites {p,a,both}` | `both` | Which downstream sites to generate. `both` writes parallel P-site and A-site codon usage + coverage plots, side by side under per-site subdirs. `p` or `a` restricts to one site (legacy single-directory layout). |
| `--codon_overlap_mode {full,any}` | `full` | `full`: read must span all 3 nt of the anchor codon. `any`: any partial overlap counts. |
| `-p, --psite_offset NT` | — | Use one fixed offset for every read length and sample. Bypasses enrichment-based selection. |
| `--offset_mode {per_sample,combined}` | `per_sample` | `per_sample`: each sample uses its own offsets; drift surfaced in `offset_drift_<align>.svg`. `combined`: pool all samples and apply one table (v0.3.x behaviour). |

#### Outputs and plotting

| Flag | Default | Description |
|---|---|---|
| `-o, --output DIR` | `analysis_results` | Base output directory. |
| `--downstream_dir NAME` | `footprint_density` | Per-sample subdirectory name for frame and codon analyses. |
| `--plot_dir NAME` | `plots_and_csv` | Subdirectory name for offset CSVs and plots. |
| `-fmt, --plot_format {png,pdf,svg}` | `png` | File format for saved plots. SVG recommended for publication-ready vector output. |
| `--x_breaks NT [NT ...]` | — | Optional custom x-axis tick marks for offset line plots. |
| `--line_plot_style {combined,separate}` | `combined` | Draw 5'/3' offsets in one panel or two. |
| `--cap_percentile FRAC` | 0.999 | Upper percentile cap for coverage-style plots. |
| `-m, --merge_density` | off | Collapse frame 1 and frame 2 into frame 0 for codon-density summaries. |
| `--order_samples NAME [NAME ...]` | — | Optional explicit sample order for plots and aggregated outputs. |

#### Read-count normalization

| Flag | Default | Description |
|---|---|---|
| `--read_counts_file PATH` | `read_counts_summary.txt` | Read-count table for RPM normalization. Auto-wired from `align/read_counts.tsv` when running through `mitoribopy all`. |
| `--read_counts_sample_col NAME` | auto | Sample column override; defaults to case-insensitive name match, then column 1. |
| `--read_counts_reads_col NAME` | auto | Read-count column override; defaults to `reads` / `read_count` / `counts`, then column 3. |
| `--read_counts_reference_col NAME` | auto | Reference column override; needed for `--rpm_norm_mode mt_mrna` if auto-detection fails. |
| `--unfiltered_read_length_range MIN MAX` | `[15, 50]` | Range for the unfiltered QC summary and heatmaps. Broaden for disome studies. |
| `--rpm_norm_mode {total,mt_mrna}` | `total` | RPM denominator: sum all rows or only mt-mRNA rows. |
| `--mrna_ref_patterns PATTERN [PATTERN ...]` | `mt_genome mt-mrna mt_mrna` | Substring patterns identifying mt-mRNA rows for `--rpm_norm_mode mt_mrna`. |

#### Optional modules

| Flag | Default | Description |
|---|---|---|
| `--structure_density` | off | Export log2 and scaled density values from footprint-density tables. |
| `--structure_density_norm_perc FRAC` | 0.99 | Upper percentile used to cap and scale structure-density values. |
| `--cor_plot` | off | Generate codon-correlation plots. |
| `--base_sample NAME` | — | Reference sample for codon-correlation comparisons. |
| `--cor_mask_method {percentile,fixed,none}` | `percentile` | Masking rule for extreme codon-correlation outliers. |
| `--cor_mask_percentile FRAC` | 0.99 | Used when `--cor_mask_method percentile`. |
| `--cor_mask_threshold FLOAT` | — | Used when `--cor_mask_method fixed`. |

`--use_rna_seq` is **deprecated** in v0.3.0 and **removed** in v0.4.0+ — use the dedicated `mitoribopy rnaseq` subcommand instead.

---

### `mitoribopy rnaseq`

Integrates a pre-computed differential-expression table (DESeq2 / Xtail / Anota2Seq) with a prior `mitoribopy rpf` run and emits TE and ΔTE tables plus diagnostic plots. Enforces a SHA256 reference-consistency gate: Ribo-seq and RNA-seq must hash to the identical transcript set.

#### DE table

| Flag | Default | Description |
|---|---|---|
| `--de-table PATH` | — | DE results table (CSV or TSV). **Required**. |
| `--de-format {auto,deseq2,xtail,anota2seq,custom}` | `auto` | DE table format; `auto` detects from column headers. |
| `--de-gene-col NAME` | — | Column name for gene IDs. Used when `--de-format custom`. |
| `--de-log2fc-col NAME` | — | Column name for log2 fold change. |
| `--de-padj-col NAME` | — | Column name for adjusted p-value. |
| `--de-basemean-col NAME` | — | Column name for basemean. |

#### Gene identifiers

| Flag | Default | Description |
|---|---|---|
| `--gene-id-convention {ensembl,refseq,hgnc,mt_prefixed,bare}` | — | **Required, no default**. Mismatched conventions silently produce zero-match runs. Examples: `hgnc` → `MT-ND1`; `bare` → `ND1`; `ensembl` → `ENSG00000198888`; `refseq` → `YP_003024026.1`; `mt_prefixed` → `MT-ND1`. |
| `--organism {h,y}` | `h` | Organism for the mt-mRNA registry. |

#### Ribo-seq inputs

| Flag | Default | Description |
|---|---|---|
| `--ribo-dir DIR` | — | Output directory of a prior `mitoribopy rpf` run. Must contain `rpf_counts.tsv` and `run_settings.json` (or `run_manifest.json`) with a recorded `reference_checksum`. |
| `--ribo-counts PATH` | `<ribo-dir>/rpf_counts.tsv` | Explicit path to the per-sample per-gene RPF counts table. |

#### Reference-consistency gate (exactly one)

| Flag | Default | Description |
|---|---|---|
| `--reference-gtf PATH` | — | Reference used by RNA-seq; `mitoribopy rnaseq` hashes this and verifies it matches the hash in the rpf run's manifest. |
| `--reference-checksum SHA256` | — | Precomputed SHA-256, useful when the reference file is not on the rnaseq host. |

#### Conditions (optional, required for replicate-based ΔTE)

| Flag | Default | Description |
|---|---|---|
| `--condition-map PATH` | — | TSV with columns `sample` and `condition` assigning Ribo-seq samples to conditions. |
| `--condition-a NAME` | — | Reference condition. |
| `--condition-b NAME` | — | Comparison condition. |

Without a condition map, `mitoribopy rnaseq` computes a point-estimate ΔTE using only the DE table's mRNA log2FC and emits rows with a `single_replicate_no_statistics` note.

#### Output

| Flag | Default | Description |
|---|---|---|
| `--output DIR` | — | Output directory for `te.tsv`, `delta_te.tsv`, and plots. |

---

### `mitoribopy all`

End-to-end orchestrator: align + rpf + (optional) rnaseq. Reads one YAML/JSON/TOML config and dispatches to the per-stage subcommands; writes a composed `run_manifest.json` covering every parameter, tool version, and reference checksum.

| Flag | Default | Description |
|---|---|---|
| `--config PATH` | — | Configuration file with `align:` / `rpf:` / `rnaseq:` sections. **Required** unless using `--print-config-template` or `--show-stage-help`. |
| `--output DIR` | — | Run root. Each stage writes under `<output>/align/`, `<output>/rpf/`, `<output>/rnaseq/`. |
| `--resume` | off | Skip stages whose sentinel output already exists: `align/read_counts.tsv`, `rpf/rpf_counts.tsv`, `rnaseq/delta_te.tsv`. |
| `--skip-align` | off | Skip the align stage even when an `align:` section exists. |
| `--skip-rpf` | off | Skip the rpf stage. |
| `--skip-rnaseq` | off | Skip the rnaseq stage even when an `rnaseq:` section exists. |
| `--manifest PATH` | `run_manifest.json` | Manifest filename (relative to `--output`). |
| `--show-stage-help STAGE` | — | Print the full help for one stage (`align`, `rpf`, or `rnaseq`) and exit. Useful because `--help` shows only orchestrator flags. |
| `--print-config-template` | — | Print a commented YAML template covering every stage and exit. |

#### Auto-wiring

When `align` and `rpf` both run, `mitoribopy all` auto-wires:

- `rpf.directory` → `<run_root>/align/bed`
- `rpf.read_counts_file` → `<run_root>/align/read_counts.tsv`

When `rpf` and `rnaseq` both run:

- `rnaseq.ribo_dir` → `<run_root>/rpf`

Each stage's own `--output` defaults to `<run_root>/<stage>/`.

#### YAML config shape

Every key under a section maps to that subcommand's CLI flag with hyphens converted to underscores:

- `align:` keys use the **hyphen** flag style (`--kit-preset` → `kit_preset`).
- `rpf:` keys use the **underscore** flag style (`--offset_type` → `offset_type`) because the rpf parser uses underscored flag names directly.
- `rnaseq:` keys use the hyphen style (`--gene-id-convention` → `gene_id_convention`).

The orchestrator handles both styles transparently. Booleans emit the bare flag (`true`) or are omitted entirely (`false`); `null` values are dropped.

---

## Output overview

For an `mitoribopy all` run with the v0.4.1 defaults (`--offset_mode per_sample`, `--analysis_sites both`):

```text
<output>/
  run_manifest.json                    # composed provenance: align + rpf + rnaseq settings,
                                       # tool versions, reference_checksum, stages_run
  align/
    mitoribopy.log
    read_counts.tsv                    # per-stage counts (input → trimmed → contam → mt → MAPQ → dedup)
    kit_resolution.tsv                 # per-sample kit + dedup decisions; the spine of v0.4 provenance
    run_settings.json                  # includes per_sample[] block with detected_kit, applied_kit,
                                       # match_rate, dedup_strategy, source per sample
    trimmed/                           # *.trimmed.fq.gz, *.cutadapt.json
    contam_filtered/                   # *.nocontam.fq.gz
    aligned/                           # *.bam, *.mapq.bam
    deduped/                           # *.dedup.bam (or hardlink of *.mapq.bam when dedup=skip)
    bed/                               # *.bed — strand-aware BED6 inputs to rpf
  rpf/
    mitoribopy.log
    rpf_counts.tsv                     # per-sample per-gene RPF counts; feeds rnaseq
    run_settings.json                  # includes reference_checksum (SHA256 of --fasta)
    plots_and_csv/
      offset_<align>.csv               # COMBINED enrichment summary (diagnostic)
      p_site_offsets_<align>.csv       # COMBINED selected offsets (diagnostic)
      offset_drift_<align>.svg         # per-sample 5'/3' offset comparison; read this FIRST
      offset_*.svg                     # heatmaps + line plots
      per_sample/
        <sample>/
          offset_<align>.csv           # per-sample enrichment summary
          p_site_offsets_<align>.csv   # per-sample selected offsets (used downstream)
    translation_profile_p/             # P-site outputs (when analysis_sites=both)
      <sample>/
        footprint_density/             # *_footprint_density.csv (P/A/E columns) + *_depth plots
        translating_frame/             # frame_usage_total.csv, frame_usage_by_transcript.csv
        codon_usage/                   # codon_usage_<transcript>.csv, codon_usage_total.csv,
                                       # a_site_codon_usage_<transcript>.csv, ...
        debug_csv/                     # raw position-by-transcript CSV dumps for debugging
    translation_profile_a/             # A-site outputs (when analysis_sites=both); same shape
      <sample>/...
    coverage_profile_plots_p/          # per-site coverage plots
      read_coverage_rpm/, p_site_coverage_rpm/, *_codon/, *_frame/
    coverage_profile_plots_a/
      ...
    structure_density/                 # if --structure_density (uses P-site by default)
    codon_correlation/                 # if --cor_plot
  rnaseq/                              # if rnaseq config supplied
    te.tsv                             # one row per (sample, gene): rpf_count, mrna_abundance, te
    delta_te.tsv                       # one row per gene: mrna_log2fc, rpf_log2fc, delta_te_log2,
                                       # padj, note
    plots/
      mrna_vs_rpf.png                  # four-quadrant log2FC scatter
      delta_te_volcano.png             # ΔTE (log2) vs -log10(padj)
    run_settings.json
```

When `--analysis_sites=p` or `=a` (single site), the per-site directories collapse back to the legacy v0.3 layout:
- `<output>/rpf/<sample>/footprint_density/` (instead of `translation_profile_p/<sample>/...`)
- `<output>/rpf/coverage_profile_plots/` (instead of `coverage_profile_plots_p/`)

### Key files to look at first

- **`align/kit_resolution.tsv`** — confirm each sample resolved to the right kit and dedup strategy. The `source` column tells you whether it was `detected`, `user_fallback`, `inferred_pretrimmed`, or set explicitly.
- **`align/read_counts.tsv`** — track per-stage drop-off. The invariants `rrna_aligned + post_rrna_filter == post_trim` and `mt_aligned + unaligned_to_mt == post_rrna_filter` must hold; if they don't, something is wrong with alignment.
- **`rpf/plots_and_csv/offset_drift_<align>.svg`** — per-sample offset comparison. Outliers are visible by eye in seconds.
- **`rpf/plots_and_csv/offset_<align>.svg`** (heatmap + line plots) — the combined enrichment around the anchor codon. A sharp peak at the canonical P-site offset (~12–15 nt from the 5' end for most mt-Ribo-seq libraries) confirms library quality.
- **`rpf/translation_profile_p/<sample>/codon_usage/codon_usage_total.csv`** — overall P-site codon occupancy. Compare against `translation_profile_a/.../codon_usage_total.csv` to see the A-site picture.

---

## Custom organisms

To run on an organism that isn't human or yeast, supply your own annotation, codon table, and start codons:

```bash
mitoribopy rpf \
  -s custom \
  -f your_reference.fa \
  --directory bed/ \
  -rpf 28 34 \
  --annotation_file examples/custom_reference/annotation_template.csv \
  --codon_tables_file examples/custom_reference/codon_tables_template.json \
  --codon_table_name custom_example \
  --start_codons ATG GTG \
  --output results/
```

Two lighter-weight strain options exist when you only need a different reference and want to keep the codon table:

- `-s vm` (vertebrate-mito codon table): supply `--annotation_file` + `-rpf`; the `vertebrate_mitochondrial` codon table is used.
- `-s ym` (yeast-mito codon table): supply `--annotation_file` + `-rpf`; the `yeast_mitochondrial` codon table is used.

Templates and a worked example:

- [examples/custom_reference/annotation_template.csv](examples/custom_reference/annotation_template.csv)
- [examples/custom_reference/codon_tables_template.json](examples/custom_reference/codon_tables_template.json)
- [examples/custom_reference/README.md](examples/custom_reference/README.md)

---

## Built-in references

MitoRiboPy ships packaged reference data for:

- Human mt-translation using the `vertebrate_mitochondrial` codon table
- Yeast mt-translation using the `yeast_mitochondrial` codon table

Built-in annotation tables are stored as CSV and built-in codon tables as JSON under [src/mitoribopy/data](src/mitoribopy/data). 27 built-in codon tables are available; pass `--codon_table_name` to select one (run `mitoribopy rpf --help` for the full list).

For bicistronic transcript regions:

- Plot titles stay consistent as `ATP8/ATP6` and `ND4L/ND4`.
- The default sequence baselines are `ATP6` and `ND4`.
- Switch them with `--atp8_atp6_baseline ATP8|ATP6` and `--nd4l_nd4_baseline ND4L|ND4`.

Legacy FASTA/BED identifiers such as `ATP86` and `ND4L4` are still recognized through built-in aliases.

---

## Examples

### A. Single-command end-to-end (recommended)

```bash
mitoribopy all --config pipeline_config.yaml --output results/ --threads 8
```

### B. Just the align stage on a directory of FASTQs

```bash
mitoribopy align \
  --kit-preset auto \
  --library-strandedness forward \
  --fastq-dir input_data/ \
  --contam-index input_data/indexes/rrna_contam \
  --mt-index input_data/indexes/mt_tx \
  --output results/align/ \
  --threads 8
```

### C. Just the rpf stage on existing BEDs (human, monosome)

```bash
mitoribopy rpf \
  -s h \
  -f references/human-mt-mRNA.fasta \
  --directory results/align/bed/ \
  -rpf 29 34 \
  --align stop \
  --offset_type 5 \
  --offset_site p \
  --offset_pick_reference p_site \
  --offset_mode per_sample \
  --analysis_sites both \
  --min_5_offset 10 \
  --max_5_offset 22 \
  --min_3_offset 10 \
  --max_3_offset 22 \
  --offset_mask_nt 5 \
  --plot_format svg \
  --read_counts_file results/align/read_counts.tsv \
  --rpm_norm_mode mt_mrna \
  --output results/rpf/ \
  -m
```

### D. Disome (collided ribosome) study

```bash
mitoribopy rpf \
  -s h \
  -f references/human-mt-mRNA.fasta \
  --directory bed/ \
  --footprint_class disome \
  --output results_disome/
```

`--footprint_class disome` widens `-rpf` to 60–90 nt (vertebrate) or 65–95 nt (yeast) and `--unfiltered_read_length_range` to 40–110 nt automatically. You can still override either with `-rpf` / `--unfiltered_read_length_range`.

### E. Custom organism

```bash
mitoribopy rpf \
  -s custom \
  -f my_reference.fa \
  --directory bed/ \
  -rpf 28 34 \
  --annotation_file my_annotation.csv \
  --codon_tables_file my_codon_tables.json \
  --codon_table_name my_table \
  --start_codons ATG GTG \
  --output results/
```

### F. RNA-seq integration for translation efficiency

```bash
mitoribopy rnaseq \
  --de-table de.tsv \
  --gene-id-convention hgnc \
  --ribo-dir results/rpf \
  --reference-gtf references/human-mt-mRNA.fasta \
  --condition-map samples.tsv \
  --condition-a control \
  --condition-b knockdown \
  --output results/rnaseq/
```

### G. Resume a partial run

```bash
mitoribopy all --config pipeline_config.yaml --output results/ --resume
```

`--resume` skips a stage when its sentinel file already exists (`align/read_counts.tsv`, `rpf/rpf_counts.tsv`, `rnaseq/delta_te.tsv`).

### H. Pre-trimmed inputs (e.g. SRA-deposited)

No special configuration needed — when adapter detection finds 0% across every kit, the resolver falls through to the `pretrimmed` kit automatically and cutadapt skips the `-a` flag. To opt in explicitly:

```yaml
align:
  kit_preset: pretrimmed
```

To force the v0.4.0 hard-fail behaviour back: pass `--no-pretrimmed-inference` (or set `allow_pretrimmed_inference: false` in YAML).

---

## Logs and provenance

- **`<output>/<stage>/mitoribopy.log`** — every stage writes a persistent log file alongside the same lines printed to the terminal.
- **Per-stage `run_settings.json`** — every stage writes its own settings JSON (resolved kit, dedup strategy, MAPQ threshold, reference checksum, tool versions, …).
- **`<output>/run_manifest.json`** — `mitoribopy all` composes per-stage settings into a top-level manifest. The `reference_checksum` from the rpf stage is promoted to the top level so a downstream `rnaseq` run can verify it without drilling into the rpf section.
- **`<output>/align/kit_resolution.tsv`** — per-sample kit + dedup decisions.

---

## Development

Run the test suite with:

```bash
PYTHONPATH=src pytest
```

Helpful tools:

```bash
PYTHONPATH=src pytest -k offset           # run a subset by keyword
PYTHONPATH=src pytest -x --tb=short       # stop at first failure, terse traceback
PYTHONPATH=src pytest tests/test_align_sample_resolve.py -v
```

Documentation lives under [docs/](docs/):

- [docs/README.md](docs/README.md) — index
- [docs/reference/cli.md](docs/reference/cli.md) — concise CLI reference
- [docs/tutorials/](docs/tutorials/) — step-by-step worked examples
- [docs/release-notes/](docs/release-notes/) — version-by-version notes
- [docs/validation/](docs/validation/) — biological validation plans
- [docs/environment/](docs/environment/) — bioconda env file + Dockerfile
- [docs/diagrams/](docs/diagrams/) — Mermaid pipeline diagrams

---

## License

MIT. See [LICENSE](LICENSE).
