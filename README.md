<p align="center">
  <img src="docs/banner.png" alt="MitoRiboPy — A Python package for mt-Ribo-seq" width="100%">
</p>

# MitoRiboPy

Mitochondrial ribosome profiling (mt-Ribo-seq) analysis, end to end.

MitoRiboPy is a Python package + CLI for analysing mt-Ribo-seq data from raw FASTQ all the way through translation-efficiency integration with paired RNA-seq. Every per-sample decision (kit, dedup, offsets) is independent, so mixed-library batches just work.

The package is built around four pipeline subcommands plus seven utility subcommands:

| Subcommand | What it does |
|---|---|
| `mitoribopy align` | FASTQ → BAM → BED6 + per-sample read counts (cutadapt + bowtie2 + umi_tools + pysam) |
| `mitoribopy rpf` | BED/BAM → offsets, translation profile, codon usage, coverage plots |
| `mitoribopy rnaseq` | Translation efficiency from paired RNA-seq + Ribo-seq. **Default flow:** pass `--rna-fastq` + `--ribo-fastq` + `--reference-fasta` and the subcommand runs trimming → bowtie2 → counting → pyDESeq2 → TE / ΔTE / plots end-to-end (requires the `[fastq]` extra: `pip install 'mitoribopy[fastq]'`). **Alternative:** pass `--de-table` from an external DESeq2 / Xtail / Anota2Seq run + a prior rpf run; this path is mutually exclusive with `--rna-fastq` and enforces a SHA256 reference-consistency gate. |
| `mitoribopy all` | End-to-end orchestrator that runs align + rpf + (optional) rnaseq from one YAML config and writes a composed `run_manifest.json` |
| `mitoribopy validate-config` | Pre-flight check: parse + canonicalise legacy keys + check paths + validate `rnaseq.mode` against supplied inputs. Exit 0 / 2. Use before launching long-running cluster jobs. |
| `mitoribopy validate-reference` | Pre-flight a custom mt-transcriptome FASTA + annotation pair (matching IDs, matching lengths, CDS divisible by 3, valid start / stop codons under the chosen codon table). Exit 0 / 2. |
| `mitoribopy validate-figures` | Mechanically validate every plot under a finished run root: check label / legend / stat-box overlap, label clipping, point counts vs source TSV, SVG text editability, PNG dpi, and metadata sidecars. Writes `figure_qc.tsv`. Exit 0 / 1 / 2 (clean / warn-only / fail; `--strict` upgrades warn → fail). |
| `mitoribopy migrate-config` | Rewrite legacy YAML keys to canonical names (e.g. `merge_density:` → `codon_density_window:`, `strain: h` → `strain: h.sapiens`). Pipe stdout to a new file or pass `--in-place`. |
| `mitoribopy summarize` | Regenerate `SUMMARY.md` and `summary_qc.tsv` from a finished run by reading `run_manifest.json` and per-stage outputs. Auto-invoked by `mitoribopy all` after every run. |
| `mitoribopy benchmark` | Time + RSS + disk-measure a `mitoribopy all` invocation. `--subsample N` reservoir-samples each FASTQ for fast tuning runs; outputs `benchmark.tsv` + `benchmark_summary.md`. |

### Which command should I use?

| Situation | Command |
|---|---|
| I have raw mt-Ribo-seq FASTQ and want everything | `mitoribopy all --config pipeline_config.yaml --output results/` |
| I only want adapter trimming + alignment | `mitoribopy align --config align_config.yaml ...` |
| I already have BED / BAM files | `mitoribopy rpf --config rpf_config.yaml ...` |
| I have full-transcriptome RNA-seq DE results + RPF counts (publication route) | `mitoribopy rnaseq --rnaseq-mode de_table --de-table de.tsv --ribo-dir runs/full/rpf` |
| I want a quick exploratory mt-only RNA/Ribo TE run | `mitoribopy rnaseq --rnaseq-mode from_fastq --rna-fastq rna/ --ribo-fastq ribo/` |
| I have a non-human / non-yeast organism | `mitoribopy rpf --strain custom --annotation_file ... --codon_tables_file ...` |
| I want to check my config without running anything | `mitoribopy validate-config pipeline_config.yaml` |
| I want to pre-flight a custom mt-transcriptome FASTA + annotation pair | `mitoribopy validate-reference --fasta custom_mt.fa --annotation custom_mt.csv` |
| I want to mechanically QC every plot in a finished run | `mitoribopy validate-figures runs/full/ --strict` |
| I want to inspect what `all` would actually do | `mitoribopy all --print-canonical-config --config ... --output ...` |
| I want to estimate cluster time / memory / disk | `mitoribopy benchmark --config ... --output bench/ --subsample 200000` |
| My config uses old key names (e.g. `merge_density:`, `strain: h`) | `mitoribopy migrate-config old.yaml > new.yaml` |
| I want to (re-)render the summary of an old run | `mitoribopy summarize runs/full/` |
| I want to safely resume a partially completed run | `mitoribopy all --config ... --output runs/full/ --resume` |

The `mitoribopy all --resume` skip is **hash-validated** against the prior run's `run_manifest.json`: edits to the config / sample sheet / reference FASTA cause the affected stage(s) to re-run rather than silently re-using stale outputs. Pass `--force-resume` (or set `MITORIBOPY_FORCE_RESUME=1`) to bypass the guard for known-safe edits.

---

## Table of contents

1. [What MitoRiboPy is for](#what-mitoribopy-is-for)
2. [Pipeline overview](#pipeline-overview)
3. [Installation](#installation)
4. [Quick start](#quick-start)
5. [Inputs you need to prepare](#inputs-you-need-to-prepare)
   - [What each module requires](#what-each-module-requires)
   - [Sample sheet (unified per-project TSV)](#sample-sheet-unified-per-project-tsv)
   - [Input files (file-by-file reference)](#input-files-file-by-file-reference)
6. [How to run — YAML vs shell wrapper](#how-to-run--yaml-vs-shell-wrapper)
7. [Strain presets and footprint classes](#strain-presets-and-footprint-classes)
8. [Subcommand reference](#subcommand-reference)
   - [`mitoribopy align`](#mitoribopy-align)
   - [`mitoribopy rpf`](#mitoribopy-rpf)
   - [`mitoribopy rnaseq`](#mitoribopy-rnaseq)
   - [`mitoribopy all`](#mitoribopy-all)
9. [What the numbers mean — RNA, RPF, TE, ΔTE](#what-the-numbers-mean--rna-rpf-te-%CE%B4te)
10. [Output overview](#output-overview)
11. [Custom organisms](#custom-organisms)
12. [Built-in references](#built-in-references)
13. [Examples](#examples)
14. [Tools](#tools)
15. [Logs and provenance](#logs-and-provenance)
16. [Development](#development)
17. [License](#license)

---

## What MitoRiboPy is for

MitoRiboPy is a focused tool for the 13 mt-mRNAs of human mitochondria (or 8 mt-mRNAs in yeast, with configurable codon tables for any other mitochondrion). It ships:

- **Per-sample adapter detection** with auto-fallback to the right kit. Mixed-kit and mixed-UMI batches resolve each sample independently. Pre-trimmed FASTQs (e.g. SRA-deposited data) are auto-detected and routed through cutadapt with no `-a` flag.
- **Per-sample offset selection** so inter-sample drift in the canonical 12–15 nt 5' P-site offset doesn't bias your downstream codon-usage tables. A combined-across-samples diagnostic is still emitted and an `offset_drift_<align>.svg` plot makes drift visible at a glance.
- **Both P-site and A-site downstream outputs** by default, side by side under per-site subdirectories. No ambiguity about which output corresponds to which site.
- **End-to-end RNA-seq + Ribo-seq → TE / ΔTE in one subcommand.** `mitoribopy rnaseq` takes raw FASTQs and a transcriptome FASTA and runs trimming, bowtie2 alignment, per-transcript counting, and pyDESeq2 itself before emitting TE, ΔTE, and a six-figure plot set. (Bringing your own pre-computed DE table from R / Python remains supported via `--de-table` and enforces a SHA256 reference-consistency gate.)
- **Strain-aware defaults**: built-in human (`-s h.sapiens`) and yeast (`-s s.cerevisiae`) annotations + codon tables, plus `custom` for any other organism with a published NCBI Genetic Code (mouse, fly, plants, fungi, ...).

What MitoRiboPy is **not**:

- Not a general-purpose nuclear Ribo-seq pipeline. The defaults, references, and dedup heuristics are calibrated for the low-complexity 13-mRNA mt universe.
- Not a general DE engine for nuclear genes. The default `mitoribopy rnaseq` flow runs **pyDESeq2 on the mt-mRNA subset only**, which is fine for exploring mt-translation efficiency end-to-end on a single library. For publication-grade DE statistics across the full transcriptome run DESeq2 / Xtail / Anota2Seq externally and pass the resulting table via `--de-table`.

---

## Pipeline overview

![Pipeline overview](docs/diagrams/01_pipeline_overview.png)

UMI handling: the UMI is extracted into the read QNAME during the cutadapt trim step (5' single-pass or 3' two-pass), so it travels through bowtie2 alignment unchanged and is available for `umi_tools dedup` after the MAPQ filter — the only stage that needs alignment coordinates AND the UMI together.

### Detailed stage diagrams

- [docs/diagrams/02_align_stage.png](docs/diagrams/02_align_stage.png) — internals of `mitoribopy align`: per-sample resolution → cutadapt + UMI → contam subtract → bowtie2 → MAPQ → dedup → BED6.
- [docs/diagrams/03_rpf_stage.png](docs/diagrams/03_rpf_stage.png) — internals of `mitoribopy rpf`: filter BED → offset enrichment + selection → translation_profile + coverage_profile_plots + optional modules.
- [docs/diagrams/04_rnaseq_stage.png](docs/diagrams/04_rnaseq_stage.png) — internals of the optional `mitoribopy rnaseq` stage: DE table + rpf_counts → SHA256 reference gate → TE → ΔTE → scatter + volcano.

Regenerate with `python docs/diagrams/render_diagrams.py` (matplotlib only; no Node / mermaid-cli required).

---

## Installation

### From PyPI (recommended)

```bash
python -m pip install mitoribopy
```

The package is published on PyPI: [pypi.org/project/mitoribopy](https://pypi.org/project/mitoribopy/). This pulls in every Python dependency (`numpy`, `pandas`, `matplotlib`, `seaborn`, `biopython`, `scipy`, `PyYAML`, `pysam`) automatically. The external bioinformatics tools (`cutadapt`, `bowtie2`, `umi_tools`, …) still need to be on `$PATH` separately — see [External tool dependencies](#external-tool-dependencies) below.

### From source (latest development version)

```bash
git clone https://github.com/Ahram-Ahn/MitoRiboPy.git
cd MitoRiboPy
python -m pip install -e .
```

Use this when you need a fix or feature that has not yet been released to PyPI. For development and tests, install the dev extras:

```bash
python -m pip install -e ".[dev]"
```

### Verify the install

```bash
mitoribopy --version
mitoribopy --help
```

If you prefer not to install at all:

```bash
PYTHONPATH=src python -m mitoribopy --help
```

### External tool dependencies

MitoRiboPy shells out to a small set of standard bioinformatics tools. All of them must be on `$PATH` for a real run:

| Tool | Used by | Required when |
|---|---|---|
| `cutadapt` | `align`, `rnaseq` (from-FASTQ) | always (length + quality filter even for pre-trimmed data) |
| `bowtie2` + `bowtie2-build` | `align`, `rnaseq` (from-FASTQ) | always |
| `umi_tools` | `align` | at least one sample's resolved kit has UMIs |
| `pysam` (Python lib) | `align`, `rpf`, `rnaseq` (from-FASTQ) | always (installed automatically via `pip`) |
| `samtools` | optional | recommended for inspecting outputs; not required |
| `pydeseq2` (Python lib) | `rnaseq` (from-FASTQ) | only when running `mitoribopy rnaseq` in from-FASTQ mode (`--rna-fastq …`); install via the `[fastq]` extra: `pip install 'mitoribopy[fastq]'`. The pre-computed-DE flow does not need it. |

The bioconda environment under [docs/environment/environment.yml](docs/environment/environment.yml) installs everything in one command:

```bash
conda env create -f docs/environment/environment.yml
conda activate mitoribopy
```

---

## Quick start

The shortest path from raw FASTQ to translation-profile + coverage outputs is one YAML file plus one command.

```bash
# 1. Start from one of two templates next to your data and fill in the paths:
#    -- examples/templates/ ships an EXHAUSTIVE template that lists every
#       available flag with its default value and a 1-line comment:
cp examples/templates/pipeline_config.example.yaml pipeline_config.yaml
#    -- OR get the curated MINIMAL template from the CLI:
mitoribopy all --print-config-template > pipeline_config.yaml

$EDITOR pipeline_config.yaml

# 2. Optional: dry-run prints the per-stage argv so you can review.
mitoribopy all --config pipeline_config.yaml --output results/ --dry-run

# 3. Run.
mitoribopy all --config pipeline_config.yaml --output results/
```

The matching shell-script templates are at [examples/templates/run_align.example.sh](examples/templates/run_align.example.sh), [examples/templates/run_rpf.example.sh](examples/templates/run_rpf.example.sh), and [examples/templates/run_pipeline.example.sh](examples/templates/run_pipeline.example.sh) — pick those if you prefer per-stage commands you can split across cluster jobs.

A minimal `pipeline_config.yaml` for a typical human mt-Ribo-seq run looks like this (annotated):

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
  strain: h.sapiens               # human mt-mRNA reference + codon table
  fasta: input_data/human-mt-mRNA.fasta
  footprint_class: monosome       # short | monosome | disome | custom
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
  codon_density_window: true
```

After the run, you'll have:

```text
results/
  align/    bed/, kit_resolution.tsv, read_counts.tsv, run_settings.json
            (intermediate trimmed/contam_filtered/aligned files are deleted as
             soon as they are consumed; pass --keep-intermediates to retain
             them; deduped/ is only created for UMI samples)
  rpf/      offset_diagnostics/{csv,plots}/, translation_profile/<sample>/...,
            coverage_profile_plots/{read_coverage_*, {p_site,a_site}_density_*},
            codon_correlation/{p_site,a_site}/, igv_tracks/<sample>/, rpf_counts.tsv
  run_manifest.json
```

See [Output overview](#output-overview) for the full directory tree.

---

## Inputs you need to prepare

This section answers "what do I need on disk before I run anything?". Read [What each module requires](#what-each-module-requires) first to scope your prep work, then see [Sample sheet](#sample-sheet-unified-per-project-tsv) for the recommended single-source-of-truth file format, and [Input files](#input-files-file-by-file-reference) for per-file specifics.

### What each module requires

Each subcommand in the pipeline consumes its own slice of the inputs below. **Required** rows must be present or the run aborts. **Optional** rows have sensible defaults but commonly need overriding.

#### `mitoribopy align` — RNase-trimmed FASTQs → BED + read counts

| Kind | Item | Required? | Notes |
|---|---|---|---|
| File | Ribo-seq FASTQs (`*.fq[.gz]` / `*.fastq[.gz]`) | **required** | Provide via `--fastq <files>` (repeatable) or `--fastq-dir <dir>`. |
| File | mt-transcriptome bowtie2 index prefix | **required** | Built once with `bowtie2-build`. Pass via `--mt-index <prefix>`. |
| File | rRNA / tRNA contaminant bowtie2 index prefix | **required** | Used to subtract contaminants. Pass via `--contam-index <prefix>`. |
| File | **Sample sheet** (`samples.tsv`) | optional | One row per Ribo-seq FASTQ; documents per-sample kit / UMI / strandedness. Strongly recommended for mixed-kit batches. See [Sample sheet](#sample-sheet-unified-per-project-tsv). |
| Option | `--kit-preset` | optional (`auto` default) | Library-prep kit; `auto` detects per FASTQ. Override per sample in the sample sheet. |
| Option | `--library-strandedness {forward,reverse,unstranded}` | optional (`forward` default) | dUTP-stranded libraries should set `reverse`. |

#### `mitoribopy rpf` — BED + reference FASTA → P-site / A-site analysis

| Kind | Item | Required? | Notes |
|---|---|---|---|
| File | Ribo-seq BEDs (or BAMs) | **required** | When run after `align`, auto-wired from `<align>/bed/`. Standalone runs use `--directory <dir>`. |
| File | mt-transcriptome FASTA | **required** | One record per mt-mRNA. Pass via `--fasta <path>`. |
| File | Annotation CSV | optional (built-in for h.sapiens / s.cerevisiae) | Required only for custom organisms. See [Custom organisms](#custom-organisms) for the full schema. |
| File | Codon-table JSON | optional (NCBI codes 1–33 bundled) | Pick one via `--codon_table_name`; supply your own only for non-NCBI codes. |
| File | Read-count table (`read_counts.tsv`) | optional (RPM normalization) | Auto-wired from `align/read_counts.tsv` when run via `mitoribopy all`. |
| Option | `-rpf <min> <max>` | optional (footprint-class default) | Read-length window. Defaults: monosome 28–34 (human), 37–41 (yeast). |
| Option | `--strain {h.sapiens, s.cerevisiae, custom}` | **required** | Drives the built-in annotation + codon table; `custom` requires the two files above. |

#### `mitoribopy rnaseq` — TE / ΔTE from RNA-seq + Ribo-seq (two flows)

**Default flow (from raw FASTQ):**

| Kind | Item | Required? | Notes |
|---|---|---|---|
| File | RNA-seq FASTQs | **required** | Via `--rna-fastq <files/dirs>` OR derived from a sample sheet. |
| File | Ribo-seq FASTQs | optional | Via `--ribo-fastq` OR sample sheet. When omitted the run short-circuits after writing the RNA DE table. |
| File | Transcriptome FASTA | **required** | Via `--reference-fasta <path>`. SHA256 recorded in the manifest. |
| File | Sample sheet OR condition map | **required (one of)** | `--sample-sheet samples.tsv` (recommended) **OR** `--condition-map conditions.tsv` (legacy two-column form). Pairs RNA and Ribo by `sample_id`, never by index. |
| Option | `--gene-id-convention {ensembl,refseq,hgnc,mt_prefixed,bare}` | **required** | Identifier scheme used in the FASTA / DE table. No default. |
| Option | `--base-sample` / `--compare-sample` | **required** | Reference vs comparison condition for the contrast (also accepted as `--condition-a` / `--condition-b`). |
| Option | `--allow-pseudo-replicates-for-demo-not-publication` | opt-in | Required to proceed when any condition has only 1 sample (publication-safe default is to fail-fast). |

**Alternative flow (bring your own DE table):**

| Kind | Item | Required? | Notes |
|---|---|---|---|
| File | DE results table (DESeq2 / Xtail / Anota2Seq) | **required** | Via `--de-table <path>`. CSV or TSV. |
| File | Prior `mitoribopy rpf` output dir | **required** | Via `--ribo-dir <dir>`. Must contain `rpf_counts.tsv` + `run_settings.json` with `reference_checksum`. |
| File | Reference FASTA / GTF | **required (one of)** | Via `--reference-gtf <path>` (gets hashed) **or** `--reference-checksum <sha256>` (when the file is not on this host). |
| File | Condition map | optional | Enables a replicate-based Ribo log2FC for ΔTE; without it, ΔTE rows carry only the mRNA log2FC. |

#### `mitoribopy all` — end-to-end orchestrator

`mitoribopy all` runs the three stages above with one shared YAML config. Required inputs are the **union** of every active stage's required inputs. The config has these top-level sections:

| Section | When to include | Notes |
|---|---|---|
| `samples:` | recommended | Top-level `samples: { table: samples.tsv }` is the canonical declaration of every input FASTQ + per-sample metadata. Auto-wires both `align` and `rnaseq`. See [Sample sheet](#sample-sheet-unified-per-project-tsv). |
| `align:` | required for align | Indexes, kit / strandedness, dedup. Omit `align.fastq` when `samples:` is set — the sheet supplies it. |
| `rpf:` | required for rpf | Strain, RPF window, FASTA, offset bounds. Auto-wires `--directory`, `--read_counts_file`, `--fasta`. |
| `rnaseq:` | required for rnaseq | Either flow's keys (the two are mutually exclusive). The sheet auto-wires `rnaseq.sample_sheet` so you do not repeat it. |

Use `mitoribopy all --print-config-template > pipeline_config.yaml` to drop a fully-commented starter into your project.

### Sample sheet (unified per-project TSV)

A single TSV declares every sample once, replacing the old pair of stage-specific tables (`--sample-overrides` and `--condition-map`). It is **the recommended way** to declare inputs for any non-trivial project: pairings between Ribo-seq and RNA-seq are by `sample_id` (never by index), per-sample kit / UMI overrides for mixed batches live in the same file, and an `exclude` column lets you drop a bad library without deleting rows.

**Required columns:** `sample_id`, `assay` (`ribo` or `rna`), `condition`, `fastq_1`.
**Optional columns:** `replicate`, `fastq_2`, `kit_preset`, `adapter`, `umi_length`, `umi_position`, `strandedness`, `dedup_strategy`, `exclude` (`true`/`false`/blank), `notes`.

Empty cells (`""`, `NA`, `None`, `-`, `null`) read as "use the default". Lines starting with `#` and blank lines are ignored. Validation is strict: a single load pass reports every row error so you can fix the sheet without iterate-and-retry.

```tsv
sample_id	assay	condition	replicate	fastq_1	fastq_2	kit_preset	umi_length	umi_position	strandedness	exclude	notes
WT_Ribo_1	ribo	WT	1	fastq/WT_Ribo_1.fq.gz		illumina_truseq_umi	8	5p	forward	false	
WT_Ribo_2	ribo	WT	2	fastq/WT_Ribo_2.fq.gz		illumina_truseq_umi	8	5p	forward	false	
KO_Ribo_1	ribo	KO	1	fastq/KO_Ribo_1.fq.gz		illumina_truseq_umi	8	5p	forward	false	
KO_Ribo_2	ribo	KO	2	fastq/KO_Ribo_2.fq.gz		illumina_truseq_umi	8	5p	forward	true	contaminated lane
WT_RNA_1	rna	WT	1	fastq/WT_RNA_1_R1.fq.gz	fastq/WT_RNA_1_R2.fq.gz	pretrimmed	0	5p	forward	false	
KO_RNA_1	rna	KO	1	fastq/KO_RNA_1_R1.fq.gz	fastq/KO_RNA_1_R2.fq.gz	pretrimmed	0	5p	forward	false	
```

How it threads through the pipeline:

| Stage | What the sheet supplies |
|---|---|
| `mitoribopy align` | The Ribo-seq FASTQ list (`assay='ribo'` rows, `exclude=false`) and a materialised `sample_overrides.tsv` carrying the per-sample kit / UMI / dedup columns. |
| `mitoribopy rnaseq` | The RNA-seq FASTQ list (`assay='rna'`), the Ribo-seq FASTQ list (re-counted from raw FASTQ in this stage's own align machinery), and the `sample_id → condition` mapping that drives the pyDESeq2 contrast. |
| `mitoribopy all` | Both of the above, auto-wired. Top-level `samples: { table: samples.tsv }` (or shorthand `samples: samples.tsv`) is the only place you need to declare inputs. |

Mutual-exclusion rules: when the sheet is set, declaring `align.fastq` / `align.fastq_dir` / `align.samples` / `align.sample_overrides` / `rnaseq.rna_fastq` / `rnaseq.ribo_fastq` / `rnaseq.condition_map` alongside it is an error — pick one input style per stage.

### Input files (file-by-file reference)

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

### Annotation CSV (custom organisms)

Built-in `h.sapiens` and `s.cerevisiae` ship complete annotation tables and need nothing here. For any other organism, supply a per-transcript CSV via `--annotation_file`. The full schema (required vs optional columns, defaults, and a worked example) lives in [Custom organisms](#custom-organisms).

### Codon-table JSON (custom organisms)

The 27 NCBI Genetic Codes are bundled. Pick one with `--codon_table_name` (full picker in [Custom organisms](#custom-organisms)). Supply your own `--codon_tables_file` only when your organism's code is not in the NCBI list.

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

For multi-sample batches, add `max_parallel_samples: N` under the `align:` section of your YAML to align samples concurrently — `--threads` is auto-divided across workers so total CPU use stays ≈ `--threads`. (`mitoribopy all` does not take `--max-parallel-samples` directly at the CLI; the flag is read from the YAML and forwarded to the align stage. The standalone `mitoribopy align` subcommand does accept it on the CLI — see [Example B](#b-just-the-align-stage-on-a-directory-of-fastqs) below.) The joint `rpf` stage stays serial (offset selection is a pooled-across-samples computation). See [Execution / concurrency](#execution--concurrency) under the `mitoribopy align` reference.

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

# Just rnaseq (default flow): raw FASTQ -> pyDESeq2 -> TE / dTE
mitoribopy rnaseq --rna-fastq rna_seq/ --ribo-fastq ribo_seq/ \
  --reference-fasta ref.fa --gene-id-convention bare \
  --condition-map samples.tsv \
  --base-sample control --compare-sample knockdown \
  --output results/rnaseq/

# `--base-sample` / `--compare-sample` are aliases for `--condition-a` /
# `--condition-b`; pick whichever spelling you prefer (the legacy form
# still works). The base condition is the reference (denominator) of
# the WT-vs-X contrast and seeds the labels on every comparison plot.

# Or rnaseq (alternative flow): existing rpf output + an external DE table
mitoribopy rnaseq --de-table de.tsv --gene-id-convention hgnc \
  --ribo-dir results/rpf --reference-gtf ref.fa --output results/rnaseq/
```

---

## Strain presets and footprint classes

### Strain (`-s` / `--strain`)

The strain preset selects the organism's mitochondrial annotation and codon table. Two organisms ship complete reference data; everything else uses `custom` and supplies its own files (see [Custom organisms](#custom-organisms)).

| Value | Organism | Codon table | Ships annotation? | Ships `-rpf` default? |
|---|---|---|:-:|:-:|
| `h.sapiens` (default) | *Homo sapiens* mt | `vertebrate_mitochondrial` (NCBI #2) | ✓ | ✓ |
| `s.cerevisiae` | *Saccharomyces cerevisiae* mt | `yeast_mitochondrial` (NCBI #3) | ✓ | ✓ |
| `custom` | Any other organism | user-supplied via `--codon_table_name` (built-in NCBI list) or `--codon_tables_file` | ✗ — pass `--annotation_file` | ✗ — pass `-rpf MIN MAX` |

`h` and `y` are also accepted as short synonyms for `h.sapiens` and `s.cerevisiae`.

### Footprint class (`--footprint_class`)

Pair `-s` with `--footprint_class` to pick sensible RPF and unfiltered-length defaults. An explicit `-rpf MIN MAX` or `--unfiltered_read_length_range MIN MAX` always wins over the footprint-class default. Built-in defaults exist for `h.sapiens` and `s.cerevisiae`; for `--strain custom` you must also pass `-rpf`.

| Value | RPF window default | `--unfiltered_read_length_range` default | Use for |
|---|---|---|---|
| `short` | h.sapiens / s.cerevisiae: 16–24 | 10–30 | Truncated RNase products. Sit just below the canonical monosome window; useful for context-dependent pausing and as a QC indicator of digest aggressiveness. |
| `monosome` (default) | h.sapiens: 28–34, s.cerevisiae: 37–41 | 15–50 | Single-ribosome footprints. The standard mt-Ribo-seq class. |
| `disome` | h.sapiens: 50–70, s.cerevisiae: 60–90 | 40–100 | Collided-ribosome footprints. eIF5A-depletion, queueing, ribosome-stalling studies. |
| `custom` | user must pass `-rpf` | unchanged | Any non-standard footprint class. |

---

## Subcommand reference

Every subcommand inherits these shared options:

| Flag | Default | Description |
|---|---|---|
| `--config PATH` | — | Configuration file (.json, .yaml, .yml, or .toml). CLI flags override values from the file. |
| `--dry-run` | off | Print planned actions and exit 0 without executing. |
| `--threads N` | 1 | Preferred thread count; exports `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `MITORIBOPY_THREADS`. When combined with `--max-parallel-samples M` (align only), each parallel worker's external tools see `max(1, N // M)` threads so the total CPU budget stays ≈ N. |
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
| `--no-pretrimmed-inference` | off | Disable the auto-fallback to `pretrimmed`. With this flag, adapter detection failure with no `--kit-preset` fallback raises an error instead of silently routing to the `pretrimmed` kit. |
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

Vendor-specific kit names (`truseq_smallrna`, `nebnext_smallrna`, `nebnext_ultra_umi`, `truseq_stranded_total`, `smarter_pico_v3`, `sequoia_express`) are also accepted as synonyms for the adapter-family preset they map to.

##### Adapter detection

| Mode | Behaviour |
|---|---|
| `auto` (default) | Scan each FASTQ; pick the matching preset per sample. Samples whose scan fails fall back to the user's `--kit-preset` / `--adapter` if supplied; otherwise fall through to `pretrimmed` when the data shows no adapter signal at all. |
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
| `--dedup-strategy {auto,umi-tools,skip}` | `auto` | Per-sample resolved. `auto` → `umi-tools` for UMI samples, `skip` otherwise. When the resolved strategy is `skip`, the orchestrator does **not** write a duplicate `deduped/<sample>.dedup.bam` — the upstream `aligned/<sample>.mapq.bam` is fed straight into BED conversion. The legacy `mark-duplicates` (picard) option was removed in v0.4.5 because coordinate-only dedup destroys codon-occupancy signal on mt-Ribo-seq libraries (see [docs/validation/taco1_ko_regression.md](docs/validation/taco1_ko_regression.md) for the empirical evidence). |
| `--umi-dedup-method {unique,percentile,cluster,adjacency,directional}` | `unique` | umi_tools `--method`. `unique` collapses only on exact coord+UMI match; other methods may over-collapse in low-complexity mt regions. |

#### Execution / concurrency

Per-sample work in `align` (cutadapt → bowtie2 → MAPQ → dedup → BAM→BED) is independent across samples — each sample writes only to its own per-sample paths. `--max-parallel-samples N` runs that work concurrently in a thread pool while the joint `mitoribopy rpf` stage stays serial (offsets are selected across all samples and aggregate `rpf_counts.tsv` requires a single pass).

| Flag | Default | Description |
|---|---|---|
| `--max-parallel-samples N` | `1` | Number of samples to align concurrently. With `--threads T`, each worker's external tools (cutadapt, bowtie2, umi_tools) get `max(1, T // N)` threads, so total CPU use stays ≈ T regardless of N. Default `1` = serial (current behaviour, fully backward-compatible). Resume-cached samples skip the pool entirely. On any per-sample failure the run is fail-fast: pending futures are cancelled and the first exception is re-raised. |

**Sizing guidance.** A reasonable starting point on a workstation or HPC node: pick `T` = total CPU budget (cores you can use), then choose `N` so each worker still gets ≥ 2 threads — i.e. `N ≤ T / 2`. Examples:

| Cores available | Suggested `--threads` | Suggested `--max-parallel-samples` | Per-tool threads |
|---|---|---|---|
| 8 | 8 | 4 | 2 |
| 16 | 16 | 4 | 4 |
| 16 | 16 | 8 | 2 |
| 32 | 32 | 8 | 4 |

I/O-bound stages (FASTQ gzip read, BAM write) can saturate disk bandwidth before they saturate CPU; if your alignment outputs live on slow networked storage, smaller `N` (more threads per worker, fewer concurrent samples) often wins over larger `N`. The joint `rpf` stage is unaffected — it runs in a single process regardless of this flag.

#### Intermediate files

| Flag | Default | Description |
|---|---|---|
| `--keep-intermediates` | off | Keep the per-step intermediate files (`trimmed/<sample>.trimmed.fq.gz`, `contam_filtered/<sample>.nocontam.fq.gz`, `aligned/<sample>.bam` pre-MAPQ). By default these are deleted as soon as the next step has consumed them, since they are large, regenerable, and not consumed by any downstream stage. Pass this flag when debugging a sample or comparing per-step intermediate counts; expect a 2–3× increase in disk footprint per sample. |

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

Plain `mitoribopy <flags>` (no subcommand) routes to `mitoribopy rpf` and prints a one-line notice; prefer the explicit subcommand form.

#### Core inputs

| Flag | Default | Description |
|---|---|---|
| `-f, --fasta REF_FASTA` | — | Reference FASTA for the transcript/annotation context. **Required**. |
| `-s, --strain {h.sapiens,s.cerevisiae,custom}` | `h.sapiens` | Reference preset; see [Strain presets](#strain-presets-and-footprint-classes). |
| `-d, --directory BED_DIR` | cwd | Directory of `.bed` and/or `.bam` files. BAMs are auto-converted to BED6 via pysam. |
| `--bam_mapq Q` | 10 | MAPQ threshold for BAM inputs before BAM→BED6 conversion. Set 0 to disable. |
| `-rpf MIN_LEN MAX_LEN` | strain/class default | Inclusive read-length filter range, e.g. `-rpf 29 34`. |
| `--footprint_class {monosome,disome,custom}` | `monosome` | Picks `-rpf` and `--unfiltered_read_length_range` defaults. |
| `--annotation_file ANNOTATION.csv` | — | Per-transcript annotation CSV. Required for `--strain custom`; the built-in `h.sapiens` and `s.cerevisiae` tables are used otherwise. See [Custom organisms](#custom-organisms) for the schema. |
| `--codon_tables_file CODON_TABLES.json` | — | Codon-table JSON override. |
| `--codon_table_name TABLE_NAME` | strain default | Codon table to load (built-in or from `--codon_tables_file`). 27 built-in tables ship; run `--help` for the full list. |
| `--start_codons CODON [CODON ...]` | strain default | Allowed start codons. Defaults: `h.sapiens` → `ATG ATA`, `s.cerevisiae` → `ATG`, `custom` → `ATG`. |
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
| `--offset_pick_reference {p_site,reported_site}` | `p_site` | `p_site`: pick in canonical P-site space, then convert to the reported space. `reported_site`: pick directly in the final `--offset_site` space (no P↔A conversion). |
| `--offset_type {5,3}` | `5` | Which read end the offset is measured from. |
| `--offset_site {p,a}` | `p` | Coordinate space for the SELECTED OFFSETS table (the values in `p_site_offsets_<align>.csv`). Does NOT control which downstream outputs are generated; use `--analysis_sites` for that. |
| `--analysis_sites {p,a,both}` | `both` | Which downstream sites to generate. `both` writes parallel P-site and A-site codon usage + coverage plots, side by side under per-site subdirs. `p` or `a` restricts to one site. |
| `--codon_overlap_mode {full,any}` | `full` | How a read counts toward the codon-level enrichment table at the anchor codon. `full` (default, right for ribo-seq): the read must cover ALL 3 nt of the anchor codon. Example with anchor codon at positions 101–103 — a 30-nt read at 100–129 counts (read covers 101, 102, 103); a read at 102–131 does NOT count (does not cover 101). `any`: any 1-nt overlap counts; the same 102–131 read above WOULD count. Use `any` only for very short reads relative to a codon (rare for ribo-seq). |
| `-p, --psite_offset NT` | — | Use one fixed offset for every read length and sample. Bypasses enrichment-based selection. |
| `--offset_mode {per_sample,combined}` | `per_sample` | `per_sample`: each sample uses its own offsets; drift surfaced in `offset_drift_<align>.svg` and audited per sample in `offset_diagnostics/csv/per_sample_offset/<sample>/offset_applied.csv`. `combined`: pool all samples and apply one offset table to every sample (use this when individual samples have very low coverage and the per-sample selection is noisy). |

#### Outputs and plotting

| Flag | Default | Description |
|---|---|---|
| `-o, --output DIR` | `analysis_results` | Base output directory. |
| `--downstream_dir NAME` | `footprint_density` | Per-sample subdirectory name for frame and codon analyses. |
| `--plot_dir NAME` | `offset_diagnostics` | Subdirectory name for offset CSVs and plots. CSVs (offset tables, per-sample offsets, summaries) land under `<plot_dir>/csv/`; PNG/SVG diagnostics under `<plot_dir>/plots/`. Per-sample offsets live at `<plot_dir>/csv/per_sample_offset/<sample>/`. |
| `-fmt, --plot_format {png,pdf,svg}` | `png` | File format for saved plots. SVG recommended for publication-ready vector output. |
| `--x_breaks NT [NT ...]` | — | Optional custom x-axis tick marks for offset line plots. |
| `--line_plot_style {combined,separate}` | `combined` | Draw 5'/3' offsets in one panel or two. |
| `--cap_percentile FRAC` | 0.999 | Upper percentile cap for coverage-style plots. |
| `-m, --codon_density_window` | off | When set, the codon-level density at each codon centre is summed with its ±1 nt neighbours (a 3-nt sliding window) before being written into the codon coverage / usage tables and plots. Smooths short-window noise around the codon centre. |
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
| `--mt_mrna_substring_patterns PATTERN [PATTERN ...]` | `mt_genome mt-mrna mt_mrna` | When `--rpm_norm_mode mt_mrna` is selected, the RPM denominator uses only rows from the read-count table whose value in the reference column (set via `--read_counts_reference_col`, or auto-detected) contains ANY of these substrings. Example: in a read-count table with rows for `rrna`, `trna`, `mt_mrna_ND1`, and `mt_mrna_COX1`, the default pattern keeps the latter two and drops the (r/t)RNA rows. |

#### Optional modules

| Flag | Default | Description |
|---|---|---|
| `--structure_density` | off | Export log2 and scaled density values from footprint-density tables. |
| `--structure_density_norm_perc FRAC` | 0.99 | Upper percentile used to cap and scale structure-density values. |
| `--cor_plot` | off | Generate publication-quality codon-correlation plots (one per non-base sample). When `--analysis_sites=both`, produces parallel outputs under `<output>/codon_correlation/p_site/` AND `<output>/codon_correlation/a_site/`; with single-site analysis only the matching folder appears. Each figure is rendered as both SVG (vector, for figures) and 300 dpi PNG (for slides), with: identity (y = x) line, OLS regression line, an `r / R² / p / slope / intercept / N` stat box, an Okabe-Ito colour-blind safe palette per Category, and leader-line labels for the 10 codons farthest from the regression line (placed greedily to avoid overlap). When `--cor_plot` is set but skipped (e.g., `--base_sample` missing), the pipeline log explicitly states the reason. |
| `--base_sample NAME` | — | Reference sample for codon-correlation comparisons. Required by `--cor_plot`. Must match an existing sample directory name. |
| `--cor_mask_method {percentile,fixed,none}` | `percentile` | Masking rule for extreme codon-correlation outliers. The masked-version plot drops codons above the threshold so the bulk of the distribution is visible. |
| `--cor_mask_percentile FRAC` | 0.99 | Used when `--cor_mask_method percentile`. |
| `--cor_mask_threshold FLOAT` | — | Used when `--cor_mask_method fixed`. |
| `--read_coverage_raw` / `--no-read_coverage_raw` | on | Toggle the `coverage_profile_plots/read_coverage_raw[_codon]/` folders. Independent of `--read_coverage_rpm`. |
| `--read_coverage_rpm` / `--no-read_coverage_rpm` | on | Toggle the `coverage_profile_plots/read_coverage_rpm[_codon]/` folders. Independent of `--read_coverage_raw`. |
| `--igv_export` / `--no-igv_export` | off | Write per-sample IGV-loadable BedGraph tracks under `<output>/igv_tracks/<sample>/<sample>_{p_site,a_site}.bedgraph`. P-site tracks default to forest green, A-site to dark orange (set via `track color=` header); load alongside the same FASTA the BEDs were aligned to. |

#### Coverage-profile output reference

`coverage_profile_plots/` always contains two kinds of figures (v0.4.4 flat layout):

- **`read_coverage_*`** (top-level, site-independent) — depth across the transcript using the FULL read footprint (broad peak). The same plot regardless of `--analysis_sites`. Variants: `_rpm` (RPM-normalised), `_raw` (raw counts), `_rpm_codon` / `_raw_codon` (codon-binned across the CDS). Each `_rpm` and `_raw` family is independently toggleable via `--read_coverage_rpm` / `--read_coverage_raw`.
- **`{p_site,a_site}_density_*`** (per-site, single-nucleotide P-site or A-site density, sibling of `read_coverage_*`) — narrow peak at the chosen ribosomal site. The legacy nested `p/` / `a/` subfolders are gone; the site is encoded in the filename prefix. Variants:
  - `{p_site,a_site}_density_rpm` / `_density_raw` — full-transcript density.
  - `{p_site,a_site}_density_rpm_codon` / `_density_raw_codon` — same density but binned per CDS codon.
  - `{p_site,a_site}_density_rpm_frame` / `_density_raw_frame` — **frame-coloured CDS-only nt plots**. Bars are coloured by reading frame (0 / 1 / 2 relative to CDS start) using a colour-blind-safe palette. Frame-0 dominance (~70–90% of CDS density on frame 0) is the canonical mt-Ribo-seq QC signature; if frames 1 and 2 are visible, the library has spurious offset selection, mis-trimmed reads, or contamination from non-ribosome-protected RNA.

For the translation-efficiency integration, use the dedicated `mitoribopy rnaseq` subcommand.

#### Synonyms

The following short / legacy spellings are accepted as synonyms for the canonical names. They emit a single `[mitoribopy] DEPRECATED:` line on first use so users know to migrate, then keep working.

| Synonym | Canonical |
|---|---|
| `-s h` / `--strain h` | `-s h.sapiens` |
| `-s y` / `--strain y` | `-s s.cerevisiae` |
| `--merge_density`, YAML `merge_density:` | `--codon_density_window`, YAML `codon_density_window:` |
| `--mrna_ref_patterns`, YAML `mrna_ref_patterns:` | `--mt_mrna_substring_patterns`, YAML `mt_mrna_substring_patterns:` |
| `--offset_pick_reference selected_site` | `--offset_pick_reference reported_site` |

---

### `mitoribopy rnaseq`

Two ways to drive this subcommand. They are mutually exclusive — passing both `--rna-fastq` and `--de-table` exits with code 2.

**Default flow — from raw FASTQ.** Pass `--rna-fastq` + `--ribo-fastq` + `--reference-fasta` and the subcommand runs the whole pipeline: SE vs PE auto-detected from filename mate tokens, adapter auto-detected via the existing `align.adapter_detect`, UMI presence inferred from per-position Shannon entropy, cutadapt + bowtie2 per sample, reads counted per transcript, pyDESeq2 fit on the RNA side, then a fall-through into the same TE / ΔTE / plot path. When any condition has only one biological sample, the run **fails fast** (publication-safe default) with a message naming the offending condition and pointing at three resolutions: supply biological replicates, switch to the `--de-table` flow with externally-computed DE, or opt into the FASTQ-record-parity pseudo-replicate fallback by passing `--allow-pseudo-replicates-for-demo-not-publication` for tutorial / smoke-test use. When the opt-in flag is used the run continues but stamps `pseudo_replicate_mode: true` in `run_settings.json` and writes an `EXPLORATORY.md` sidecar that calls out which outputs are NOT defensible (padj, "significant" markers, dispersion estimates). The reference-consistency hard-gate is skipped in this flow — instead the FASTA SHA256 is recorded in `run_settings.json` under `from_fastq.reference_checksum`. Requires the `[fastq]` optional-dependency extra:

```bash
pip install 'mitoribopy[fastq]'
```

**Alternative flow — bring your own DE table.** Pass `--de-table` from a prior external DESeq2 / Xtail / Anota2Seq run together with a prior `mitoribopy rpf` run (`--ribo-dir`). This is the right path for publication-grade DE statistics — run DE on the full transcriptome in R / Python and feed the result here. Enforces a SHA256 reference-consistency gate: Ribo-seq and RNA-seq must hash to the identical transcript set, otherwise the run aborts with a `MISMATCH` banner.

Both flows produce the same `te.tsv` / `delta_te.tsv` / `plots/` shape. The plot set splits into three tiers by availability:

- **Always emitted (both flows):** `mrna_vs_rpf`, `delta_te_volcano`, `ma`, `de_volcano_mrna`, `te_bar_by_condition`, `te_heatmap` — six comparison plots driven only by the DE table + Ribo counts.
- **Conditional (both flows, when `--condition-map` + `--base-sample` + `--compare-sample` are all set):** `te_compare_scatter`, `te_log2fc_bar` — these need replicate-per-condition information to compute per-condition TE means, so they are skipped in `--de-table` runs that omit a condition map.
- **Default flow only:** `sample_pca` (PCA needs the full sample × gene counts matrix that only the default flow has on hand), and `de_volcano_rpf` (driven by an `rpf_de_table.tsv` from a second pyDESeq2 fit on the Ribo-seq subset; skipped with a stderr WARNING when the Ribo subset has fewer than two condition levels).

Every PNG is rendered at **300 dpi** under a shared publication style (Okabe-Ito colour-blind-safe palette, sans-serif fonts, white-bbox gene labels with leader lines, stat boxes on the volcanos, replicate dots overlaid on the bar plot, condition strip on the heatmap, quadrant captions on `mrna_vs_rpf`) and is co-emitted as an **editable-text SVG sidecar** with the same stem (`mrna_vs_rpf.png` + `mrna_vs_rpf.svg`, etc.). The default flow additionally writes intermediate counts matrices, a generated `de_table.tsv`, and the `rpf_de_table.tsv` mentioned above.

#### Inputs (default flow, from raw FASTQ)

| Flag | Default | Description |
|---|---|---|
| `--rna-fastq PATH [PATH ...]` | — | RNA-seq FASTQ files or directories. **Required for the default flow.** Mutually exclusive with `--de-table`. |
| `--ribo-fastq PATH [PATH ...]` | — | Ribo-seq FASTQ files or directories. Optional; when omitted the run short-circuits after writing `de_table.tsv` (manifest mode `from-fastq-rna-only`). |
| `--reference-fasta PATH` | — | Transcriptome FASTA. **Required for the default flow.** SHA256 recorded under `from_fastq.reference_checksum` in `run_settings.json`. |
| `--bowtie2-index PREFIX` | — | Pre-built bowtie2 index prefix; skips `bowtie2-build`. `--reference-fasta` is still required for hashing. |
| `--workdir DIR` | `<output>/work` | Scratch directory for trim / index / per-sample BAM artefacts. |
| `--align-threads N` | `4` | Threads passed to cutadapt and bowtie2. (`--threads` separately drives BLAS / pyDESeq2 thread caps.) |
| `--allow-pseudo-replicates-for-demo-not-publication` | off | **Opt-IN** to the FASTQ-record-parity split for any condition with n=1. Default behaviour without this flag is a hard error (exit 2) — pyDESeq2 cannot fit dispersion on n=1 designs, and pseudo-replicates produce padj that is **NOT** publication-grade. Use only for demos / tutorials. When set, `run_settings.json` records `pseudo_replicate_mode: true` and an `EXPLORATORY.md` sidecar is written. The legacy `--no-auto-pseudo-replicate` flag is accepted as a deprecated no-op. |

In the default flow the `--condition-map` table must list every input sample (RNA + Ribo). When pseudo-replicate mode is opted into, the auto-generated `rep1` / `rep2` sample names are written to `<output>/condition_map.augmented.tsv` and the existing TE / ΔTE step picks them up automatically. PE + UMI is currently `NotImplementedError` — preprocess UMIs into the read name first, or use the alternative `--de-table` flow with your own pre-computed DE results.

#### Alternative flow inputs (bring your own DE table)

| Flag | Default | Description |
|---|---|---|
| `--de-table PATH` | — | DE results table (CSV or TSV). **Required to enter the alternative flow.** Mutually exclusive with `--rna-fastq`. |
| `--de-format {auto,deseq2,xtail,anota2seq,custom}` | `auto` | DE table format; `auto` detects from column headers. |
| `--de-gene-col NAME` | — | Column name for gene IDs. Used when `--de-format custom`. |
| `--de-log2fc-col NAME` | — | Column name for log2 fold change. |
| `--de-padj-col NAME` | — | Column name for adjusted p-value. |
| `--de-basemean-col NAME` | — | Column name for basemean. |

#### Gene identifiers

| Flag | Default | Description |
|---|---|---|
| `--gene-id-convention {ensembl,refseq,hgnc,mt_prefixed,bare}` | — | **Required, no default**. Mismatched conventions silently produce zero-match runs. Examples: `hgnc` → `MT-ND1`; `bare` → `ND1`; `ensembl` → `ENSG00000198888`; `refseq` → `YP_003024026.1`; `mt_prefixed` → `MT-ND1`. |
| `--organism {h.sapiens,s.cerevisiae}` | `h.sapiens` | Organism for the mt-mRNA registry. Short / spelled-out forms (`h`, `y`, `human`, `yeast`) are accepted as synonyms. |

#### Ribo-seq inputs (alternative flow only)

| Flag | Default | Description |
|---|---|---|
| `--ribo-dir DIR` | — | Output directory of a prior `mitoribopy rpf` run. Must contain `rpf_counts.tsv` and `run_settings.json` (or `run_manifest.json`) with a recorded `reference_checksum`. Required when using `--de-table`. |
| `--ribo-counts PATH` | `<ribo-dir>/rpf_counts.tsv` | Explicit path to the per-sample per-gene RPF counts table. |

#### Reference-consistency gate (alternative flow only; exactly one)

| Flag | Default | Description |
|---|---|---|
| `--reference-gtf PATH` | — | Reference used by RNA-seq; `mitoribopy rnaseq` hashes this and verifies it matches the hash in the rpf run's manifest. |
| `--reference-checksum SHA256` | — | Precomputed SHA-256, useful when the reference file is not on the rnaseq host. |

#### Conditions

| Flag | Default | Description |
|---|---|---|
| `--condition-map PATH` | — | TSV with columns `sample` and `condition`. **Required in the default flow** (drives the pyDESeq2 contrast). In the `--de-table` flow it is optional and enables a replicate-based Ribo log2FC for ΔTE. |
| `--base-sample NAME` | — | **Reference condition** (denominator of the contrast; the "baseline" used as the x-axis of every TE / DE comparison plot). Required in the default flow. Alias for `--condition-a` — pick whichever spelling you prefer; the YAML key matches `base_sample` for parity with the `rpf` config block. When both `--base-sample` and `--condition-a` are set they must agree. |
| `--compare-sample NAME` | — | **Comparison condition** (numerator of the contrast). Required in the default flow. Alias for `--condition-b`; same equality rule applies when both are set. |
| `--condition-a NAME` | — | Legacy spelling of `--base-sample` (still accepted; identical semantics). |
| `--condition-b NAME` | — | Legacy spelling of `--compare-sample`. |

Without a condition map in the `--de-table` flow, `mitoribopy rnaseq` computes a point-estimate ΔTE using only the DE table's mRNA log2FC and emits rows with a `single_replicate_no_statistics` note.

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
| `--print-canonical-config` | — | Load `--config`, apply every auto-wiring + sample-sheet expansion that `mitoribopy all` would normally apply, then print the resulting **canonical config** to stdout (YAML if PyYAML is available, JSON otherwise) and exit. The same blob is embedded in `run_manifest.json` under `config_canonical` on real runs, so this is the right way to preview "what is `mitoribopy all` actually going to run with?" before committing to a full pipeline. |

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

## What the numbers mean — RNA, RPF, TE, ΔTE

Translation efficiency is one of the most-asked-for and most-misread quantities in Ribo-seq. This section pins each quantity to a specific output file and a specific equation so reviewers don't have to guess what `te.tsv` means versus a coverage-ratio plot.

### The four quantities

| # | Quantity | Definition | What it is NOT |
|---|---|---|---|
| 1 | **RNA abundance** | Per-sample, per-transcript RNA-seq counts (or a normalised abundance estimate downstream of pyDESeq2). | Not a coverage profile, not a Ribo-seq quantity. |
| 2 | **RPF abundance** | Per-sample, per-transcript ribosome footprint counts after offset-aware filtering. | Not raw aligned-read counts (the read-length auto-filter prunes noise bins) and not a per-position density. |
| 3 | **TE (translation efficiency)** | Per-sample, per-gene RPF abundance normalised by RNA abundance. See equation below. | Not an inferential statistic — TE is a point estimate per (sample, gene). |
| 4 | **ΔTE (delta-TE)** | Per-gene log2 change in TE between two conditions. See equation below. | Not the ratio of two TEs computed sample-by-sample; it is a contrast on the log scale. |

### TE equation (`te.tsv`)

For sample *s* and gene *g*:

```
TE(s, g) = (RPF_count(s, g) + δ) / (mRNA_abundance(g) + δ)
```

where δ is the package's pseudocount (default 0.5; see [`src/mitoribopy/rnaseq/te.py`](MitoRiboPy-refactor/src/mitoribopy/rnaseq/te.py)). The pseudocount goes on **both** numerator and denominator to avoid div-by-zero and log(0). Per-gene mRNA abundance is the DE table's `baseMean` column (DESeq2 / Xtail / Anota2Seq) or a pyDESeq2-derived estimate in the from-FASTQ flow. Genes missing from the DE table are skipped — TE is undefined without an mRNA denominator.

Output: `<output>/rnaseq/te.tsv` with columns `sample, gene, rpf_count, mrna_abundance, te`. One row per (sample, gene) that has both an RPF count and a mRNA abundance.

### ΔTE equation (`delta_te.tsv`)

For gene *g*, comparing condition *B* to condition *A* (with *A* = `--base-sample` / `--condition-a`):

```
ΔTE_log2(g) = log2(RPF_FC) - log2(mRNA_FC)

         where  RPF_FC  = mean(RPF_B) / mean(RPF_A)   (with pseudocount on each mean)
                mRNA_FC = exp2(mRNA_log2fc(g))         (from the DE table)
```

ΔTE is **not** computed as `TE(B) / TE(A)` per sample. It is a contrast on log-fold-changes: each side's denominator (RNA abundance) and numerator (RPF abundance) are aggregated separately before the ratio. This is the standard formulation used by Xtail, Anota2Seq, and downstream summary statistics, and it is what removes the spurious sample-level pseudocount asymmetry you'd get from a per-sample TE ratio.

Output: `<output>/rnaseq/delta_te.tsv` with columns `gene, mrna_log2fc, rpf_log2fc, delta_te_log2, padj, note`. The `note` column carries any qualifier MitoRiboPy attached to the row:

| `note` value | Meaning |
|---|---|
| `""` (empty) | Both RPF and mRNA log2FC available; ΔTE is computed from replicate-based means. |
| `single_replicate_no_statistics` | No condition map / contrast was supplied (or only one replicate per condition); ΔTE row carries only the mRNA log2FC and `rpf_log2fc=None`. |
| `insufficient_ribo_replicates` | Condition map present but the gene has zero Ribo counts in one of the two conditions; `rpf_log2fc=None` and ΔTE is `None`. |
| `missing_from_de_table` | Gene is in the Ribo counts but absent from the DE table; everything except `gene` is `None`. |

### What is NOT in te.tsv / delta_te.tsv

These two files are **gene-level summary tables**. If you want any of the following, look elsewhere:

| You want | Look at |
|---|---|
| Per-position P-site / A-site density | `<output>/translation_profile/<sample>/footprint_density/<transcript>_footprint_density.csv` |
| Per-codon usage by sample | `<output>/translation_profile/<sample>/codon_usage/` |
| Periodicity / 3-nt phasing diagnostic | `<output>/qc/periodicity_metagene.png` (start- and stop-aligned) and `qc/frame_summary.tsv` |
| The descriptive coverage-normalised RPF metric (the legacy "rna-seq ratio" the v0.2.x module emitted) | **Removed in this refactor.** Use `te.tsv` for sample-level TE, or compute your own ratio from the per-position density CSVs above. |

### Pseudo-replicate runs

When a from-FASTQ flow run is launched with `--allow-pseudo-replicates-for-demo-not-publication`, both `te.tsv` and `delta_te.tsv` are still written, but `run_settings.json` carries `pseudo_replicate_mode: true` and an `EXPLORATORY.md` sidecar lists every column you should NOT cite (padj, "significant" markers, dispersion estimates). The TE and ΔTE point estimates remain readable as exploratory; only the inferential statistics are unsafe.

---

## Output overview

For a `mitoribopy all` run with the defaults (`--offset_mode per_sample`, `--analysis_sites both`):

```text
<output>/
  run_manifest.json                    # composed provenance: align + rpf + rnaseq settings,
                                       # tool versions, reference_checksum, stages_run
  align/
    mitoribopy.log
    read_counts.tsv                    # per-stage counts (input → trimmed → contam → mt → MAPQ → dedup)
    kit_resolution.tsv                 # per-sample kit + dedup decisions; the provenance spine
    run_settings.json                  # includes per_sample[] block with detected_kit, applied_kit,
                                       # match_rate, dedup_strategy, source per sample
    trimmed/                           # only *.cutadapt.json (per-sample);
                                       # the trimmed *.fq.gz is deleted after
                                       # contam_filter runs unless
                                       # --keep-intermediates is passed
    contam_filtered/                   # empty unless --keep-intermediates
    aligned/                           # *.mapq.bam (kept; pre-mapq *.bam is
                                       # deleted unless --keep-intermediates)
    deduped/                           # only created when at least one sample
                                       # uses umi-tools; for dedup=skip the
                                       # orchestrator wires mapq.bam straight
                                       # into BED conversion so no duplicate
                                       # dedup.bam is written
    bed/                               # *.bed -- strand-aware BED6 inputs to rpf
  rpf/
    mitoribopy.log
    rpf_counts.tsv                     # per-sample per-gene RPF counts; feeds rnaseq
    run_settings.json                  # includes reference_checksum (SHA256 of --fasta)
    offset_diagnostics/                # renamed from plots_and_csv/ in v0.4.4
      csv/                             # all CSV diagnostics live here
        offset_<align>.csv             # COMBINED enrichment summary
        p_site_offsets_<align>.csv     # COMBINED selected offsets
        read_length_summary.csv        # filtered RPF window per sample
        per_sample_offset/             # renamed from per_sample/ in v0.4.4
          <sample>/
            offset_<align>.csv         # per-sample enrichment summary
            p_site_offsets_<align>.csv # per-sample selected offsets
            offset_applied.csv         # exact offsets row applied to <sample>
                                       # (audit so reviewers can verify
                                       # per-sample selection from disk)
      plots/                           # all PNG/SVG diagnostics live here
        offset_drift_<align>.svg       # per-sample 5'/3' offset comparison; read this FIRST
        offset_*.svg                   # heatmaps + line plots
        unfiltered_heatmap_*.svg
        <sample>_read_length_distribution.svg
    translation_profile/               # FLAT layout (v0.4.4): one folder per sample,
      <sample>/                        # site is encoded in filename prefix.
        footprint_density/             # <transcript>_footprint_density.csv (cols:
                                       #   Position, Nucleotide, A_site, P_site)
                                       # + <transcript>_p_site_depth.png and/or
                                       #   <transcript>_a_site_depth.png
        translating_frame/             # p_site_frame_usage_total.csv, p_site_*_by_transcript.csv,
                                       # a_site_frame_usage_total.csv, a_site_*_by_transcript.csv,
                                       # plus matching _plot.png files
        codon_usage/                   # p_site_codon_usage_<transcript>.csv (P-site, stop-masked),
                                       # p_site_codon_usage_total.csv,
                                       # a_site_codon_usage_<transcript>.csv (A-site, no mask),
                                       # a_site_codon_usage_total.csv, plus matching _plot.png files
    coverage_profile_plots/            # FLAT layout (v0.4.4): site is the filename prefix.
      read_coverage_rpm/               # SITE-INDEPENDENT (gated by --read_coverage_rpm)
      read_coverage_raw/               # SITE-INDEPENDENT (gated by --read_coverage_raw)
      read_coverage_rpm_codon/
      read_coverage_raw_codon/
      p_site_density_rpm/              # P-site nt-resolution density (RPM)
      p_site_density_raw/              # P-site nt-resolution density (raw)
      p_site_density_rpm_codon/        # codon-binned across CDS, RPM
      p_site_density_raw_codon/        # codon-binned across CDS, raw
      p_site_density_rpm_frame/        # CDS-only nt plot, frame-coloured (0/1/2 vs CDS start);
                                       # frame-0 dominance is the canonical mt-Ribo-seq QC.
      p_site_density_raw_frame/        # same but raw counts
      a_site_density_*/                # A-site mirror folders (only when analysis_sites in {a,both})
    structure_density/                 # if --structure_density (uses P-site by default)
    codon_correlation/                 # if --cor_plot
      p_site/<base>_vs_<sample>_{all,masked}.{csv,svg,png}    # P-site correlation
      a_site/<base>_vs_<sample>_{all,masked}.{csv,svg,png}    # A-site correlation
    igv_tracks/                        # if --igv_export
      <sample>/
        <sample>_p_site.bedgraph       # P-site footprint density (forest green track)
        <sample>_a_site.bedgraph       # A-site footprint density (dark orange track)
  rnaseq/                              # if rnaseq config supplied
    te.tsv                             # one row per (sample, gene): rpf_count, mrna_abundance, te
    delta_te.tsv                       # one row per gene: mrna_log2fc, rpf_log2fc, delta_te_log2,
                                       # padj, note
    plots/                               # every entry below also has a sibling .svg
                                         # at 300 dpi (PNG) + editable-text SVG (vector,
                                         # Illustrator-friendly). Okabe-Ito colour-blind-
                                         # safe palette across all plots; titles wear
                                         # `<base> vs <compare>` from --base-sample /
                                         # --compare-sample so reviewers know the
                                         # contrast at a glance.
      mrna_vs_rpf.png   (+ .svg)         # four-quadrant log2FC scatter (mRNA vs RPF) with
                                         # faint quadrant captions (co-regulated up / buffered
                                         # up / co-regulated down / buffered down) and points
                                         # coloured by quadrant; identity y=x line drawn
      delta_te_volcano.png   (+ .svg)    # ΔTE (log2) vs -log10(padj) with stat box
                                         # (n_TE_up / n_TE_down / n_n.s.) and threshold
                                         # guides at padj=0.05 + |ΔTE|=1
      ma.png   (+ .svg)                  # log10(baseMean) vs log2FoldChange; sig up = vermillion,
                                         # sig down = blue, n.s. = grey; sig genes labelled
                                         # with white-bbox + leader line
      de_volcano_mrna.png   (+ .svg)     # mRNA DE volcano (log2FC vs -log10 padj),
                                         # vermillion = sig up, blue = sig down, grey = n.s.,
                                         # threshold guides at padj=0.05 and |L2FC|=1;
                                         # stat box reports counts; smart label placer
                                         # spreads labels around dense clusters
      de_volcano_rpf.png   (+ .svg)      # Ribo-seq DE volcano (default flow only;
                                         # second pyDESeq2 fit on the Ribo subset; skipped
                                         # with stderr WARNING when Ribo subset has fewer
                                         # than two condition levels)
      te_bar_by_condition.png   (+ .svg) # log2(TE) per gene, bars grouped by condition + SE
                                         # error bars; every replicate also drawn as a
                                         # black dot jittered along x within the bar's
                                         # footprint (so the sample size behind each mean
                                         # is visible by eye)
      te_heatmap.png   (+ .svg)          # gene × sample log2(TE) heatmap (RdBu_r, centred
                                         # at 0); coloured condition strip above the columns
                                         # spells the assignment in white-bold; cell-value
                                         # annotations get auto-contrast (white on saturated,
                                         # black on faint)
      te_compare_scatter.png   (+ .svg)  # per-gene mean log2(TE) in <base> (x) vs <compare>
                                         # (y) with identity y=x line; per-replicate cloud
                                         # behind gene-mean dots; gene-means coloured by
                                         # direction; Pearson r in stat box. Emitted only
                                         # when --condition-map + --base-sample +
                                         # --compare-sample are all set.
      te_log2fc_bar.png   (+ .svg)       # sorted bar of log2(TE_compare / TE_base) per gene
                                         # with the numeric log2FC printed at each bar's
                                         # tip; bars above zero (vermillion) mean TE up
                                         # in the compare condition. Same conditions as
                                         # te_compare_scatter.
      sample_pca.png   (+ .svg)          # PC1 vs PC2 from log1p counts; condition = colour,
                                         # assay = marker shape; variance % on axes
                                         # (default flow only — needs the full counts matrix)
    run_settings.json
    # --- default flow (--rna-fastq) only -----------------------------
    de_table.tsv                       # mRNA pyDESeq2 result in canonical DESeq2 schema
    rpf_de_table.tsv                   # Ribo-seq pyDESeq2 result (same schema, drives
                                       # de_volcano_rpf.png; absent when the Ribo subset
                                       # has fewer than two condition levels)
    rna_counts.tsv                     # wide gene × sample RNA-seq counts
    rpf_counts.tsv                     # long-format Ribo-seq counts (sample\tgene\tcount)
    rpf_counts_matrix.tsv              # wide gene × sample Ribo-seq counts
    condition_map.augmented.tsv        # original entries + auto-generated rep1 / rep2 names
```

In v0.4.4 the legacy `p/` / `a/` subfolders under `translation_profile/` and `coverage_profile_plots/` were dropped. The site is encoded in filename prefixes (`p_site_*` / `a_site_*`) so a single per-sample folder hosts both sites without ambiguity, and downstream tooling never has to branch on `--analysis_sites`.

### What should I inspect first?

After every run, walk these files in order. The first three (`SUMMARY.md`,
`warnings.tsv`, `summary_qc.tsv`) catch ~90% of the gotchas without having
to dig into per-stage TSVs.

1. `SUMMARY.md` — one-page human-readable view of what ran and what was produced.
2. `warnings.tsv` — every structured warning (`stage / sample_id / severity / code / message / suggested_action`). An empty file (header only) means nothing was flagged.
3. `summary_qc.tsv` — per-sample QC across all stages.
4. `align/kit_resolution.tsv` — confirm each sample resolved to the right kit / UMI / dedup strategy.
5. `align/read_counts.tsv` — per-stage read funnel (trim → contam → mt → MAPQ → dedup).
6. `rpf/offset_diagnostics/plots/offset_drift_*.svg` — per-sample offset drift in a single glance.
7. `rpf/offset_diagnostics/csv/per_sample_offset/<sample>/offset_applied.csv` — exact offset table actually applied.
8. `rpf/rpf_counts.tsv` (+ `rpf/rpf_counts.metadata.json` sidecar) — per-(sample, gene) RPF count matrix.
9. `rnaseq/te.tsv` and `rnaseq/delta_te.tsv` — translation efficiency and ΔTE tables (when the rnaseq stage ran).
10. `figure_qc.tsv` (after `mitoribopy validate-figures runs/full/`) — mechanical pass/warn/fail per plot (label overlap, point counts, SVG editability).

### When NOT to trust the output

The pipeline finishes "successfully" (exit 0) under a number of conditions
that should still make a reviewer pause. Check `warnings.tsv` for any of
the codes below; each row's `suggested_action` tells you what to do.

- **Adapter detection low confidence** — `align/kit_resolution.tsv` has `source=user_fallback` or a small `confidence_margin`. The detection picked something but wasn't sure; spot-check the assigned adapter against the FASTQ.
- **High contaminant fraction** — `read_counts.tsv` shows >50% of post-trim reads filtered out as rRNA / contaminant. The library may not be the libraryyou think it is.
- **Low mt-mRNA alignment fraction** — `mt_aligned / post_rrna_filter < 0.05`. The reference, the strain preset, or the input may be off.
- **No clear offset peak** — `rpf/offset_diagnostics/plots/offset_*.svg` shows a flat enrichment landscape. The selected offset was a fallback, not a real peak.
- **Weak frame-0 dominance** — coverage plots in `rpf/coverage_profile_plots/p_site_density_rpm_frame/` show ≪70% frame-0 density. Either contamination or a wrong offset.
- **Sample-specific offset fallback** — a sample's `offset_applied.csv` row shows the global fallback offset rather than the per-sample optimum. Either too few reads or a flat landscape.
- **Pseudo-replicate mode** — `rnaseq/EXPLORATORY.md` exists. The DE statistics are exploratory only — re-run with biological replicates for publication.
- **Inferred UMI without declaration** — `kit_resolution.tsv` shows `umi_source=inferred`. Confirm with the wet-lab record before trusting dedup.
- **Gene-ID match rate below threshold** — `warnings.tsv` has a `GENE_ID_MATCH_RATE_LOW` row. Your DE table and rpf reference may not match.
- **`figure_qc.tsv` has any `fail` row** — overlapping labels, clipped text, or a point-count mismatch with the source TSV. A `--strict` `validate-figures` run will exit 2 in this case.

### Publication route (recommended for figures + tables in a paper)

For a publication-grade analysis, run the pipeline on the full
transcriptome externally (DESeq2 / Xtail / Anota2Seq) and let MitoRiboPy
handle the mt-translation-efficiency layer:

1. `mitoribopy align` + `mitoribopy rpf` to get RPF counts.
2. Run RNA-seq DE externally on the FULL transcriptome with your tool of choice.
3. `mitoribopy rnaseq --rnaseq-mode de_table --de-table de.tsv --ribo-dir runs/full/rpf/` to compute TE / ΔTE on top.
4. `mitoribopy validate-figures runs/full/ --strict` to mechanically QC every plot.
5. Bundle `run_manifest.json`, `summary_qc.tsv`, `figure_qc.tsv`, and the SVG sidecars into the paper's supplement / repository — together they record every parameter and pass each plot through a defined QC contract.

### Key files to look at first

- **`align/kit_resolution.tsv`** — confirm each sample resolved to the right kit and dedup strategy. The `source` column tells you whether it was `detected`, `user_fallback`, `inferred_pretrimmed`, or set explicitly.
- **`align/read_counts.tsv`** — track per-stage drop-off. The invariants `rrna_aligned + post_rrna_filter == post_trim` and `mt_aligned + unaligned_to_mt == post_rrna_filter` must hold; if they don't, something is wrong with alignment.
- **`rpf/offset_diagnostics/plots/offset_drift_<align>.svg`** — per-sample offset comparison. Outliers are visible by eye in seconds.
- **`rpf/offset_diagnostics/plots/offset_<align>.svg`** (heatmap + line plots) — the combined enrichment around the anchor codon. A sharp peak at the canonical P-site offset (~12–15 nt from the 5' end for most mt-Ribo-seq libraries) confirms library quality.
- **`rpf/offset_diagnostics/csv/per_sample_offset/<sample>/offset_applied.csv`** — exact offsets row applied to `<sample>` downstream. New in v0.4.4 so a reviewer can confirm from disk that per-sample offset selection was honoured by the translation-profile and coverage-profile steps.
- **`rpf/translation_profile/<sample>/codon_usage/p_site_codon_usage_total.csv`** — overall P-site codon occupancy. Compare against `a_site_codon_usage_total.csv` (same folder) to see the A-site picture.
- **`rpf/coverage_profile_plots/p_site_density_rpm_frame/<transcript>_*.svg`** — frame-coloured single-nucleotide P-site density across the CDS only. Bars are coloured by reading frame (0 / 1 / 2 relative to CDS start) using a colour-blind safe palette. **Frame-0 dominance is the canonical mt-Ribo-seq QC signature**: a healthy library shows ~70–90% of the CDS density on frame 0 with much smaller frames 1 and 2; flat or jittery frames suggest contamination, poor offset selection, or a low-complexity region. The `p_site_density_raw_frame/` companion has the same plots in raw counts (use this when comparing against published per-codon counts; otherwise prefer the RPM version for cross-sample comparison).
- **`rpf/codon_correlation/{p_site,a_site}/<base>_vs_<other>_*.svg|png`** — publication-quality scatter of codon usage between samples, one folder per requested site (only when `--cor_plot` is set with `--base_sample`). Includes identity (y = x) line, OLS regression, an `r / R² / p / slope / intercept / N` stat box, and leader-line labels for the 10 codons farthest from the regression line.
- **`rpf/igv_tracks/<sample>/<sample>_{p_site,a_site}.bedgraph`** — IGV-loadable footprint density tracks (only when `--igv_export` is set). Open them alongside the FASTA the BEDs were aligned to; P-site is forest green, A-site dark orange.

---

## Custom organisms

To analyse an organism other than human or yeast, run with `--strain custom` and supply three things: a **reference FASTA** of mt-mRNA transcripts, a **per-transcript annotation CSV**, and a **codon-table choice**. Start codons can be left at the `[ATG]` default or overridden with `--start_codons`.

### Picking a codon table

MitoRiboPy bundles every NCBI Genetic Code as a named table. The names match NCBI's organism-group labels; pick the one for your organism's mitochondrial code (or nuclear code for ciliates / *Candida* / similar).

| `--codon_table_name` | NCBI # | Use for |
|---|:---:|---|
| `standard` | 1 | Most plant mitochondria and plastids; many fungal nuclear genomes |
| `vertebrate_mitochondrial` | 2 | Mouse, rat, zebrafish, Xenopus, chicken — any non-human vertebrate mt |
| `yeast_mitochondrial` | 3 | *S. cerevisiae* and close relatives (matches `s.cerevisiae` preset) |
| `mold_mitochondrial` | 4 | *Neurospora*, *Aspergillus*, *Trichoderma*, mycoplasmas |
| `invertebrate_mitochondrial` | 5 | *Drosophila*, mosquito, *C. elegans* and other invertebrates |
| `echinoderm_mitochondrial` | 9 | Sea urchin, starfish |
| `ascidian_mitochondrial` | 13 | Sea squirts (tunicates) |
| `alternative_flatworm_mitochondrial` | 14 | Flatworms |
| `trematode_mitochondrial` | 21 | Trematodes |

Run `mitoribopy rpf --help` to see the full list of all 27 bundled tables (NCBI #1–34, with a few historical numbers omitted by NCBI). Sources: [NCBI Taxonomy Genetic Codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) — the JSON tables under [`src/mitoribopy/data/codon_tables.json`](src/mitoribopy/data/codon_tables.json) are a faithful transcription of those NCBI tables.

If the table you need isn't in the built-in list (e.g. an unusual reassignment in a non-model organism), pass `--codon_tables_file your_tables.json --codon_table_name your_name` to load your own JSON. Each entry maps the 64 codons to single-letter amino-acid codes (use `*` for stop):

```json
{
  "your_organism_mitochondrial": {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TAA": "*", "TAG": "*", "TGA": "W",
    "...": "..."
  }
}
```

### Annotation CSV format

`--annotation_file` takes a CSV with one row per mt-mRNA transcript. The loader requires four columns and accepts four optional columns; CDS bounds (`start_codon`, `stop_codon`, `l_cds`) are computed from the UTR lengths and do NOT need to be supplied.

| Column | Required? | Type | Meaning |
|---|:---:|---|---|
| `transcript` | yes | string | Display name used in plots and CSVs (e.g. `ND1`, `COX1`). Must be unique per row. |
| `l_tr` | yes | int | Total transcript length in nt (`l_utr5 + l_cds + l_utr3`). |
| `l_utr5` | yes | int | Length of the 5' UTR in nt (0 if none). |
| `l_utr3` | yes | int | Length of the 3' UTR in nt (0 if none). |
| `l_cds` | optional | int | Length of the CDS in nt. Computed as `l_tr − l_utr5 − l_utr3` when omitted. |
| `sequence_name` | optional | string | The exact FASTA header / BED chrom field for this transcript. Defaults to the `transcript` column when blank, but **must match your FASTA header** if those names differ. |
| `display_name` | optional | string | Human-readable label used in plot titles. Defaults to `transcript` when blank. |
| `sequence_aliases` | optional | string | Semicolon-separated list of legacy IDs that should also map to this transcript (e.g. `ATP86;ATP8_6` if older BEDs use those names). Leave blank when there are none. |

Minimal example for a hypothetical 3-transcript mouse mt-mRNA reference:

```csv
transcript,sequence_name,l_tr,l_utr5,l_utr3
ND1,mouse_ND1,957,0,0
COX1,mouse_COX1,1545,0,0
ATP6,mouse_ATP6,681,0,0
```

The `sequence_name` field MUST match the FASTA header you pass to `-f` AND the BED `chrom` field produced by `mitoribopy align` (so the bowtie2 indexes have to be built from the same FASTA). [`examples/custom_reference/annotation_template.csv`](examples/custom_reference/annotation_template.csv) is a complete worked example.

### Putting it together: a mouse run

```bash
# Build the bowtie2 indexes from your mouse mt-mRNA FASTA + an rRNA decoy.
bowtie2-build mouse-mt-mRNA.fasta indexes/mt
bowtie2-build mouse_rrna.fa       indexes/rrna

# Run the full pipeline.
mitoribopy all --config mouse_pipeline.yaml --output mouse_results/
```

```yaml
# mouse_pipeline.yaml
align:
  kit_preset: auto
  fastq: input_data/seq
  contam_index: indexes/rrna
  mt_index: indexes/mt

rpf:
  strain: custom
  fasta: mouse-mt-mRNA.fasta
  annotation_file: mouse_annotation.csv      # CSV described above
  codon_table_name: vertebrate_mitochondrial # NCBI #2
  rpf: [28, 34]                              # mouse mt-monosome window
  align: stop
  offset_type: "5"
  offset_site: p
```

For other organisms the recipe is identical — only the codon table name and the annotation CSV change.

### Bicistronic transcript pairs

The two overlapping mt-mRNA pairs (`ATP8/ATP6` and `ND4L/ND4`) are kept as paired display names in plot titles. Choose which member's coordinates seed the merged sequence with `--atp8_atp6_baseline ATP8|ATP6` and `--nd4l_nd4_baseline ND4L|ND4` (defaults: `ATP6` and `ND4`). Legacy FASTA / BED identifiers `ATP86` and `ND4L4` are also recognised through the built-in alias map.

---

## Built-in references

MitoRiboPy ships packaged reference data for two organisms:

- *Homo sapiens* mt-translation using the `vertebrate_mitochondrial` codon table (`-s h.sapiens`)
- *Saccharomyces cerevisiae* mt-translation using the `yeast_mitochondrial` codon table (`-s s.cerevisiae`)

Built-in annotation tables are stored as CSV and the bundled NCBI Genetic Codes as JSON under [src/mitoribopy/data](src/mitoribopy/data). 27 codon tables are available out of the box; the picker for non-human, non-yeast organisms lives in [Custom organisms](#custom-organisms).

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

For batches with many samples, add `--max-parallel-samples N` to align several samples concurrently. `--threads` is divided across workers, so total CPU use stays ≈ 8 below — `4 workers × 2 threads/tool`:

```bash
mitoribopy align \
  --kit-preset auto \
  --library-strandedness forward \
  --fastq-dir input_data/ \
  --contam-index input_data/indexes/rrna_contam \
  --mt-index input_data/indexes/mt_tx \
  --output results/align/ \
  --threads 8 \
  --max-parallel-samples 4
```

The joint `mitoribopy rpf` stage is unaffected by this flag; offset selection and the aggregated `rpf_counts.tsv` always require a single pass over all samples together. See the [Execution / concurrency](#execution--concurrency) reference for sizing guidance.

### C. Just the rpf stage on existing BEDs (human, monosome)

```bash
mitoribopy rpf \
  -s h.sapiens \
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
  --codon_density_window
```

### D. Disome (collided ribosome) study

```bash
mitoribopy rpf \
  -s h.sapiens \
  -f references/human-mt-mRNA.fasta \
  --directory bed/ \
  --footprint_class disome \
  --output results_disome/
```

`--footprint_class disome` widens `-rpf` to 50–70 nt (h.sapiens) or 60–90 nt (s.cerevisiae) and `--unfiltered_read_length_range` to 40–100 nt automatically. You can still override either with `-rpf` / `--unfiltered_read_length_range`.

### E. Short truncated RNase products

```bash
mitoribopy rpf \
  -s h.sapiens \
  -f references/human-mt-mRNA.fasta \
  --directory bed/ \
  --footprint_class short \
  --output results_short/
```

`--footprint_class short` targets the 16–24 nt RNase truncation tail and widens the unfiltered QC window to 10–30 nt so the entire short tail is visible in the read-length heatmap.

### F. Custom organism (mouse mt example)

```bash
mitoribopy rpf \
  -s custom \
  -f references/mouse-mt-mRNA.fasta \
  --directory bed/ \
  -rpf 28 34 \
  --annotation_file references/mouse_annotation.csv \
  --codon_table_name vertebrate_mitochondrial \
  --output results_mouse/
```

See [Custom organisms](#custom-organisms) for the annotation CSV schema, the codon-table picker for other organisms, and how to supply your own JSON when the organism's code is not in the built-in NCBI list.

### G. RNA-seq integration for translation efficiency

**Default flow — raw FASTQ → pyDESeq2 → TE / ΔTE in one shot (`--rna-fastq`).**

Requires the `[fastq]` extra (one-off):

```bash
pip install 'mitoribopy[fastq]'
```

Then:

```bash
mitoribopy rnaseq \
  --rna-fastq input_data/rna_seq/        \
  --ribo-fastq input_data/ribo_seq/      \
  --reference-fasta references/human-mt-mRNA.fasta \
  --gene-id-convention bare              \
  --condition-map samples.tsv            \
  --base-sample control                  \
  --compare-sample knockdown             \
  --output results/rnaseq/               \
  --align-threads 8
```

`--base-sample` is the **reference** condition for the contrast (denominator of log2FC); `--compare-sample` is the **comparison** condition (numerator). Both are required in the default flow. They are aliases for the legacy `--condition-a` / `--condition-b` flags — pick whichever spelling reads better; passing both an alias and the legacy form simultaneously is allowed only when the values agree. The base condition seeds the labels on every WT-vs-X comparison plot (`mrna_vs_rpf`, `delta_te_volcano`, `de_volcano_mrna`, `de_volcano_rpf`, `te_compare_scatter`, `te_log2fc_bar`).

`--rna-fastq` accepts files OR directories; SE vs PE is auto-detected from filename mate tokens (`_R1/_R2`, `_1/_2`, `.1./.2.`, `_R1_001/_R2_001`, `_read1/_read2`). When any condition has only one biological sample the run **fails fast** by default (publication-safe). To opt into the FASTQ-record-parity pseudo-replicate fallback for demo / tutorial use, pass `--allow-pseudo-replicates-for-demo-not-publication`; the run then continues but stamps `pseudo_replicate_mode: true` in `run_settings.json` and writes an `EXPLORATORY.md` sidecar that names the outputs (padj, "significant" markers, dispersion estimates) that are not statistically defensible.

In the default flow, pyDESeq2 is fit **twice** — once on the RNA subset (writes `de_table.tsv`) and once on the Ribo subset (writes `rpf_de_table.tsv`). The Ribo fit is what powers `de_volcano_rpf.png` and lets reviewers see at a glance whether a gene's TE shift is driven by a footprint-level change, an mRNA-level change, or both. If the Ribo subset has fewer than two condition levels (e.g. you passed `--ribo-fastq` for only one condition) the second fit is skipped with a stderr WARNING and the RPF volcano is omitted.

**Alternative flow — bring your own DE table (`--de-table`).**

When you already ran DESeq2 / Xtail / Anota2Seq externally (which is the recommended path for publication-grade DE statistics over the full transcriptome), feed the result via `--de-table` and skip the alignment / pyDESeq2 work inside the subcommand:

```bash
mitoribopy rnaseq \
  --de-table de.tsv \
  --gene-id-convention hgnc \
  --ribo-dir results/rpf \
  --reference-gtf references/human-mt-mRNA.fasta \
  --condition-map samples.tsv \
  --base-sample control \
  --compare-sample knockdown \
  --output results/rnaseq/
```

The two flows are mutually exclusive — passing both `--rna-fastq` and `--de-table` exits with code 2.

Templates: the shell-script template at [examples/templates/run_rnaseq.example.sh](examples/templates/run_rnaseq.example.sh) and the YAML template at [examples/templates/rnaseq_config.example.yaml](examples/templates/rnaseq_config.example.yaml) list every flag for both flows with its default value, with the default flow on top.

### H. Resume a partial run

```bash
mitoribopy all --config pipeline_config.yaml --output results/ --resume
```

`--resume` works at two levels:

1. **Stage level.** Each stage is skipped when its final sentinel file already exists: `align/read_counts.tsv`, `rpf/rpf_counts.tsv`, `rnaseq/delta_te.tsv`.
2. **Per-sample inside the align stage.** Every sample that finishes successfully writes a small JSON marker at `<output>/align/.sample_done/<sample>.json` containing its read-count row. On a subsequent `--resume` run (when `read_counts.tsv` is missing because the previous run crashed mid-batch), samples whose marker is present are reloaded from JSON instead of re-run; the markers from completed samples are merged with the freshly-processed ones to produce the aggregated `read_counts.tsv`. This means a 50-sample run that died at sample 30 picks up at sample 31 the next time you pass `--resume`.

The per-sample markers are written atomically (`.tmp` then rename) so a kill mid-write cannot leave a half-flushed file that would silently feed the wrong row into the aggregated table; corrupt or schema-drifted markers are treated as absent and the sample re-runs cleanly.

### I. Pre-trimmed inputs (e.g. SRA-deposited)

The `pretrimmed` preset name describes the **input state** ("the FASTQ has already been adapter-clipped"), not an action this pipeline takes. SRA-deposited FASTQs typically arrive in this state — the adapter and (when present) the UMI are already stripped before upload.

No special configuration is needed — when adapter detection finds 0% across every kit, the resolver falls through to the `pretrimmed` kit automatically and cutadapt skips the `-a` flag (length and quality filtering still run). To opt in explicitly:

```yaml
align:
  kit_preset: pretrimmed
```

To opt out of the auto-inference and force adapter detection failure to raise instead, pass `--no-pretrimmed-inference` (or set `allow_pretrimmed_inference: false` in YAML).

> **Wording note:** `pretrimmed` describes the input. The per-sample read-count column `post_trim` describes the output stage (read count **after** the cutadapt step has run, regardless of whether anything was actually clipped).

### J. Per-sample UMI / kit overrides (mixed-UMI batches)

When samples in the same batch have different UMI lengths or positions (e.g. one library preps with the NEBNext UMI 5' kit, another with QIAseq miRNA 3' UMI, a third already-trimmed from SRA), put a per-sample list under `align.samples:`:

```yaml
align:
  kit_preset: auto                 # global default for samples not listed below
  fastq: input_data/seq            # FASTQ directory
  contam_index: input_data/indexes/rrna_contam
  mt_index: input_data/indexes/mt_tx
  samples:
    - name: sampleA                # FASTQ basename (sampleA.fq.gz -> sampleA)
      kit_preset: illumina_truseq_umi
      umi_length: 8
      umi_position: 5p
    - name: sampleB
      kit_preset: qiaseq_mirna
      umi_length: 12
      umi_position: 3p
    - name: sampleC                # SRA-deposited, already adapter-clipped
      kit_preset: pretrimmed
      umi_length: 0
```

`mitoribopy all` materializes this block as a sidecar file at `<output>/align/sample_overrides.tsv` and passes it to `mitoribopy align` via `--sample-overrides`. Any field left unset in a sample entry falls through to the global default for that field, so a sample can override only its UMI length without restating the rest of the kit. Per-sample overrides are reported in `kit_resolution.tsv` with a `source` column starting with `per_sample_override:` so the provenance file tells the truth.

You can also write the TSV by hand and skip the YAML samples list:

```
sample      kit_preset            adapter  umi_length  umi_position  dedup_strategy
sampleA     illumina_truseq_umi              8           5p
sampleB     qiaseq_mirna                     12          3p
sampleC     pretrimmed                       0
```

then call `mitoribopy align --sample-overrides path/to/overrides.tsv …`.

---

## Tools

Standalone helper scripts under `mitoribopy.tools.*` — useful for shrinking inputs to a quick smoke-test size before a full pipeline run.

### Subsample a BED file

```bash
python -m mitoribopy.tools.subsample \
    --input  results/align/bed/sample.bed \
    --output /tmp/sample_subsampled.bed \
    --n 50000 \
    --seed 42
```

Reservoir-samples `--n` BED rows (Algorithm R, deterministic with `--seed`). Use the subsampled BED as `--directory` input to `mitoribopy rpf` for a fast end-to-end test.

### Subsample a FASTQ file

```bash
python -m mitoribopy.tools.subsample_fastq \
    --input  raw/sample.fastq.gz \
    --output /tmp/sample_subsampled.fastq.gz \
    --n 200000 \
    --seed 42
```

Reservoir-samples `--n` FASTQ records (4 lines each); gzip is auto-detected on either end via the `.gz` suffix. Plain `.fastq` works too. The first record header must start with `@` or the tool fails fast.

---

## Logs and provenance

- **`<output>/<stage>/mitoribopy.log`** — every stage writes a persistent log file alongside the same lines printed to the terminal.
- **Per-stage `run_settings.json`** — every stage writes its own settings JSON (resolved kit, dedup strategy, MAPQ threshold, reference checksum, tool versions, …).
- **`<output>/run_manifest.json`** — `mitoribopy all` composes per-stage settings into a top-level manifest. Schema version 1.0.0 carries: `schema_version`, `mitoribopy_version`, `git_commit` (best-effort, `null` outside a repo), `command` (the original argv), `config_source` + `config_source_sha256` (hash of the YAML the user wrote), `config_canonical` (the merged + auto-wired config that actually drove the run), `sample_sheet` + `sample_sheet_sha256` when applicable, a `stages: { align: {status, runtime_seconds}, rpf: {...}, rnaseq: {...} }` map (status is `completed` / `skipped` / `not_configured`; skipped stages carry a `reason`), the per-stage `run_settings.json` payloads, a flat `tools: {python, mitoribopy, cutadapt, ...}` map lifted from those payloads, and a `warnings` placeholder. The `reference_checksum` from the rpf stage is promoted to the top level so a downstream `rnaseq` run can verify it without drilling into the rpf section. Pin to a manifest layout by reading `schema_version` first.
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
