# Inputs you need to prepare

This page answers "what do I need on disk before I run anything?".

* Read [What each module requires](#what-each-module-requires) first
  to scope your prep work.
* See [Sample sheet](#sample-sheet-unified-per-project-tsv) for the
  recommended single-source-of-truth file format.
* See [Input files](#input-files-file-by-file-reference) for
  per-file specifics.

For column-by-column schema details on every output the pipeline
writes (so you can see what the inputs map to), see
[`docs/reference/output_schema.md`](reference/output_schema.md).

---

## What each module requires

Each subcommand in the pipeline consumes its own slice of the inputs
below. **Required** rows must be present or the run aborts.
**Optional** rows have sensible defaults but commonly need
overriding.

### `mitoribopy align` — RNase-trimmed FASTQs → BED + read counts

| Kind | Item | Required? | Notes |
|---|---|---|---|
| File | Ribo-seq FASTQs (`*.fq[.gz]` / `*.fastq[.gz]`) | **required** | Provide via `--fastq <files>` (repeatable) or `--fastq-dir <dir>`. |
| File | mt-transcriptome bowtie2 index prefix | **required** | Built once with `bowtie2-build`. Pass via `--mt-index <prefix>`. |
| File | rRNA / tRNA contaminant bowtie2 index prefix | **required** | Used to subtract contaminants. Pass via `--contam-index <prefix>`. |
| File | **Sample sheet** (`samples.tsv`) | optional | One row per Ribo-seq FASTQ; documents per-sample kit / UMI / strandedness. Strongly recommended for mixed-kit batches. See [Sample sheet](#sample-sheet-unified-per-project-tsv). |
| Option | `--adapter <SEQ>` | optional | 3' adapter sequence. Auto-detection runs by default; pin this when detection cannot identify the library. Override per sample in the sample sheet. |
| Option | `--pretrimmed` | optional (`false` default) | Declare already-trimmed FASTQs (cutadapt skips `-a`). Mutually exclusive with `--adapter`. |
| Option | `--library-strandedness {forward,reverse,unstranded}` | optional (`forward` default) | dUTP-stranded libraries should set `reverse`. |

### `mitoribopy rpf` — BED + reference FASTA → P-site / A-site analysis

| Kind | Item | Required? | Notes |
|---|---|---|---|
| File | Ribo-seq BEDs (or BAMs) | **required** | When run after `align`, auto-wired from `<align>/bed/`. Standalone runs use `--directory <dir>`. |
| File | mt-transcriptome FASTA | **required** | One record per mt-mRNA. Pass via `--fasta <path>`. |
| File | Annotation CSV | optional (built-in for h.sapiens / s.cerevisiae) | Required only for custom organisms. See [`docs/custom_organisms.md`](custom_organisms.md) for the full schema. |
| File | Codon-table JSON | optional (NCBI codes 1–33 bundled) | Pick one via `--codon_table_name`; supply your own only for non-NCBI codes. |
| File | Read-count table (`read_counts.tsv`) | optional (RPM normalization) | Auto-wired from `align/read_counts.tsv` when run via `mitoribopy all`. |
| Option | `-rpf <min> <max>` | optional (footprint-class default) | Read-length window. Defaults: monosome 28–34 (human), 37–41 (yeast). |
| Option | `--strain {h.sapiens, s.cerevisiae, custom}` | **required** | Drives the built-in annotation + codon table; `custom` requires the two files above. |

### `mitoribopy rnaseq` — TE / ΔTE from RNA-seq + Ribo-seq (two flows)

For the publication-grade boundary between the two flows (and the
`--strict` gates that protect it), see
[`docs/rnaseq_te.md`](rnaseq_te.md).

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

### `mitoribopy all` — end-to-end orchestrator

`mitoribopy all` runs the three stages above with one shared YAML
config. Required inputs are the **union** of every active stage's
required inputs. The config has these top-level sections:

| Section | When to include | Notes |
|---|---|---|
| `samples:` | recommended | Top-level `samples: { table: samples.tsv }` is the canonical declaration of every input FASTQ + per-sample metadata. Auto-wires both `align` and `rnaseq`. See [Sample sheet](#sample-sheet-unified-per-project-tsv). |
| `align:` | required for align | Indexes, kit / strandedness, dedup. Omit `align.fastq` when `samples:` is set — the sheet supplies it. |
| `rpf:` | required for rpf | Strain, RPF window, FASTA, offset bounds. Auto-wires `--directory`, `--read_counts_file`, `--fasta`. |
| `rnaseq:` | required for rnaseq | Either flow's keys (the two are mutually exclusive). The sheet auto-wires `rnaseq.sample_sheet` so you do not repeat it. |

Use `mitoribopy all --print-config-template > pipeline_config.yaml`
to drop a fully-commented starter into your project.

---

## Sample sheet (unified per-project TSV)

A single TSV declares every sample once, replacing the old pair of
stage-specific tables (`--sample-overrides` and `--condition-map`).
It is **the recommended way** to declare inputs for any non-trivial
project: pairings between Ribo-seq and RNA-seq are by `sample_id`
(never by index), per-sample kit / UMI overrides for mixed batches
live in the same file, and an `exclude` column lets you drop a bad
library without deleting rows.

**Required columns:** `sample_id`, `assay` (`ribo` or `rna`),
`condition`, `fastq_1`.

**Optional columns:** `replicate`, `fastq_2`, `adapter`,
`pretrimmed`, `umi_length`, `umi_position` (`5p` / `3p` / `both`),
`umi_length_5p`, `umi_length_3p` (per-end lengths for
`umi_position=both`), `strandedness`, `dedup_strategy`, `exclude`
(`true`/`false`/blank), `notes`.

The full schema reference (cell-value validation rules, conflicts
with per-stage inputs, the canonical token map for
`dedup_strategy`) lives in
[`docs/reference/sample_sheet_schema.md`](reference/sample_sheet_schema.md).

Empty cells (`""`, `NA`, `None`, `-`, `null`) read as "use the
default". Lines starting with `#` and blank lines are ignored.
Validation is strict: a single load pass reports every row error so
you can fix the sheet without iterate-and-retry.

```tsv
sample_id	assay	condition	replicate	fastq_1	fastq_2	adapter	pretrimmed	umi_length	umi_position	strandedness	exclude	notes
WT_Ribo_1	ribo	WT	1	fastq/WT_Ribo_1.fq.gz		AGATCGGAAGAGCACACGTCTGAACTCCAGTCA	false	8	5p	forward	false	
WT_Ribo_2	ribo	WT	2	fastq/WT_Ribo_2.fq.gz		AGATCGGAAGAGCACACGTCTGAACTCCAGTCA	false	8	5p	forward	false	
KO_Ribo_1	ribo	KO	1	fastq/KO_Ribo_1.fq.gz		AGATCGGAAGAGCACACGTCTGAACTCCAGTCA	false	8	5p	forward	false	
KO_Ribo_2	ribo	KO	2	fastq/KO_Ribo_2.fq.gz		AGATCGGAAGAGCACACGTCTGAACTCCAGTCA	false	8	5p	forward	true	contaminated lane
WT_RNA_1	rna	WT	1	fastq/WT_RNA_1_R1.fq.gz	fastq/WT_RNA_1_R2.fq.gz		true	0		forward	false	
KO_RNA_1	rna	KO	1	fastq/KO_RNA_1_R1.fq.gz	fastq/KO_RNA_1_R2.fq.gz		true	0		forward	false	
```

> Auto-detection (the default) makes the `adapter` column optional in
> the common case. It still travels in `kit_resolution.tsv` after the
> run as `detected_kit` / `applied_kit` so reviewers can see which
> adapter family the detector matched.

How it threads through the pipeline:

| Stage | What the sheet supplies |
|---|---|
| `mitoribopy align` | The Ribo-seq FASTQ list (`assay='ribo'` rows, `exclude=false`) and a materialised `sample_overrides.tsv` carrying the per-sample kit / UMI / dedup columns. |
| `mitoribopy rnaseq` | The RNA-seq FASTQ list (`assay='rna'`), the Ribo-seq FASTQ list (re-counted from raw FASTQ in this stage's own align machinery), and the `sample_id → condition` mapping that drives the pyDESeq2 contrast. |
| `mitoribopy all` | Both of the above, auto-wired. Top-level `samples: { table: samples.tsv }` (or shorthand `samples: samples.tsv`) is the only place you need to declare inputs. |

Mutual-exclusion rules: when the sheet is set, declaring
`align.fastq` / `align.fastq_dir` / `align.samples` /
`align.sample_overrides` / `rnaseq.rna_fastq` / `rnaseq.ribo_fastq` /
`rnaseq.condition_map` alongside it is an error — pick one input
style per stage.

---

## Input files (file-by-file reference)

### FASTQ (primary input to `align`)

Accepted file extensions: `*.fq`, `*.fq.gz`, `*.fastq`,
`*.fastq.gz`. Both gzipped and uncompressed are auto-detected.

Two ways to point the pipeline at your FASTQs:

1. **Directory** (recommended): pass a directory containing every
   input FASTQ.
   - CLI: `--fastq-dir input_data/`
   - YAML: `align.fastq: input_data/` (a single string is treated as
     a directory)
2. **Explicit list**: name each FASTQ.
   - CLI: `--fastq sampleA.fq.gz --fastq sampleB.fq.gz` (repeatable)
   - YAML: `align.fastq: [sampleA.fq.gz, sampleB.fq.gz]` (a list is
     treated as explicit paths)

Sample names are derived from the FASTQ filename with the extension
stripped (`WT_R1.fq.gz` → `WT_R1`). The same name flows through every
per-sample table (`read_counts.tsv`, `kit_resolution.tsv`, downstream
profile and codon usage subdirs).

### BED (input to `rpf` if you already have aligned BEDs)

Expected columns:

1. `chrom`
2. `start`
3. `end`

Additional BED columns are tolerated. Coordinates are 0-based,
end-exclusive intervals (standard BED).

When the `align` stage runs first (or you use `mitoribopy all`),
`mitoribopy rpf` consumes `<align>/bed/` automatically — you never
need to handle BED files by hand.

### BAM (alternative input to `rpf`)

`mitoribopy rpf --directory <dir>` accepts BAM files mixed with BED.
BAMs are auto-converted to BED6 under `<output>/bam_converted/` via
pysam. The `--bam_mapq` flag (default 10) filters BAM reads on MAPQ
before conversion to suppress NUMT cross-talk; set to 0 to disable.

### Reference FASTA

One FASTA record per mt-mRNA. Headers must match the `sequence_name`
column of the annotation CSV (or any of its `sequence_aliases`).
The built-in human and yeast annotations cover the canonical
mt-mRNAs and ship under `src/mitoribopy/data/`.

For total-genome FASTAs (rare in mt-Ribo-seq), use
`--annotation_file` to map FASTA records onto your own annotation
rows.

### Annotation CSV (custom organisms)

Built-in `h.sapiens` and `s.cerevisiae` ship complete annotation
tables and need nothing here. For any other organism, supply a
per-transcript CSV via `--annotation_file`. The full schema lives in
[`docs/custom_organisms.md`](custom_organisms.md).

### Codon-table JSON (custom organisms)

The 27 NCBI Genetic Codes are bundled. Pick one with
`--codon_table_name` (full picker in
[`docs/custom_organisms.md`](custom_organisms.md)). Supply your own
`--codon_tables_file` only when your organism's code is not in the
NCBI list.

### Read-count table (optional, for RPM normalization)

`.csv`, `.tsv`, and `.txt` accepted; delimiter is auto-detected.
Column matching is flexible and case-insensitive, with positional
fallback:

- column 1: sample name
- column 2: reference (used when `--rpm_norm_mode mt_mrna`)
- column 3: read count

When `mitoribopy all` runs `align` first, the read-count table is
auto-wired from `<run_root>/align/read_counts.tsv`.
