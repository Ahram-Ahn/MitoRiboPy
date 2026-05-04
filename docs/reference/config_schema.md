# Pipeline config schema (YAML / JSON / TOML)

The orchestrator (`mitoribopy all`) and each per-stage subcommand
accept a config file via `--config <path>`. The format is auto-
detected from the suffix: `.yaml` / `.yml` (PyYAML), `.json`, or
`.toml` (Python 3.11+ stdlib `tomllib`, or the `tomli` fallback).

This document is the canonical reference for the **orchestrator**
config, which has these top-level sections:

```yaml
samples: { table: samples.tsv }     # RECOMMENDED single source of truth

execution:
  ...

periodicity:
  ...

align:
  ...

rpf:
  ...

rnaseq:
  ...
```

For per-stage flat configs (no `align: / rpf: / rnaseq:` wrapper), see
the templates under `examples/templates/<stage>_config.example.yaml`.

---

## `samples:` (RECOMMENDED)

Two equivalent forms:

```yaml
samples:
  table: samples.tsv
```

```yaml
samples: samples.tsv
```

When set, the orchestrator:

* loads and validates the sheet via :func:`mitoribopy.sample_sheet.load_sample_sheet`;
* derives `align.fastq` (one entry per `assay='ribo'` row) and
  `align.sample_overrides` (a materialised TSV under
  `<output>/align/sample_overrides.tsv`) when any per-sample adapter,
  UMI, pretrimmed, strandedness, or dedup override is set;
* threads the sheet path through to `rnaseq.sample_sheet` when the
  rnaseq section is also set.

Conflicts with `align.fastq`, `align.fastq_dir`, `align.samples`,
`align.sample_overrides`, `rnaseq.rna_fastq`, `rnaseq.ribo_fastq`,
or `rnaseq.condition_map` are rejected at config load time. See
`sample_sheet_schema.md` for the column reference.

---

## `align:` section

Every key corresponds 1:1 to a `mitoribopy align` CLI flag (with
underscores instead of dashes). Run `mitoribopy align --help` (or
`mitoribopy all --show-stage-help align`) for the full flag list.

Selected highlights:

| Key | Meaning |
|---|---|
| `adapter` | 3' adapter sequence (string). Auto-detection runs by default; pin `adapter:` only when detection cannot identify the library. Mutually exclusive with `pretrimmed`. |
| `pretrimmed` | `true` / `false` (default `false`). When `true`, declares already-trimmed FASTQs (cutadapt skips `-a`). Mutually exclusive with `adapter`. |
| `adapter_detection` | `auto` / `off` / `strict` |
| `umi_length`, `umi_position` | per-run defaults; per-sample overrides go in the sample sheet. `umi_position` accepts `5p`, `3p`, or `both` |
| `umi_length_5p`, `umi_length_3p` | per-end UMI lengths for `umi_position: both` (e.g. xGen Duplex, Twist). Required when `umi_position=both`; `umi_length` MUST equal their sum and is auto-derived when omitted |
| `contam_index`, `mt_index` | bowtie2 index prefixes (sidecar `.1.bt2` must exist) |
| `min_length`, `max_length`, `quality`, `mapq` | length / quality / MAPQ filters |
| `dedup_strategy` | `auto` / `umi_coordinate` / `skip` (`umi-tools` and `umi_tools` are accepted as legacy aliases and rewritten) |
| `max_parallel_samples` | per-stage worker count (see HPC docs) |

Inputs (FASTQs) come from either `align.fastq:` (string = directory
shortcut, list = explicit paths) OR the top-level `samples:` sheet —
not both.

---

## `rpf:` section

Every key corresponds 1:1 to a `mitoribopy rpf` CLI flag. YAML keys
are always snake_case (`offset_type`, `min_5_offset`, …); the rpf CLI
flag is the same name but in hyphen-case (`--offset-type`,
`--min-5-offset`). The pre-v0.6.0 underscore CLI form
(`--offset_type`) is still accepted for one transition cycle but is
not surfaced in `--help` or generated docs. The orchestrator
serialises every section to hyphenated argv as of v0.6.0.

Selected highlights:

| Key | Meaning |
|---|---|
| `strain` | `h.sapiens` / `s.cerevisiae` / `custom` (legacy `h` / `y` accepted) |
| `footprint_class` | `short` / `monosome` / `disome` |
| `rpf` | `[min, max]` read-length window override |
| `fasta` | reference FASTA (mt-transcriptome) |
| `directory` | BED input dir (auto-wired to `<run>/align/bed/` by `mitoribopy all`) |
| `align` | `start` / `stop` (codon to anchor offsets at) |
| `offset_type` | `5` / `3` (which read end to measure offsets from) |
| `offset_site` | `p` / `a` (P-site or A-site coordinate space for the SELECTED OFFSETS table) |
| `offset_pick_reference` | `p_site` / `reported_site` (legacy `selected_site`) |
| `offset_mode` | `per_sample` (default) / `combined` |
| `analysis_sites` | `both` / `p` / `a` (which downstream outputs to produce) |
| `codon_density_window` | smooth codon density with ±1 nt window (legacy `merge_density`) |

---

## `execution:` section

The optional top-level `execution:` block defines the run-level
resource budget. The orchestrator cascades these values into stage
flags unless a stage pins its own value.

| Key | Meaning |
|---|---|
| `threads` | total CPU-thread budget; also set by `mitoribopy all --threads` when YAML does not pin it |
| `memory_gb` | memory budget used by the resource planner |
| `parallel_samples` | align-stage sample workers; `auto` lets the planner choose |
| `single_sample_mode` | force serial per-sample align work |
| `min_threads_per_sample` | floor used by the automatic parallel-sample calculation |
| `estimated_memory_per_sample_gb` | memory estimate used by the automatic parallel-sample calculation |
| `scheduler` | free-form scheduler label recorded in `resource_plan.json` |

## `periodicity:` section

The optional top-level `periodicity:` block tunes the Fourier
periodicity QC emitted by `mitoribopy rpf`. Explicit
`rpf.periodicity_*` keys win when both are set.

| Key | Meaning |
|---|---|
| `enabled` | `true` / `false`; disable only for debugging or very fast smoke runs |
| `fourier_window_nt` | Fourier metagene window size in nt; default 99 |
| `metagene_nt` | start/stop metagene plot window in nt; default 300 |
| `metagene_normalize` | `per_gene_unit_mean` (default) / `none` |
| `fourier_bootstrap_n` | bootstrap iterations for period-3 CI; default 200 |
| `fourier_permutations_n` | circular-shift permutations for the empirical null; default 200 |
| `fourier_ci_alpha` | two-sided CI alpha; default 0.10 |
| `fourier_random_seed` | RNG seed recorded in `periodicity.metadata.json`; default 42 |
| `no_fourier_stats` | `true` skips CI and permutation p-values |

### Standalone `mitoribopy periodicity` flags

The `rpf` stage emits the metagene Fourier QC bundle under
`<output>/rpf/qc/`. Built-in defaults: window 99 nt = 33 codons, drop
5 codons after AUG (initiation peak), drop 1 codon before stop
(termination peak), min mean coverage 0.1, min total counts 30. To
re-score periodicity from a saved site table without re-running
offset selection, use the standalone subcommand:

```bash
mitoribopy periodicity \
  --site-table runs/full/rpf/qc/site_table.tsv \
  --output     runs/full/rpf/qc/standalone_periodicity \
  --site p \
  --fourier-window-nt 99 \
  --drop-codons-after-start 5 \
  --drop-codons-before-stop 1 \
  --min-mean-coverage 0.1 \
  --min-total-counts 30
```

The Fourier QC bundle outputs (`fourier_spectrum_combined.tsv`,
`fourier_period3_score_combined.tsv`, plus the per-`(sample,
read_length)` figures under `fourier_spectrum/<sample>/`) are
**always emitted**. Three figures are written per `(sample,
read_length)`: `*_combined.png` (canonical mt-mRNAs aggregated into a
single trace per panel), `*_ATP86.png` (junction-bracketed
ATP8/ATP6 bicistronic analysis: top panel = ATP8 frame, bottom panel
= ATP6 frame, windows overlapping at the bicistronic junction at
nt 177-202), and `*_ND4L4.png` (junction-bracketed ND4L/ND4 analysis;
windows flank the 4-nt overlap junction).

Required columns in the input site table: `sample`, `gene`,
`transcript_id`, `read_length`, `site_type` (`p` or `a`), `site_pos`
(transcript-coordinate, 0-based), `cds_start`, `cds_end`. Optional
`count` column is honoured for weighted counting (defaults to 1 per
row when absent).

The older frame-fraction QC bundle (`qc_summary.tsv`,
`frame_counts_*.tsv`, `gene_periodicity.tsv`, the frame heatmaps,
`--phase-score`, and the threshold flags) is not emitted by the
current package; see [periodicity.md](periodicity.md) for the current
Fourier contract.

---

## `rnaseq:` section

Two mutually exclusive flows, selected by `rnaseq_mode`:

| `rnaseq_mode` | Required keys | Notes |
|---|---|---|
| `de_table` | `de_table`, `gene_id_convention`, one of `reference_gtf` / `reference_checksum`, plus `ribo_dir` (auto-wired in `mitoribopy all`) | **Publication-grade.** Consumes a pre-computed external DE table. SHA256 reference-consistency gate. |
| `from_fastq` | `rna_fastq` (or `sample_sheet`), `reference_fasta`, `gene_id_convention`, `condition_a` / `condition_b` | **Exploratory.** Runs cutadapt + bowtie2 + pyDESeq2 on the mt-mRNA subset. |
| `none` | (no inputs) | Stage section present but inert. |

If `rnaseq_mode` is omitted, the orchestrator infers it from supplied
inputs and emits a stderr banner when from-FASTQ is selected.

Selected keys:

| Key | Meaning |
|---|---|
| `rnaseq_mode` | `de_table` / `from_fastq` / `none` (alias: `mode:`) |
| `rna_fastq` | dir or list of RNA FASTQs (from-FASTQ flow) |
| `ribo_fastq` | optional; in `mitoribopy all` from-FASTQ defaults to **reusing** the rpf stage's `rpf_counts.tsv` instead of re-aligning |
| `recount_ribo_fastq` | `true` / `false` (default `false`); set `true` to force a second alignment pass over Ribo FASTQs in the rnaseq stage |
| `de_table` | path to external DE table (DESeq2 / Xtail / Anota2Seq; auto-detected) |
| `reference_gtf` / `reference_checksum` | exactly one is required in the de_table flow |
| `gene_id_convention` | `hgnc` / `ensembl` / `refseq` / `mt_prefixed` / `bare`. Required (no default) — silent zero-match runs are otherwise the failure mode. |
| `condition_a` / `condition_b` (or `base_sample` / `compare_sample`) | DESeq2 contrast |

---

## Auto-wiring inside `mitoribopy all`

The orchestrator sets a number of cross-stage defaults so the user
doesn't have to repeat paths. Visible via
`mitoribopy all --print-canonical-config --config ... --output ...`:

| Auto-wired key | Source |
|---|---|
| `align.output` | `<output>/align/` |
| `rpf.output` | `<output>/rpf/` |
| `rpf.directory` | `<output>/align/bed/` (when align is in the run) |
| `rpf.read_counts_file` | `<output>/align/read_counts.tsv` |
| `rnaseq.output` | `<output>/rnaseq/` |
| `rnaseq.ribo_dir` | `<output>/rpf/` (de_table flow) |
| `rnaseq.reference_fasta` | `rpf.fasta` (from-FASTQ flow, when not set explicitly) |
| `rnaseq.upstream_rpf_counts` | `<output>/rpf/rpf_counts.tsv` (from-FASTQ flow with rpf stage in the run, unless `recount_ribo_fastq: true`) |

---

## Validation

Use `mitoribopy validate-config <path>` to pre-flight a config:

* parses the YAML;
* canonicalises legacy keys (`merge_density` → `codon_density_window`,
  `strain: h` → `strain: h.sapiens`, etc.) and reports each rewrite;
* checks file paths;
* enforces sample-sheet vs per-stage conflict rules;
* resolves `rnaseq.mode` against supplied inputs.

`--strict` makes any legacy-key rewrite a hard failure (publication
mode); `--no-path-checks` skips on-disk existence (CI-side validation
on a different host).

## Migration

Use `mitoribopy migrate-config <old.yaml>` to rewrite legacy keys to
their canonical names. Outputs canonical YAML to stdout; the per-rewrite
change log goes to stderr. Pass `--in-place` to overwrite the file
(with a `.bak` backup).

## Provenance

Every `mitoribopy all` run writes the canonical (post-migration,
post-auto-wire) config under `run_manifest.json` → `config_canonical`,
along with the SHA256 of the source file (`config_source_sha256`). The
hash is what the `--resume` guard compares against on subsequent runs.
