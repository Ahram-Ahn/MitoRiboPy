# Unified sample sheet — schema reference

The sample sheet is the **canonical project spine** for MitoRiboPy
runs. A single tab-delimited file declares every sample's identity,
assay, condition, FASTQ paths, and (optionally) per-sample kit / UMI /
strandedness overrides. When the sample sheet is set, the orchestrator
derives every per-stage input from it and **rejects** simultaneous use
of conflicting per-stage flags.

The loader is :func:`mitoribopy.sample_sheet.load_sample_sheet`. The
schema is enforced at load time; invalid sheets fail with a single
error message that lists every problem found.

## Required columns

| Column | Type | Description |
|---|---|---|
| `sample_id` | string | Unique sample identifier. Doubles as the FASTQ basename for the per-sample resolver in `mitoribopy align`. |
| `assay` | `ribo` \| `rna` | Lower-cased on load. Only these two values are accepted. |
| `condition` | string | Free-form condition label. Drives the DESeq2 contrast in the rnaseq stage. Cannot be empty. |
| `fastq_1` | path | Path to R1 FASTQ (or the SE FASTQ). Relative paths are resolved relative to the working directory. |

## Optional columns

| Column | Type | Description |
|---|---|---|
| `replicate` | string | Free-form replicate label (`1`, `2`, `A`, …). Informational only — pairing between assays is by `sample_id`, not by row order. |
| `fastq_2` | path | R2 FASTQ for paired-end reads. |
| `adapter` | string | Per-sample 3' adapter sequence override. Auto-detection runs by default; pin this when detection cannot identify the library. Mutually exclusive with `pretrimmed`. |
| `pretrimmed` | `true` \| `false` | When `true`, declares the FASTQ as already adapter-trimmed (cutadapt skips `-a` and only enforces length + quality). Mutually exclusive with `adapter`. |
| `umi_length` | int ≥ 0 | Per-sample UMI length override. **Recommended for publication** (otherwise the rnaseq from-FASTQ path may infer it from R1 entropy and emit a `UMI_INFERRED_NO_DECLARATION` warning). For `umi_position=both`, this MUST equal `umi_length_5p + umi_length_3p`. |
| `umi_position` | `5p` \| `3p` \| `both` | Per-sample UMI position. Required when `umi_length > 0` to avoid ambiguity. `both` is dual-end UMI (xGen Duplex, Twist) — supply `umi_length_5p` and `umi_length_3p`. |
| `umi_length_5p` | int ≥ 0 | Per-end 5' UMI length when `umi_position=both`. Required (must be > 0) for that mode; ignored otherwise. |
| `umi_length_3p` | int ≥ 0 | Per-end 3' UMI length when `umi_position=both`. Required (must be > 0) for that mode; ignored otherwise. |
| `strandedness` | `forward` \| `reverse` \| `unstranded` | Per-sample library strandedness. |
| `dedup_strategy` | `auto` \| `umi_coordinate` \| `skip` | Per-sample dedup choice. `auto` resolves to `umi_coordinate` (UMI-aware coordinate dedup via `umi_tools`) when `umi_length > 0`, else `skip`. The legacy aliases `umi-tools` and `umi_tools` are accepted and rewritten to canonical `umi_coordinate` at parse time so downstream `canonical_config.yaml` / `run_manifest.json` record one stable token. |
| `exclude` | `true` \| `false` \| (empty) | When `true`, the row is dropped from `active()` views without deleting the row. Lets a bad library be quarantined without altering the sheet structure. |
| `notes` | string | Free-form. Carried through to the manifest's `sample_sheet` entry but otherwise ignored. |

Empty cells, `NA`, `None`, `null`, and `-` are all read as `None`.
Lookups on string fields are case-sensitive. `assay` is canonicalised
to lowercase.

## Validation rules

The loader reports every error in one pass; you don't have to fix and
retry one row at a time.

* The header row must contain every required column. Unknown columns
  are rejected (catches typos like `samply_id`).
* `sample_id` must be unique across the whole sheet, including
  excluded rows.
* `assay` must be exactly `ribo` or `rna` (after lowercase).
* `umi_length` must parse as a non-negative integer.
* `strandedness` and `dedup_strategy`, when set, must use the
  literal values declared by `mitoribopy.align._types`.
* `exclude` must be a boolean-ish token.
* The data section must contain at least one row.

## Pairing semantics

RNA-seq and Ribo-seq rows are paired by `sample_id`, NOT by row
order. To use a different pairing, give the paired rows the same
`sample_id` even when they are different assays — the loader's
duplicate-id check catches the case where you accidentally duplicate
within an assay.

## Conflicts with per-stage inputs

When the orchestrator's top-level `samples:` block is set, the
following per-stage keys must NOT also be set (the unified sheet
supersedes them):

| Stage | Forbidden when `samples:` is set |
|---|---|
| align | `align.fastq`, `align.fastq_dir`, `align.samples`, `align.sample_overrides` |
| rnaseq | `rnaseq.rna_fastq`, `rnaseq.ribo_fastq`, `rnaseq.condition_map` |

`mitoribopy validate-config` and `mitoribopy all` both reject these
combinations with a single error message naming the offending key.
The standalone `mitoribopy rnaseq --sample-sheet` subcommand does the
same for the per-flag equivalents.

## Example

```
sample_id	assay	condition	replicate	fastq_1	fastq_2	adapter	pretrimmed	umi_length	umi_position	strandedness	dedup_strategy	exclude	notes
WT_Ribo_1	ribo	WT	1	ribo/WT_Ribo_1.fq.gz		AGATCGGAAGAGCACACGTCTGAACTCCAGTCA	false	8	5p	forward	umi_coordinate	false	first replicate
WT_Ribo_2	ribo	WT	2	ribo/WT_Ribo_2.fq.gz		AGATCGGAAGAGCACACGTCTGAACTCCAGTCA	false	8	5p	forward	umi_coordinate	false
KO_Ribo_1	ribo	KO	1	ribo/KO_Ribo_1.fq.gz		AGATCGGAAGAGCACACGTCTGAACTCCAGTCA	false	8	5p	forward	umi_coordinate	false
KO_Ribo_2	ribo	KO	2	ribo/KO_Ribo_2.fq.gz		AGATCGGAAGAGCACACGTCTGAACTCCAGTCA	false	8	5p	forward	umi_coordinate	true	failed library
WT_RNA_1	rna	WT	1	rna/WT_R1.fq.gz	rna/WT_R2.fq.gz		true	0		forward	skip	false
KO_RNA_1	rna	KO	1	rna/KO_R1.fq.gz	rna/KO_R2.fq.gz		true	0		forward	skip	false
```

Note: starting in v0.7.1 the user-facing `kit_preset` column was removed
in favour of the explicit `adapter` and `pretrimmed` columns. Internal
kit names (`illumina_truseq`, `qiaseq_mirna`, …) are still reported in
`kit_resolution.tsv` (`detected_kit` and `applied_kit`) so reviewers can
see which adapter family the detector matched.

## Cross-references

* `docs/reference/cli.md` — full CLI flag reference.
* `docs/validation/umi_dedup_validation.md` — UMI source labels and the
  conservative inference policy.
* `mitoribopy.sample_sheet` module docstring — canonical spec with
  loader behaviour.
