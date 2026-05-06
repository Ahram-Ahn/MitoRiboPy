# Output schema reference

Machine-precise reference for every file MitoRiboPy writes under a
run root. The intent is that a downstream script — or a manuscript
reviewer — can consume any output without reading source code.

The single source of truth for the file inventory and version
versioning is the ``_KNOWN_OUTPUTS`` registry in
[`src/mitoribopy/io/outputs_index.py`](../../src/mitoribopy/io/outputs_index.py)
and the ``OUTPUT_SCHEMA_VERSIONS`` map in
[`src/mitoribopy/io/schema_versions.py`](../../src/mitoribopy/io/schema_versions.py).
[`tests/test_output_schema_doc_parity.py`](../../tests/test_output_schema_doc_parity.py)
fails the suite when any registry entry is missing a row in this doc,
so the table cannot silently drift behind the runtime.

---

## Conventions

* **Schema version.** Every TSV in the table below carries a leading
  ``# schema_version: X.Y`` comment line. Standard pandas / csv
  readers ignore ``#``-prefixed lines when the caller passes
  ``comment="#"``; the column header is the first non-comment line.
  MINOR bumps (``1.0`` → ``1.1``) only append columns at the end, so
  fixed-position readers stay forwards-compatible across MINOR.
  MAJOR bumps (``1.x`` → ``2.0``) rename / remove / reorder columns.
* **Counts vs density vs RPM vs log2.** The package is explicit about
  what each numeric column carries:
  * **raw count** — integer, never normalised.
  * **density** — count divided by the gene's window length, in reads
    per nt or reads per codon (the file's column comment says which).
  * **RPM** — reads per million library-size after the relevant
    filter; the per-file note states which filter denominator is used.
  * **log2** — log2-transformed; pseudocount where applicable is
    documented in [`docs/rnaseq_te.md`](../rnaseq_te.md) (default
    δ = 0.5 in TE / ΔTE).
* **Coordinate space.** Genomic coordinates in BED outputs are
  zero-based, half-open (BED6 convention). Transcript-relative
  coordinates inside annotation rows (``start_codon``, ``stop_codon``)
  are zero-based first-nt-of-trinucleotide. P-site / A-site offsets
  are documented per-file: P-site is the canonical read 5'-end-to-P
  shift; A-site is `P + 3`.
* **Missing values.** Empty cell, ``NA``, ``None``, ``-`` and
  ``null`` are equivalent on read; on write the package emits the
  empty string for unknown / not-applicable, ``NaN`` for genuine
  numerical NaN, and an explicit string note otherwise (see
  ``te.tsv``'s ``note`` column for examples).
* **Provenance triplet.** Every run writes
  ``run_manifest.json`` (overall), ``canonical_config.yaml`` (the
  fully resolved config), and ``outputs_index.tsv`` (one row per
  emitted output with stage, recommended use, and schema version).
  These three files together let a reviewer reproduce the run from
  the same `mitoribopy --version`.

---

## align stage outputs

### `align/read_counts.tsv` (schema 1.0)

Per-sample read funnel from raw FASTQ through trim, contaminant
filter, mt-transcriptome alignment, MAPQ filter, and dedup. Rows are
written in sample-processing order; the orchestrator sorts by sample
name before writing so the file is deterministic.

| Column | Type | Units | Meaning |
|---|---|---|---|
| `sample` | string | — | Sample id (matches the unified sample sheet's `sample_id`). |
| `total_reads` | int | reads | Raw FASTQ read count. |
| `post_trim` | int | reads | Reads remaining after the cutadapt step. For a `pretrimmed` library this only reflects the length / quality filter; the name does **not** mean "reads were adapter-trimmed". |
| `rrna_aligned` | int | reads | Reads removed by the contaminant subtract step (typically rRNA / NUMT). |
| `post_rrna_filter` | int | reads | Reads entering mt-transcriptome alignment. Invariant: `rrna_aligned + post_rrna_filter == post_trim`. |
| `mt_aligned` | int | reads | Reads aligned to the mt-transcriptome. Invariant: `mt_aligned + unaligned_to_mt == post_rrna_filter`. |
| `unaligned_to_mt` | int | reads | `post_rrna_filter` minus `mt_aligned`. |
| `mt_aligned_after_mapq` | int | reads | After MAPQ filter (default ≥ 10). Always `≤ mt_aligned`. |
| `mt_aligned_after_dedup` | int | reads | After dedup (UMI-coordinate or skip). Always `≤ mt_aligned_after_mapq`. |

**Producer:** `mitoribopy align`. **Bumped by:** never since 1.0.
**Consumed by:** `summary_qc.tsv`, `SUMMARY.md`, every plotting code
path that needs library-depth normalisation.

### `align/kit_resolution.tsv` (schema 1.2)

Per-sample resolution of library kit, adapter sequence, UMI position
and length, and dedup strategy. Carries the detector's match-rate
columns so a reviewer can spot a sample that auto-resolved by a tiny
margin. The resolved-kit name is the *output* of detection — there is
no user-facing `--kit-preset` since v0.7.1; users supply `--adapter`
or `--pretrimmed` (or set `adapter:` / `pretrimmed:` in the
sample sheet).

Key columns: `sample`, `input_fastq`, `applied_kit`, `adapter`,
`umi_length`, `umi_position`, `dedup_strategy`, `library_strandedness`,
`expected_output_bed`, `detected_kit`, `best_adapter`,
`best_match_rate`, `second_best_match_rate`, `confidence_margin`,
`source` (`detected` / `inferred_pretrimmed` / `user_fallback` /
`per_sample_override:*`).

**Producer:** `mitoribopy align`. **Recent additions (1.1 → 1.2):**
`best_adapter`, `best_match_rate`, `second_best_match_rate`,
`confidence_margin`.

### `align/run_settings.json`

Resolved align-stage configuration plus the resolved external-tool
versions (`cutadapt`, `bowtie2`, `samtools`, `umi_tools`). Two
top-level sections matter for downstream reproducibility:
``args_resolved`` (flag values that actually drove the run, after
sample-sheet expansion + per-sample overrides) and ``tool_versions``
(SHA-stamped where the tool exposes one).

### `align/umi_qc.tsv` (schema 1.0)

Per-sample UMI / dedup audit.

| Column | Type | Meaning |
|---|---|---|
| `sample_id` | string | |
| `umi_present` | `true` / `false` | Whether a non-empty UMI was extracted. |
| `umi_length` | int | UMI length actually used (post-resolution). |
| `umi_position` | `5p` / `3p` / `both` / empty | Position the UMI was extracted from. |
| `n_reads_pre_dedup` | int | Reads entering dedup. |
| `n_reads_post_dedup` | int | Reads remaining after dedup. |
| `duplicate_fraction` | float in [0, 1] | `1 - n_reads_post_dedup / n_reads_pre_dedup`. |
| `dedup_strategy` | `auto` / `umi_coordinate` / `skip` | The resolved strategy (canonical token). |
| `dedup_method` | string | The actual `umi_tools dedup --method` value when `dedup_strategy=umi_coordinate`; empty on `skip`. |
| `warning_code` | string | One of `no_umi`, `umi_short`, `directional_on_short_umi`, `high_dup_fraction`, `umi_dedup_skipped_with_umi_present`, `dedup_skipped_no_umi`, or empty on a clean sample. |

**Producer:** `mitoribopy align`. The codes match the
`warnings.tsv` / `warning_codes.md` registry.

### `align/bed/`

Per-sample BED6 of mt-transcriptome alignments after dedup and MAPQ
filter. One file per sample, named `<sample>.bed`. BED6 columns:
`chrom`, `start` (0-based, inclusive), `end` (0-based, exclusive),
`name` (read id), `score` (MAPQ), `strand` (`+` / `-`).
Direct input to `mitoribopy rpf`.

---

## rpf stage outputs

### `rpf/rpf_counts.tsv` (schema 1.0)

Per `(sample, gene)` RPF count after the read-length filter and offset
assignment. Long format.

| Column | Type | Meaning |
|---|---|---|
| `sample` | string | Sample id. |
| `gene` | string | Transcript / gene name (as reported in the annotation CSV). |
| `count` | int | Reads assigned to that gene under the resolved P-site offsets. |

**Producer:** `mitoribopy rpf` (default downstream consumer of
`align/bed/`). **Note:** counts are post-offset, post-length-filter
— `mt_aligned_after_dedup` from `read_counts.tsv` is an upper bound.

### `rpf/rpf_counts.metadata.json`

Sidecar JSON recording the read-length window, offset selection mode
(`per_sample` / `combined`), the per-(sample, length) offsets actually
applied, the codon-edge masking nt counts, and the annotation CSV
checksum. Diff this against the rnaseq stage's reference checksum
(``rnaseq/run_settings.json``) to confirm the two stages used the
same reference.

### `rpf/offset_diagnostics/`

Subdirectory carrying the offset selection diagnostic plots
(per-sample P-site offsets per read length, combined / per-sample
overlays, masked-region annotations). The CSV table at
`rpf/offset_diagnostics/csv/offset_start.csv` is the per-(sample,
length) enrichment table the offset selector consumed; per-sample
copies live under `csv/per_sample_offset/<sample>/`.

### `rpf/offset_applied.csv`

The final per-sample, per-length offset table actually applied
downstream. One row per (sample, read_length, offset_type). Columns:
`sample`, `read_length`, `offset_type` (`5p` / `3p`), `offset`,
`source` (`combined` / `per_sample` / `manual`).

### `rpf/run_settings.json`

Resolved rpf-stage configuration (strain, fasta, annotation CSV
path, codon table, RPF window, offset mode + bounds, codon-edge
masking, periodicity knobs) plus the resolved tool versions and the
reference FASTA + annotation SHA256 checksums. The checksums are the
input to the rnaseq stage's reference-consistency gate.

### `rpf/qc/fourier_spectrum_combined.tsv` (schema 1.0)

One row per `(sample, read_length, gene_set, region, period_nt)`.
Full column table and the biological reading guide live in the
dedicated [periodicity reference](periodicity.md#output-schema). Key
columns: `period_nt` (float, grid 2.0..10.0 step 0.05) and
`amplitude` (float in [0, 1]; 1.0 = pure sinusoid at that period,
0.0 = no projection).

### `rpf/qc/fourier_period3_score_combined.tsv` (schema 1.1)

One row per `(sample, read_length, gene_set, region)`. The headline
QC artefact for periodicity. Carries the global and local 3-nt
spectral ratios, bootstrap CIs (1.1+), circular-shift permutation
p-values (1.1+), and the `snr_call` tier. See the
[periodicity reference's "Output schema" section](periodicity.md#output-schema)
for the full column-by-column breakdown plus the
[when-not-to-overinterpret caveats](periodicity.md#when-not-to-overinterpret).

### `rpf/qc/periodicity.metadata.json`

Sidecar recording every knob that drove the Fourier QC: window width,
codon-skip nt counts, gene_set definitions, DFT period grid, and
(1.1+) the bootstrap / permutation parameters (`n_bootstrap`,
`n_permutations`, `ci_alpha`, `random_seed`, `ci_method`,
`null_method`). Reproducing the CI / p-value bounds from this
sidecar alone is the contract.

### `rpf/qc/strand_sanity.tsv`

Per-sample minus-strand fraction averaged across mt-transcripts.
Columns: `sample`, `n_total`, `n_minus`, `minus_fraction`. The mt-transcriptome
FASTA is sense-oriented, so a non-zero `minus_fraction` signals a
wrong `--library-strandedness` upstream (most often `forward` set when
the library is `reverse`, or vice versa).

### `rpf/qc/strand_sanity_per_transcript.tsv` (schema 1.0)

Per-(sample, transcript) minus-strand fraction. Catches a localised
antisense bleed-through on one mt-mRNA — typically an unmasked NUMT or
a misannotated reference contig — that the all-transcripts mean in
`strand_sanity.tsv` would hide. Columns: `sample`, `transcript`,
`n_total`, `n_minus`, `minus_fraction`.

### `rpf/codon_correlation/{p_site,a_site}/codon_correlation_metrics.tsv`

Per-codon, per-(sample, transcript) correlation metrics with
normalised codon density, log2 fold-change vs the base sample, robust
residual, and label decisions (`auto_label` / `forced_label`). The
units of `density` are RPM-normalised reads-per-codon when
`--cor-metric log2_density_rpm` (the default for publication) is set.
The companion sidecars at
`rpf/codon_correlation/{p_site,a_site}/codon_correlation.metadata.json`
record which metric, regression, support filter, and label policy
produced each table — the P-site and A-site sidecars are emitted
independently so the metadata can diverge if the two sites were
configured separately.

---

## rnaseq stage outputs

### `rnaseq/te.tsv` (schema 2.0)

One row per (sample, gene) translation efficiency.

| Column | Type | Meaning |
|---|---|---|
| `sample_id` | string | Sample id. |
| `condition` | string | Condition label (matches the unified sheet's `condition`). |
| `assay` | `ribo` / `rna` | Which side of the pair the row describes. |
| `gene` | string | Gene name in the convention selected by `gene_id_convention` (`hgnc` / `ensembl` / `refseq` / `bare`). |
| `rpf_count` | int | Raw RPF count for this (sample, gene). |
| `rna_abundance` | float | mRNA abundance (DESeq2 baseMean for `de_table` mode; pyDESeq2-derived in `from_fastq`). |
| `te` | float | `(rpf_count + δ) / (rna_abundance + δ)`, δ = 0.5. |
| `log2_te` | float | `log2(te)`. NaN where `te` is 0 / undefined. |
| `note` | string | Free-form, e.g. `missing_from_de_table` when the gene is in the RPF counts but absent from the DE table. |

**Producer:** `mitoribopy rnaseq`. See
[`docs/rnaseq_te.md`](../rnaseq_te.md) for the publication-grade
boundary between the `de_table` and `from_fastq` modes.

### `rnaseq/delta_te.tsv` (schema 2.0)

One row per gene; the contrast `compare_condition − base_condition`.

| Column | Type | Meaning |
|---|---|---|
| `gene` | string | Gene name. |
| `base_condition`, `compare_condition` | string | Names of the two conditions in the contrast. |
| `mrna_log2fc` | float | DESeq2 log2 fold-change on the mRNA side. |
| `rpf_log2fc` | float | log2 fold-change on the RPF side. |
| `delta_te_log2` | float | `rpf_log2fc - mrna_log2fc`. |
| `padj_mrna` | float | DE table's adjusted p-value for the mRNA contrast. |
| `padj_rpf` | float | Adjusted p-value for the RPF contrast (set in `from_fastq` mode; NaN in `de_table` mode). |
| `padj_delta_te` | float | Adjusted p-value for the ΔTE contrast (when computable; NaN when method is exploratory). |
| `method` | string | The label written by the rnaseq pipeline — see [`src/mitoribopy/rnaseq/te.py`](../../src/mitoribopy/rnaseq/te.py) for the canonical list (`deseq_external`, `pydeseq2_in_tree`, etc.). |
| `note` | string | Free-form caveats (e.g. `pseudo_replicate_mode`, `mt_only_subset`). |

### `rnaseq/rna_counts.tsv` (schema 1.0)

`from_fastq` mode only. One row per (sample, gene) RNA count after
mt-transcriptome alignment. Columns: `sample`, `gene`, `count`.
Mirror of `rpf/rpf_counts.tsv` on the RNA side.

### `rnaseq/rpf_counts_matrix.tsv` (schema 1.0)

`from_fastq` mode only, and only emitted when the rnaseq stage is
not reusing an upstream `rpf_counts.tsv` (i.e. when the rnaseq
subcommand re-aligned the Ribo FASTQs itself). Wide matrix: one row
per gene, one column per sample, integer counts. Companion to the
long-form `rnaseq/rpf_counts.tsv`. The wide form exists so a
downstream R / pandas pivot is unnecessary; the long form is the
canonical input to the in-tree pyDESeq2 fit.

### `rnaseq/EXPLORATORY.md`

Sidecar written when the rnaseq run was launched in pseudo-replicate
mode or in a path that does not back a publication claim. Lists the
columns / claims that must NOT be cited (typically `padj_*`,
"significant gene" markers). **Absence of this file means
publication-grade.**

### `rnaseq/run_settings.json`

Resolved rnaseq-stage configuration, the resolved DE table format
and column map, the reference checksums consumed by the consistency
gate, and the resolved `rnaseq_mode`.

---

## Orchestrator (run-root) outputs

These files sit at the run root, not under a stage directory.

### `run_manifest.json`

The top-level provenance file. Its JSON schema is bundled at
[`src/mitoribopy/data/run_manifest.schema.json`](../../src/mitoribopy/data/run_manifest.schema.json)
and is exercised by [`tests/test_run_manifest_schema.py`](../../tests/test_run_manifest_schema.py).
Top-level keys:

* `schema_version` — the manifest's own schema version (string,
  `MAJOR.MINOR.PATCH`).
* `mitoribopy_version`, `git_commit` — installed package version
  and (when run from a checkout) the commit SHA.
* `created_at` — ISO-8601 timestamp.
* `command_line` — the argv that started the run.
* `config_source`, `canonical_config_path` — input config path and
  the resolved-config sidecar at `canonical_config.yaml`.
* `sample_sheet` — `{path, sha256, n_rows, n_active}` block.
* `tool_versions` — map of external tool name to resolved version.
* `output_schemas` — mirror of `OUTPUT_SCHEMA_VERSIONS` so a
  consumer can verify the TSVs they're about to read are the
  versions they expect.
* `warnings` — every `warnings.tsv` row, embedded for one-stop
  parsing.
* `outputs` — every `outputs_index.tsv` row, ditto.
* `resource_plan` — mirror of `resource_plan.json`.

### `resource_plan.json`

Resolved execution plan written before any stage runs (so it survives
a mid-pipeline crash). Records `threads`, `parallel_samples`,
`single_sample_mode`, `min_threads_per_sample`,
`estimated_memory_per_sample_gb`, `scheduler` (e.g. `auto`, `lsf`,
`slurm`), and a `reason` string explaining how the plan was derived
from the user input + machine probes. Mirror of the `resource_plan`
block in `run_manifest.json`.

### `canonical_config.yaml`

The fully-resolved config that drove the run, after legacy-key
canonicalisation, sample-sheet expansion (`samples:` block →
per-stage `align.fastq` + `align.sample_overrides` + (when rnaseq is
explicit) `rnaseq.sample_sheet`), execution-block cascade, and
rnaseq-mode resolution. Diff this against the input config to see
exactly what the orchestrator decided. Loadable as input to a future
`mitoribopy all --config canonical_config.yaml` for a re-run with the
same effective settings.

### `warnings.tsv` (schema 2.0)

Structured warnings collected across all stages. Columns:
`stage`, `sample_id`, `severity`, `code`, `message`, `suggested_action`.
Codes are stable; the registry plus per-code remediation lives in
[`docs/reference/warning_codes.md`](warning_codes.md).
``severity ∈ {info, warn, error}``; under `mitoribopy all --strict`
any `error` row promotes the run to a non-zero exit.

### `summary_qc.tsv` (schema 1.0)

Per-sample QC roll-up. One row per active sample.

Columns: `sample_id`, `assay`, `condition`, `kit_applied`,
`adapter_match_rate`, `post_trim_reads`, `contam_fraction`,
`mt_mrna_fraction`, `mapq_kept_fraction`, `dedup_removed_fraction`,
`umi_source`, `dominant_rpf_length`, `frame0_fraction`,
`offset_mode`, `offset_confidence`, `qc_status` (`pass` / `warn`),
`qc_notes` (comma-joined reasons for any warn).

The thresholds that drive `qc_status` live in
[`src/mitoribopy/summary/qc.py`](../../src/mitoribopy/summary/qc.py)
(`QC_THRESHOLDS`); they are conservative defaults pinned by the
unit-test suite.

### `SUMMARY.md`

Human-readable run summary. Links to every major output, includes
the periodicity statistical-confidence table, the warnings count by
severity, and a links-to-the-figures block. Generated from
`summary_qc.tsv` + `warnings.tsv` + the figure registry.

### `outputs_index.tsv` (schema 1.0)

One row per emitted output file. Columns: `output_type`, `stage`,
`path`, `description`, `recommended_for` (`downstream-scripting` /
`reviewer-spot-check`), `schema_version`. The
``schema_version`` cell is empty for outputs without a registered
schema (sidecar JSONs, plot directories, etc.).

The registry that drives the file is
[`src/mitoribopy/io/outputs_index.py`](../../src/mitoribopy/io/outputs_index.py).
This output_schema doc is the human mirror of the same registry.

### `figure_qc.tsv` (schema 1.0)

Per-plot mechanical QC produced by `mitoribopy validate-figures`.
Columns: `plot_path`, `stage`, `plot_type`, `status` (`pass` / `warn` /
`fail`), `n_points_expected`, `n_points_drawn`, `n_labels`,
`label_overlap_count`, `label_outside_axes_count`, `legend_overlap`,
`stat_box_overlap`, `clipped_text`, `min_font_size`,
`svg_text_editable`, `png_dpi`, `has_source_data`, `warnings`.
Missing metadata sidecars are not warnings by themselves; sidecar-backed
columns remain empty when no sidecar exists. Under
`validate-figures --strict`, any `warn` row is
promoted to `fail` and the subcommand exits non-zero.

### `progress.jsonl`

Newline-delimited JSON event stream, one event per line, emitted by
the progress system during the run. Each line carries `event_type`
(e.g. `run_start`, `stage_start`, `stage_end`, `progress`,
`run_end`), a UTC `timestamp`, and event-specific fields. Useful for
HPC schedulers that want to surface live progress without parsing
stage logs.

---

## See also

* [`docs/reference/sample_sheet_schema.md`](sample_sheet_schema.md)
  — the input-side counterpart to this output catalogue.
* [`docs/reference/periodicity.md`](periodicity.md) — full reference
  for the periodicity QC bundle (windows, DFT method, statistical
  hardening, when-not-to-overinterpret).
* [`docs/reference/warning_codes.md`](warning_codes.md) — every
  `warnings.tsv` code with severity, stage, summary, and remediation.
* [`docs/rnaseq_te.md`](../rnaseq_te.md) — the publication-boundary
  reference for `te.tsv` / `delta_te.tsv` and the strict-mode gates
  that protect them.
* [`docs/reference/cli.md`](cli.md) — every flag in every
  subcommand, regenerated from the live argparse parsers.
