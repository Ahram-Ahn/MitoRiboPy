# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Forward-looking entries only. Implemented changes are recorded under a
versioned release header. Long-running design notes live in
[`docs/developer/roadmap.md`](docs/developer/roadmap.md).

### Planned
- Extend `examples/smoke/` with richer numerical assertions once the
  synthetic reference and read simulator carry a stable biological
  signal. The fixture itself is shipped as an opt-in wiring test in
  v0.7.0+.

## [0.7.1] - 2026-05-02

### Changed (BREAKING)
- **Adapter-trimming UX simplified — kit-preset input removed.** The
  user-facing `--kit-preset` flag, `align.kit_preset:` config key, and
  `kit_preset` columns (sample sheet + sample-overrides TSV +
  `align.samples:` YAML block) are gone. Auto-detection of the 3'
  adapter is now the only path; pass `--adapter <SEQ>` (or
  `adapter: <SEQ>` in YAML) when detection cannot identify the library,
  and the new `--pretrimmed` flag (or `pretrimmed: true` in YAML)
  declares already-trimmed FASTQs explicitly. Detection still
  identifies named adapter families internally and reports the kit name
  in `kit_resolution.tsv` (`detected_kit` and `applied_kit` columns) so
  provenance is preserved — the asymmetry is intentional: kit names are
  output for reporting, not input from users.

  **Migration:**
  - `--kit-preset illumina_smallrna` → `--adapter TGGAATTCTCGGGTGCCAAGG`
  - `--kit-preset illumina_truseq` → `--adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`
  - `--kit-preset illumina_truseq_umi` → same `--adapter` plus `--umi-length 8 --umi-position 5p`
  - `--kit-preset qiaseq_mirna` → `--adapter AACTGTAGGCACCATCAAT --umi-length 12 --umi-position 3p`
  - `--kit-preset pretrimmed` → `--pretrimmed`
  - `--kit-preset auto` → omit (this is the default)
  - YAML / sample-sheet `kit_preset:` entries: replace with `adapter:` /
    `pretrimmed:` per the same mapping. `mitoribopy migrate-config`
    drops the key with a logged warning so you can review case-by-case.
  - The `--adapter-detection off` mode now requires `--adapter` or
    `--pretrimmed` (previously required `--kit-preset`).
  - `KIT_PRESET_ALIASES` and `resolve_kit_alias()` are removed (legacy
    vendor names like `truseq_smallrna` / `nebnext_smallrna` no longer
    needed canonicalisation since the user input vector is gone).

  Internal sentinels (`pretrimmed`, `custom`) and the named adapter-
  family presets (`illumina_smallrna`, `illumina_truseq`,
  `illumina_truseq_umi`, `qiaseq_mirna`) remain in
  `mitoribopy.align._types.KIT_PRESETS` because the detector still
  emits them as the `applied_kit` / `detected_kit` labels.

- **`SampleOverride` dataclass:** `kit_preset` field replaced with
  `pretrimmed: bool | None`. `mitoribopy.align._types.SampleOverride`
  imports that referenced `kit_preset` need updating.

- **`SampleRow` dataclass:** `kit_preset` field replaced with
  `pretrimmed: bool | None`. Sample sheets carrying a `kit_preset`
  column raise `SampleSheetError` with a migration message.

- **`run_settings.json`:** `kit_preset_default` field removed; replaced
  with `pretrimmed_default: bool` alongside the existing
  `adapter_default`.

## [0.7.0] - 2026-05-02

### Publication-readiness — consolidated release

This is the publication-grade snapshot of MitoRiboPy. It supersedes
the in-flight v0.8.0 / v0.9.0 / v0.9.1 entries that previously sat
between v0.6.2 and HEAD: those bumps were the cumulative iteration
notes for the periodicity rewrite, the statistical hardening, and the
follow-through wiring. Re-issued here as a single v0.7.0 release
because there is no public artefact for v0.8.x / v0.9.x and a coherent
SemVer trajectory v0.6.2 → v0.7.0 → v1.0.0 (API freeze) is preferable
to four-version churn for users who pin against PyPI.

The biological invariant is preserved: TACO1-KO polyproline stalling
in `MT-CO1` survives the metagene-Fourier rewrite, the new metagene
normalisation default, the dedup defaults, and the offset-selection
machinery (see [`docs/validation/taco1_ko_regression.md`](docs/validation/taco1_ko_regression.md)).

### Changed (BREAKING)
- **Periodicity QC bundle replaced with an aggregate-then-DFT metagene
  Fourier analysis.** The frame-fraction QC bundle (`qc_summary.tsv` /
  `qc_summary.md`, `frame_counts_by_sample_length.tsv`,
  `frame_counts_by_gene.tsv`, `gene_periodicity.tsv`, `phase_score`)
  and the four named frame-heatmap plots
  (`frame_fraction_heatmap.svg`, `frame_by_length_heatmap.{png,svg}`,
  `read_length_periodicity_barplot.svg`,
  `gene_phase_score_dotplot.svg`) are removed. The Fourier-spectrum
  bundle is rewritten end-to-end: per-gene normalise (divide by window
  mean) → element-wise mean across qualifying transcripts →
  mean-centre + linear detrend → Hann window → direct DFT at exactly
  period 3.0 (not bin-snapped). Default window: **99 nt = 33 codons**
  (multiple of 3, no period-3 leakage); skip 5 codons after AUG
  (initiation peak) and 1 codon before stop (termination peak). The
  `utr3` region is dropped — human mt-mRNA 3' UTRs are too short for a
  meaningful 99-nt window. New outputs: `fourier_spectrum_combined.tsv`
  and `fourier_period3_score_combined.tsv` (one row per `[sample,
  read_length, gene_set, region]`; `gene_set ∈ {combined, ATP86,
  ND4L4}`; `region ∈ {orf_start, orf_stop}`). Three figures per
  `(sample, read_length)`: `*_combined.png` (single trace per panel
  from canonical mt-mRNAs), `*_ATP86.png` (junction-bracketed
  ATP8/ATP6 bicistronic — top: ATP8 frame, bottom: ATP6 frame, windows
  overlapping at the bicistronic junction nt 177-202), `*_ND4L4.png`
  (junction-bracketed ND4L/ND4 — windows flank the 4-nt junction).
  Each panel shows the `spectral_ratio_3nt` and `snr_call` as
  in-figure annotations.

  **Migration**: scripts that read `qc_summary.tsv` for the per-sample
  QC verdict should switch to `fourier_period3_score_combined.tsv` and
  filter by `gene_set == "combined"`, `region == "orf_start"` (or
  `"orf_stop"`); the headline column is `spectral_ratio_3nt`. Scripts
  that read `frame_counts_by_sample_length.tsv` to pick a best read
  length should use `mitoribopy summarize`'s output (which now reads
  from the Fourier score table). The `--phase-score`,
  `--periodicity-phase-score`, `--periodicity-exclude-start-codons`,
  `--periodicity-exclude-stop-codons`, `--periodicity-min-reads-per-length`,
  `--good-frame-fraction`, `--warn-frame-fraction`,
  `--min-reads-per-length`, `--min-reads-per-gene`,
  `--exclude-start-codons`, `--exclude-stop-codons` CLI flags are
  removed. The standalone `mitoribopy periodicity` subcommand exposes
  `--fourier-window-nt`, `--drop-codons-after-start`,
  `--drop-codons-before-stop`, `--min-mean-coverage`,
  `--min-total-counts`, `--no-plots`, plus the new statistical-
  hardening flags listed below. The annotation `sequence_aliases`
  field is consulted by the chrom-to-transcript lookup so fused-FASTA
  datasets resolve to BOTH constituent transcripts (ATP8 + ATP6 /
  ND4L + ND4) for the bicistronic-figure panels.

- **Metagene aggregation defaults to per-gene unit-mean normalisation.**
  `compute_metagene` now divides each transcript's per-position density
  by its own mean before averaging, so a single high-expression
  transcript no longer dominates the metagene shape (the Fourier path
  was already per-gene-unit-mean; the metagene TSV / plot were not).
  The legacy depth-weighted sum is preserved behind
  `compute_metagene(..., normalize="none")`.
  `metagene_{start,stop}.tsv` and the legacy alias
  `periodicity_{start,stop}.tsv` gain new columns `normalize` and
  `n_transcripts` plus a `# normalize: ...` header comment so the
  flavour is unambiguous; the legacy alias also carries a
  `# DEPRECATED: ... will be removed in v1.0.0` header. Plot y-axis
  labels read `"P-site density (per-gene unit-mean)"` (vs. the legacy
  `"P-site reads"` only when the user opts back into `normalize="none"`).
  Exposed at the pipeline level via `mitoribopy rpf
  --periodicity-metagene-normalize {per_gene_unit_mean,none}` and the
  matching `metagene_normalize` key in the `mitoribopy all` YAML
  `periodicity:` block.

- **`from_fastq` rnaseq mode no longer surfaces Wald p-values.** The
  in-tree pyDESeq2 path runs on the 13-mt-mRNA subset only; the
  dispersion-shrinkage estimator is unreliable at that gene count, so
  `de_table.tsv` and `rpf_de_table.tsv` written by the from-FASTQ
  flow now have `padj = NA` for every row. `log2FoldChange` and
  `baseMean` are preserved (point estimates remain meaningful). The
  publication-grade flow (`rnaseq_mode=de_table` with an external
  full-transcriptome DE table) is unchanged. New `nullify_padj`
  parameter on `mitoribopy.rnaseq.de_analysis.deseq2_to_de_table` lets
  a downstream caller opt into the same behaviour. Stderr prints a
  banner explaining the choice when the from-FASTQ path runs.

- **Native dual-end UMI support.** `umi_position` accepts `both`
  (alongside `5p` / `3p`); `--umi-length-5p` and `--umi-length-3p`
  flags + matching sample-sheet columns wire dual-end UMI libraries
  (xGen Duplex, Twist) through cutadapt's two-pass trim. Pass 1
  extracts the 5' UMI; pass 2 extracts the 3' UMI and APPENDS to the
  QNAME so `umi_tools` dedups on the concatenated `<5pUMI><3pUMI>`
  token. `recommend_umi_method` uses the combined length when picking
  `directional` vs `unique`.

### Added — statistical hardening
- **Bootstrap CI + circular-shift permutation null for the metagene
  Fourier QC.** `fourier_period3_score_combined.tsv` (schema bumped
  to `1.1`) gains `amp_3nt_ci_{low,high}`,
  `spectral_ratio_3nt(_local)_ci_{low,high}`, `permutation_p`,
  `permutation_p_local`, plus audit fields `n_bootstrap`,
  `n_permutations`, `ci_alpha`, `ci_method`, `null_method`.
  Defaults: 200 bootstrap iterations + 200 circular-shift permutations
  + 90 % percentile CI + RNG seed 42 (recorded in
  `periodicity.metadata.json`). Vectorised via a precomputed DFT basis
  matrix so the additional cost is on the order of seconds. New CLI
  knobs on `mitoribopy periodicity`: `--bootstrap-n`,
  `--permutations-n`, `--ci-alpha`, `--random-seed`, `--no-stats`.
  Same knobs reachable from `mitoribopy rpf` via
  `--periodicity-fourier-{bootstrap-n,permutations-n,ci-alpha,random-seed}`
  and `--periodicity-no-fourier-stats`, and from `mitoribopy all` via
  the YAML `periodicity:` block (`fourier_bootstrap_n`,
  `fourier_permutations_n`, `fourier_ci_alpha`, `fourier_random_seed`,
  `no_fourier_stats`). Below `MIN_GENES_FOR_BOOTSTRAP_CI = 3`
  qualifying tracks the CI is skipped (NaN columns + `ci_method ==
  "skipped_too_few_genes"`) rather than emitted at misleadingly tight
  precision.

### Added — provenance + audit
- **`run_manifest.json` JSON Schema** at
  `src/mitoribopy/data/run_manifest.schema.json` (draft 2020-12). Two
  public helpers in `mitoribopy.data`: `load_run_manifest_schema()`
  and `validate_run_manifest(manifest)`. The validator uses
  `jsonschema` when available (added to the `[dev]` and new `[schema]`
  extras) and falls back to a minimal required-key + type check
  otherwise so the package gains no new hard dependency.
  `mitoribopy all` validates its own manifest at the end of every
  run; a schema drift records a `W_MANIFEST_SCHEMA_DRIFT` warning
  (mirrored to `warnings.tsv` and embedded in the manifest under
  `warnings`) rather than failing the run, so the user keeps every
  other output. `tests/test_run_manifest_schema.py` validates a
  representative manifest in CI.
- **Per-(sample, transcript) strand-sanity audit.**
  `<output>/rpf/qc/strand_sanity_per_transcript.tsv` (schema `1.0`) is
  written alongside the per-sample summary so a localised antisense
  bleed-through on one mt-mRNA (an unmasked NUMT, a misannotated
  reference contig) surfaces at a glance instead of being averaged
  away. Public helper:
  `mitoribopy.analysis.periodicity.compute_per_transcript_strand_sanity`.
  Both `strand_sanity.tsv` and `strand_sanity_per_transcript.tsv` are
  registered in `outputs_index.tsv`.
- **`E_OFFSET_FALLBACK_USED` warning** records (mirrored to
  `warnings.tsv` and the manifest) when a per-sample offset table was
  supplied but the requested sample is missing from it — the silent
  fallback to the combined-across-samples offsets that previously
  biased downstream codon-usage tables for that sample is now visible
  to a reviewer. Per-(stage, sample) deduplication keeps the log
  compact.

### Added — packaging + reproducibility
- **Bundled smoke fixture under [`examples/smoke/`](examples/smoke/).**
  Three synthetic mt-mRNAs (`MT-CO1`, `MT-ND1`, `MT-ATP6`), two
  samples (one WT, one KO), no UMIs, no contaminants. A user installs
  `mitoribopy>=0.7.0` and runs `python generate_smoke_fastqs.py`
  followed by `mitoribopy all --config pipeline_config.smoke.yaml
  --output results/` — the whole pipeline (cutadapt → bowtie2 → BED →
  offsets → translation profile → coverage → metagene Fourier QC)
  finishes in 10–30 s and asserts every file in
  [`expected_outputs.txt`](examples/smoke/expected_outputs.txt) exists +
  is non-empty. New pytest marker `smoke`; `pytest -m smoke` re-runs
  the same flow in CI when the external tools are available, skipped
  cleanly otherwise. README's *Installation* section gains a
  "30-second smoke test" callout.
- **`mitoribopy all --print-config-template --profile {minimal,
  publication,exhaustive}`** unifies the three template flavours behind
  a single flag. `minimal` (default; backwards-compatible) is the
  curated short template. `publication` appends a documented overlay
  declaring `strict: true`, `rnaseq_mode: de_table` (commented),
  and explicit Fourier statistical-hardening defaults
  (`fourier_bootstrap_n: 200`, `fourier_permutations_n: 200`,
  `fourier_ci_alpha: 0.10`, `metagene_normalize: per_gene_unit_mean`),
  so a methods-paper run sees the right gates without hunting for
  them. `exhaustive` prints the full annotated example from
  [`examples/templates/pipeline_config.example.yaml`](examples/templates/pipeline_config.example.yaml).
  Unknown profiles fall back to `minimal` with a warning. Tests in
  [`tests/test_print_config_template_profiles.py`](tests/test_print_config_template_profiles.py).
- **Theil-Sen slope CI surfaced through the codon-correlation outputs.**
  The CI was always computed (`scipy.stats.theilslopes` returns
  `lo_slope` / `hi_slope` by default) but was dropped before reaching
  the user. v0.7.0 routes `RobustSlope_CI_low` / `RobustSlope_CI_high`
  into the codon-correlation summary CSV, embeds them in the per-plot
  metadata sidecar (under `slope`, `slope_ci_low`, `slope_ci_high`,
  `slope_ci_alpha=0.05`), and renders them inline in the figure as
  `slope = X.XXX  95% CI [lo, hi]`. Tests in
  [`tests/test_codon_correlation_slope_ci.py`](tests/test_codon_correlation_slope_ci.py).
- **`pipeline/manifest.py`** extracted from `cli/all_.py` (which had
  grown past 2 350 LoC). The new module owns `MANIFEST_SCHEMA_VERSION`,
  the manifest builder (`write_manifest`), and the helpers
  (`sha256_of`, `read_stage_settings`, `git_commit`, `yaml_dump`,
  `build_stages_block`, `lift_tool_versions`,
  `W_MANIFEST_SCHEMA_DRIFT_CODE`). `cli/all_.py` re-exports the
  public surface so downstream code that imports from the original
  location keeps working; the CLI module is now ~300 lines smaller
  (2 353 → 2 050 LoC). Tests in
  [`tests/test_pipeline_manifest_module.py`](tests/test_pipeline_manifest_module.py).

### Added — UX surfaces
- **`SUMMARY.md` "Periodicity statistical confidence" section.**
  `summary/render.py` reads the score TSV and renders a per-sample
  table with `spectral_ratio_3nt`, its 90 % bootstrap CI, the
  circular-shift permutation p, and the `snr_call` tier. A reviewer no
  longer has to open the TSV to see the QC verdict + uncertainty.
- **Fourier figure annotations show CI + permutation_p** inline next
  to `spectral_ratio_3nt` (e.g. `8.42x  CI [6.10, 10.55]  p<0.001`).
  The per-figure metadata sidecar embeds the same fields so
  `validate-figures` can reason about them without re-reading the
  score TSV.
- **Per-plot `.metadata.json` sidecars** for codon-correlation
  scatters, translation-profile footprint-density depth plots, and
  codon-usage bar plots. Eliminates `FIGURE_QC_WARN` in
  `mitoribopy validate-figures` for those three families.
- **Per-frame split coverage plot.**
  `<output>/rpf/coverage_profile_plots/{p_site,a_site}_density_{rpm,raw}_frame_split/`
  directories complement the existing `_frame` overlay folders. Each
  figure stacks three sub-rows per sample (frame 0, +1, +2) sharing
  the y-axis, so low-frame signal is no longer hidden under tall
  same-position bars from another frame.
- **Stage-level config validation.** `mitoribopy validate-config`
  rejects unknown keys inside `align`, `rpf`, `rnaseq`, `execution`,
  and `periodicity` sections. Typos surface a `did you mean …`
  suggestion via `difflib.get_close_matches`. `--strict` upgrades the
  warning to a hard error.
- **`CITATION.cff`** at the repository root for software citation. The
  README's *Citation* section reproduces the suggested manuscript
  block.
- **Test coverage.** New test files: `test_fourier_stats.py`
  (12 tests), `test_metagene_normalization.py` (5 tests),
  `test_offset_fallback_warning.py` (5 tests),
  `test_per_transcript_strand_sanity.py` (4 tests),
  `test_from_fastq_padj_guard.py` (3 tests),
  `test_run_manifest_schema.py` (9 tests),
  `test_fourier_edge_cases.py` (11 tests for empty BED / zero-coverage
  / short-CDS / single-period basis), and
  `test_version_consistency.py` (cross-source version pin check).

### Removed (BREAKING)
- **`fft_period3_power.tsv` and the `--fft-period3` /
  `--periodicity-fft-period3` flags.** Replaced by the metagene
  Fourier bundle.
- **`mark-duplicates` dedup strategy** (already removed in v0.4.5;
  re-confirmed here as the publication-safe default — coordinate-only
  dedup destroys codon-occupancy signal on low-complexity mt-Ribo-seq
  libraries; see `docs/validation/taco1_ko_regression.md`).

### Fixed
- **`build_qc_summary` `best_read_length_dominant_fraction` bug.** The
  field was computed via a `A and B or C` ternary that also called
  `Series.get(...)` like a dict; when the dominant frame differed
  from the expected frame the value collapsed to the expected-frame
  fraction instead of the actual maximum. Replaced with a direct
  `max(f0, f1, f2)` from the chosen row.
- **`build_qc_summary` depth-aware best-length pick.** The candidate
  pool is now restricted to rows clearing `min_reads_per_length`
  whenever any qualify, with fallback to all rows otherwise — a 5-read
  length with lucky frame-0 dominance no longer outranks a 5,000-read
  length at 0.75.
- **`exclude_start_codons` / `exclude_stop_codons` plumbed through the
  pipeline path** as parameters to `compute_frame_summary`,
  `compute_frame_summary_by_length`, `build_gene_periodicity`,
  `run_periodicity_qc`, and `run_periodicity_qc_bundle`.

### Notes
- `MANIFEST_SCHEMA_VERSION` is at `1.3.0`. The `output_schemas` map
  in the manifest records `fourier_period3_score_combined.tsv: "1.1"`
  and `strand_sanity_per_transcript.tsv: "1.0"`. Downstream consumers
  that gate on schema versions should accept the new minor version.
- The `[schema]` optional-dependency extra installs `jsonschema` for
  full draft 2020-12 validation of `run_manifest.json` without pulling
  in the full dev tooling: `pip install 'mitoribopy[schema]'`.
- `CLAUDE.md`'s back-compat-shim removal note now correctly cites
  v1.0.0 (was the long-passed v0.4.0).

## [0.6.2] - 2026-05-01

### Publication-readiness — fourth-edit cleanup

Bug-fix + reproducibility release on top of v0.6.0. Every change is
software-engineering hygiene; no biological logic change.

### Added
- **Top-level `execution:` config block.** Resource-planning lives in
  one place at the run root. Keys: `threads`, `memory_gb`,
  `parallel_samples`, `single_sample_mode`, `min_threads_per_sample`,
  `estimated_memory_per_sample_gb`. The values cascade into `align`
  and `rnaseq` (`align_threads` / `max_parallel_samples` /
  `single_sample_mode` / `memory_gb`) when those stage keys are
  unset; explicit stage values still win. The `mitoribopy all
  --threads N` CLI flag now populates `execution.threads` so the
  global budget is honoured by every stage.
- **`<run_root>/resource_plan.json`** is written by `mitoribopy all`
  for every real run. It records the resolved thread budget,
  per-sample worker count, single-sample mode, scheduler, and the
  reason auto-mode picked those values. Embedded into
  `run_manifest.json` under `resource_plan`. Already-existing
  `<run_root>/align/resource_plan.json` continues to record the
  align-stage view.
- **Periodicity QC bundle extension.** `analysis.periodicity` now
  emits a top-level `qc_summary.tsv` (best read length, expected /
  dominant frame fraction, entropy bias, per-call totals,
  `overall_qc_call`) and a `qc_summary.md` companion under
  `<run_root>/rpf/periodicity/`. Per-length frame tables also gain a
  `expected_frame_enrichment` column. Optional `--phase-score` and
  `--fft-period3` add ribotricer-style consistency and FFT period-3
  power columns. Default thresholds:
  `good_frame_fraction=0.60`, `warn_frame_fraction=0.50`,
  `min_reads_per_length=1000`. QC labels: `good`, `warn`, `poor`,
  `low_depth`.
- **Stable warning-/error-code registry.** `docs/reference/warning_codes.md`
  enumerates every code emitted via `mitoribopy.io.warnings_log`
  (severity, stage, sample-id requirement, remediation). The
  `WARNING_CODES` Python registry in
  `src/mitoribopy/io/warning_codes.py` is the single source of
  truth; each `record(...)` call validates `code` against it (in
  test mode) so a typo cannot ship a new alias silently.
- **Release checklist.** `docs/developer/release_checklist.md` is the
  mandatory gate for every PyPI cut: version-source agreement,
  template validation, toy-data run, sdist/wheel build, TestPyPI
  smoke, PyPI upload, Zenodo archive, manuscript-version pin.
- **Architecture-history doc.** `docs/developer/architecture_history.md`
  collects the refactor-era design notes that previously lived in
  module docstrings (Phase 6 v0.3.0, etc.).

### Changed
- **`dedup_strategy` canonical value is `umi_coordinate`.** The legacy
  `umi-tools` / `umi_tools` strings are accepted and rewritten to
  `umi_coordinate` by `migrate-config`; the resolved value in
  `canonical_config.yaml` and the run manifest is always
  `umi_coordinate`. CLI help (`mitoribopy align --help`) lists the
  three canonical choices `auto | umi_coordinate | skip`. The
  underlying implementation is still UMI-aware coordinate-collapse
  via `umi_tools`; only the public name changed to describe the
  statistical operation rather than the tool.
- **`mitoribopy all --threads N` propagates.** The top-level CLI flag
  now writes into `execution.threads` and cascades into every stage
  whose own `threads:` (or `align_threads:`) is unset. The previous
  template comment instructing users to copy `threads:` into every
  stage section is removed. Explicit stage settings still take
  precedence.
- **YAML `allow_pseudo_replicates:` rewritten by `migrate-config`.**
  The legacy template key now serialises to the real CLI flag
  `--allow-pseudo-replicates-for-demo-not-publication` instead of
  the never-existed `--allow-pseudo-replicates`. The migration is
  silent for `false` (the safe default) and emits a clearly named
  warning code when set to `true`.
- **Public module docstrings.** Phase / refactor / "v0.3.0" tags
  removed from every module docstring under `src/mitoribopy/`. The
  history is preserved in the new
  `docs/developer/architecture_history.md`.

### Fixed
- `examples/templates/pipeline_config.example.yaml` rnaseq section
  no longer carries the misnamed `allow_pseudo_replicates: false`
  comment without a migration note. The template is now strict-mode
  clean.
- `_dict_to_argv` now serialises `dedup_strategy: umi_coordinate`
  and `allow_pseudo_replicates: true` to the canonical CLI flags.

## [0.6.0] - 2026-05-01

### Publication-readiness freeze (refactor-4)

Manuscript-target release. The blocking issues identified in the
2026-05-01 publication-readiness audit are addressed; no new
biological features.

### Removed
- **Pre-subcommand fallback removed.** Invoking `mitoribopy <flags>`
  without an explicit subcommand previously routed to `mitoribopy rpf`
  with a one-line deprecation warning that incorrectly claimed the
  fallback would be removed in v0.4.0. The fallback is now removed
  for real; the call exits with code 2 and a clear message that
  points at the right subcommand. Always use `mitoribopy rpf …`,
  `mitoribopy align …`, etc.

### Changed
- **`mitoribopy rpf` parser: hyphen-style flags are now canonical.**
  Every flag that was previously underscore-only (`--offset_type`,
  `--min_5_offset`, `--codon_density_window`, `--mt_mrna_substring_patterns`,
  …) now has a hyphenated public form (`--offset-type`,
  `--min-5-offset`, `--codon-density-window`,
  `--mt-mrna-substring-patterns`, …). Underscore aliases continue to
  parse for one transition cycle but are no longer surfaced in
  `--help` or in the `mitoribopy all` dry-run plan. The parser
  also reports as `mitoribopy rpf` rather than the bare `mitoribopy`,
  and its examples and `--config` metavar describe the modern
  YAML / JSON / TOML inputs.
- **`mitoribopy all` orchestrator emits hyphenated argv for every
  stage.** The previous rpf-specific underscore special-case in
  `_dict_to_argv` is now legacy-only; `align`, `rpf`, and `rnaseq`
  all receive the same hyphenated flag style. The underscore mode
  is still selectable but is not used by any first-party call site.

### Docs
- Quick-start, CLI reference, and config schema now lead with the
  canonical hyphenated rpf flags. The `--merge_density` /
  `--mrna_ref_patterns` / `selected_site` aliases are still listed
  under "Synonyms" with their canonical replacements.
- Tutorials and the `examples/templates/run_*.example.sh` scripts
  are rewritten to use hyphenated flags.

## [0.5.1] - 2026-04-28

### Added
- **`mitoribopy rnaseq` gains four WT-vs-X comparison plots.** The
  default plot set under `<output>/rnaseq/plots/` now includes
  `de_volcano_mrna.png` (mRNA DE volcano with red / blue / grey
  significance colouring + threshold guides at `padj<0.05` and
  `|L2FC|>=1`), `de_volcano_rpf.png` (Ribo-seq DE volcano; default
  flow only), `te_compare_scatter.png` (per-gene mean log2(TE) in
  the base condition vs the compare condition with the identity
  line), and `te_log2fc_bar.png` (sorted bar of
  `log2(TE_compare/TE_base)` per gene). Plot titles include the
  contrast label (`<base> vs <compare>`) when both conditions are
  set. The `te_compare_scatter` and `te_log2fc_bar` plots also
  emit in the `--de-table` (alternative) flow whenever a
  condition map plus base / compare conditions are provided.
- **Second pyDESeq2 fit on the Ribo-seq subset (default flow).**
  The default flow now runs pyDESeq2 twice — once on the RNA
  subset (writes the existing `de_table.tsv`) and once on the
  Ribo-seq subset (writes the new `rpf_de_table.tsv`). The Ribo
  fit feeds `de_volcano_rpf.png`. When the Ribo subset has fewer
  than two condition levels (e.g. only one `--ribo-fastq`
  condition), the second fit is skipped with a stderr WARNING and
  the RPF volcano is omitted; the rest of the run is unaffected.
- **`--base-sample` / `--compare-sample` CLI aliases for
  `--condition-a` / `--condition-b`** (with matching YAML keys
  `base_sample:` / `compare_sample:`). The new spelling mirrors
  the rpf config's `base_sample` key so the contrast reference
  reads the same across both stages. The legacy `--condition-a` /
  `--condition-b` flags keep working unchanged. Passing both an
  alias and the legacy form simultaneously is allowed only when
  the values agree (a mismatch exits with code 2 and a clear
  error). `run_settings.json` records both `condition_a` /
  `condition_b` and `base_sample` / `compare_sample` so downstream
  tooling can pick either field.

### Changed
- `run_deseq2(...)` now takes an `assay` keyword (default
  `"rna"`); pass `assay="ribo"` to fit the Ribo-seq subset. The
  default behaviour is unchanged for existing callers.
- **Publication-style refresh of every rnaseq plot.** Every
  figure now renders under a shared rc-context: Okabe-Ito
  colour-blind-safe palette (up = vermillion, down = blue,
  n.s. = light grey, consistent across volcanos / MA / TE),
  sans-serif fonts (Arial / Helvetica fallback), top + right
  spines hidden, minor ticks visible, 300 dpi PNG output with
  an editable-text SVG sidecar (`svg.fonttype = none`) for
  every figure. Gene labels on the volcanos and scatters wear
  white-bbox tags with thin leader lines; a small greedy
  8-direction label placer (no `adjustText` dependency)
  prefers above-the-point positions and avoids overlap. Stat
  boxes on the volcanos report `genes plotted / up / down /
  n.s.`; the TE compare scatter adds Pearson r. The
  `te_bar_by_condition` plot now overlays every replicate as
  a black dot jittered along x within its bar's footprint.
  The `te_heatmap` gains a colour strip above the columns
  spelling the condition assignment, with auto-contrast
  cell annotations. The `mrna_vs_rpf` four-quadrant scatter
  carries faint "co-regulated up / buffered up / co-regulated
  down / buffered down" captions. The `te_log2fc_bar`
  prints the numeric log2FC at each bar's tip.
- README, release notes, and the four template files
  (`rnaseq_config.example.yaml`, `pipeline_config.example.yaml`,
  `run_rnaseq.example.sh`, `run_pipeline.example.sh`) were
  updated to use the new `--base-sample` / `--compare-sample`
  spelling and document the new plots and `rpf_de_table.tsv`.

## [0.5.0] - 2026-04-27

### Added
- **`mitoribopy rnaseq` is now from-FASTQ-by-default.** The default
  flow takes raw RNA-seq + Ribo-seq FASTQs and a transcriptome FASTA
  and runs trimming, bowtie2 alignment, per-transcript counting, and
  pyDESeq2 itself before emitting TE / ΔTE / plots. The original
  "bring your own DE table" flow stays fully supported as the
  alternative flow (`--de-table`). The two flows are mutually
  exclusive (passing both `--rna-fastq` and `--de-table` exits with
  code 2). Running `mitoribopy rnaseq` with no flags now reports the
  default-flow required args in the missing-flag error message and
  emits a HINT line pointing at the alternative flow. SE vs PE is auto-detected from
  filename mate tokens (`_R1/_R2`, `_1/_2`, `.1./.2.`,
  `_R1_001/_R2_001`, `_read1/_read2`). Adapter is auto-detected via
  the existing `align.adapter_detect`; UMI presence is inferred from
  per-position Shannon entropy on the first/last 16 nt of reads
  (conservative — returns length=0 when uncertain). The bowtie2 index
  is content-addressed (cached at
  `<workdir>/bt2_cache/transcriptome_<sha12_of_fasta>`) so repeated
  runs against the same FASTA skip the rebuild. The reference-
  consistency hard-gate is skipped in Mode B; the FASTA SHA256 is
  recorded under `from_fastq.reference_checksum` in `run_settings.json`
  instead. The pre-computed-DE flow (Mode A) is unchanged.
- **`[fastq]` optional-dependency extra** (`pydeseq2 >= 0.4`). Soft-
  imported so the existing pre-computed-DE flow keeps working without
  it. Install with `pip install 'mitoribopy[fastq]'`.
- **Auto-pseudo-replicate fallback for n=1 conditions in the default flow.**
  pyDESeq2 needs at least two samples per condition to estimate
  dispersion. When a condition has exactly one sample, its FASTQ is
  automatically stream-split by record parity (record N → rep1 if
  even, rep2 if odd) so pyDESeq2 sees n=2. A loud stderr WARNING is
  emitted per split. The augmented condition map (original entries +
  rep1 / rep2 entries) lands at
  `<output>/condition_map.augmented.tsv`. Disable with
  `--no-auto-pseudo-replicate`.
- **Four new diagnostic plots** added to the `rnaseq` output (both
  flows for the first three; default flow only for the fourth):
  `plots/ma.png` (DESeq2 MA plot), `plots/te_bar_by_condition.png`
  (log2(TE) per gene grouped by condition with SE error bars),
  `plots/te_heatmap.png` (gene × sample log2(TE) heatmap, RdBu
  centred at 0), `plots/sample_pca.png` (PC1 vs PC2 from log1p
  counts, SVD-based — no scikit-learn dep).
- **New CLI flags on `mitoribopy rnaseq`** for the default flow:
  `--rna-fastq PATH [PATH ...]`, `--ribo-fastq PATH [PATH ...]`,
  `--reference-fasta PATH`, `--bowtie2-index PREFIX`,
  `--workdir DIR`, `--align-threads N`,
  `--no-auto-pseudo-replicate`.
- **New modules under `mitoribopy.rnaseq`:** `fastq_pairing`
  (`FastqSample`, `enumerate_fastqs`, `detect_samples`),
  `umi_detect` (`UmiDetectionResult`, `detect_umi`), `alignment`
  (`SampleAlignmentResult`, `align_sample`, `count_per_transcript`,
  `write_counts_matrix`, `write_long_counts`,
  `build_bowtie2_index`), `de_analysis` (`build_sample_sheet`,
  `run_deseq2`, `deseq2_to_de_table`, `write_de_table_tsv`),
  `split_replicates` (`stream_split_by_parity`,
  `split_sample_into_pseudo_replicates`).
- **New example templates**:
  [`examples/templates/run_rnaseq.example.sh`](examples/templates/run_rnaseq.example.sh)
  and
  [`examples/templates/rnaseq_config.example.yaml`](examples/templates/rnaseq_config.example.yaml)
  exhaustively cover both modes; the existing
  [`run_pipeline.example.sh`](examples/templates/run_pipeline.example.sh)
  gains an `RNASEQ_MODE=a|b` toggle, and
  [`pipeline_config.example.yaml`](examples/templates/pipeline_config.example.yaml)
  documents the Mode A / Mode B keys side-by-side.

### Changed
- **`mitoribopy rnaseq` reframed around the from-FASTQ default flow.**
  The CLI help banner, the argparse argument-group ordering (default-
  flow inputs first, alternative `--de-table` flow at the bottom), the
  README rnaseq sections, the tutorial, the CLI reference, and the
  example templates all now lead with the default flow. The label
  "Mode A / Mode B" used in interim drafts has been retired in favour
  of "default flow" / "alternative flow (--de-table)".
- **Missing-flags error message uses the default flow.** Running
  `mitoribopy rnaseq` with no flags now lists the default-flow
  required args (`--rna-fastq`, `--reference-fasta`, …) and emits a
  HINT line about the alternative flow. Previously listed the
  `--de-table` flow's required args. Tests updated accordingly.
- **Example templates default to the from-FASTQ flow.**
  `run_rnaseq.example.sh` defaults to `FLOW="default"`;
  `run_pipeline.example.sh` defaults to `RNASEQ_FLOW="default"`.
  Pass `FLOW="de-table"` / `RNASEQ_FLOW="de-table"` to run the
  alternative flow.

### Known gaps
- **PE + UMI in the default flow is currently `NotImplementedError`.**
  Preprocess UMIs into the read name first (e.g. via `umi_tools
  extract`), or pass `--de-table` with your own pre-computed DE
  results.
- **Genome / splice-aware aligner in the default flow is out of
  scope.** The bowtie2 index is built from a transcriptome FASTA,
  matching the package's mt-mRNA focus.

### Added (from the pre-0.5.0 `[Unreleased]` work; rolled into 0.5.0)
- **Per-step timing in `mitoribopy rpf`.** Each of the 7 pipeline
  steps (initialize → load read counts → unfiltered read-length QC →
  filter BED → offset enrichment → offset selection → downstream
  modules) is now wrapped in a stopwatch in `pipeline.runner` and
  emits a one-line duration immediately after its existing
  `[PIPELINE] Step K/7 OK: …` banner, e.g.
  `[PIPELINE] step 4 filter BED: 8.7s`. After all steps complete (or
  the pipeline bails out early), a per-step timing summary table is
  appended to both the console and `mitoribopy.log` plus an overall
  wall-clock line — symmetric with the per-sample timing already in
  `mitoribopy align`. New `mitoribopy.progress.render_step_timeline`
  helper backs the table. The internal `_emit_step_ok` / `Step K/7
  OK: …` banner format is unchanged, so existing consumers that
  grep the log for step completion continue to work.
- **Per-stage timing in `mitoribopy align`.** Each per-sample stage
  (cutadapt → bowtie2 contam → mt-align → MAPQ → dedup → BAM→BED) now
  emits one compact line carrying its wall-clock duration, e.g.
  `[ALIGN] sampleA: trim          12.3s — kit=auto, kept 941k/1.0M reads`.
  At the end of the run a per-stage summary table is appended to both
  the console and `mitoribopy.log`, showing total/mean/max across all
  samples plus an overall wall-clock line — useful for spotting which
  step dominates a run and for sizing `--threads`/`--max-parallel-samples`.
  Resume-only runs (every sample loaded from cache) skip the table since
  no fresh stage timings exist to summarize. Backed by a new
  `mitoribopy.progress` module (`Stopwatch`, `StageTimings`,
  `format_duration`, `render_summary_lines`) plus an optional tqdm-based
  `SampleCounter` that activates when tqdm is installed and the run is
  in parallel mode (graceful no-op fallback otherwise).
- **`--max-parallel-samples N` for `mitoribopy align`.** Run multiple
  samples concurrently in a `ThreadPoolExecutor`. Default `1` (serial,
  fully backward-compatible). When combined with `--threads T`, each
  worker's external tools (cutadapt, bowtie2, umi_tools) get
  `max(1, T // N)` threads, so total CPU use stays ≈ T regardless of
  N — set `--threads` to your total CPU budget and `--max-parallel-samples`
  to how many samples should run side-by-side. Resume-cached samples
  short-circuit the pool. On any per-sample failure the run is
  fail-fast: pending futures are cancelled and the first exception is
  re-raised, matching the previous serial behaviour. The joint
  `mitoribopy rpf` stage is unaffected — offset selection and
  `rpf_counts.tsv` aggregate across all samples and remain serial.
  Also wired through `mitoribopy all` via the `align.max_parallel_samples`
  YAML key.

## [0.4.5] - 2026-04-25

### Added
- **`--rpf_min_count_frac FRAC`** auto-filter for offset selection. The
  default `0.20` (on by default) keeps only read lengths whose total
  count across all samples is >= 20% of the most-enriched length's
  count. Pair with a wide `--rpf` range (e.g. `--rpf 27 36` for human)
  to let the data decide which lengths carry real signal. Set to `0`
  to disable. The kept and dropped lengths are reported once via
  `[BED] read-length auto-filter:`.

### Changed
- **Compact `[COVERAGE]` log.** The per-transcript line inside
  `run_coverage_profile_plots` is replaced by one summary line per
  requested site plus one for read coverage (~4 `[COVERAGE]` lines
  total instead of one per transcript per site).

### Changed (cont'd)
- **Templates moved to `examples/templates/`.** The shell-script
  templates (`run_align.example.sh`, `run_rpf.example.sh`,
  `run_pipeline.example.sh`) and the YAML templates
  (`align_config.example.yaml`, `rpf_config.example.yaml`,
  `pipeline_config.example.yaml`) now sit side-by-side under
  `examples/templates/` instead of being split between the repo root
  and `examples/run_scripts/`. Each file carries a
  `Compatible with: MitoRiboPy 0.4.5+` header line. README and the
  `01_end_to_end_fastq.md` tutorial point at the new path.

### Removed
- **`--dedup-strategy mark-duplicates` and the confirmation flag.**
  The picard `MarkDuplicates` option (coordinate-only dedup) is
  removed. It destroys codon-occupancy signal on low-complexity
  mt-Ribo-seq libraries, and the v0.4.4 confirmation gate did not
  protect users who copy-pasted the flag from an internet recipe.
  `--dedup-strategy` choices are now `{auto, umi-tools, skip}` and
  the long `--i-understand-mark-duplicates-destroys-mt-ribo-seq-signal`
  flag no longer exists. The `picard` external-tool dependency is
  also removed from the install instructions and the startup tool-
  check. See `docs/validation/taco1_ko_regression.md` for the
  empirical evidence that motivated the removal.

## [0.4.4] and earlier (collapsed)

Detailed entries for v0.4.4 → v0.1.0 are preserved in git history (`git log -- CHANGELOG.md`) and the per-version pages under [docs/release-notes/](docs/release-notes/). Highlights:

- **v0.4.4** — flat per-sample translation_profile / coverage_profile_plots layouts (site encoded in filename prefix); per-site `--cor_plot` outputs; `offset_diagnostics/` rename; per-sample `offset_applied.csv` audit.
- **v0.4.2** — per-sample UMI / kit overrides via `align.samples:` YAML list; `short` footprint class.
- **v0.4.1 / v0.4.0** — `mitoribopy all` orchestrator with `run_manifest.json`; auto-wiring of `align/` outputs into `rpf/`; `--print-canonical-config`.
- **v0.3.0** — `mitoribopy align` subcommand (cutadapt + bowtie2 + umi_tools wrappers); kit-preset registry; per-sample adapter detection.
- **v0.2.0** — `mitoribopy rpf` and `mitoribopy rnaseq` subcommands; legacy `_legacy` modules removed.
- **v0.1.0** — initial package release with `mitoribopy` CLI entry point and P/A-site offset selection.
