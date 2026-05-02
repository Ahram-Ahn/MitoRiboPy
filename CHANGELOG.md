# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed (BREAKING)
- **Periodicity QC bundle replaced with a Wakigawa-faithful metagene Fourier analysis.** The frame-fraction QC bundle (`qc_summary.tsv`/`.md`, `frame_counts_by_sample_length.tsv`, `frame_counts_by_gene.tsv`, `gene_periodicity.tsv`, `phase_score`) and the four named frame-heatmap plots (`frame_fraction_heatmap.svg`, `frame_by_length_heatmap.{png,svg}`, `read_length_periodicity_barplot.svg`, `gene_phase_score_dotplot.svg`) are removed. The Fourier-spectrum bundle is rewritten end-to-end to follow Wakigawa et al. 2025 (bioRxiv 2025.05.03.652009): per-gene normalise (divide by window mean) → element-wise mean across qualifying transcripts → mean-centre → Hann window → direct DFT at exactly period 3.0 (not bin-snapped). New default window: **99 nt = 33 codons** (multiple of 3, no period-3 leakage); skip 5 codons after AUG (initiation peak) and 1 codon before stop (termination peak). The `utr3` region was dropped — human mt-mRNA 3' UTRs are too short for a meaningful 99-nt window, and the user-facing artefact is cleaner without a sparse third panel. New outputs: `fourier_spectrum_combined.tsv` and `fourier_period3_score_combined.tsv` (one row per `[sample, read_length, gene_set, region]`; `gene_set ∈ {combined, ATP86, ND4L4}`; `region ∈ {orf_start, orf_stop}`). The score table includes `spectral_ratio_3nt` and an `snr_call` tier (`excellent ≥ 10×`, `healthy ≥ 5×`, `modest ≥ 2×`, `broken < 2×`). **Three figures per (sample, read_length)**: `*_combined.png` (single trace per panel from canonical mt-mRNAs), `*_ATP86.png` (junction-bracketed ATP8/ATP6 bicistronic — top: ATP8 frame, bottom: ATP6 frame, windows overlapping at the bicistronic junction nt 177-202), `*_ND4L4.png` (junction-bracketed ND4L/ND4 — windows flank the 4-nt junction). Each panel shows the spectral_ratio_3nt and snr_call as in-figure annotations.

  **Migration**: scripts that read `qc_summary.tsv` for the per-sample QC verdict should switch to `fourier_period3_score_combined.tsv` and filter by `gene_set == "combined"`, `region == "orf_start"` (or `"orf_stop"`); the headline column is `spectral_ratio_3nt`. Scripts that read `frame_counts_by_sample_length.tsv` to pick a best read length should use `mitoribopy summarize`'s output (which now reads from the Fourier score table). The `--phase-score`, `--periodicity-phase-score`, `--periodicity-exclude-start-codons`, `--periodicity-exclude-stop-codons`, `--periodicity-min-reads-per-length`, `--good-frame-fraction`, `--warn-frame-fraction`, `--min-reads-per-length`, `--min-reads-per-gene`, `--exclude-start-codons`, `--exclude-stop-codons` CLI flags are removed. The standalone `mitoribopy periodicity` subcommand now exposes `--fourier-window-nt`, `--drop-codons-after-start`, `--drop-codons-before-stop`, `--min-mean-coverage`, `--min-total-counts`, `--no-plots` instead. The annotation `sequence_aliases` field (already populated by `_apply_bicistronic_defaults` for ATP86 / ND4L4) is now consulted by the chrom-to-transcript lookup, so fused-FASTA datasets resolve to BOTH constituent transcripts (ATP8 + ATP6 / ND4L + ND4) for the bicistronic-figure panels.

  Statistical hardening (bootstrap CI over genes, circular-shift permutation null) is marked TODO in `src/mitoribopy/analysis/fourier_spectrum.py` — the `spectral_ratio_3nt` value does not change, but the `amp_3nt_ci_low/high` and `permutation_p` columns are not yet populated. Implement when reviewers ask.
- **Native dual-end UMI support.** `umi_position` accepts `both` (alongside `5p` / `3p`); `--umi-length-5p` and `--umi-length-3p` flags + matching sample-sheet columns wire dual-end UMI libraries (xGen Duplex, Twist) through cutadapt's two-pass trim. Pass 1 extracts the 5' UMI; pass 2 extracts the 3' UMI and APPENDS to the QNAME so umi_tools dedups on the concatenated `<5pUMI><3pUMI>` token. `recommend_umi_method` now uses the combined length when picking `directional` vs `unique`.
- **Per-plot `.metadata.json` sidecars** for codon-correlation scatters, translation-profile footprint-density depth plots, and codon-usage bar plots. Eliminates `FIGURE_QC_WARN` in `mitoribopy validate-figures` for those three families. Folder-level `codon_correlation.metadata.json` is preserved (downstream tooling references it) — the per-plot sidecar coexists with it.

### Removed (BREAKING)
- **`fft_period3_power.tsv` and the `--fft-period3` / `--periodicity-fft-period3` flags.** The single-scalar `power[k=1/3] / median(background)` computation collapsed too much information. Replaced by the Wakigawa metagene Fourier bundle (see the BREAKING entry above). Migration: scripts that read `fft_period3_power.tsv` should switch to `fourier_period3_score_combined.tsv` (`spectral_ratio_3nt` column with the Wakigawa-style `snr_call` tier).

- **Per-frame split coverage plot.** New
  `<output>/rpf/coverage_profile_plots/{p_site,a_site}_density_{rpm,raw}_frame_split/`
  directories complement the existing `_frame` overlay folders. Each
  figure stacks three sub-rows per sample (frame 0, +1, +2) sharing the
  y-axis, so low-frame signal at the same CDS position is no longer
  hidden under tall same-position bars from another frame. Useful for
  fused overlapping ORFs (e.g. human `ATP86` where the ATP6 ORF is
  +2 nt offset from ATP8). Each folder ships a
  `coverage_plot.metadata.json` sidecar with the same frame formula
  and palette as the overlay variant. Implemented as
  `_plot_frame_split` in
  `src/mitoribopy/plotting/coverage_profile_plots.py`.
- **Human mt-mRNA overlap-pair annotation.** `gene_periodicity.tsv`
  gains an `is_overlap_pair` column (default on) flagging known
  overlap regions: canonical names (`MT-ATP8`, `MT-ATP6`, `MT-ND4L`,
  `MT-ND4`), no-hyphen short forms (`MTATP8`, …), AND the fused-ORF
  spellings used by some FASTAs (`ATP86` = ATP8+ATP6,
  `ND4L4` = ND4L+ND4 — including the bundled
  `input_data/human-mt-mRNA.fasta`). Treat `is_overlap_pair=true` rows
  as inherently ambiguous-frame even when their per-frame fractions
  look clean. Public surface:
  `mitoribopy.analysis.periodicity_qc.HUMAN_MT_OVERLAPPING_GENES` and
  `is_known_overlap_gene(name)`.
- **`exclude_start_codons` / `exclude_stop_codons` plumbed through the
  pipeline path.** Added as parameters to `compute_frame_summary`,
  `compute_frame_summary_by_length`, `build_gene_periodicity`,
  `run_periodicity_qc`, and `run_periodicity_qc_bundle`. Defaults stay
  at 0 to preserve historical pooled numbers; the standalone
  `mitoribopy periodicity` CLI applies the spec defaults of 6 / 3
  uniformly across the per-length and gene-level tables.
- **`compute_phase_score` plumbing.** `run_periodicity_qc` now accepts
  `compute_phase_score=True` and forwards it to
  `run_periodicity_qc_bundle`. Previously the orchestrator hardcoded
  it off and there was no override path.
- **`periodicity.metadata.json` records `exclude_start_codons`,
  `exclude_stop_codons`, and `phase_score_enabled`** alongside the
  existing thresholds, so a reviewer can re-derive every periodicity
  number from disk without re-running the CLI.

### Fixed
- **`build_qc_summary` `best_read_length_dominant_fraction` bug.** The
  field was computed via a tortured `A and B or C` ternary that also
  called `Series.get(...)` like a dict; in the case where the dominant
  frame differed from the expected frame the value collapsed to the
  expected-frame fraction instead of the actual maximum. Replaced with
  a direct `max(f0, f1, f2)` from the chosen row.
- **`build_qc_summary` depth-aware best-length pick.** Previously the
  best read length was selected purely by `expected_frame_fraction`,
  so a 5-read length with lucky frame-0 dominance could outrank a
  5,000-read length at 0.75. Now the candidate pool is restricted to
  rows clearing `min_reads_per_length` whenever any qualify; falls
  back to all rows otherwise.

### Tests
- New `tests/test_periodicity_refinements.py` covers the depth-aware
  best-length pick, the dominant-fraction bug fix, the codon-edge
  exclusion plumbing, and the overlap-pair helper (including the
  fused-ORF spellings).

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

