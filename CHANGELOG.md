# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

## [0.4.4] - 2026-04-25

### Added
- **Per-site codon correlation.** `--cor_plot` with `--analysis_sites both`
  now produces parallel outputs under `codon_correlation/p_site/` AND
  `codon_correlation/a_site/`. The function reads the matching site's
  `codon_usage_total.csv` from the new flat translation-profile layout
  and the existing P-site stop-codon override is preserved.
- **Granular read-coverage toggles.** New `--read_coverage_raw` /
  `--no-read_coverage_raw` and `--read_coverage_rpm` /
  `--no-read_coverage_rpm` flags (default both on) control the
  `coverage_profile_plots/read_coverage_*` folders independently.
- **IGV-compatible BedGraph export.** New `--igv_export` flag writes
  per-sample BedGraph tracks at
  `igv_tracks/<sample>/<sample>_{p_site,a_site}.bedgraph`. P-site
  tracks default to forest green, A-site to dark orange via the
  BedGraph `track color=` header. No new dependency.
- **FASTQ subsampler.** `python -m mitoribopy.tools.subsample_fastq`
  reservoir-samples 4-line FASTQ records (Algorithm R, deterministic
  with `--seed`); gzip auto-detected via the `.gz` suffix on either end.
  Sibling of the existing BED subsampler.
- **Per-sample offset audit.** New
  `offset_diagnostics/csv/per_sample_offset/<sample>/offset_applied.csv`
  records the exact offsets row each sample applied downstream, so a
  reviewer can confirm `--offset_mode per_sample` was honoured by the
  translation-profile and coverage-profile steps.

### Changed
- **Flat output layout — `p/` / `a/` subfolders dropped.**
  - `translation_profile/<sample>/codon_usage/` now holds
    `p_site_codon_usage_*.csv` and `a_site_codon_usage_*.csv` side by
    side. The legacy `translation_profile/{p,a}/<sample>/...` nesting
    is gone.
  - `coverage_profile_plots/p_site_density_*` and `a_site_density_*`
    are siblings of `read_coverage_*`, replacing the legacy
    `coverage_profile_plots/{p,a}/site_density_*` nesting.
  - The footprint-density CSV (already site-agnostic, columns
    `Position, Nucleotide, A_site, P_site`) is written once per
    transcript instead of once per site.
- **Default `plot_dir` renamed** from `plots_and_csv` to
  `offset_diagnostics`. CSVs are now under
  `offset_diagnostics/csv/`, plots under `offset_diagnostics/plots/`.
  `per_sample/` is renamed to `per_sample_offset/` and lives under
  `offset_diagnostics/csv/`. The `--plot_dir` flag still accepts a
  custom name.
- **`run_translation_profile_analysis()` signature**: prefer the new
  `requested_sites=["p", "a"]` keyword over `site_override=`. The
  legacy keyword is honoured for one cycle.
- **`run_coverage_profile_plots()` signature**: prefer the granular
  `write_read_coverage_raw=` / `write_read_coverage_rpm=` flags over
  the coarse `write_read_coverage=`. Legacy keyword still works.
- **Example run scripts** (`run_align.example.sh`,
  `run_rpf.example.sh`, `run_pipeline.example.sh`) moved from the repo
  root to `examples/run_scripts/`. Example YAML files stay at the repo
  root.

### Fixed
- **`cor_plot: true` silent skip bug.** Three of the four skip branches
  in `run_downstream_modules` previously appended `"codon correlation"`
  to `skipped_modules` without emitting a log line, so users with
  `cor_plot: true` and a valid `base_sample` could see no output and no
  error. Every skip path now logs an explicit `[PIPELINE]` reason.

## [0.4.2] - 2026-04-25

### Added
- **Per-sample UMI / kit overrides** via `align.samples:` YAML list (or
  `--sample-overrides TSV`). Each sample can independently set
  `kit_preset`, `adapter`, `umi_length`, `umi_position`, `dedup_strategy`.
  Overrides are reported in `kit_resolution.tsv` with a
  `per_sample_override:*` source label.
- **`short` footprint class** (16-24 nt RNase truncation products);
  unfiltered window 10-30 nt. Same defaults across `h.sapiens` and
  `s.cerevisiae`.
- **Per-sample resume markers** at `<output>/.sample_done/<sample>.json`.
  `mitoribopy align --resume` (auto-set by `mitoribopy all --resume`
  when `read_counts.tsv` is missing) reloads completed samples from
  marker JSON instead of re-running. Markers are written atomically.
- **`--keep-intermediates` flag** to retain the trimmed FASTQ /
  contam-filtered FASTQ / pre-MAPQ BAM that are otherwise deleted as
  soon as the next step consumes them.
- **PyPI install** path: `pip install mitoribopy`.
- **Custom organisms section in README** documenting the codon-table
  picker (NCBI Genetic Codes by organism group), annotation CSV
  schema, and a worked mouse mt example.
- **Comprehensive `pipeline_config.example.yaml`** at the repo root
  showing every flag with its default and a 1-line comment.
- **Four new pipeline diagrams** under `docs/diagrams/` (PNG, >=1920 px
  wide): pipeline overview + per-stage detail. Generated with
  `python docs/diagrams/render_diagrams.py` (matplotlib; no Node
  required).

### Changed
- **Strain values renamed** to species names: `h` -> `h.sapiens`
  (default), `y` -> `s.cerevisiae`. `vm` and `ym` removed; use
  `--strain custom` and pick a built-in NCBI Genetic Code with
  `--codon_table_name`. Old short forms `h` and `y` are accepted as
  deprecated synonyms with a one-line `[mitoribopy] DEPRECATED:`
  warning.
- **Output layout unified** under per-stage roots:
  `translation_profile/<sample>/{p,a}/...` and
  `coverage_profile_plots/{read_coverage_*, {p,a}/site_density_*}/`.
  The `translation_profile_p/`/`translation_profile_a/` and
  `coverage_profile_plots_p/`/`coverage_profile_plots_a/` split
  directories are gone. `read_coverage_*` plots are written once per
  run instead of duplicated per site.
- **`p_site_coverage_*` subdirs renamed to `site_density_*`** under the
  per-site (`p/` or `a/`) parent.
- **`footprint_density.csv` columns** are now
  `Position, Nucleotide, A_site, P_site`. `E_site` (always
  `P_site - 3`) and `*_selected_depth` (an outlier-capped duplicate)
  removed; A-site comes before P-site.
- **`disome` footprint windows rebalanced** to current literature:
  `h.sapiens` 50-70 (was 60-90), `s.cerevisiae` 60-90 (was 65-95).
  Disome unfiltered window tightened to 40-100 nt.
- **`dedup_strategy=skip` no longer writes a duplicate BAM** under
  `deduped/`. The orchestrator wires `aligned/<sample>.mapq.bam`
  straight into BED conversion for no-UMI samples.
- **Codon-correlation scatter is now publication-ready**: identity
  line, OLS regression, `r / R^2 / p / slope / intercept / N` stat
  box, leader-line labels for top-10 residual codons (no overlap),
  Okabe-Ito colour-blind-safe palette, SVG + 300 dpi PNG outputs.
- **Renamed CLI flags** with deprecated synonyms: `--merge_density` ->
  `--codon_density_window`; `--mrna_ref_patterns` ->
  `--mt_mrna_substring_patterns`; `--offset_pick_reference selected_site`
  -> `reported_site`. All old names emit a one-line
  `[mitoribopy] DEPRECATED:` notice on first use.
- **README rewritten as a standalone reference** (no version-history
  language). Synonyms for deprecated flags listed in one neutral
  table.
- **Tutorial 01 rewritten** against the current API (output tree,
  YAML config, codon-usage paths).

### Removed
- **`debug_csv/`** subdirectory under translation-profile output --
  unread by anything downstream.
- **Strain values `vm` and `ym`** (use `--strain custom`).
- **`src/mitoribopy/utils/`** empty subpackage scaffold.
- **`pipeline_config.example.json`** -- replaced by
  `pipeline_config.example.yaml` which is YAML, supports comments, and
  matches the actual three-section shape that `mitoribopy all --config`
  expects.
- **Mermaid diagram sources** (`docs/diagrams/*.mmd`) and the
  `generate_mitoribopy_diagram.py` / `render_mitoribopy_diagram.sh`
  helpers -- replaced by `render_diagrams.py` + checked-in PNGs.

### Fixed
- **Two-pass 3' UMI cutadapt** (QIAseq miRNA and similar) now uses
  `min_length + umi_length` for pass 1, so the post-pass-2 insert is
  guaranteed to meet the user-requested minimum length floor. Without
  this fix, a 12 nt UMI with the default `--min-length 15` could
  produce 3 nt inserts.
- **`mitoribopy align --config`** silently ignored every value in the
  loaded YAML / JSON / TOML file. The flag was registered on the
  parser (so it appeared in `--help`) but the values were never
  merged into the argparse namespace. Now applies the loaded values
  as parser defaults BEFORE the real parse, so explicit CLI flags
  still win on conflict (matches the documented contract and matches
  how `mitoribopy rpf --config` already worked). Adds
  `align_config.example.yaml` and `rpf_config.example.yaml` at the
  repo root as exhaustive flat-shape templates for the standalone
  subcommands.

## [0.4.1] - 2026-04-24

### Added
- **`pretrimmed` kit preset** for FASTQs that have already been
  adapter-trimmed (SRA-deposited data, output of a prior trim step).
  cutadapt still runs for length and quality filtering but skips the
  ``-a`` flag entirely. Auto-inferred when adapter detection finds 0%
  match across every kit (the strongest signal that no adapter is
  there); also selectable explicitly via ``--kit-preset pretrimmed``.
  Mixed batches with both raw and pre-trimmed FASTQs in the same run
  resolve each sample independently.
- **Detection tuning flags**: ``--adapter-detect-reads`` (default 5000),
  ``--adapter-detect-min-rate`` (default 0.30), ``--adapter-detect-min-len``
  (default 12), ``--adapter-detect-pretrimmed-threshold`` (default 0.05),
  and ``--no-pretrimmed-inference``. All exposed via the YAML ``align``
  section as well.
- ``DetectionResult.pretrimmed`` field carrying the universal-absence
  signal so callers can route to the right fallback.

### Changed
- **``KIT_PRESETS`` consolidated by adapter family** rather than vendor.
  v0.4.0 had a long list (``truseq_smallrna``, ``nebnext_smallrna``,
  ``nebnext_ultra_umi``, ``qiaseq_mirna``, ``truseq_stranded_total``,
  ``smarter_pico_v3``, ``sequoia_express``, …); v0.4.1 ships a short
  canonical set:
  - ``illumina_smallrna`` &mdash; TruSeq Small RNA adapter, no UMI.
  - ``illumina_truseq`` &mdash; TruSeq Read 1 adapter, no UMI. Covers
    NEBNext Multiplex Small RNA, TruSeq Stranded Total RNA Gold,
    Takara SMARTer Stranded Total v3 Pico, Bio-Rad SEQuoia Express
    Standard, and any other prep on the standard Illumina R1 adapter.
  - ``illumina_truseq_umi`` &mdash; TruSeq R1 + 8 nt 5' UMI. Covers
    NEBNext Ultra II UMI, SEQuoia Complete UMI.
  - ``qiaseq_mirna``, ``custom``, ``auto``, ``pretrimmed``.
  Old vendor-specific names remain accepted as aliases via
  ``KIT_PRESET_ALIASES`` and ``resolve_kit_alias``; existing YAML
  configs and CLI invocations keep working without changes.
- ``ResolvedKit.adapter`` is now ``str | None`` (was ``str``). Only the
  ``pretrimmed`` kit holds ``None``; cutadapt's ``_build_pass1_command``
  skips the ``-a`` flag when it sees ``None``.
- Adapter-detection failure with no fallback now falls through to the
  ``pretrimmed`` kit by default instead of hard-erroring. This is the
  fix for the ``"adapter auto-detection found no known kit"`` failure
  on SRA-deposited inputs reported in v0.4.0. Pass
  ``--no-pretrimmed-inference`` to keep the v0.4.0 strict behaviour.

## [0.4.0] - 2026-04-24

### Added
- **Per-sample adapter detection**: `mitoribopy align` now scans every
  input FASTQ independently and resolves the kit per sample. Mixed-kit
  runs (e.g. one TruSeq + one NEBNext UMI sample in the same batch) are
  fully supported. Detection logic lives in
  `mitoribopy.align.sample_resolve` and produces a
  `kit_resolution.tsv` next to `read_counts.tsv` recording each
  sample's detected kit, applied kit, match rate, and dedup decision.
- **Per-sample dedup resolution**: `--dedup-strategy auto` (the
  default) now resolves to `umi-tools` for each UMI-bearing sample and
  `skip` for each UMI-less sample within the same run. The required
  external tools (umi_tools, picard) are checked as the union of every
  sample's resolved strategy.
- **`--kit-preset auto`** (the new default) tells the CLI to rely
  entirely on per-sample auto-detection. An explicit preset
  (truseq_smallrna, nebnext_smallrna, nebnext_ultra_umi, qiaseq_mirna,
  custom) becomes a per-sample fallback used only when detection cannot
  identify a known adapter.
- **Polymorphic `align.fastq` YAML key**: `align.fastq` now accepts
  either a single directory string (treated as a `--fastq-dir`
  shortcut) or an explicit list of paths. Saves YAML users from listing
  every input file by hand.
- **Per-sample offset selection (`--offset_mode per_sample`, default)**:
  per-sample 5'/3' offsets drive each sample's downstream translation
  profile and coverage plots, so inter-sample offset drift no longer
  biases the combined output. The combined-across-samples offset table
  is still written as a diagnostic, and an `offset_drift_<align>.svg`
  plot surfaces per-sample drift at a glance. Pass
  `--offset_mode combined` to recover the v0.3.x global behaviour.
- **`--analysis_sites {p,a,both}`** (default `both`) explicitly controls
  which sites get downstream codon-usage and coverage outputs. With
  `both`, parallel `translation_profile_p/` and `translation_profile_a/`
  (and `coverage_profile_plots_p/`, `coverage_profile_plots_a/`)
  directories are produced so users can compare P-site and A-site
  results side by side without rerunning the pipeline.
- New module `mitoribopy.align.sample_resolve` with `SampleResolution`,
  `resolve_sample_resolutions`, `required_dedup_tools`, and
  `write_kit_resolution_tsv`.
- New helper `mitoribopy.analysis.build_per_sample_summaries` exposed
  for downstream tooling that wants the per-sample enrichment table
  without re-running the full pipeline.

### Changed
- **`mitoribopy align`** writes a per-sample provenance table
  (`kit_resolution.tsv`) and embeds the same data under
  `run_settings.json -> per_sample`. The legacy top-level keys
  `kit_preset`, `adapter`, `umi_length`, `umi_position` are replaced by
  `kit_preset_default`, `adapter_default`, and the per-sample list.
  `dedup_strategy` at the top level reports `mixed` when the run mixes
  UMI and non-UMI samples.
- **`--offset_site`** help text clarifies that the flag controls only
  the offset *selection* coordinate space, not which downstream outputs
  are produced. To pick which downstream outputs to generate, use
  `--analysis_sites`.
- **Translation profile and coverage plots** now write per-site
  subdirectories when `--analysis_sites=both` (the default). Single-site
  runs (`--analysis_sites=p` or `=a`) keep the legacy top-level layout
  for backward compatibility.
- **YAML config template** updated by `mitoribopy all
  --print-config-template`: shows `kit_preset: auto` as the new
  default, the polymorphic `fastq:` shortcut, the new
  `offset_mode: per_sample` and `analysis_sites: both` keys, and a
  per-sample resolution explainer.

### Removed
- The legacy single-FASTQ-only adapter sanity check (`_apply_adapter_detection`)
  in `mitoribopy.cli.align`. The per-sample resolver replaces it.

## [0.3.0] - 2026-04-23
### Added
- **`mitoribopy align --adapter-detection {auto|off|strict}`** (default
  `auto`): scans the head of the first input FASTQ for known adapter
  signatures from `KIT_PRESETS` and either auto-selects the matching
  preset (when `--kit-preset` is left at `custom`) or hard-fails on
  mismatch (`strict` mode). Catches the silent failure where a wrong
  `--kit-preset` drops ~99% of reads as "too long" without any error.
  New module `mitoribopy.align.adapter_detect` exposes `detect_adapter`
  and `DetectionResult` for programmatic use. The chosen mode is
  recorded under `adapter_detection_mode` in `run_settings.json`.
- **`mitoribopy all --print-config-template`**: prints a fully-commented
  YAML covering every stage with sensible defaults and inline preset
  choices. New users can pipe it into a file to start a new project in
  one step.
- **`--footprint_class {monosome,disome,custom}`** (default `monosome`):
  biological footprint preset that injects sensible defaults for `-rpf`
  and `--unfiltered_read_length_range`. Monosome defaults are h/vm
  28-34, y/ym 37-41, unfiltered 15-50; disome widens to h/vm 60-90,
  y/ym 65-95, unfiltered 40-110. Explicit user values always win.
- **New strain presets `vm` (any vertebrate mito) and `ym` (any
  fungus with yeast-mito codon code)**: both inherit the codon table
  for their anchor strain but require user-supplied `--annotation_file`
  and an explicit `-rpf` range.
- **Frame-colored P-site plots at CDS nucleotide resolution**: color
  each position by `(nt_pos - cds_start) % 3` using the Okabe-Ito
  colorblind-safe palette. Frame-0 dominance is the canonical
  mt-Ribo-seq QC signature; the three-color partition makes it visible
  at a glance. Written under `p_site_coverage_{rpm,raw}_frame/`.
- **RPF-window band on the per-sample read-length distribution plot**:
  `plot_read_length_distribution` takes an optional `rpf_range` kwarg
  and shades the RPF window so the reader sees which lengths entered
  downstream analysis.
- **Defensive guards** for predictable failure modes:
  malformed BED rows (`end <= start`) are dropped with a per-file
  WARNING in both filtered and unfiltered paths; FASTA records with no
  matching annotation row now emit a loud WARNING listing up to 10
  missing records; `--threads <= 0` raises `SystemExit` with a clear
  message.
- **Logging-style convention documented** at the top of
  `mitoribopy.console`: every console line routes through
  `log_info` / `log_warning` / `log_error` / `log_progress` with a
  `[COMPONENT]` prefix; bare `print(...)` is reserved for pre-logger
  argparse validation. Audited and cleaned the last two offenders
  (`tools/subsample.py`, `tools/read_composition.py`).
- **`PipelineContext` convention**: every piece of derived state lives
  on the context, not on the parsed `argparse.Namespace`.
  `run_coverage_profile_plots` now takes `total_mrna_map` as an
  explicit keyword argument.
- **Regression test**: read-length CSV summary and the value that feeds
  the distribution plot must always agree, row-for-row. Locks in the
  invariant against the originally-reported ~3 nt shift.
- **New subcommand CLI** under `src/mitoribopy/cli/`:
  - `mitoribopy align`: FASTQ -> BAM + BED6 + per-sample read counts.
  - `mitoribopy rpf`: Ribo-seq analysis pipeline from BED/BAM inputs.
  - `mitoribopy rnaseq`: DE-table + Ribo-seq integration; TE /
    &Delta;TE with a SHA256 reference-consistency gate.
  - `mitoribopy all`: align + rpf + (optional) rnaseq orchestrator
    with a composed `run_manifest.json`.
  - Top-level `mitoribopy --help` lists all subcommands; shared
    options: `--config`, `--dry-run`, `--threads N`,
    `--log-level {DEBUG,INFO,WARNING,ERROR}`.
- **Config files** in JSON / YAML / TOML (auto-detected by suffix).
  YAML uses PyYAML; TOML uses stdlib `tomllib` on Py 3.11+ or the
  optional `tomli` fallback via the `toml-py310` extra.
- New dependencies: `PyYAML>=6.0`, `pysam>=0.22`.

### Changed
- **Orchestrator `_dict_to_argv` is now flag-style-aware** so
  `mitoribopy all` correctly runs rpf. align / rnaseq sections keep
  hyphen-converted flags (`--kit-preset`); rpf keeps underscored
  flags (`--offset_type`, `--min_5_offset`, ...). Fixes the previous
  end-to-end failure where the orchestrator emitted `--offset-type 5`
  and rpf rejected it.

### Deprecated
- `mitoribopy rpf --use_rna_seq` and its four companion flags
  (`--rna_seq_dir`, `--rna_order`, `--rna_out_dir`,
  `--do_rna_ribo_ratio`) remain available but print a runtime
  DEPRECATION notice; scheduled for removal in v0.4.0. Use the
  dedicated `mitoribopy rnaseq` subcommand for DE-based TE /
  &Delta;TE analysis with a SHA256 reference-consistency gate.

### Phase 3 and earlier (all below) landed during the v0.3.0 development cycle
- New subcommand CLI under `src/mitoribopy/cli/`:
  - `mitoribopy align`: FASTQ -> BAM + BED6 + per-sample read counts. See "Phase 3 (`mitoribopy align`)" below.
  - `mitoribopy rpf`: Ribo-seq analysis pipeline from BED inputs (delegates to the existing `mitoribopy.pipeline.runner`).
  - `mitoribopy rnaseq` (Phase 5 stub): DE-table + Ribo-seq integration; flag surface defined.
  - `mitoribopy all` (Phase 6 stub): align + rpf + (optional) rnaseq orchestrator; flag surface defined.
  - Top-level `mitoribopy --help` lists all subcommands; `mitoribopy <sub> --help` shows per-subcommand flags.
- Every subcommand supports shared options: `--config`, `--dry-run`, `--threads N`, and `--log-level {DEBUG,INFO,WARNING,ERROR}`.
- `mitoribopy.config.load_user_config` now auto-detects JSON, YAML (`.yaml`/`.yml`), and TOML (`.toml`) config files. YAML support uses `PyYAML`; TOML uses stdlib `tomllib` on Python 3.11+ and optionally `tomli` on Python 3.10 (install via the `toml-py310` extra).
- `--threads N` exports `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, and `MITORIBOPY_THREADS` so BLAS/OpenMP libraries and external tools invoked by `align` honor the request.
- `--dry-run` prints each subcommand's planned steps and exits 0.
- New dependencies: `PyYAML>=6.0`, `pysam>=0.22`.
- New optional extra: `mitoribopy[toml-py310]` pulling in `tomli` for Python 3.10 users who want TOML configs.

#### Phase 3 (`mitoribopy align`)
- New `src/mitoribopy/align/` package implementing the FASTQ preprocessing pipeline for mt-Ribo-seq:
  - `tool_check.py`: PATH + version verification that reports every missing required tool in one error with a bioconda install command. Exposes `bowtie2_strand_flag(forward|reverse|unstranded)`.
  - `trim.py`: cutadapt wrapper with kit-aware adapter + UMI handling. Kit presets: `truseq_smallrna`, `nebnext_smallrna`, `nebnext_ultra_umi`, `qiaseq_mirna`, `custom` (default; forces explicit `--adapter`). Two-pass cutadapt for 3' UMI kits (QIAseq).
  - `contam.py`: bowtie2 end-to-end `--very-sensitive -L 18 --no-unal` for rRNA/tRNA subtraction; refuses to silently skip when `--contam-index` is missing.
  - `align.py`: bowtie2 end-to-end mt-transcriptome alignment (Path A). `pysam.sort` + `pysam.index` for sort/index so samtools is not a PATH dep.
  - `dedup.py`: Approach D - `auto` resolves to `umi-tools` when UMIs are present or `skip` otherwise (the required default for low-complexity mt-Ribo-seq). Opt-in `mark-duplicates` gated behind `--i-understand-mark-duplicates-destroys-mt-ribo-seq-signal` because coordinate-only dedup destroys codon-occupancy signal.
  - `bam_utils.py`: pysam-based MAPQ filter, BAM -> BED6 conversion (strand-aware), and `count_mapped_reads`.
  - `read_counts.py`: per-sample provenance TSV with columns `sample`, `total_reads`, `post_trim`, `rrna_aligned`, `post_rrna_filter`, `mt_aligned`, `unaligned_to_mt`, `mt_aligned_after_mapq`, `mt_aligned_after_dedup`.
- `mitoribopy align` CLI flags:
  - `--kit-preset`, `--adapter`, `--umi-length`, `--umi-position` (library-prep).
  - `--library-strandedness {forward,reverse,unstranded}` (default `forward`, enforced by bowtie2 `--norc`). Resolves the ND5 / ND6 antisense overlap on Path A by construction.
  - `--min-length` / `--max-length` defaults 15 / 45 nt (mt-RPF range).
  - `--contam-index`, `--mt-index` (bowtie2 index prefixes; both required; no silent defaults).
  - `--mapq` default 10 (NUMT cross-talk suppression).
  - `--dedup-strategy {auto,umi-tools,skip,mark-duplicates}` (default `auto`).
  - `--umi-dedup-method` default `unique` (strictest; avoids directional over-collapse).
  - `--seed` default 42.
  - `--fastq-dir` / `--fastq` (repeatable) inputs, `--output` run root.
- Output layout under `<output>/`: `trimmed/`, `contam_filtered/`, `aligned/`, `deduped/`, `bed/`, `read_counts.tsv`, `run_settings.json` (Phase 6 provenance record with kit, adapter, strandedness, dedup strategy, and tool versions).
- `docs/environment/environment.yml`: pinned bioconda env (cutadapt>=4.0, bowtie2>=2.4, samtools>=1.15, umi_tools>=1.1, fastqc, bedtools, pysam, PyYAML, dev tools).
- `docs/environment/Dockerfile`: miniconda3-based image; `ENTRYPOINT mitoribopy`.
- `tests/test_align_integration.py` (@pytest.mark.requires_tools): end-to-end run on a tiny synthetic library, auto-skipped when cutadapt / bowtie2 / bowtie2-build / samtools are absent via `conftest.py::pytest_collection_modifyitems`.

#### Phase 4 (BAM input for `mitoribopy rpf`)
- `mitoribopy rpf` now accepts BAM files alongside BED files under `--directory`. Each `<sample>.bam` is converted to BED6 in `<output>/bam_converted/<sample>.bed` via pysam and then flows through the unchanged BED pipeline. No samtools / bedtools PATH dependency (pysam bundles htslib).
- New CLI flag `--bam_mapq Q` on the `rpf` subcommand (default `10`; matches `mitoribopy align --mapq`). Set `0` to disable the pre-conversion MAPQ filter. Filtering defaults to `10` to preserve NUMT-suppression semantics when a raw alignment BAM is passed in from a non-MitoRiboPy source.
- `src/mitoribopy/io/bam_reader.py` (new module):
  - `convert_bam_to_bed(bam_in, bed_out, mapq_threshold)` - single-file BAM -> BED6 conversion, optionally MAPQ-filtered.
  - `prepare_bam_inputs(input_dir, converted_dir, mapq_threshold)` - batch conversion for every `*.bam` in `input_dir`; skips with a warning when a `<sample>.bed` already exists (explicit BED wins on name conflict).
  - Module docstring documents the exact BAM -> BED6 field mapping (`chrom`, `start`, `end`, `name`, `score`, `strand`) so downstream consumers can audit it.
- `src/mitoribopy/io/bed_reader.process_bed_files` gains an optional `bam_mapq` kwarg and now accepts a mixed `.bed` + `.bam` input directory. Sample-name ordering is deterministic regardless of source type.
- `src/mitoribopy/config/runtime.DEFAULT_CONFIG` gains a `bam_mapq: 10` entry so JSON / YAML / TOML configs can set the threshold declaratively.
- `tests/test_bam_reader.py` (10 tests): conversion happy path, MAPQ filter, mixed BED + BAM ingestion, name-conflict precedence, `bam_mapq` plumbing through `process_bed_files`, empty-directory handling.

#### Phase 5 (`mitoribopy rnaseq`)
- New `mitoribopy rnaseq` subcommand integrates a pre-computed differential-expression table (DESeq2 / Xtail / Anota2Seq auto-detected, or user-supplied column mapping via `--de-format custom`) with a prior `mitoribopy rpf` run and emits TE + ΔTE tables and plots.
- SHA256 reference-consistency gate: `rpf` now writes `<output>/run_settings.json` with a `reference_checksum` (SHA-256 of the `-f / --fasta` file). `rnaseq` verifies the caller-supplied `--reference-gtf` (hashed locally) or `--reference-checksum` against that recorded hash and HARD-FAILS on mismatch. This prevents the silently-invalid TE comparisons that happen when Ribo-seq and RNA-seq were aligned to different references.
- Required CLI flags (no defaults): `--gene-id-convention {ensembl,refseq,hgnc,mt_prefixed,bare}`, `--de-table`, `--ribo-dir`, `--output`, and exactly one of `--reference-gtf` / `--reference-checksum`. Optional: `--de-format`, `--de-gene-col / --de-log2fc-col / --de-padj-col / --de-basemean-col` for custom DE schemas, `--condition-map / --condition-a / --condition-b` for replicate-based ΔTE, `--organism {h,y}`.
- `src/mitoribopy/rnaseq/` modules:
  - `_types.py`: `DeFormat`, `GeneIdConvention`, `DeColumnMap`, `DeTable`, `TeRow`, `DTeRow`, plus the `DE_COLUMN_ALIASES` registry pinning the DESeq2 / Xtail / Anota2Seq column names.
  - `gene_ids.py`: curated `HUMAN_MT_MRNAS` (13 genes, Ensembl + RefSeq + HGNC + bare) and `YEAST_MT_MRNAS` (8 genes). `match_mt_mrnas(de_ids, convention, organism)` returns matched + missing lists; a partial-match WARNING names every missing gene so the user can see which IDs got lost.
  - `de_loader.py`: `detect_de_format(header)` keyed on the log2FC column name; `load_de_table(path, column_map)` auto-detects CSV/TSV delimiter and canonicalizes NA / Inf values.
  - `reference_gate.py`: `compute_reference_checksum(path)` (stdlib sha256), `verify_reference_consistency(ribo_dir, reference_path|reference_checksum)`, `ReferenceMismatchError` on any disagreement or missing recorded hash.
  - `counts.py`: `load_ribo_counts(path)` reads the new `<rpf>/rpf_counts.tsv` into `{gene: {sample: count}}`.
  - `te.py`: `compute_te` with a 0.5 pseudocount (Laplace-style); `compute_delta_te` uses the DE table's mRNA log2FC and computes an internal Ribo-seq log2FC from replicate means when `--condition-map` is provided, else emits per-gene rows with a `single_replicate_no_statistics` note.
  - `plots.py`: `plot_mrna_vs_rpf_scatter` (four-quadrant log2FC scatter with per-gene labels) and `plot_delta_te_volcano` (ΔTE_log2 vs −log10(padj)).
- `rpf` now emits `<output>/rpf_counts.tsv` (columns: `sample`, `gene`, `count`) at the end of the BED-filter step so `rnaseq` has deterministic input. Emitted regardless of whether any downstream modules were requested.
- `--use_rna_seq` flag on `mitoribopy rpf` is now marked DEPRECATED in the help text and emits a stderr deprecation warning at invocation. Scheduled for removal in v0.4.0. Users should migrate to `mitoribopy rnaseq`.
- `tests/test_rnaseq.py` (24 tests) + `tests/test_rnaseq_cli.py` (5 end-to-end tests): gene ID matching across all five conventions and both organisms, DE format auto-detection for every supported schema + `custom` fallback, reference-consistency match/mismatch/no-digest paths, ribo-counts loader, TE + ΔTE math (single-replicate note, replicate-based log2FC recovery to within 0.05, missing-gene-from-DE handling), full CLI orchestrator producing `te.tsv`, `delta_te.tsv`, scatter + volcano plots, and the reference-mismatch HARD FAIL path.

#### Phase 6 (`mitoribopy all`)
- New `mitoribopy all` subcommand wires `align -> rpf -> (optional) rnaseq` behind a single shared config file (YAML / JSON / TOML). Each stage's flags live under a top-level section (`align:`, `rpf:`, `rnaseq:`) matching the subcommand's own CLI; the orchestrator serializes the section into an argv list and calls the subcommand's `run()` entry point.
- Auto-wiring: `<output>/align/`, `<output>/rpf/`, `<output>/rnaseq/` stage directories are populated automatically; `rpf`'s `--directory` defaults to `<output>/align/bed/`; `rnaseq`'s `--ribo-dir` defaults to `<output>/rpf/` so downstream stages pick up upstream outputs without extra configuration.
- `--resume` skips stages whose sentinel output already exists (`align/read_counts.tsv`, `rpf/rpf_counts.tsv`, `rnaseq/delta_te.tsv`).
- `--skip-align / --skip-rpf / --skip-rnaseq` toggle individual stages even when the config has that section.
- Writes `<output>/run_manifest.json` at the end with the full composition: per-stage `run_settings.json` payloads, `stages_run`, `stages_skipped`, and a promoted top-level `reference_checksum` from the rpf stage so future rnaseq runs reading this manifest can find it without drilling.
- The `rnaseq` section is auto-activated only when it has a `de_table` entry; otherwise it is skipped with a log line.
- `tests/test_all_cli.py` (11 tests): dict-to-argv serializer covers scalars, bools, lists, and `None`; help text lists every toggle; dry-run prints per-stage argv; required-args aggregation; end-to-end orchestration (mocked subcommands) records three stages in the manifest and promotes `reference_checksum`; `--resume` skips align when `read_counts.tsv` is present; non-zero exit from any stage propagates upward; rnaseq section without `de_table` is skipped; `--skip-rpf` honored.

### Explicitly excluded
- Dispersion-based statistics (DESeq2 shrinkage, Xtail modeling) over only the 13 mt-mRNAs. `rnaseq` consumes pre-computed DE output rather than running its own DE; the small gene universe violates the assumptions those methods need. Users should run DE on the full transcriptome and hand the resulting table to `rnaseq`.

### Changed
- `mitoribopy <flags>` (no subcommand) still works in v0.3.x but now emits a stderr `DEPRECATION` warning and routes to `mitoribopy rpf <flags>`. This fallback will be removed in v0.4.0.
- The `src/mitoribopy/cli.py` module has been replaced by a `src/mitoribopy/cli/` package. The public names `mitoribopy.cli.main`, `mitoribopy.cli._normalize_args`, and `mitoribopy.cli.run_pipeline_cli` are preserved.

### Breaking
- Removed the root-level `main.py` compatibility shim. Use the installed `mitoribopy` entrypoint or `python -m mitoribopy` instead.
- Removed vestigial Phase-I configuration scaffolding:
  - Module `mitoribopy.config.defaults` (with `DEFAULT_PACKAGE_NAME`, `DEFAULT_IMPORT_NAME`).
  - Module `mitoribopy.config.loader` (with `load_phase_one_config`).
  - Module `mitoribopy.config.models` (with `PhaseOneConfig`).
  - Re-exports of `PhaseOneConfig` and `load_phase_one_config` from `mitoribopy.config`.
- Removed the hidden legacy CLI flags `-v`, `--varna`, and `--varna_norm_perc`. Use `--structure_density` and `--structure_density_norm_perc`.
- Removed the legacy JSON config key alias handling in `mitoribopy.config.load_user_config`. The keys `varna` and `varna_norm_perc` are no longer silently remapped; use `structure_density` and `structure_density_norm_perc` instead. Unknown keys are still reported as a warning.

### Changed
- Updated `pipeline_config.example.json` to use the `structure_density` key.
- Cleaned up stale `main.py` references in docstrings and comments.

## [0.2.0] - 2026-04-07
### Added
- Packaged JSON and CSV reference loading for codon tables and annotations.
- Custom organism support through `--annotation_file`, `--codon_tables_file`, `--codon_table_name`, and `--start_codons`.
- Example custom reference templates under `examples/custom_reference/`.
- Persistent run logging in `<output>/mitoribopy.log`.
- Additional tests for reference-data loading, bicistronic aliases, visualization behavior, and CLI coverage.

### Changed
- Replaced hardcoded reference tables with data-driven loaders under `mitoribopy.data.reference_data`.
- Split human and yeast bicistronic CDS annotations into separate logical transcript rows while preserving stable display labels.
- Renamed analysis and plotting modules to `translation_profile_analysis`, `coverage_profile_plots`, and `structure_density_export`.
- Reworked the README into a publishable package-oriented format.
- Updated package version metadata to `0.2.0`.

### Removed
- Runtime dependency on legacy `_legacy` pipeline modules.
- On-disk filtered BED output as part of normal pipeline execution.

## [0.1.0] - 2026-03-30
### Added
- First package-oriented release of MitoRiboPy.
- `mitoribopy` CLI entrypoint and package modules under `src/mitoribopy`.
- P-site/A-site offset selection support with `--offset_site` and `--offset_pick_reference`.
- mt-mRNA-only RPM normalization mode with configurable reference pattern matching.
- Human NT3/NT5 validation runs and migration documentation in `docs/reports/`.
