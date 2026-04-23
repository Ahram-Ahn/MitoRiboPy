# MitoRiboPy v0.3.0 — Phase 0 CLI & Public-API Surface

_Complete enumeration of user-visible CLI flags, importable symbols, and file paths that subsequent phases will touch._

## 1. CLI flags (baseline `mitoribopy` at v0.2.0)

Exhaustively extracted from `src/mitoribopy/pipeline/runner.py::build_parser`. Every entry lists: flag(s), argparse type/action, default source, config-file key, and which phase will modify or wrap it.

### Core Inputs

| Flag | Type / action | Default | Config key | Touched in |
|---|---|---|---|---|
| `--config <path>` | path | `None` | — | Phase 2 (add to every subcommand) |
| `-f`, `--fasta <path>` | path (**required**) | — | — | Phase 2 (rpf + all) |
| `-s`, `--strain {y,h,custom}` | choice | `y` | `strain` | Phase 2 |
| `-d`, `--directory <path>` | path | cwd | `directory` | Phase 2 |
| `-rpf MIN MAX` | 2×int | `None` (→ strain preset) | `rpf` | Phase 2 |
| `--annotation_file <csv>` | path | `None` | `annotation_file` | Phase 2 |
| `--codon_tables_file <json>` | path | `None` | `codon_tables_file` | Phase 2 |
| `--codon_table_name <name>` | str | `None` | `codon_table_name` | Phase 2 |
| `--start_codons <CODON ...>` | list[str] | strain preset | `start_codons` | Phase 2 |
| `--atp8_atp6_baseline {ATP6,ATP8}` | choice | `ATP6` | `atp8_atp6_baseline` | Phase 2 |
| `--nd4l_nd4_baseline {ND4,ND4L}` | choice | `ND4` | `nd4l_nd4_baseline` | Phase 2 |

### Offset Enrichment and Selection

| Flag | Type / action | Default | Config key | Touched in |
|---|---|---|---|---|
| `-a`, `--align {start,stop}` | choice | `start` | `align` | Phase 2 |
| `-r`, `--range <nt>` | int | `20` | `range` | Phase 2 |
| `--min_offset <nt>` | int | `11` | `min_offset` | Phase 2 |
| `--max_offset <nt>` | int | `20` | `max_offset` | Phase 2 |
| `--min_5_offset <nt>` | int | `None` | `min_5_offset` | Phase 2 |
| `--max_5_offset <nt>` | int | `None` | `max_5_offset` | Phase 2 |
| `--min_3_offset <nt>` | int | `None` | `min_3_offset` | Phase 2 |
| `--max_3_offset <nt>` | int | `None` | `max_3_offset` | Phase 2 |
| `--offset_mask_nt <nt>` | int | `5` | `offset_mask_nt` | Phase 2 |
| `--offset_pick_reference {selected_site,p_site}` | choice | `p_site` | `offset_pick_reference` | Phase 2 |
| `--offset_type {5,3}` | choice | `5` | `offset_type` | Phase 2 |
| `--offset_site {p,a}` | choice | `p` | `offset_site` | Phase 2 |
| `--codon_overlap_mode {full,any}` | choice | `full` | `codon_overlap_mode` | Phase 2 |
| `-p`, `--psite_offset <nt>` | int | `None` | `psite_offset` | Phase 2 |

### Outputs and Plotting

| Flag | Type / action | Default | Config key | Touched in |
|---|---|---|---|---|
| `-o`, `--output <dir>` | path | `analysis_results` | `output` | Phase 2 |
| `--downstream_dir <name>` | str | `footprint_density` | `downstream_dir` | Phase 2 |
| `--plot_dir <name>` | str | `plots_and_csv` | `plot_dir` | Phase 2 |
| `-fmt`, `--plot_format {png,pdf,svg}` | choice | `png` | `plot_format` | Phase 2 |
| `--x_breaks <nt ...>` | list[int] | `None` | `x_breaks` | Phase 2 |
| `--line_plot_style {combined,separate}` | choice | `combined` | `line_plot_style` | Phase 2 |
| `--cap_percentile <float>` | float | `0.999` | `cap_percentile` | Phase 2 |
| `-m`, `--merge_density` | store_true | `False` | `merge_density` | Phase 2 |
| `--order_samples <name ...>` | list[str] | `None` | `order_samples` | Phase 2 |

### Read-Count Normalization

| Flag | Type / action | Default | Config key | Touched in |
|---|---|---|---|---|
| `--read_counts_file <path>` | path | `read_counts_summary.txt` | `read_counts_file` | Phase 2 |
| `--read_counts_sample_col <name>` | str | `None` | `read_counts_sample_col` | Phase 2 |
| `--read_counts_reads_col <name>` | str | `None` | `read_counts_reads_col` | Phase 2 |
| `--rpm_norm_mode {total,mt_mrna}` | choice | `total` | `rpm_norm_mode` | Phase 2 |
| `--read_counts_reference_col <name>` | str | `None` | `read_counts_reference_col` | Phase 2 |
| `--mrna_ref_patterns <pat ...>` | list[str] | `["mt_genome","mt-mrna","mt_mrna"]` | `mrna_ref_patterns` | Phase 2 |

### Optional Modules

| Flag | Type / action | Default | Config key | Touched in |
|---|---|---|---|---|
| `--structure_density` | store_true | `False` | `structure_density` | Phase 2 |
| `--structure_density_norm_perc <float>` | float | `0.99` | `structure_density_norm_perc` | Phase 2 |
| `-v`, `--varna` | **HIDDEN** alias → `structure_density` | SUPPRESS | `varna` (legacy) | Phase 1 (decide) |
| `--varna_norm_perc` | **HIDDEN** alias → `structure_density_norm_perc` | SUPPRESS | `varna_norm_perc` (legacy) | Phase 1 (decide) |
| `--cor_plot` | store_true | `False` | `cor_plot` | Phase 2 |
| `--base_sample <name>` | str | `None` | `base_sample` | Phase 2 |
| `--cor_mask_method {percentile,fixed,none}` | choice | `percentile` | `cor_mask_method` | Phase 2 |
| `--cor_mask_percentile <float>` | float | `0.99` | `cor_mask_percentile` | Phase 2 |
| `--cor_mask_threshold <float>` | float | `None` | `cor_mask_threshold` | Phase 2 |
| `--use_rna_seq` | store_true | `False` | `use_rna_seq` | Phase 5 (deprecate) |
| `--rna_seq_dir <dir>` | path | `None` | `rna_seq_dir` | Phase 5 (deprecate) |
| `--rna_order <name ...>` | list[str] | `None` | `rna_order` | Phase 5 (deprecate) |
| `--rna_out_dir <name>` | str | `rna_seq_results` | `rna_out_dir` | Phase 5 (deprecate) |
| `--do_rna_ribo_ratio` | store_true | `False` | `do_rna_ribo_ratio` | Phase 5 (deprecate) |

### Top-level only

| Flag | Behavior |
|---|---|
| `--version`, `-V` | Print `MitoRiboPy <__version__>` (handled in `cli.py`) |
| `--help`, `-h` | Built-in argparse help |
| (positional) `run --` | Stripped by `_normalize_args` (legacy `mitoribopy run -- <args>` form) |

### Validation rules enforced in `parse_pipeline_args`

- `min_offset <= max_offset`, same for 5' and 3' variants
- `offset_mask_nt >= 0`
- `range > 0`
- `0 < cap_percentile <= 1`
- `0 < structure_density_norm_perc <= 1`
- `0 < cor_mask_percentile < 1` (when `cor_mask_method == percentile`)
- `--use_rna_seq` requires `--rna_seq_dir` and `--rna_order`
- `--strain custom` requires `--annotation_file` and `-rpf MIN MAX` and (`--codon_tables_file` or `--codon_table_name`)
- `min_5_offset/max_5_offset/min_3_offset/max_3_offset` fall back to `min_offset/max_offset` if `None`

## 2. Public importable symbols

Public = re-exported from a package `__init__.py` or documented in README.

| Import path | Kind | Touched in |
|---|---|---|
| `mitoribopy.__version__` | attr | — |
| `mitoribopy.cli.main` | function | Phase 2 (delegates to new subcommand dispatcher) |
| `mitoribopy.cli._normalize_args` | function (test-only) | Phase 2 |
| `mitoribopy.pipeline.PipelineContext` | dataclass | — |
| `mitoribopy.pipeline.build_parser` | function | Phase 2 |
| `mitoribopy.pipeline.parse_pipeline_args` | function | Phase 2 |
| `mitoribopy.pipeline.run_pipeline` | function | Phase 2 (wrapped by `cli/rpf.py`) |
| `mitoribopy.pipeline.run_pipeline_cli` | function | Phase 2 |
| `mitoribopy.pipeline.main` | function | Phase 2 |
| `mitoribopy.config.DEFAULT_CONFIG` | dict | Phase 1 (keys added/removed) |
| `mitoribopy.config.load_user_config` | function | Phase 1 (legacy key map) |
| `mitoribopy.config.resolve_rpf_range` | function | — |
| `mitoribopy.config.PhaseOneConfig` | dataclass | Phase 1 (possible removal) |
| `mitoribopy.config.load_phase_one_config` | function | Phase 1 (possible removal) |
| `mitoribopy.data.available_codon_table_names` | function | — |
| `mitoribopy.data.annotation_sequence_candidates` | function | — |
| `mitoribopy.data.build_sequence_display_map` | function | — |
| `mitoribopy.data.load_annotation_table` | function | — |
| `mitoribopy.data.load_codon_table` | function | — |
| `mitoribopy.data.load_codon_tables` | function | — |
| `mitoribopy.data.resolve_sequence_name` | function | — |
| `mitoribopy.data.resolve_start_codons` | function | — |
| `mitoribopy.data.transcript_display_title` | function | — |
| `mitoribopy.data.standard_codon_table` | dict | — |
| `mitoribopy.data.yeast_mitochondrial_codon_table` | dict | — |
| `mitoribopy.data.human_mitochondrial_codon_table` | dict | — |
| `mitoribopy.data.yeast_annotation_df` | DataFrame | — |
| `mitoribopy.data.human_annotation_df` | DataFrame | — |
| `mitoribopy.analysis.run_codon_correlation` | function | — |
| `mitoribopy.analysis.run_translation_profile_analysis` | function | — |
| `mitoribopy.analysis.run_rna_seq_analysis` | function | Phase 5 (replaced by `mitoribopy rnaseq`) |
| `mitoribopy.analysis.compute_offsets` | function | — |
| `mitoribopy.analysis.create_csv_for_offset_enrichment` | function | — |
| `mitoribopy.analysis.plot_offset_enrichment` | function | — |
| `mitoribopy.analysis.determine_p_site_offsets` | function | — |
| `mitoribopy.io.compute_total_counts` | function | — |
| `mitoribopy.io.process_bed_files` | function | Phase 4 (BAM input branch upstream) |
| `mitoribopy.io.compute_unfiltered_read_length_summary` | function | — |
| `mitoribopy.io.plot_unfiltered_read_length_heatmap` | function | — |
| `mitoribopy.plotting.apply_publication_style` | function | — |
| `mitoribopy.plotting.plot_codon_usage_dataframe` | function | — |
| `mitoribopy.plotting.plot_offset_enrichment` | function | — |
| `mitoribopy.plotting.plot_frame_usage_by_transcript` | function | — |
| `mitoribopy.plotting.plot_frame_usage_total` | function | — |
| `mitoribopy.plotting.plot_read_length_distribution` | function | — |
| `mitoribopy.plotting.plot_site_depth_profile` | function | — |
| `mitoribopy.plotting.plot_unfiltered_read_length_heatmap` | function | — |
| `mitoribopy.plotting.run_coverage_profile_plots` | function | — |
| `mitoribopy.plotting.run_structure_density_export` | function | — |
| `mitoribopy.tools.summarize_and_plot` | function | — |
| `mitoribopy.tools.reservoir_sample_lines` | function | — |

## 3. Config-file keys recognized by `load_user_config`

Source: `src/mitoribopy/config/runtime.py::DEFAULT_CONFIG`.

All `DEFAULT_CONFIG` keys above are accepted. Legacy keys remapped on load:
- `varna` → `structure_density`
- `varna_norm_perc` → `structure_density_norm_perc`

Unknown keys are logged as warnings and ignored.

## 4. File paths that subsequent phases will modify

### Phase 1 (breaking — `_legacy` and root wrappers)
- `main.py` (root) — decision pending (shim vs. remove)
- `pipeline_config.example.json` — update `"varna"` key → `"structure_density"`
- `src/mitoribopy/config/runtime.py` — decision on `legacy_key_map` retention
- `src/mitoribopy/pipeline/runner.py` — decision on hidden `--varna` / `--varna_norm_perc` retention
- `pyproject.toml` — version bump to `0.3.0`
- `CHANGELOG.md` — add `[0.3.0] - Unreleased` section
- `src/mitoribopy/__init__.py` — fallback version string
- `src/mitoribopy/config/defaults.py`, `config/loader.py`, `config/models.py` — vestigial; removal pending Ahram confirmation

_(Phase 1 is largely a **no-op on `_legacy` itself** since it was already removed in v0.2.0.)_

### Phase 2 (subcommand CLI)
- `src/mitoribopy/cli.py` — converted into `src/mitoribopy/cli/` package (per ToT approval)
- `src/mitoribopy/cli/__init__.py` (new)
- `src/mitoribopy/cli/__main__.py` or main entrypoint (new)
- `src/mitoribopy/cli/align.py` (new)
- `src/mitoribopy/cli/rpf.py` (new; wraps `pipeline.runner`)
- `src/mitoribopy/cli/rnaseq.py` (new)
- `src/mitoribopy/cli/all_.py` (new)
- `src/mitoribopy/cli/common.py` (new; shared `--config/--dry-run/--threads/--log-level`)
- `src/mitoribopy/__main__.py` — updated target
- `pyproject.toml` — `mitoribopy` entry point may stay `mitoribopy.cli:main`
- `tests/test_cli.py` — extended for subcommand dispatch + backward-compat fallback

### Phase 3 (`mitoribopy align`)
- `src/mitoribopy/align/__init__.py` (new)
- `src/mitoribopy/align/tool_check.py` (new)
- `src/mitoribopy/align/trim.py` (new)
- `src/mitoribopy/align/contam.py` (new)
- `src/mitoribopy/align/align.py` (new)
- `src/mitoribopy/align/dedup.py` (new)
- `src/mitoribopy/align/bam_utils.py` (new)
- `src/mitoribopy/align/read_counts.py` (new)
- `docs/environment/environment.yml` (new)
- `docs/environment/Dockerfile` (new)
- `tests/test_align_tool_check.py` (new)
- `tests/test_align_trim.py` (new; mocked subprocess)
- `tests/test_align_integration.py` (new; `@pytest.mark.requires_tools`)
- `tests/fixtures/mini_fastq/` + `tests/fixtures/mini_mtdna/` (new)
- `pyproject.toml` — add optional `pysam` extra or core dependency (decision at Phase 3.6)

### Phase 4 (BAM input for `rpf`)
- `src/mitoribopy/io/bed_reader.py` — gain BAM input path
- `src/mitoribopy/io/bam_reader.py` (new) or extension of `bed_reader.py`
- `src/mitoribopy/cli/rpf.py` — accept BAM
- `tests/test_bam_input.py` (new)

### Phase 5 (`mitoribopy rnaseq`)
- `src/mitoribopy/rnaseq/__init__.py` (new)
- `src/mitoribopy/rnaseq/de_loader.py` (new; DESeq2/Xtail/Anota2Seq auto-detect)
- `src/mitoribopy/rnaseq/te.py` (new; TE + ΔTE computation)
- `src/mitoribopy/rnaseq/reference_gate.py` (new; SHA256 reference-match check)
- `src/mitoribopy/rnaseq/plots.py` (new)
- `src/mitoribopy/cli/rnaseq.py` (new)
- `src/mitoribopy/analysis/rna_seq.py` — emit deprecation warning pointing to `mitoribopy rnaseq`
- `src/mitoribopy/pipeline/runner.py` — add `--use_rna_seq` deprecation banner
- `tests/test_rnaseq_*` (new)

### Phase 6 (`mitoribopy all`)
- `src/mitoribopy/cli/all_.py`
- `src/mitoribopy/pipeline/manifest.py` (new; `run_manifest.json`)
- `tests/test_pipeline_manifest.py` (new)

### Phase 7 (docs and validation)
- `README.md`
- `docs/tutorials/01_end_to_end_fastq.md` (new)
- `docs/tutorials/02_rnaseq_integration.md` (new)
- `docs/reference/cli.md` (new; auto-generated)
- `docs/diagrams/mitoribopy_package_flow.mmd` — update
- `docs/validation/taco1_ko_regression.md` (new)

### Phase 8 (release candidate)
- `pyproject.toml` — verify version
- `CHANGELOG.md` — finalize `[0.3.0]`
- git tag `v0.3.0-rc1`
