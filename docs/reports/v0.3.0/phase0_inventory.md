# MitoRiboPy v0.3.0 — Phase 0 File Inventory

_Baseline snapshot of the repository at the start of the v0.3.0 refactor._

- Branch: `claude/mitoribopy-v0.3.0-refactor-rUWcy`
- Baseline package version: `0.2.0` (per `pyproject.toml`)
- Baseline HEAD: `7babcf5` "Refine docs and add release notes"
- Python: `>= 3.10` (CI tests 3.10 / 3.11 / 3.12, run locally on 3.11.15)
- Build backend: `setuptools>=68` via `pyproject.toml` / `[tool.setuptools.packages.find] where=["src"]`
- Entry point: `mitoribopy = mitoribopy.cli:main`
- Baseline tests: `29 passed` (see `phase0_baseline_tests.txt`)
- Baseline build: sdist + wheel succeed (see `phase0_baseline_build.txt`)

Legend:
- **Public API**: imported by external users or documented in `README.md`.
- **Runtime**: executed on a normal `mitoribopy` run.
- **Test-only**: loaded only from `tests/`.
- **Dev/meta**: build metadata, CI, docs, packaging.

## Repository root

| Path | Lines | Purpose | Public API? | Runtime? | Test-only? | Notes |
|---|---:|---|---|---|---|---|
| `pyproject.toml` | 64 | PEP 517/518 build + entry points + pytest config | — | — | — | Dev/meta. Version bump target. |
| `LICENSE` | — | MIT license | — | — | — | Dev/meta. |
| `README.md` | 335 | User-facing docs and CLI usage examples | — | — | — | Dev/meta. |
| `CHANGELOG.md` | 36 | Keep-a-Changelog history | — | — | — | Dev/meta. Needs `[0.3.0] — Unreleased` section. |
| `main.py` | 19 | Root-level thin shim: `from mitoribopy.pipeline.runner import main` | yes (shim) | yes | — | Documented in README; decision pending for v0.3.0. |
| `pipeline_config.example.json` | 24 | Example JSON config for `--config` | yes | — | — | Still has legacy `"varna": false`; needs migration to `structure_density`. |
| `.gitignore` | — | Git ignore rules | — | — | — | Dev/meta. |
| `.github/` | — | GitHub Actions CI + release workflows | — | — | — | Dev/meta. |

## `src/mitoribopy/` — package

### Top-level modules

| Path | Lines | Purpose | Public API? | Runtime? | Test-only? | Notes |
|---|---:|---|---|---|---|---|
| `__init__.py` | 8 | Exposes `__version__` via `importlib.metadata` | yes | yes | — | |
| `__main__.py` | 8 | `python -m mitoribopy` entrypoint → `cli.main` | yes | yes | — | |
| `cli.py` | 30 | Entrypoint, strips `run --`, handles `--version`, delegates to `pipeline.runner.run_pipeline_cli` | yes | yes | — | Will be replaced by `cli/` package in Phase 2. |
| `console.py` | 132 | Logging helpers (`log_message`, `log_warning`, `configure_file_logging`) | yes (internal) | yes | — | |

### `analysis/` — core Ribo-seq analysis

| Path | Lines | Purpose | Public API? | Runtime? | Test-only? | Notes |
|---|---:|---|---|---|---|---|
| `__init__.py` | 21 | Re-exports `run_codon_correlation`, `run_translation_profile_analysis`, `run_rna_seq_analysis`, `compute_offsets`, `create_csv_for_offset_enrichment`, `plot_offset_enrichment`, `determine_p_site_offsets` | yes | yes | — | |
| `codon_correlation.py` | 160 | Codon-correlation plots and masking | yes | yes (optional) | — | |
| `offset_enrichment.py` | 208 | Offset enrichment CSV and P/A-site enrichment calculations | yes | yes | — | **Biology-critical.** |
| `offset_selection.py` | 164 | `determine_p_site_offsets`, end-specific 5'/3' bound selection, P/A-site conversion | yes | yes | — | **Biology-critical.** Frame-0 reference + P/A-site shift logic. |
| `rna_seq.py` | 288 | Current `--use_rna_seq` integration (to be superseded by `mitoribopy rnaseq` in Phase 5) | yes | yes (opt-in) | — | |
| `translation_profile_analysis.py` | 728 | `run_translation_profile_analysis`: frame usage, codon usage, footprint density | yes | yes | — | **Biology-critical.** Largest module. |

### `config/` — configuration

| Path | Lines | Purpose | Public API? | Runtime? | Test-only? | Notes |
|---|---:|---|---|---|---|---|
| `__init__.py` | 13 | Re-exports `DEFAULT_CONFIG`, `load_user_config`, `resolve_rpf_range`, `PhaseOneConfig`, `load_phase_one_config` | yes | yes | — | |
| `defaults.py` | 9 | `DEFAULT_PACKAGE_NAME`, `DEFAULT_IMPORT_NAME` constants | yes (internal) | no | no | Vestigial Phase-I scaffold; candidate for removal. |
| `loader.py` | 11 | `load_phase_one_config()` returning `PhaseOneConfig` | yes | no | no | Vestigial Phase-I scaffold; candidate for removal. |
| `models.py` | 14 | `PhaseOneConfig` dataclass | yes | no | no | Vestigial Phase-I scaffold; candidate for removal. |
| `runtime.py` | 115 | `DEFAULT_CONFIG`, `load_user_config`, `resolve_rpf_range`, `legacy_key_map` (`varna` → `structure_density`) | yes | yes | — | Legacy map still active. |

### `data/` — packaged reference tables

| Path | Lines | Purpose | Public API? | Runtime? | Test-only? | Notes |
|---|---:|---|---|---|---|---|
| `__init__.py` | 36 | Re-exports loaders + in-memory codon / annotation tables | yes | yes | — | |
| `codon_tables.json` | — | Packaged codon tables (vertebrate_mitochondrial, yeast_mitochondrial, standard, custom_example) | yes (data) | yes | — | **Biology-critical.** |
| `codon_tables.py` | 19 | Convenience loaders returning dicts | yes | yes | — | |
| `human_annotation.csv` | — | Human mt-mRNA annotation | yes (data) | yes | — | **Biology-critical.** |
| `yeast_annotation.csv` | — | Yeast mt-mRNA annotation | yes (data) | yes | — | **Biology-critical.** |
| `reference_data.py` | 389 | Annotation loader, bicistronic baseline resolution, sequence-alias resolution | yes | yes | — | **Biology-critical.** |
| `transcript_annotations.py` | 16 | Legacy-style DataFrame accessors (`human_annotation_df`, `yeast_annotation_df`) | yes | yes | — | |

### `io/` — readers/writers

| Path | Lines | Purpose | Public API? | Runtime? | Test-only? | Notes |
|---|---:|---|---|---|---|---|
| `__init__.py` | 15 | Re-exports `process_bed_files`, `compute_total_counts`, `compute_unfiltered_read_length_summary`, `plot_unfiltered_read_length_heatmap` | yes | yes | — | |
| `bed_reader.py` | 226 | BED parsing, RPF length filtering, unfiltered read-length QC | yes | yes | — | |
| `read_counts.py` | 164 | Read-count table parsing (CSV/TSV/TXT) for RPM normalization | yes | yes | — | |

### `pipeline/` — orchestration

| Path | Lines | Purpose | Public API? | Runtime? | Test-only? | Notes |
|---|---:|---|---|---|---|---|
| `__init__.py` | 13 | Re-exports `PipelineContext`, `build_parser`, `parse_pipeline_args`, `run_pipeline`, `run_pipeline_cli`, `main` | yes | yes | — | |
| `context.py` | 28 | `PipelineContext` dataclass | yes | yes | — | |
| `runner.py` | 596 | `argparse` CLI (all flags), `parse_pipeline_args`, `run_pipeline_cli`, `main` | yes | yes | — | Moves into `cli/rpf.py` + this file in Phase 2. |
| `steps.py` | 386 | Pipeline step helpers consumed by `runner.run_pipeline` | yes (internal) | yes | — | |

### `plotting/` — figures

| Path | Lines | Purpose | Public API? | Runtime? | Test-only? | Notes |
|---|---:|---|---|---|---|---|
| `__init__.py` | 29 | Re-exports all plotting entrypoints | yes | yes | — | |
| `coverage_profile_plots.py` | 319 | Coverage profile plots | yes | yes | — | |
| `structure_density_export.py` | 97 | Structure-density export (log2 + scaled) | yes | yes (opt-in) | — | |
| `style.py` | 26 | Publication-style matplotlib defaults | yes | yes | — | |
| `translation_profile_plots.py` | 115 | Codon usage, frame usage, site-depth profile plots | yes | yes | — | |
| `visualization.py` | 389 | Offset enrichment plots, read-length distribution, heatmaps | yes | yes | — | |

### `tools/` — auxiliary utilities

| Path | Lines | Purpose | Public API? | Runtime? | Test-only? | Notes |
|---|---:|---|---|---|---|---|
| `__init__.py` | 6 | Re-exports `summarize_and_plot`, `reservoir_sample_lines` | yes | no | yes (partial) | Not wired into the normal pipeline run. |
| `read_composition.py` | 246 | Read-composition summary + plots | yes | no | yes (partial) | |
| `subsample.py` | 70 | Reservoir sampler for large BED/FASTQ | yes | no | yes | |

### `utils/`

| Path | Lines | Purpose | Public API? | Runtime? | Test-only? | Notes |
|---|---:|---|---|---|---|---|
| `__init__.py` | 2 | Empty placeholder | — | — | — | Scaffold only. |

## `tests/`

| Path | Lines | Purpose | Notes |
|---|---:|---|---|
| `conftest.py` | 11 | Adds `src/` to `sys.path` for local source-tree imports | |
| `test_cli.py` | 81 | CLI entry, normalization, parser help structure, custom-strain validation | Asserts `--varna` not in help. |
| `test_codon_tables.py` | 32 | Mito-specific codon reassignments present | **Biology gate.** |
| `test_config_runtime.py` | 41 | `DEFAULT_CONFIG`, `load_user_config` unknown-key handling, `resolve_rpf_range` | |
| `test_integration_pipeline.py` | 379 | End-to-end smoke + offset selection P/A-site shift + mask/bounds | **Biology gate.** |
| `test_read_counts.py` | 37 | Flexible column matching for read-count tables | |
| `test_reference_data.py` | 53 | Built-in annotation and codon-table loading, bicistronic aliases | **Biology gate.** |
| `test_subsample.py` | 39 | Deterministic reservoir sampler | |
| `test_visualization.py` | 43 | Offset enrichment plot axis behavior | |

## `docs/`

| Path | Purpose |
|---|---|
| `docs/README.md` | Docs index |
| `docs/PACKAGE_REFACTOR_SCHEME.md` | Original v0.1/0.2 refactor plan |
| `docs/release-notes/README.md` | Release-note index |
| `docs/release-notes/v0.1.0.md` | 0.1.0 release notes |
| `docs/release-notes/v0.2.0.md` | 0.2.0 release notes |
| `docs/reports/*.txt` | Historical migration reports |
| `docs/diagrams/mitoribopy_package_flow.mmd` | Existing architecture diagram (v0.2.0 era) |
| `docs/diagrams/generate_mitoribopy_diagram.py` | Diagram generator |
| `docs/diagrams/render_mitoribopy_diagram.sh` | Diagram render helper |
| `docs/reports/v0.3.0/` | **New** — Phase 0 artifacts live here. |

## `examples/`

| Path | Purpose |
|---|---|
| `examples/custom_reference/README.md` | Custom-organism how-to |
| `examples/custom_reference/annotation_template.csv` | Annotation CSV template |
| `examples/custom_reference/codon_tables_template.json` | Codon-table JSON template |

## Findings that change the v0.3.0 plan

1. **`src/mitoribopy/_legacy/` already removed** in v0.2.0 per `CHANGELOG.md`. Phase 1 task 1.1 item "delete `_legacy/` entire directory" is **no-op** — to confirm, `grep -rn _legacy src/ tests/` returns no runtime hits; only doc mentions remain.
2. **Root compatibility wrappers `offset_analysis.py`, `bed_processing.py`, `codon_data.py`, `pipeline_config.py` do NOT exist** at the repo root in this baseline. Phase 1 task 1.1 reduces to the `main.py` decision.
3. **`main.py` at repo root** is a 19-line shim delegating to `mitoribopy.pipeline.runner.main`. README references `python main.py --help`. Decision pending: keep as shim through v0.3.x and remove in v0.4.0, or remove now.
4. **`legacy_key_map` and hidden `--varna` / `--varna_norm_perc` CLI flags are still present** in `src/mitoribopy/config/runtime.py` and `src/mitoribopy/pipeline/runner.py`. `pipeline_config.example.json` still uses the legacy `"varna": false` key. Phase 1 task 1.4 applies only if the rename has fully propagated; **it has not** — the example config needs updating as part of this phase.
5. **`config/defaults.py`, `config/loader.py`, `config/models.py` appear to be vestigial Phase-I scaffolding** (only consumed by `config/__init__.py` re-exports, not used in the runtime pipeline). Candidates for removal; flagged for Ahram before deletion.
6. **`tools/` subpackage is not exercised by the main pipeline** — only `test_subsample.py` consumes it. Safe to leave as a supported utility surface during the refactor.
