# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.0] - Unreleased
### Added
- New subcommand CLI under `src/mitoribopy/cli/`:
  - `mitoribopy align` (Phase 3 stub): FASTQ preprocessing pipeline; flag surface defined, implementation in a later Phase-3 commit.
  - `mitoribopy rpf`: Ribo-seq analysis pipeline from BED inputs (delegates to the existing `mitoribopy.pipeline.runner`).
  - `mitoribopy rnaseq` (Phase 5 stub): DE-table + Ribo-seq integration; flag surface defined.
  - `mitoribopy all` (Phase 6 stub): align + rpf + (optional) rnaseq orchestrator; flag surface defined.
  - Top-level `mitoribopy --help` lists all subcommands; `mitoribopy <sub> --help` shows per-subcommand flags.
- Every subcommand supports shared options: `--config`, `--dry-run`, `--threads N`, and `--log-level {DEBUG,INFO,WARNING,ERROR}`.
- `mitoribopy.config.load_user_config` now auto-detects JSON, YAML (`.yaml`/`.yml`), and TOML (`.toml`) config files. YAML support uses `PyYAML`; TOML uses stdlib `tomllib` on Python 3.11+ and optionally `tomli` on Python 3.10 (install via the `toml-py310` extra).
- `--threads N` exports `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, and `MITORIBOPY_THREADS` so BLAS/OpenMP libraries and external tools invoked by later `align` phases honor the request.
- `--dry-run` prints each subcommand's planned steps and exits 0.
- New dependency: `PyYAML>=6.0`.
- New optional extra: `mitoribopy[toml-py310]` pulling in `tomli` for Python 3.10 users who want TOML configs.

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
