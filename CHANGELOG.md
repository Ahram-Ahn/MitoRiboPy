# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
