# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- GitHub Actions CI workflow to run `pytest` on pushes and pull requests targeting `main`.
- Integration tests with tiny fixture data for offset-selection logic and end-to-end CLI smoke coverage.
- Release metadata polish in `pyproject.toml` (authors, classifiers, URLs, keywords, license metadata).
- Structured release notes directory under `docs/release-notes/`.

## [0.1.0] - 2026-03-30
### Added
- First package-oriented release of MitoRiboPy.
- `mitoribopy` CLI entrypoint and package modules under `src/mitoribopy`.
- P-site/A-site offset selection support with `--offset_site` and `--offset_pick_reference`.
- mt-mRNA-only RPM normalization mode with configurable reference pattern matching.
- Human NT3/NT5 validation runs and migration documentation in `docs/reports/`.
