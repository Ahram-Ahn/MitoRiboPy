# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.0] - Unreleased
### Added
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
