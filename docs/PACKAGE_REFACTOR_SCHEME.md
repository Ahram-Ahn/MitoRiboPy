# Python Package Refactor Scheme

This document is a historical design note from the package migration. It is
kept for project history, but parts of the original plan now differ from the
current `0.2.0` package layout.

## Status Note (2026-04-07)

- The package now runs through `mitoribopy.pipeline.runner` and
  `mitoribopy.pipeline.steps` without a runtime dependency on `_legacy`.
- The former `inframe_analysis` workflow now lives under
  `mitoribopy.analysis.translation_profile_analysis`.
- The former IGV-style plotting module is now
  `mitoribopy.plotting.coverage_profile_plots`.
- The former VARNA export module is now
  `mitoribopy.plotting.structure_density_export`.
- Built-in codon tables and annotations are packaged as JSON and CSV data files
  under `src/mitoribopy/data/`.

## Current Status (2026-03-29)
- Phase I: completed (package scaffold + CLI forwarding).
- Phase II slice 1: completed for offset analysis modules with compatibility wrappers.
- Phase II slice 2: completed for BED/read-count I/O modules with compatibility wrappers.
- Phase II slice 3: completed for config and codon/annotation data modules with compatibility wrappers.
- Phase II slice 4: completed for in-frame/RNA/IGV/VARNA/correlation modules and tool migration wrappers.
- Phase III: completed (readability pass on package modules, consistent codon-axis ordering, and reduced script side effects).

## Project Goal
Turn this mitochondrial ribosome profiling codebase into a public, installable, well-documented Python package that is:
- Easy to read and modify
- Scientifically transparent
- Stable for users
- Clear in naming and structure

---

## 1) Package Identity and Scope

### Proposed package names
- Distribution name (PyPI): `MitoRiboPy`
- Import name: `mitoribopy`
- CLI command: `mitoribopy`

These names are short, descriptive, and easy to remember.

### Primary user workflows
1. Run full mitochondrial Ribo-seq pipeline from BED + FASTA + read-count files.
2. Run selected modules only (offsets, in-frame analysis, IGV plots, VARNA export, codon correlation).
3. Use reproducible config files and command-line overrides.

---

## 2) Target Directory Layout (Final State)

```text
project-root/
  pyproject.toml
  README.md
  LICENSE
  CHANGELOG.md
  .gitignore

  src/
    mitoribopy/
      __init__.py
      __main__.py
      cli.py

      config/
        __init__.py
        defaults.py
        models.py
        loader.py
        validators.py

      data/
        __init__.py
        codon_tables.py
        transcript_annotations.py

      io/
        __init__.py
        bed_reader.py
        fasta_reader.py
        read_counts.py
        filesystem.py

      analysis/
        __init__.py
        offset_enrichment.py
        offset_selection.py
        site_projection.py
        footprint_density.py
        frame_usage.py
        codon_usage.py
        codon_correlation.py
        rna_seq.py

      plotting/
        __init__.py
        style.py
        offset_plots.py
        footprint_plots.py
        igv_plots.py
        varna_export.py
        correlation_plots.py

      pipeline/
        __init__.py
        context.py
        runner.py
        steps.py

      tools/
        __init__.py
        subsample.py

      utils/
        __init__.py
        logging.py
        typing.py
        text.py

  tests/
    unit/
    integration/
    regression/

  examples/
    configs/
      human_stop_psite.json
      human_stop_asite.json
      yeast_start_psite.json
    commands/
      README.md

  docs/
    PACKAGE_REFACTOR_SCHEME.md
    CONTRIBUTING.md
    API.md
    SCIENTIFIC_ASSUMPTIONS.md
```

---

## 3) Old -> New File Mapping

### Entry and config
- `main.py` -> `src/mitoribopy/cli.py` + `src/mitoribopy/pipeline/runner.py`
- `pipeline_config.py` -> `src/mitoribopy/config/defaults.py` + `models.py` + `loader.py`
- `pipeline_config.example.json` -> `examples/configs/*.json`

### Core analysis
- `offset_analysis.py` -> `analysis/offset_enrichment.py` + `analysis/offset_selection.py`
- `inframe_analysis.py` -> `analysis/footprint_density.py` + `analysis/frame_usage.py` + `analysis/codon_usage.py`
- `codon_correlation.py` -> `analysis/codon_correlation.py`
- `rna_seq_analysis.py` -> `analysis/rna_seq.py`

### IO/data
- `bed_processing.py` -> `io/bed_reader.py` + `io/read_counts.py`
- `codon_data.py` -> `data/codon_tables.py` + `data/transcript_annotations.py`

### Plotting
- `plot_style.py` -> `plotting/style.py`
- `igv_style_plot.py` -> `plotting/igv_plots.py`
- `varna_plot.py` -> `plotting/varna_export.py`

### Tools
- `subsample_bed.py` -> `tools/subsample.py` (plus CLI subcommand)
- `aligned_ratio.py` -> optional `tools/read_composition.py` or move to `examples/legacy/`

---

## 4) Naming Rules (Human-readable, Toggle-friendly)

### File names
- Use `snake_case`, and include domain intent:
  - Good: `offset_selection.py`, `transcript_annotations.py`
  - Avoid vague names: `analysis2.py`, `utils_misc.py`

### Function names
- Verb + object + context when helpful:
  - `compute_offset_enrichment_table`
  - `select_offsets_with_tiebreak`
  - `build_footprint_density_table`

### Variable names
- Explicit and scientific:
  - `read_length_nt` (not `rl`)
  - `offset_5p_nt` / `offset_3p_nt`
  - `normalization_mode`
  - `selected_site`

### Toggle semantics
- All toggles should be enumerated and documented in one place.
- Avoid hidden toggles in function internals.
- Prefer config dataclass fields with clear defaults.

---

## 5) Script/Function Explanation Standard

Every public module and function must include:

1. **Module docstring**
- What this module does
- Key assumptions
- Input/output files or data structures
- Small usage example

2. **Function docstring (Google-style)**
- Purpose
- Args (type + meaning)
- Returns
- Raises
- Notes/assumptions
- Example call

3. **Inline comments**
- Explain *why* for non-obvious logic
- Avoid noisy comments for obvious lines

4. **Type hints**
- Required for public functions and dataclasses

---

## 6) README Requirements (Public-ready)

`README.md` should include all of these sections:
1. Project overview
2. Scientific scope and assumptions
3. Installation (pip and editable)
4. Quick start (human and yeast examples)
5. Full CLI usage with option explanations
6. Input file formats (BED, FASTA, read-count file)
7. Output directory structure and key files
8. P-site vs A-site semantics
9. Offset selection strategy and tie-break behavior
10. Troubleshooting and common mistakes
11. Reproducibility (versioning/config/seed)
12. Developer guide and contribution steps

---

## 7) Phased Execution Plan

### Phase 0: Baseline freeze and regression references
- Keep current behavior as baseline.
- Create small deterministic test fixtures (subsampled NT3/NT5).
- Record expected key outputs:
  - selected offsets table
  - codon totals (row counts, key stats)
  - frame usage totals

**Exit criteria**
- Baseline regression checks are reproducible.

### Phase 1: Packaging skeleton (no behavior changes)
- Add `pyproject.toml`.
- Create `src/mitoribopy` package structure.
- Add console entrypoint `mitoribopy`.
- Keep compatibility wrappers so old commands still work.

**Exit criteria**
- `pip install -e .` works.
- `mitoribopy --help` works.

### Phase 2: Module split and renaming (behavior-preserving)
- Move code into target module locations.
- Keep old module files as temporary thin wrappers.
- Add deprecation warnings for old import paths.

**Exit criteria**
- Existing workflows produce same outputs within tolerance.
- Import paths are stable in new package layout.

### Phase 3: Code cleanup and readability pass
- Rename unclear variables.
- Reduce oversized functions into focused units.
- Add/upgrade docstrings and targeted comments.

**Exit criteria**
- Public functions fully documented.
- No giant “do everything” functions remain.

### Phase 4: Documentation completion
- Rewrite `README.md` into complete public documentation.
- Add `docs/API.md`, `docs/SCIENTIFIC_ASSUMPTIONS.md`, `docs/CONTRIBUTING.md`.

**Exit criteria**
- New user can run pipeline from docs only.

### Phase 5: Tests and quality tooling
- Add unit tests for offset logic and site transforms.
- Add integration tests for mini datasets.
- Add formatter/lint/type checks in CI.

**Exit criteria**
- CI passes on each commit.
- Core scientific outputs are regression-protected.

---

## 8) Quality Gates for Public Release

Before first package release:
- [ ] Package installs cleanly (`pip install .`)
- [ ] CLI and config docs complete
- [ ] Human and yeast example runs validated
- [ ] P-site/A-site and offset-picking semantics documented and tested
- [ ] Minimum test coverage on critical logic (offset selection, site transforms, codon aggregation)
- [ ] Version tag and changelog entry prepared

---

## 9) Suggested First Release Strategy

- Version `0.1.0`: package skeleton + migrated modules + complete docs + baseline tests.
- Version `0.2.0`: refined API, subcommands, enhanced plotting/export options.
- Version `1.0.0`: stable public API and tested workflow guarantees.

---

## 10) Immediate Next Step

Start with **Phase 1 + Phase 0 fixtures** in one focused change set:
1. Add package skeleton + CLI entrypoint.
2. Add regression mini fixtures from NT3/NT5.
3. Verify old and new entrypoints give matching outputs on fixture set.

This gives a safe foundation to perform larger cleanup without scientific drift.
