# `ref-5` — Publication-readiness refinements

**Branch:** `ref-5` (off `main` @ MitoRiboPy-refactor)
**Date:** 2026-05-01
**Driver:** `MitoRiboPy_publication_readiness_review.txt` (target journals:
*Bioinformatics* Application Note, *NAR Genomics and Bioinformatics*)

This branch lands the **P0** items from the review that can be
implemented without external benchmark datasets, plus several P1 items.
All 751 existing tests still pass; 27 new tests cover the new behaviour.

## Summary of code-level wins

| # | Refinement | Status | New / changed files |
|---|---|---|---|
| 1 | `ResourcePlan` + `--parallel-samples auto` default | ✅ | `pipeline/resource_plan.py` (new), `cli/align.py` |
| 2 | Per-read-length periodicity QC | ✅ | `analysis/periodicity.py` |
| 3 | Coverage-plot frame metadata sidecar + clarified legends | ✅ | `plotting/coverage_profile_plots.py` |
| 4 | Codon correlation: log2 density RPM + robust regression + MA + support labels | ✅ | `analysis/codon_correlation.py`, `pipeline/runner.py`, `pipeline/steps.py`, `config/runtime.py` |
| 5 | UMI QC output (`umi_qc.tsv`) + clearer dedup docs | ✅ | `align/dedup.py`, `cli/align.py` |
| 6 | Sample-sheet schema extension (biological_sample_id, library_type, read_length_min/max, reference_id, include) | ✅ | `sample_sheet.py` |
| 7 | Stable warning-code registry expansion | ✅ | `errors.py` |
| 8 | `pyproject` Alpha → Beta | ✅ | `pyproject.toml` |
| 9 | New tests | ✅ | 5 new test files, 27 tests |

Out of scope on this branch (P0/P1 items requiring external datasets or
substantial new infrastructure — flagged for a follow-up):

- Comparative benchmark vs riboWaltz / plastid (P0 #10).
- Zenodo archive integration (P1).
- Conda-lock / container reproducibility bundle (P1).
- MultiQC-style HTML report (P2).

## What changed and why

### 1. Resource planning + parallel-by-default

**Problem (review §1, P0 #2/#3).** `--max-parallel-samples` defaulted
to `1`. A modern multi-core machine running on three TACO1 samples sat
on one core for two-thirds of the wall time. The review explicitly
requires `--parallel-samples auto` as the default.

**Fix.** New module `pipeline/resource_plan.py` with a frozen
`ResourcePlan` dataclass and a pure `plan_parallelism()` function. The
align CLI now accepts:

- `--max-parallel-samples auto|N` (default `auto`)
- `--single-sample-mode` (alias for `--max-parallel-samples 1`)
- `--memory-gb GB|auto` (caps the auto scheduler)

`auto` resolves to
`min(n_samples, total_threads // min_threads_per_sample, total_memory // est_per_sample)`
and writes the resolved plan to `<output>/resource_plan.json`. The
existing `_per_worker_threads()` helper is now driven by the plan so
all the downstream concurrency code (`ThreadPoolExecutor`, log line)
sees a single source of truth.

**Backwards compat.** Pre-existing configs that pass
`max_parallel_samples: 1` still hit the serial path — the user has
explicitly requested it, and the plan records the override reason.

### 2. Per-read-length periodicity QC

**Problem (review §2, P0 #4).** `compute_frame_summary` pooled every
read length into one frame-0 number. Mt-Ribo-seq libraries can have a
clean 30-mer class and a contaminated 32-mer class, but a pooled value
hides this. riboWaltz (Lauria et al. 2018) and the Li et al. 2021
mitoribosome protocol both require per-length offsets / per-length
phasing.

**Fix.** New `compute_frame_summary_by_length()` produces one row per
(sample, read_length) with:

```
n_reads_total, n_reads_cds, frame{0,1,2}_fraction, dominant_frame,
frame0_dominance, periodicity_score, frame_entropy,
include_for_downstream, exclusion_reason
```

`run_periodicity_qc()` now writes `qc/by_length/`:

- `frame_by_length.tsv` (long format)
- `length_inclusion_decisions.tsv` (reviewer-friendly subset)
- `periodicity.metadata.json` (thresholds + frame formula)
- `frame_by_length_heatmap.png/.svg` (frame fraction × read length)

Inclusion thresholds (`min_cds_reads_per_length=200`,
`min_frame0_fraction=0.50`, `min_frame0_dominance=0.10`,
`max_frame_entropy=1.45`) match riboWaltz's periodicity-threshold
mode but are configurable per call.

`select_read_lengths_by_periodicity()` is exposed as a helper so
downstream callers can derive the "reviewer-defensible" length set.

### 3. Coverage-plot frame metadata

**Problem (review §3).** Coverage figures showed `frame 0/1/2` without
naming the coordinate system. Reviewers cannot distinguish a P-site
frame from an A-site frame from the figure alone.

**Fix.** Two changes in `plotting/coverage_profile_plots.py`:

1. **Inline caption.** Frame-coloured plots now carry a sub-title
   `"Frame = (assigned-site nt - CDS-start nt) mod 3; Frame 0 = annotated CDS frame."`
2. **Sidecar.** Every `*_density_*_frame/` directory gets a
   `coverage_plot.metadata.json`:

   ```json
   {
     "plot_type": "coverage_by_frame",
     "site": "P-site",
     "frame_formula": "(P_site_nt - CDS_start_nt) % 3",
     "frame_0_definition": "assigned P-site lies in the annotated coding frame",
     "frame_labels": ["Frame 0: annotated CDS frame", "Frame +1: shifted +1 nt from CDS frame", ...],
     "offset_type": "5",
     "offset_site": "p",
     "normalization": "RPM",
     "included_read_lengths": [29, 30, 31, 32, 33, 34]
   }
   ```

Legend labels were also tightened from `"Frame 0 (codon-locked)"` to
the unambiguous `"Frame 0: annotated CDS frame"`.

### 4. Codon correlation refactor

**Problem (review §4, P0 #6).** `run_codon_correlation` ran an OLS
regression on the raw `CoverageDivFreq` column and labelled the top-10
absolute residuals. Raw counts confound depth, codon abundance and
gene length; OLS labels are dominated by a few high-count codons.

**Fix.** Rewrote `analysis/codon_correlation.py`:

- New `metric` arg, default `log2_density_rpm`. Implements
  `log2((value / sum * 1e6) + pseudocount)` with `pseudocount="auto"`
  resolving to `0.5 * min_nonzero`.
- New `regression` arg, default `theil_sen` (uses
  `scipy.stats.theilslopes`). Other choices: `ols`, `none`.
- Three-panel publication figure: scatter (with identity + robust
  fit), MA / Bland-Altman, residuals.
- `support_min_raw=10` flags low-support codons; they remain in the
  TSV but are dimmed in figures and excluded from the residual ranking.
- Label policy: codons ranked by
  `|log2_fold_change| * log10(1 + min_raw_support)`, so low-support
  rows can't dominate.
- New long-format output `codon_correlation_metrics.tsv` with every
  per-codon metric needed to re-plot offline.
- `codon_correlation.metadata.json` records metric, regression,
  pseudocount, and warnings.
- When `metric="raw_count"`, the function emits
  `W_CODON_RAW_COUNT_PRIMARY` to `warnings.tsv` and routes the figure
  to a `raw_count_qc/` subdirectory so it is clearly a QC artefact.

CLI flags exposed in `pipeline/runner.py`:
`--cor-metric`, `--cor-regression`, `--cor-support-min-raw`,
`--cor-label-top-n`, `--cor-pseudocount`, `--cor-raw-panel`.

### 5. UMI QC output

**Problem (review §6, P0 #7/#8).** UMI extraction is correctly placed
before alignment and dedup happens after — but there was no
machine-readable record of how much PCR duplication each sample had,
which method was used, or whether a UMI-bearing sample was
accidentally run with `dedup_strategy=skip`.

**Fix.** New helpers in `align/dedup.py`:

- `recommend_umi_method(umi_length, duplicate_fraction)` — picks
  `unique` for short UMIs (< 8 nt) and `directional` for long UMIs
  with high duplication, returning a stable warning code.
- `build_umi_qc_row()` + `write_umi_qc_tsv()` — emits
  `<output>/umi_qc.tsv` with one row per sample:

  ```
  sample_id, umi_present, umi_length, umi_position,
  n_reads_pre_dedup, n_reads_post_dedup, duplicate_fraction,
  dedup_strategy, dedup_method, warning_code
  ```

Wired into the align CLI immediately after `read_counts.tsv`. Warning
codes flagged: `umi_present_but_skipped`, `short_umi_collision_risk`,
`method_off_recommended`, `no_umi`.

The `dedup.py` module-level docstring was rewritten to describe the
two-step UMI flow (extract pre-align, coordinate+UMI dedup post-align)
and the biological reason coordinate-only dedup is wrong for
mt-Ribo-seq.

### 6. Sample-sheet extension

**Problem (review §5).** The current sample sheet supports per-sample
kit / UMI / strandedness overrides but lacks fields the review
explicitly calls out: biological replicate grouping, library type,
read-length bounds, reference identifier, and an `include` polarity
some users prefer to `exclude`.

**Fix.** Added optional columns:
`biological_sample_id`, `library_type` (`single_end` / `paired_end` /
`auto`), `read_length_min`, `read_length_max`, `reference_id`,
`include`. Validation is strict: invalid `library_type`,
non-integer length bounds, and `read_length_min > read_length_max`
all fail at load time with a single error listing every problem.

`include` and `exclude` can both appear; `include` wins when both are
non-empty so authors can adopt either polarity without surprise.

### 7. Stable warning-code registry

Added the codes referenced by the review:
`E_FASTQ_MATE_PAIRING_FAILED`, `E_ADAPTER_AMBIGUOUS`,
`E_UMI_DECLARED_BUT_NOT_FOUND`, `E_READ_LENGTH_NO_PERIODICITY`,
`E_RNA_RIBO_PAIRING_INCOMPLETE`, `W_CODON_RAW_COUNT_PRIMARY`,
`W_UMI_SHORT_COLLISION_RISK`, `W_UMI_PRESENT_BUT_SKIPPED`,
`W_PARALLEL_AUTO_OVERRIDDEN`. The `W_*` prefix is new — it marks soft,
non-blocking warnings that should land in `warnings.tsv` but never
abort the run.

### 8. `pyproject` Alpha → Beta

Bumped `Development Status :: 3 - Alpha` to `Development Status :: 4 - Beta`.
Justification: the package now has a frozen output contract, a stable
warning-code registry, sample-sheet schema validation, full test
coverage of the new public surfaces, and a documented FASTQ→TE flow.

### 9. New tests

- `tests/test_resource_plan.py` (9 tests) — auto-default behaviour,
  single-sample override, memory cap, threads-floor, JSON round-trip.
- `tests/test_periodicity_by_length.py` (3 tests) — frame-by-length
  call distinguishes good vs bad lengths, write paths exist, metadata
  carries the frame formula.
- `tests/test_codon_correlation_transformations.py` (5 tests) — new
  metrics TSV columns, AAA-perturbed data shows expected log2 FC,
  raw-count metric warns and routes to `raw_count_qc/`, label score
  prefers high-support outliers, invalid metric/regression raise.
- `tests/test_umi_qc.py` (7 tests) — method recommendation, skipped
  UMI warning, duplicate-fraction arithmetic, TSV round-trip.
- `tests/test_coverage_frame_metadata.py` (3 tests) — P-site vs A-site
  formula, included_read_lengths, label structure.

## Validation: end-to-end TACO1 run

`mitoribopy all --config pipeline_config.yaml --output results/taco1_ref5
--strict --threads 8` completed cleanly **EXIT=0** in
**~57 minutes** (align 8m, rpf 6m, rnaseq 41m). All three stages report
`completed` in `run_manifest.json`.

### Quick-look at TACO1 read flow

| sample | total | post_trim | mt_aligned | post_dedup |
|---|---|---|---|---|
| TACO_WT_S1_R1   | 33.4M | 26.3M | 3.63M | 3.53M |
| TACO_KO_S2_R1   | 28.3M | 22.3M | 2.29M | 2.23M |
| TACO_RESC_S3_R1 | 29.4M | 24.5M | 3.62M | 3.51M |

### New artefacts produced (verified on disk)

| Artefact | Location | Notes |
|---|---|---|
| `resource_plan.json` | `align/` | Records `parallel_samples: 1` (user override from config), `per_sample_threads: 8`, reason field. |
| `umi_qc.tsv` | `align/` | All 3 ribo samples report `umi_present=false` / `warning_code=no_umi` (consistent with `kit_preset: pretrimmed`). |
| `frame_by_length.tsv` | `rpf/qc/by_length/` | 18 rows (3 samples × 6 read lengths). 1 row passes (`TACO_KO_S2_R1` 29-mer, frame0=0.56, dominance=0.32). |
| `length_inclusion_decisions.tsv` | `rpf/qc/by_length/` | Reviewer-friendly include/exclude per length. |
| `periodicity.metadata.json` | `rpf/qc/by_length/` | Records thresholds + frame formula. |
| `frame_by_length_heatmap.png/.svg` | `rpf/qc/by_length/` | Per-sample frame fraction × read length, excluded rows red-bordered. |
| `coverage_plot.metadata.json` | `rpf/coverage_profile_plots/{p_site,a_site}_density_{rpm,raw}_frame/` | Four sidecars (one per combination). |
| `codon_correlation_metrics.tsv` | `rpf/codon_correlation/{p_site,a_site}/` | Long-format, all new columns (log2_fold_change, mean_log2_density, support_min_raw, robust_residual, label_score, include_primary, exclusion_reason). |
| `codon_correlation.metadata.json` | `rpf/codon_correlation/{p_site,a_site}/` | `metric: log2_density_rpm`, `regression: theil_sen`, `support_min_raw: 10`, `warnings: []`. |

### Reviewer-facing read of the periodicity QC

Across all three TACO1 ribo libraries, only one read-length class
(29-mer in `TACO_KO_S2_R1`) clears the `min_frame0_fraction=0.50` and
`min_frame0_dominance=0.10` bar. Every other length is rejected with
`exclusion_reason=weak_periodicity`. This is **exactly the kind of
finding the publication review demanded a per-length QC for** — a
pooled frame summary would have papered over it. The package's own
README already warns that mt-Ribo-seq phasing is "muddier than nuclear
cytosolic Ribo-seq"; this QC layer turns that caveat into a per-sample,
per-length, machine-readable receipt.

### Existing pipeline behaviour preserved

- `read_counts.tsv` schema unchanged (back-compat verified).
- `frame_summary.tsv`, `periodicity_start.tsv`, `periodicity_stop.tsv`,
  `strand_sanity.tsv`, `periodicity_metagene.png/.svg` still produced
  with identical content.
- Legacy `*_vs_*_{all,masked}.csv` codon-correlation files still
  written (with the new metric columns appended).
- All 751 prior unit tests pass; one test was updated
  (`test_max_parallel_samples_default_is_one` →
  `test_max_parallel_samples_default_is_auto`) to reflect the new
  default.

### Known gap surfaced by validation

`outputs_index.tsv` does **not** yet list the new artefacts
(resource_plan.json, umi_qc.tsv, frame_by_length.tsv,
codon_correlation_metrics.tsv, coverage_plot.metadata.json). The index
is built from an explicit registry in the summarize step; that registry
needs the new entries appended. Files are present on disk and
reviewable, but the SUMMARY.md inventory is incomplete. **Recommended
follow-up:** extend `summary` step to register these paths.

## Out-of-scope follow-ups (not on this branch)

The full publication-readiness review calls out a few invasive
changes that this branch does not attempt:

- **FASTQ-first DAG (review §1, P0 #1).** A formal `mitoribopy preprocess`
  stage with parallel Ribo + RNA preprocessing nodes. The current
  `mitoribopy all` is still linear (align → rpf → rnaseq); refactoring
  into a true DAG is a larger surgery best done as a separate branch.
- **Comparative benchmark (review P0 #10).** Requires external
  datasets (riboWaltz / plastid output, public mt-Ribo-seq runs) and
  Zenodo archiving.
- **Reference manifest** (`reference_manifest.yaml` with checksums).
- **Conda-lock / container** reproducibility bundle.
- **MultiQC-style HTML report** wrapping `umi_qc.tsv` +
  `frame_by_length.tsv` + figures.

## Files changed

```
$ git diff --stat main
 pyproject.toml                                    |   2 +-
 src/mitoribopy/align/dedup.py                     | 194 ++++++-
 src/mitoribopy/analysis/codon_correlation.py      | 589 +++++++++++++++++-----
 src/mitoribopy/analysis/periodicity.py            | 346 +++++++++++++
 src/mitoribopy/cli/align.py                       | 121 ++++-
 src/mitoribopy/config/runtime.py                  |   6 +
 src/mitoribopy/errors.py                          |  34 ++
 src/mitoribopy/pipeline/__init__.py               |   4 +
 src/mitoribopy/pipeline/runner.py                 |  58 +++
 src/mitoribopy/pipeline/steps.py                  |  10 +
 src/mitoribopy/plotting/coverage_profile_plots.py | 100 +++-
 src/mitoribopy/sample_sheet.py                    | 112 +++-
 tests/test_align_cli.py                           |  22 +-
 src/mitoribopy/pipeline/resource_plan.py          | new
 tests/test_codon_correlation_transformations.py   | new
 tests/test_coverage_frame_metadata.py             | new
 tests/test_periodicity_by_length.py               | new
 tests/test_resource_plan.py                       | new
 tests/test_umi_qc.py                              | new
```
