# MitoRiboPy roadmap

Forward-looking notes that should NOT live in
[`CHANGELOG.md`](../../CHANGELOG.md). The changelog is a record of
what shipped; this file is a record of what is being considered.

Reviewers and downstream tools read the changelog to decide what to
upgrade to. They should not have to wade through speculative entries.
Move items here as soon as they slip out of the active release cycle.

---

## Periodicity QC — statistical hardening

The metagene Fourier bundle shipped in v0.8.0 reports a point estimate
(`spectral_ratio_3nt`) and a four-tier verdict (`snr_call`). Two
companion columns are tracked but not yet emitted:

* `amp_3nt_ci_low` / `amp_3nt_ci_high` — bootstrap CI over genes:
  resample the per-gene tracks with replacement, recompute the
  metagene + DFT, take the percentile interval.
* `permutation_p` — circular-shift permutation null: shift each
  per-gene track by a random offset, recompute the ratio, derive the
  empirical p-value.

Implementation lives next to the score-table builder in
[`src/mitoribopy/analysis/fourier_spectrum.py`](../../src/mitoribopy/analysis/fourier_spectrum.py).
Both columns are added together with the implementation rather than
landing as empty columns first — empty columns in a publication-facing
TSV look unfinished even when documented.

Tests to add together with the implementation:

* synthetic period-3 signal — ratio confidently in the `excellent`
  tier; CI tight; permutation p ≈ 0.
* uniform random noise — ratio in the `broken` tier; permutation p
  uniform.
* circular-shift null sanity — shifted version of a real metagene
  reproduces the same null distribution as the synthetic-noise case.
* small-gene-set behavior — `n_genes < 3` should report `NaN` CI and
  flag a warning rather than emit a misleadingly tight interval.

---

## Templates — single source of truth

Currently `examples/templates/pipeline_config.example.yaml` is the
authoritative exhaustive template, and `mitoribopy all
--print-config-template` prints a curated minimal one. There are two
things to write next to each other.

Plan: extend the orchestrator with `--profile`:

```
mitoribopy all --print-config-template --profile minimal
mitoribopy all --print-config-template --profile publication
mitoribopy all --print-config-template --profile exhaustive
```

Then generate the published YAML templates from these profiles in
docs CI so a flag rename automatically updates the example.

---

## Smoke-test fixture

`branch-mrp/smoke_test/` exists locally as a tiny end-to-end fixture,
but it is not bundled with the package. A future release should ship
`examples/smoke/` containing:

* tiny references (`human_mt_tiny.fa`, `contaminant_tiny.fa`) plus an
  index-build script,
* tiny FASTQs (`WT_Ribo_1.fq.gz`, `KO_Ribo_1.fq.gz`),
* `samples.tsv`, `pipeline_config.smoke.yaml`,
* `expected_outputs_manifest.tsv` listing the files that must exist
  after a clean run, with sizes and content hashes.

Hooks: `pytest -m smoke` runs it on a fresh CI runner; the README
gains a "30-second smoke test" line under Installation.

---

## Manifest schema validation

`run_manifest.json` is the reproducibility artifact every run writes.
There is no JSON schema yet — adding one (and a `tests/test_manifest_schema.py`
that loads a real manifest from `branch-mrp/results/` and validates it)
would prevent a silent field rename from drifting reviewer-facing
provenance.

Minimum fields the schema must require:
`mitoribopy_version`, `manifest_schema_version`, `command`,
`config_source_sha256`, `config_canonical`, `stages`, `input_hashes`,
`output_hashes`, `tool_versions`, `resource_plan`, `warnings_summary`.

---

## Output-contract tests

Per-command tests that assert the exact set of output files produced
and the exact set of columns in each TSV. Three files:

* `tests/test_align_outputs_contract.py`
* `tests/test_rpf_outputs_contract.py`
* `tests/test_rnaseq_outputs_contract.py`

These should run against the smoke fixture above so they don't depend
on a real TACO1-KO run.
