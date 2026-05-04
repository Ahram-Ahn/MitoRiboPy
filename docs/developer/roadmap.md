# MitoRiboPy roadmap

Forward-looking notes that should NOT live in
[`CHANGELOG.md`](../../CHANGELOG.md). The changelog records shipped
user-visible behavior; this file tracks remaining maintenance work.

The following items are intentionally limited to work that is not
already implemented in `src/mitoribopy/`.

---

## Benchmark Reference Rows

`mitoribopy benchmark` is implemented and writes `benchmark.tsv` plus
`benchmark_summary.md`. The missing piece is a measured reference table
for the datasets named in [`docs/benchmarking.md`](../benchmarking.md):
smoke synthetic, public human 1M, public yeast 1M, TACO1 WT/KO, and one
non-subsampled lab-realistic run.

Before publishing numbers:

* run each case on the chosen reference hardware;
* record the exact MitoRiboPy version, commit, external tool versions,
  thread count, and whether FASTQs were subsampled;
* keep measured results out of the docs until every table cell can be
  backed by a real run artifact.

## Template Single Source

The CLI already supports:

```bash
mitoribopy all --print-config-template --profile minimal
mitoribopy all --print-config-template --profile publication
mitoribopy all --print-config-template --profile exhaustive
```

The remaining work is to stop maintaining the long YAML examples by
hand. Generate `examples/templates/*.yaml` from the same profile data
used by the CLI, then add a docs check that fails when a flag is added
without updating the generated examples.

## Wheel-Friendly Examples

`--profile exhaustive` reads
`examples/templates/pipeline_config.example.yaml` when the checkout
layout is available and falls back to the minimal template when a wheel
installation does not include `examples/`. Decide whether the exhaustive
template should be packaged as package data, or whether the CLI should
generate the exhaustive profile entirely from parser metadata.

## Benchmark Sample-Sheet Subsampling

`mitoribopy benchmark --subsample N` rewrites explicit FASTQ paths in
`align.fastq`, `align.fastq_dir`, `rnaseq.rna_fastq`, and
`rnaseq.ribo_fastq`. It does not yet rewrite FASTQ paths inside a
top-level `samples:` sheet. Add sample-sheet-aware subsampling before
recommending `--subsample` as the default tuning path for all
orchestrator configs.
