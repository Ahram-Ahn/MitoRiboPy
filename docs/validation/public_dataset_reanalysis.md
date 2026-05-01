# Public-dataset reanalysis recipe

This document is the reproducible smoke recipe for re-running
MitoRiboPy on a small public mt-Ribo-seq dataset. The goal is twofold:

1. give a new user a known-input → known-output exercise that touches
   every stage end to end;
2. provide a reference transcript for the manuscript's "package-level
   standardisation" claims — when the paper says "we reanalysed
   dataset X with workflow Y", this is the recipe Y.

For purely synthetic invariants (no external download, runs in the
test suite), see [`synthetic_mini.md`](synthetic_mini.md).

---

## Pick a dataset

Recommended starter datasets (small, well-documented, mt-relevant):

| Dataset | Organism | Description | Source |
|---|---|---|---|
| GSE161975 | *H. sapiens* | TACO1 KO regression set (already exercised by the in-repo regression test) | NCBI GEO |
| GSE110381 | *S. cerevisiae* | mt-Ribo-seq under heat shock | NCBI GEO |
| GSE127794 | *H. sapiens* | small-RNA mt-RPF series | NCBI GEO |

Use `prefetch` + `fasterq-dump` (or `pyega3` for ENA-staged data) to
download FASTQs to a working directory.

```bash
mkdir -p data/raw
cd data/raw
prefetch SRR12345678
fasterq-dump --threads 4 SRR12345678
gzip *.fastq
```

---

## Build the references

The package needs:

* a bowtie2 index of the mt-transcriptome (the `mt_index` config key);
* a bowtie2 index of contaminant rRNA / tRNA (the `contam_index` key);
* the same transcriptome FASTA you used for `mt_index` (the rpf
  stage's `fasta` key), so the reference checksums match across stages.

```bash
mkdir -p ref
cd ref

# mt-transcriptome (curated copy ships in the repo for h.sapiens; for
# other organisms use Ensembl's transcript FASTA filtered to the MT
# contig).
bowtie2-build human-mt-mRNA.fasta mt_tx

# Contam (rRNA + tRNA, optional but recommended).
bowtie2-build rRNA_tRNA.fasta rrna_contam
```

For a non-human / non-yeast organism, use the custom-organism path
documented in `docs/tutorials/04_custom_organism.md`: pass
`--annotation_file` and `--codon_tables_file` to `mitoribopy rpf`.

---

## Pre-flight the config

```bash
mitoribopy validate-config pipeline_config.yaml
```

`validate-config` will:

* parse the YAML;
* canonicalise legacy keys (any `merge_density:` → `codon_density_window:`,
  `strain: h` → `strain: h.sapiens`, etc.) and report each rewrite;
* check on-disk paths (FASTQs, references, indexes);
* resolve `rnaseq.mode` against the supplied inputs;
* exit `0` on success, `2` on any structural error.

Strict pre-publication runs add `--strict`, which makes legacy-key
rewrites a hard error so the YAML is committed as fully canonical.

---

## Tune with `benchmark`

Before launching a 6-hour cluster job on the full dataset, run a
subsampled benchmark on one or two samples:

```bash
mitoribopy benchmark \
    --config pipeline_config.yaml \
    --output benchmarks/200k \
    --subsample 200000 \
    --threads 8
cat benchmarks/200k/benchmark_summary.md
```

The summary lists per-stage wall time, peak RSS, and the on-disk
footprint — useful for setting `--mem` and `--time` on the actual
cluster job.

---

## Run end to end

```bash
mitoribopy all \
    --config pipeline_config.yaml \
    --output runs/full
```

The orchestrator writes the per-stage outputs under
`runs/full/{align,rpf,rnaseq}/`, plus the canonical run-level files:

* `runs/full/run_manifest.json` — every input hash, tool version,
  and per-stage runtime in one place;
* `runs/full/SUMMARY.md` — human-readable run overview;
* `runs/full/summary_qc.tsv` — per-sample QC roll-up
  (`qc_status` ∈ {`pass`, `warn`}, `qc_notes` lists threshold trips);
* `runs/full/warnings.tsv` — every structured warning emitted during
  the run (e.g. `UMI_INFERRED_NO_DECLARATION` flags a sample whose
  UMI length was inferred from R1 entropy without a sample-sheet
  declaration).

For RNA-seq integration with publication-grade DE, run DESeq2 /
Xtail / Anota2Seq externally on the full transcriptome and feed the
result back via `rnaseq.mode: de_table`:

```yaml
rnaseq:
  rnaseq_mode: de_table
  de_table: external/deseq2_results.tsv
  gene_id_convention: hgnc
  reference_gtf: ref/Homo_sapiens.GRCh38.gtf  # SHA256-gated against rpf
```

The orchestrator hashes `reference_gtf` and refuses to run if the
checksum does not match the one recorded by the `rpf` stage.

---

## Reproduce later

The on-disk artefacts above are sufficient to re-render summaries
without re-executing any stage:

```bash
mitoribopy summarize runs/full
```

…and to re-run the pipeline with confidence that nothing drifted:

```bash
mitoribopy all --config pipeline_config.yaml --output runs/full --resume
```

The hash-validated resume re-checks every input against the prior
manifest; it refuses to skip stages when any of `config_source_sha256`,
`sample_sheet_sha256`, `reference_checksum`, `mitoribopy_version`, or
`schema_version` has drifted, so a stale partial run cannot silently
mix with a new config.

---

## What to put in the manuscript

For a methods-paper reanalysis the pertinent provenance is:

* MitoRiboPy version + git commit (in `run_manifest.json` →
  `mitoribopy_version`, `git_commit`);
* canonical config blob (in `run_manifest.json` → `config_canonical`);
* reference checksum (in `run_manifest.json` → `reference_checksum`);
* tool versions (in `run_manifest.json` → `tools`);
* per-stage runtime + status (in `run_manifest.json` → `stages`);
* output schema versions (in `run_manifest.json` → `output_schemas`).

Quoting these directly from the manifest gives the same provenance
metadata the package itself uses for resume validation.
