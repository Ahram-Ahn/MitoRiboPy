<p align="center">
  <img src="docs/banner.png" alt="MitoRiboPy — A Python package for mt-Ribo-seq" width="100%">
</p>

# MitoRiboPy

Mitochondrial ribosome profiling (mt-Ribo-seq) analysis, end to end.

MitoRiboPy is a Python package + CLI for analysing mt-Ribo-seq data from raw FASTQ all the way through translation-efficiency integration with paired RNA-seq. Every per-sample decision (kit, dedup, offsets) is independent, so mixed-library batches just work.

The package is built around four pipeline subcommands plus six utilities. Pass `--help` on any of them for the full flag list, or see [docs/reference/cli.md](docs/reference/cli.md).

| Subcommand | What it does | Typical use |
|---|---|---|
| `mitoribopy align` | FASTQ → BAM → BED6 (cutadapt + bowtie2 + umi_tools + pysam) | "I only want trimming + alignment" |
| `mitoribopy rpf` | BED/BAM → offsets, translation profile, codon usage, coverage plots, metagene Fourier QC | "I already have aligned BEDs" |
| `mitoribopy rnaseq` | Translation efficiency from paired RNA-seq + Ribo-seq. `--rna-fastq` runs trimming → bowtie2 → counting → pyDESeq2 → TE/ΔTE end-to-end (needs `pip install 'mitoribopy[fastq]'`); `--de-table` accepts an external DE table + a prior rpf run (SHA256-gated). | "TE / ΔTE" |
| `mitoribopy all` | End-to-end orchestrator (align + rpf + optional rnaseq from one YAML); writes a composed `run_manifest.json`. `--resume` is hash-validated against the prior manifest. | "I have raw FASTQ and want everything" |
| `mitoribopy validate-config` | Parse + canonicalise + check paths + validate `rnaseq.mode`. Exit 0 / 2. | Before long cluster jobs |
| `mitoribopy validate-reference` | Pre-flight a custom mt-transcriptome FASTA + annotation pair. Exit 0 / 2. | Custom organisms |
| `mitoribopy validate-figures` | Mechanically QC every plot under a finished run; writes `figure_qc.tsv`. Exit 0 / 1 / 2 (`--strict` upgrades warn → fail). | After a run |
| `mitoribopy periodicity` | Standalone metagene Fourier QC on a saved site table. | Re-tune Fourier window without re-running offsets |
| `mitoribopy summarize` | Regenerate `SUMMARY.md` from `run_manifest.json`. Auto-invoked by `all`. | Re-render an old run's summary |
| `mitoribopy benchmark` | Time + RSS + disk for a `mitoribopy all` invocation. `--subsample N` reservoir-samples each FASTQ. | Cluster sizing |

---

## Table of contents

1. [What MitoRiboPy is for](#what-mitoribopy-is-for)
2. [Pipeline overview](#pipeline-overview)
3. [Installation](#installation)
4. [Quick start](#quick-start)
5. [Inputs you need to prepare](#inputs-you-need-to-prepare)
6. [How to run — YAML vs shell wrapper](#how-to-run--yaml-vs-shell-wrapper)
7. [Strain presets, custom organisms, and footprint classes](#strain-presets-custom-organisms-and-footprint-classes)
8. [Subcommand reference](#subcommand-reference)
9. [Further reading](#further-reading)
10. [What the numbers mean — RNA, RPF, TE, ΔTE](#what-the-numbers-mean--rna-rpf-te-%CE%B4te)
11. [Output overview](#output-overview)
12. [Examples](#examples)
13. [Tools](#tools)
14. [Logs and provenance](#logs-and-provenance)
15. [Development](#development)
16. [Citation](#citation)
17. [Known limitations](#known-limitations)
18. [License](#license)

## What MitoRiboPy is for

MitoRiboPy is a focused tool for the 13 mt-mRNAs of human mitochondria (or 8 mt-mRNAs in yeast, with configurable codon tables for any other mitochondrion). It ships:

- **Per-sample adapter detection** with auto-fallback to the right kit. Mixed-kit and mixed-UMI batches resolve each sample independently. Pre-trimmed FASTQs (e.g. SRA-deposited data) are auto-detected and routed through cutadapt with no `-a` flag.
- **Per-sample offset selection** so inter-sample drift in the canonical 12–15 nt 5' P-site offset doesn't bias your downstream codon-usage tables. A combined-across-samples diagnostic is still emitted and an `offset_drift_<align>.svg` plot makes drift visible at a glance.
- **Both P-site and A-site downstream outputs** by default, side by side under per-site subdirectories. No ambiguity about which output corresponds to which site.
- **End-to-end RNA-seq + Ribo-seq → TE / ΔTE in one subcommand.** `mitoribopy rnaseq` takes raw FASTQs and a transcriptome FASTA and runs trimming, bowtie2 alignment, per-transcript counting, and pyDESeq2 itself before emitting TE, ΔTE, and a six-figure plot set. (Bringing your own pre-computed DE table from R / Python remains supported via `--de-table` and enforces a SHA256 reference-consistency gate.)
- **Strain-aware defaults**: built-in human (`-s h.sapiens`) and yeast (`-s s.cerevisiae`) annotations + codon tables, plus `custom` for any other organism with a published NCBI Genetic Code (mouse, fly, plants, fungi, ...).

What MitoRiboPy is **not**:

- Not a general-purpose nuclear Ribo-seq pipeline. The defaults, references, and dedup heuristics are calibrated for the low-complexity 13-mRNA mt universe.
- Not a general DE engine for nuclear genes. The default `mitoribopy rnaseq` flow runs **pyDESeq2 on the mt-mRNA subset only**, which is fine for exploring mt-translation efficiency end-to-end on a single library. For publication-grade DE statistics across the full transcriptome run DESeq2 / Xtail / Anota2Seq externally and pass the resulting table via `--de-table`.

---

## Pipeline overview

![Pipeline overview](docs/diagrams/01_pipeline_overview.png)

UMI handling: the UMI is extracted into the read QNAME during the cutadapt trim step (5' single-pass or 3' two-pass), so it travels through bowtie2 alignment unchanged and is available for `umi_tools dedup` after the MAPQ filter — the only stage that needs alignment coordinates AND the UMI together.

### Detailed stage diagrams

- [docs/diagrams/02_align_stage.png](docs/diagrams/02_align_stage.png) — internals of `mitoribopy align`: per-sample resolution → cutadapt + UMI → contam subtract → bowtie2 → MAPQ → dedup → BED6.
- [docs/diagrams/03_rpf_stage.png](docs/diagrams/03_rpf_stage.png) — internals of `mitoribopy rpf`: filter BED → offset enrichment + selection → translation_profile + coverage_profile_plots + optional modules.
- [docs/diagrams/04_rnaseq_stage.png](docs/diagrams/04_rnaseq_stage.png) — internals of the optional `mitoribopy rnaseq` stage: DE table + rpf_counts → SHA256 reference gate → TE → ΔTE → scatter + volcano.

Regenerate with `python docs/diagrams/render_diagrams.py` (matplotlib only; no Node / mermaid-cli required).

---

## Installation

The README and `examples/templates/` describe the current **v0.7.1** interface. Verify with `mitoribopy --version`. See [CHANGELOG.md](CHANGELOG.md) for the consolidated v0.7.0 release notes — this is the publication-readiness release: aggregate-then-DFT metagene Fourier QC with bootstrap CI + circular-shift permutation null, per-gene unit-mean metagene aggregation by default (legacy depth-weighted sum still available behind `normalize="none"`), nulled Wald p-values in the from-FASTQ rnaseq mode, JSON Schema for `run_manifest.json`, per-(sample, transcript) strand-sanity audit, and a "Periodicity statistical confidence" table in `SUMMARY.md`.

### From source (recommended)

```bash
git clone https://github.com/Ahram-Ahn/MitoRiboPy.git
cd MitoRiboPy
git checkout v0.7.1          # current published version; omit for HEAD
python -m pip install -e .
mitoribopy --version          # MUST print 0.7.1 or later
```

This pulls every Python dependency (`numpy`, `pandas`, `matplotlib`, `seaborn`, `biopython`, `scipy`, `PyYAML`, `pysam`) automatically. The external bioinformatics tools (`cutadapt`, `bowtie2`, `umi_tools`, …) still need to be on `$PATH` separately — see [External tool dependencies](#external-tool-dependencies) below.

For development and tests, add the dev extras:

```bash
python -m pip install -e ".[dev]"
```

### From PyPI

```bash
python -m pip install 'mitoribopy>=0.7.1'
```

The package is published on PyPI: [pypi.org/project/mitoribopy](https://pypi.org/project/mitoribopy/). Pin the lower bound (`>=0.7.1`) so a stale PyPI cache cannot install a pre-publication-freeze build.

### Verify the install

```bash
mitoribopy --version
mitoribopy --help
```

If you prefer not to install at all:

```bash
PYTHONPATH=src python -m mitoribopy --help
```

### 30-second smoke test

The repo ships a tiny end-to-end fixture under [`examples/smoke/`](examples/smoke/) so a fresh install can be verified in one command. Three synthetic mt-mRNAs, two samples, no UMIs, no contaminants — designed to exercise every stage (cutadapt → bowtie2 → BED → offsets → translation profile → coverage → metagene Fourier QC) without external data:

```bash
cd examples/smoke
python generate_smoke_fastqs.py    # writes *.fastq.gz + bowtie2 index
mitoribopy all --config pipeline_config.smoke.yaml --output results/
```

Expected wall-clock: 10–30 s on a 2024 laptop. Every file in [`examples/smoke/expected_outputs.txt`](examples/smoke/expected_outputs.txt) must exist + be non-empty after the run; the same assertion runs in CI under `pytest -m smoke` when the external tools are present.

### External tool dependencies

MitoRiboPy shells out to a small set of standard bioinformatics tools. All of them must be on `$PATH` for a real run:

| Tool | Used by | Required when |
|---|---|---|
| `cutadapt` | `align`, `rnaseq` (from-FASTQ) | always (length + quality filter even for pre-trimmed data) |
| `bowtie2` + `bowtie2-build` | `align`, `rnaseq` (from-FASTQ) | always |
| `umi_tools` | `align` | at least one sample's resolved kit has UMIs |
| `pysam` (Python lib) | `align`, `rpf`, `rnaseq` (from-FASTQ) | always (installed automatically via `pip`) |
| `samtools` | optional | recommended for inspecting outputs; not required |
| `pydeseq2` (Python lib) | `rnaseq` (from-FASTQ) | only when running `mitoribopy rnaseq` in from-FASTQ mode (`--rna-fastq …`); install via the `[fastq]` extra: `pip install 'mitoribopy[fastq]'`. The pre-computed-DE flow does not need it. |

The bioconda environment under [docs/environment/environment.yml](docs/environment/environment.yml) installs everything in one command:

```bash
conda env create -f docs/environment/environment.yml
conda activate mitoribopy
```

---

## Quick start

The shortest path from raw FASTQ to translation-profile + coverage outputs is one YAML file plus one command.

```bash
# 1. Start from one of two templates next to your data and fill in the paths:
#    -- examples/templates/ ships an EXHAUSTIVE template that lists every
#       available flag with its default value and a 1-line comment:
cp examples/templates/pipeline_config.example.yaml pipeline_config.yaml
#    -- OR get the curated MINIMAL template from the CLI:
mitoribopy all --print-config-template --profile minimal > pipeline_config.yaml
#       Other profiles: --profile publication, --profile exhaustive

$EDITOR pipeline_config.yaml

# 2. Optional: dry-run prints the per-stage argv so you can review.
mitoribopy all --config pipeline_config.yaml --output results/ --dry-run

# 3. Run.
mitoribopy all --config pipeline_config.yaml --output results/
```

The matching shell-script templates are at [examples/templates/run_align.example.sh](examples/templates/run_align.example.sh), [examples/templates/run_rpf.example.sh](examples/templates/run_rpf.example.sh), and [examples/templates/run_pipeline.example.sh](examples/templates/run_pipeline.example.sh) — pick those if you prefer per-stage commands you can split across cluster jobs.

### Publication-safe recipe

Use this recipe for any run that backs a manuscript, preprint, or shared dataset. **One switch — `mitoribopy all --strict` — turns on every publication-readiness gate** (config preflight, align strict-publication-mode, post-run figure QC, warning promotion). Each step also leaves a self-auditing artifact you can drop into a methods section.

```bash
# 0. Pin the manuscript version
python -m pip install 'mitoribopy>=0.7.1'
mitoribopy --version

# 1. Optional pre-flights (also run automatically by --strict below)
mitoribopy validate-config pipeline_config.yaml --strict
mitoribopy validate-reference \
    --fasta references/human_mt_transcriptome.fa \
    --strain h.sapiens

# 2. Dry-run to inspect the per-stage commands the orchestrator will issue
mitoribopy all --config pipeline_config.yaml --output results/full_run --dry-run

# 3. Run, publication-safe, with structured progress events for audit logs
mitoribopy all \
    --config pipeline_config.yaml \
    --output results/full_run \
    --strict \
    --progress jsonl \
    --progress-file results/full_run/progress.jsonl
```

That single `--strict` invocation:

* runs `validate-config --strict` up-front and aborts before any stage if a deprecated key, unknown key, or missing input is found,
* forwards `--strict-publication-mode` into the align stage so non-default policies that would invalidate a publication run fail-fast,
* writes `<output>/canonical_config.yaml` (the fully-resolved config the run actually executed — auto-wiring + sample-sheet expansion + rnaseq-mode resolution applied) so a reviewer can diff it against your input config,
* runs `validate-figures --strict` after the pipeline finishes, promoting warn-only QC findings to fail in `figure_qc.tsv`,
* still emits `SUMMARY.md`, `outputs_index.tsv`, `warnings.tsv`, and `progress.jsonl` regardless of the strictness gates.

For TE / ΔTE, the publication route is `--rnaseq-mode de_table`: run a full-transcriptome DESeq2 / Xtail / Anota2Seq externally and feed the table back. The in-tree from-FASTQ path runs pyDESeq2 on the **mt-mRNA subset only** and is exploratory; n=1 designs fail-fast unless you pass `--allow-pseudo-replicates-for-demo-not-publication`. See [`docs/rnaseq_te.md`](docs/rnaseq_te.md) for the full publication-boundary reference (mode comparison, strict-mode gates, output files, when to use which).

```yaml
# Publication TE route (place under `rnaseq:` in your pipeline_config.yaml).
rnaseq:
  rnaseq_mode: de_table
  de_table: external_full_transcriptome_deseq2.tsv
  reference_gtf: references/gencode_or_refseq.gtf   # SHA-gated against rpf
  gene_id_convention: hgnc                          # or ensembl / refseq / mt_prefixed / bare
  base_sample: WT
  compare_sample: KO
```

A minimal `pipeline_config.yaml` for a typical human mt-Ribo-seq run looks like this (annotated). **The recommended idiom is a top-level `samples:` block** — a single TSV that names every Ribo-seq and (optional) RNA-seq FASTQ once, and is auto-wired into both the align and rnaseq stages. `align.fastq:` is still accepted as a standalone shortcut for a single-stage run.

```yaml
# Top-level: unified per-project sample sheet (recommended). The same
# TSV drives both align (Ribo-seq rows) and rnaseq (RNA-seq rows). See
# `mitoribopy --help` and docs/reference/sample_sheet_schema.md for
# the column reference. Pairings between Ribo-seq and RNA-seq are by
# sample_id, never by index.
samples:
  table: samples.tsv

align:
  # Adapter handling. Auto-detection is the default since v0.7.1 — the
  # pipeline scans the head of every FASTQ, picks the matching adapter
  # family, and reports the result in kit_resolution.tsv. Pin `adapter:`
  # only when detection cannot identify the library; set `pretrimmed:`
  # for already-trimmed FASTQs.
  # adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
  # pretrimmed: false
  adapter_detection: auto         # auto | off | strict
  library_strandedness: forward
  contam_index: input_data/indexes/rrna_contam
  mt_index: input_data/indexes/mt_tx
  mapq: 10
  min_length: 15
  max_length: 45
  dedup_strategy: auto            # umi_coordinate per sample if UMI, else skip

rpf:
  strain: h.sapiens               # human mt-mRNA reference + codon table
  fasta: input_data/human-mt-mRNA.fasta
  footprint_class: monosome       # short | monosome | disome | custom
                                  # monosome's biological default window:
                                  #   h.sapiens 28-34 nt; s.cerevisiae 37-41 nt
  # rpf: [29, 34]                 # USER OVERRIDE — only if read-length QC
                                  # confirms a tighter window carries the
                                  # periodic signal. Omit to use the
                                  # footprint_class default.
  align: stop                     # anchor offsets at the stop codon
  offset_type: "5"                # offsets reported from the read 5' end
  offset_site: p                  # coordinate space of the SELECTED OFFSETS
                                  # table (controls offset_applied.csv only)
  offset_pick_reference: p_site
  offset_mode: per_sample         # per-sample offsets drive downstream
  analysis_sites: both            # which downstream P-site/A-site outputs
                                  # to generate (independent of offset_site)
  min_5_offset: 10
  max_5_offset: 22
  min_3_offset: 10
  max_3_offset: 22
  offset_mask_nt: 5
  plot_format: svg
  codon_density_window: true
```

> `offset_site` controls the coordinate system of the selected-offset table. `analysis_sites` controls which downstream P-site / A-site outputs are produced. They are independent — for example, you can pick offsets in P-site space and still emit both P-site and A-site coverage profiles.

After the run, you'll have:

```text
results/
  align/    bed/, kit_resolution.tsv, read_counts.tsv, run_settings.json
            (intermediate trimmed/contam_filtered/aligned files are deleted as
             soon as they are consumed; pass --keep-intermediates to retain
             them; deduped/ is only created for UMI samples)
  rpf/      offset_diagnostics/{csv,plots}/, translation_profile/<sample>/...,
            coverage_profile_plots/{read_coverage_*, {p_site,a_site}_density_*},
            codon_correlation/{p_site,a_site}/, igv_tracks/<sample>/, rpf_counts.tsv
  run_manifest.json
```

See [Output overview](#output-overview) for the full directory tree.

---

## Inputs you need to prepare

What every subcommand needs on disk before you run anything (the per-stage required vs optional matrix), the unified-sample-sheet schema, and the file-by-file reference for FASTQs / BEDs / FASTAs / annotations live in [`docs/inputs.md`](docs/inputs.md). Pair with [`docs/reference/sample_sheet_schema.md`](docs/reference/sample_sheet_schema.md) for the cell-value validation rules.

---

## How to run — YAML vs shell wrapper

### Recommended: YAML config

```bash
mitoribopy all --config pipeline_config.yaml --output results/ --threads 8
```

This is the canonical invocation. The YAML is self-documenting, version-controllable, and loads exactly the same way the CLI reads it programmatically.

For multi-sample batches, add `max_parallel_samples: N` under the `align:` section of your YAML to align samples concurrently — `--threads` is auto-divided across workers so total CPU use stays ≈ `--threads`. (`mitoribopy all` does not take `--max-parallel-samples` directly at the CLI; the flag is read from the YAML and forwarded to the align stage. The standalone `mitoribopy align` subcommand does accept it on the CLI — see [Example B](#b-just-the-align-stage-on-a-directory-of-fastqs) below.) The joint `rpf` stage stays serial (offset selection is a pooled-across-samples computation). See [Execution / concurrency](#execution--concurrency) under the `mitoribopy align` reference.

### Alternative: bash wrapper for batch / cluster jobs

When every flag should be visible in the job script (e.g. for a cluster scheduler that captures stdout/stderr per task), wrap the YAML invocation in a thin shell wrapper:

```bash
#!/usr/bin/env bash
set -uo pipefail

ENV_BIN=/path/to/conda/envs/mitoribopy/bin
export PATH="$ENV_BIN:$PATH"

ROOT=/path/to/project
cd "$ROOT"

OUT=results/full_run
mkdir -p "$OUT"

mitoribopy all \
  --config pipeline_config.yaml \
  --output "$OUT" \
  --threads 8
RC=$?

echo
echo "================ kit_resolution.tsv ================"
cat "$OUT/align/kit_resolution.tsv"
echo
echo "================ read_counts.tsv ================"
cat "$OUT/align/read_counts.tsv"

echo "$RC" > "$OUT.exitcode"
exit $RC
```

### Direct subcommand invocation

You can run any single stage directly without going through `mitoribopy all`. This is useful when you only have BED inputs (skip `align`), or when you want to iterate on `rpf` parameters without re-running alignment.

```bash
# Just align (auto-detection picks the adapter; pass --adapter <SEQ>
# only when detection cannot identify your library, or --pretrimmed
# for already-trimmed FASTQs).
mitoribopy align --fastq-dir fastqs/ \
  --contam-index idx/rrna --mt-index idx/mt --output results/align/

# Just rpf, against an existing BED dir
mitoribopy rpf -s h -f ref.fa --directory bed/ -rpf 29 34 --output results/rpf/

# Just rnaseq (default flow): raw FASTQ -> pyDESeq2 -> TE / dTE
mitoribopy rnaseq --rna-fastq rna_seq/ --ribo-fastq ribo_seq/ \
  --reference-fasta ref.fa --gene-id-convention bare \
  --condition-map samples.tsv \
  --base-sample control --compare-sample knockdown \
  --output results/rnaseq/

# `--base-sample` / `--compare-sample` are aliases for `--condition-a` /
# `--condition-b`; pick whichever spelling you prefer (the legacy form
# still works). The base condition is the reference (denominator) of
# the WT-vs-X contrast and seeds the labels on every comparison plot.

# Or rnaseq (alternative flow): existing rpf output + an external DE table
mitoribopy rnaseq --de-table de.tsv --gene-id-convention hgnc \
  --ribo-dir results/rpf --reference-gtf ref.fa --output results/rnaseq/
```

---

## Strain presets, custom organisms, and footprint classes

Built-in `h.sapiens` and `s.cerevisiae` strain presets ship complete annotation tables and codon tables. For any other organism — mouse / rat / *Drosophila* / *C. elegans* / etc. — pass `--strain custom` plus a per-transcript annotation CSV and a codon-table choice from the 27 NCBI Genetic Codes bundled with the package. The full reference (strain matrix, footprint-class window defaults, codon-table picker, annotation-CSV schema, the worked mouse example, and the bicistronic-pair handling) lives in [`docs/custom_organisms.md`](docs/custom_organisms.md).

---

## Subcommand reference

Every subcommand inherits these shared options (full per-flag detail in [`docs/reference/cli.md`](docs/reference/cli.md), regenerated from the live argparse parsers):

| Flag | Default | Description |
|---|---|---|
| `--config PATH` | — | Configuration file (.json, .yaml, .yml, or .toml). CLI flags override values from the file. |
| `--dry-run` | off | Print planned actions and exit 0 without executing. |
| `--threads N` | 1 | Preferred thread count; exports `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `MITORIBOPY_THREADS`. When combined with `--max-parallel-samples M` (align only), each parallel worker's external tools see `max(1, T // M)` threads so the total CPU budget stays ≈ T. |
| `--log-level {DEBUG,INFO,WARNING,ERROR}` | `INFO` | Python logging level for console output. |

The publication-grade highlights for each subcommand:

* **`mitoribopy align`** — FASTQ → BAM → BED6 (cutadapt + bowtie2 + umi_tools + pysam). Adapter auto-detection is the default; pin with `--adapter <SEQ>` or `--pretrimmed`. Dedup canonical token is `umi_coordinate` (legacy `umi-tools` / `umi_tools` accepted as aliases). Required: `--contam-index`, `--mt-index`, `--output`, FASTQ inputs. The detector reports the matched kit family in `kit_resolution.tsv` for provenance.
* **`mitoribopy rpf`** — Ribo-seq analysis from BED / BAM. Required: `-f FASTA`, `-d BED_DIR`, `-o OUTPUT_DIR`, `-s {h.sapiens,s.cerevisiae,custom}`. End-specific offsets (`--min_5_offset` / `--max_5_offset` / `--min_3_offset` / `--max_3_offset`) are preferred over the shared bounds. Periodicity QC ships under `qc/`; see [`docs/reference/periodicity.md`](docs/reference/periodicity.md). For `--strain custom`, see [`docs/custom_organisms.md`](docs/custom_organisms.md).
* **`mitoribopy rnaseq`** — Two mutually exclusive flows: `de_table` (publication-grade — external full-transcriptome DESeq2 / Xtail / Anota2Seq) and `from_fastq` (exploratory mt-mRNA-only pyDESeq2). `--strict` refuses `from_fastq` unless `allow_exploratory_from_fastq_in_strict: true` is set. Full publication-boundary reference: [`docs/rnaseq_te.md`](docs/rnaseq_te.md).
* **`mitoribopy all`** — End-to-end orchestrator. Reads one YAML/JSON/TOML config and dispatches to the per-stage subcommands; auto-wires `rpf.directory`, `rpf.read_counts_file`, `rnaseq.ribo_dir`. Use `--print-config-template --profile {minimal,publication,exhaustive}` to bootstrap a config; pair with `--strict` for publication-safe runs.

YAML config shape: every key under a section maps to that subcommand's CLI flag, hyphens converted to underscores. `align:` and `rnaseq:` use hyphen style (`--adapter-detection` → `adapter_detection`); `rpf:` uses underscore style (`--offset_type` → `offset_type`). Booleans emit the bare flag (`true`) or are omitted (`false`); `null` is dropped.

---

## Further reading

The README intentionally stays short. The detail lives under [`docs/`](docs/) — every page is a single, focused topic with cross-links into the reference set:

* **Inputs & sample sheet** — [`docs/inputs.md`](docs/inputs.md) (per-stage required vs optional matrix; sample-sheet TSV; file-by-file reference) plus [`docs/reference/sample_sheet_schema.md`](docs/reference/sample_sheet_schema.md) (cell-value validation rules, the canonical token map).
* **Custom organisms** — [`docs/custom_organisms.md`](docs/custom_organisms.md) (strain matrix, footprint-class window defaults, codon-table picker, annotation-CSV schema, mouse worked example, bicistronic-pair handling).
* **TE / ΔTE numerics** — [`docs/te_numerics.md`](docs/te_numerics.md) (TE / ΔTE equations, per-row `note` taxonomy, what is and isn't in the gene-level summary tables).
* **RNA-seq publication boundaries** — [`docs/rnaseq_te.md`](docs/rnaseq_te.md) (when to use `de_table` vs `from_fastq`, the strict-mode gates, output-file matrix).
* **Periodicity QC** — [`docs/reference/periodicity.md`](docs/reference/periodicity.md) (metagene Fourier method, statistical hardening, when-not-to-overinterpret caveats).
* **Output schema** — [`docs/reference/output_schema.md`](docs/reference/output_schema.md) (column-by-column reference for every TSV / CSV / JSON the package writes; units; coordinate spaces).
* **Warning / error codes** — [`docs/reference/warning_codes.md`](docs/reference/warning_codes.md) (every `warnings.tsv` code with severity, stage, summary, remediation, and publication impact).
* **CLI reference** — [`docs/reference/cli.md`](docs/reference/cli.md) (every flag in every subcommand, regenerated from the live argparse parsers; CI fails when it drifts).
* **Validation & regression** — [`docs/validation/`](docs/validation/) (TACO1-KO biological gate, synthetic-mini integration test, public-dataset reanalysis, UMI dedup, RNA-seq TE).
* **Tutorials** — [`docs/tutorials/`](docs/tutorials/) (end-to-end FASTQ flow, RNA-seq integration, HPC cluster run).
* **Benchmarking** — [`docs/benchmarking.md`](docs/benchmarking.md) (table schema, command, reference cases).
* **Release & developer** — [`docs/developer/release_checklist.md`](docs/developer/release_checklist.md), [`docs/developer/architecture_history.md`](docs/developer/architecture_history.md), [`docs/developer/roadmap.md`](docs/developer/roadmap.md).
* **Smoke fixture** — [`examples/smoke/README.md`](examples/smoke/README.md) (opt-in via `pytest -m smoke`).

---

## What the numbers mean — RNA, RPF, TE, ΔTE

Translation efficiency is one of the most-asked-for and most-misread quantities in Ribo-seq, so MitoRiboPy pins each quantity to a specific output file and a specific equation. The full reference — TE / ΔTE equations, the per-row `note` taxonomy, what is and isn't in the gene-level summary tables, and pseudo-replicate caveats — lives in [`docs/te_numerics.md`](docs/te_numerics.md). The metagene Fourier 3-nt periodicity QC bundle has its own dedicated reference at [`docs/reference/periodicity.md`](docs/reference/periodicity.md).

---

## Output overview

A column-by-column reference for every TSV / CSV / JSON the pipeline writes — including units, coordinate spaces, and which downstream consumer reads each file — lives in [`docs/reference/output_schema.md`](docs/reference/output_schema.md). The tree below is the at-a-glance shape; pair it with that doc when you need to understand a specific column.

For a `mitoribopy all` run with the defaults (`--offset_mode per_sample`, `--analysis_sites both`):

```text
<output>/
  run_manifest.json                    # composed provenance: settings, tool versions,
                                       # reference_checksum, stages_run
  SUMMARY.md, summary_qc.tsv,          # auto-rendered after every run
  warnings.tsv, figure_qc.tsv,
  outputs_index.tsv

  align/
    read_counts.tsv                    # per-stage funnel (input → trim → contam → mt → MAPQ → dedup)
    kit_resolution.tsv                 # per-sample kit + UMI + dedup decisions
    run_settings.json
    aligned/<sample>.mapq.bam          # MAPQ-filtered BAMs (kept)
    deduped/<sample>.dedup.bam         # only when at least one sample resolves to dedup_strategy=umi_coordinate (i.e. invokes the umi_tools binary)
    bed/<sample>.bed                   # strand-aware BED6 inputs to rpf
    trimmed/, contam_filtered/         # only --keep-intermediates

  rpf/
    rpf_counts.tsv                     # per-(sample, gene) counts; feeds rnaseq
    run_settings.json                  # includes reference_checksum
    offset_diagnostics/
      csv/                             # offset_<align>.csv (combined),
                                       # per_sample_offset/<sample>/offset_applied.csv
      plots/                           # offset_drift_<align>.svg, heatmaps, line plots
    translation_profile/<sample>/
      footprint_density/, translating_frame/, codon_usage/   # site is in the filename prefix
    coverage_profile_plots/
      read_coverage_{rpm,raw}[_codon]/         # full-footprint depth (site-independent)
      {p_site,a_site}_density_{rpm,raw}/       # single-nt density at the chosen site
      {p_site,a_site}_density_*_frame/         # frame-coloured CDS overlay (frame-0 dominance = QC)
      {p_site,a_site}_density_*_frame_split/   # 3 stacked sub-rows per sample (frame 0, +1, +2)
    qc/                                # metagene Fourier QC bundle (always emitted)
      fourier_spectrum_combined.tsv          # per-(sample, length, gene_set, region) amplitude curve
      fourier_period3_score_combined.tsv     # spectral_ratio_3nt + snr_call tier (excellent/healthy/modest/broken)
      periodicity.metadata.json
      fourier_spectrum/<sample>/<sample>_<length>nt_{combined,ATP86,ND4L4}.{png,svg}
      metagene_{start,stop}.tsv              # P-site metagene tables
      metagene_{start,stop}_p_site.svg       # 3-nt periodicity plots
      strand_sanity.tsv                      # per-sample minus-strand fraction (should be 0)
    codon_correlation/{p_site,a_site}/       # if --cor_plot
    igv_tracks/<sample>/<sample>_{p_site,a_site}.bedgraph     # if --igv_export
    structure_density/                       # if --structure_density

  rnaseq/                              # if rnaseq config supplied
    te.tsv, delta_te.tsv               # headline tables
    plots/                             # 6 always-emitted publication plots; up to 3 more
                                       # (mrna_vs_rpf, delta_te_volcano, ma, de_volcano_mrna,
                                       # te_bar_by_condition, te_heatmap, [te_compare_scatter,
                                       # te_log2fc_bar, sample_pca]) at 300 dpi PNG + editable SVG
    de_table.tsv, rpf_de_table.tsv,    # default flow only (from raw FASTQ)
    rna_counts.tsv, rpf_counts.tsv,
    rpf_counts_matrix.tsv,
    condition_map.augmented.tsv
```

### What to inspect first

After every run, walk these files in order. The first three (`SUMMARY.md`, `warnings.tsv`, `summary_qc.tsv`) catch ~90 % of gotchas without digging into per-stage TSVs.

1. **`SUMMARY.md`** — one-page human-readable view of what ran and what was produced.
2. **`warnings.tsv`** — every structured warning (`stage / sample_id / severity / code / message / suggested_action`). Header-only = nothing flagged.
3. **`summary_qc.tsv`** — per-sample QC roll-up across all stages.
4. **`align/kit_resolution.tsv`** — per-sample kit / UMI / dedup decisions; the `source` column tells you whether it was `detected`, `user_fallback`, `inferred_pretrimmed`, or set explicitly.
5. **`align/read_counts.tsv`** — per-stage drop-off; invariants `rrna_aligned + post_rrna_filter == post_trim` and `mt_aligned + unaligned_to_mt == post_rrna_filter` must hold.
6. **`rpf/offset_diagnostics/plots/offset_drift_<align>.svg`** — per-sample offset drift, visible by eye.
7. **`rpf/offset_diagnostics/csv/per_sample_offset/<sample>/offset_applied.csv`** — exact offset row applied downstream.
8. **`rpf/qc/fourier_period3_score_combined.tsv`** — metagene Fourier periodicity verdict per `(sample, read_length, gene_set, region)`. Look at `gene_set=combined` rows: `snr_call ∈ {excellent, healthy}` is publication-grade; `modest` is borderline; `broken` means offset assignment is suspect for that read length. **Read this before trusting any downstream codon-level table.**
9. **`rpf/qc/fourier_spectrum/<sample>/*.png`** — three figures per (sample, length): `combined`, `ATP86`, `ND4L4`. Each panel reports its `spectral_ratio_3nt` and `snr_call` in-figure.
10. **`rpf/coverage_profile_plots/p_site_density_rpm_frame/<transcript>_*.svg`** — frame-coloured CDS density. Frame-0 dominance (~70-90 %) is the canonical mt-Ribo-seq QC signature; flat or jittery frames suggest contamination or poor offset selection. The `_frame_split/` companion stacks each frame in its own sub-row when the overlay's tallest frame would hide low-frame signal.
11. **`rpf/rpf_counts.tsv`** + sidecar — per-(sample, gene) RPF count matrix.
12. **`rnaseq/te.tsv` + `rnaseq/delta_te.tsv`** — translation efficiency and ΔTE (when the rnaseq stage ran).
13. **`figure_qc.tsv`** (after `mitoribopy validate-figures runs/full/`) — mechanical pass/warn/fail per plot.

### Common reasons NOT to trust the output

The pipeline can finish exit-0 under several conditions that should still make a reviewer pause. Check `warnings.tsv` first; the `suggested_action` column tells you what to do.

- **Adapter detection low confidence** (`source=user_fallback` in `kit_resolution.tsv`).
- **High contaminant fraction** (>50 % of post-trim reads filtered as rRNA / contaminant).
- **Low mt-mRNA alignment fraction** (`mt_aligned / post_rrna_filter < 0.05`).
- **No clear offset peak** in `offset_*.svg` — selected offset was a fallback, not a real peak.
- **Periodicity `broken` or `no_signal`** in `fourier_period3_score_combined.tsv` — offsets did not produce codon phasing for that read length; cross-check `metagene_start_p_site.svg` for the same length.
- **Pseudo-replicate mode** — `rnaseq/EXPLORATORY.md` exists. DE statistics are exploratory only.
- **Gene-ID match rate below threshold** — `warnings.tsv` row `GENE_ID_MATCH_RATE_LOW`. DE table and rpf reference do not match.
- **`figure_qc.tsv` has any `fail` row** — overlapping labels, clipped text, or point-count mismatch.

### Publication route

1. `mitoribopy align` + `mitoribopy rpf` for RPF counts + Fourier QC.
2. Run RNA-seq DE externally on the FULL transcriptome with your tool of choice (DESeq2 / Xtail / Anota2Seq) for publication-grade statistics.
3. `mitoribopy rnaseq --de-table de.tsv --ribo-dir runs/full/rpf/` to compute TE / ΔTE on top.
4. `mitoribopy validate-figures runs/full/ --strict` to mechanically QC every plot.
5. Bundle `run_manifest.json`, `summary_qc.tsv`, `figure_qc.tsv`, and the SVG sidecars into the paper's supplement.

---

## Examples

Three representative invocations. For more (footprint classes, custom organisms, mixed-UMI batches, resume semantics), see the templates under [examples/templates/](examples/templates/) and the tutorials under [docs/tutorials/](docs/tutorials/).

### A. End-to-end on raw FASTQ (recommended)

```bash
mitoribopy all --config pipeline_config.yaml --output results/ --threads 8
```

`--print-config-template` prints a commented YAML covering every stage. `--resume` is hash-validated — config / sample sheet / FASTA edits force the affected stage(s) to re-run. For mixed-kit / mixed-UMI batches, declare a per-sample list under `align.samples:` (overrides materialised to `<output>/align/sample_overrides.tsv`). For batches with many samples, set `align.max_parallel_samples` to align N samples concurrently while the joint `rpf` stage stays serial.

### B. Just the rpf stage on existing BEDs

```bash
mitoribopy rpf \
  -s h.sapiens \
  -f references/human-mt-mRNA.fasta \
  --directory results/align/bed/ \
  -rpf 29 34 \
  --offset_mode per_sample --analysis_sites both \
  --rpm_norm_mode mt_mrna \
  --read_counts_file results/align/read_counts.tsv \
  --output results/rpf/ \
  --plot_format svg
```

For disome studies, replace the explicit `-rpf 29 34` with `--footprint_class disome` (auto-widens to 50-70 nt for human, 60-90 nt for yeast). For custom organisms, use `-s custom` plus `--annotation_file` and `--codon_table_name`.

### C. RNA-seq integration for translation efficiency

```bash
# Default flow: raw FASTQ -> pyDESeq2 -> TE / ΔTE in one shot
# (one-off setup: pip install 'mitoribopy[fastq]')
mitoribopy rnaseq \
  --rna-fastq input_data/rna_seq/  --ribo-fastq input_data/ribo_seq/ \
  --reference-fasta references/human-mt-mRNA.fasta \
  --gene-id-convention bare \
  --condition-map samples.tsv --base-sample control --compare-sample knockdown \
  --output results/rnaseq/  --align-threads 8
```

For publication-grade DE statistics, run DE on the full transcriptome externally (DESeq2 / Xtail / Anota2Seq) and feed the result via `--de-table`. The two flows are mutually exclusive (passing both exits 2). Templates: [examples/templates/run_rnaseq.example.sh](examples/templates/run_rnaseq.example.sh) and [examples/templates/rnaseq_config.example.yaml](examples/templates/rnaseq_config.example.yaml).

---

## Tools

Standalone helper scripts under `mitoribopy.tools.*` — useful for shrinking inputs to a quick smoke-test size before a full pipeline run.

### Subsample a BED file

```bash
python -m mitoribopy.tools.subsample \
    --input  results/align/bed/sample.bed \
    --output /tmp/sample_subsampled.bed \
    --n 50000 \
    --seed 42
```

Reservoir-samples `--n` BED rows (Algorithm R, deterministic with `--seed`). Use the subsampled BED as `--directory` input to `mitoribopy rpf` for a fast end-to-end test.

### Subsample a FASTQ file

```bash
python -m mitoribopy.tools.subsample_fastq \
    --input  raw/sample.fastq.gz \
    --output /tmp/sample_subsampled.fastq.gz \
    --n 200000 \
    --seed 42
```

Reservoir-samples `--n` FASTQ records (4 lines each); gzip is auto-detected on either end via the `.gz` suffix. Plain `.fastq` works too. The first record header must start with `@` or the tool fails fast.

---

## Logs and provenance

- **`<output>/<stage>/mitoribopy.log`** — every stage writes a persistent log file alongside the same lines printed to the terminal.
- **Per-stage `run_settings.json`** — every stage writes its own settings JSON (resolved kit, dedup strategy, MAPQ threshold, reference checksum, tool versions, …).
- **`<output>/run_manifest.json`** — `mitoribopy all` composes per-stage settings into a top-level manifest. Schema version 1.0.0 carries: `schema_version`, `mitoribopy_version`, `git_commit` (best-effort, `null` outside a repo), `command` (the original argv), `config_source` + `config_source_sha256` (hash of the YAML the user wrote), `config_canonical` (the merged + auto-wired config that actually drove the run), `sample_sheet` + `sample_sheet_sha256` when applicable, a `stages: { align: {status, runtime_seconds}, rpf: {...}, rnaseq: {...} }` map (status is `completed` / `skipped` / `not_configured`; skipped stages carry a `reason`), the per-stage `run_settings.json` payloads, a flat `tools: {python, mitoribopy, cutadapt, ...}` map lifted from those payloads, and a `warnings` placeholder. The `reference_checksum` from the rpf stage is promoted to the top level so a downstream `rnaseq` run can verify it without drilling into the rpf section. Pin to a manifest layout by reading `schema_version` first.
- **`<output>/align/kit_resolution.tsv`** — per-sample kit + dedup decisions.

---

## Development

Run the test suite with:

```bash
PYTHONPATH=src pytest
```

Helpful tools:

```bash
PYTHONPATH=src pytest -k offset           # run a subset by keyword
PYTHONPATH=src pytest -x --tb=short       # stop at first failure, terse traceback
PYTHONPATH=src pytest tests/test_align_sample_resolve.py -v
```

Documentation lives under [docs/](docs/):

- [docs/README.md](docs/README.md) — index
- [docs/reference/cli.md](docs/reference/cli.md) — concise CLI reference
- [docs/tutorials/](docs/tutorials/) — step-by-step worked examples
- [docs/release-notes/](docs/release-notes/) — version-by-version notes
- [docs/validation/](docs/validation/) — biological validation plans
- [docs/environment/](docs/environment/) — bioconda env file + Dockerfile
- [docs/diagrams/](docs/diagrams/) — Mermaid pipeline diagrams

---

## Citation

For exact reproducibility, cite the GitHub release tag and the Zenodo DOI together — not just the package name. PyPI's `mitoribopy` slot can move forward; a tagged release does not.

Suggested manuscript block:

```
Software: MitoRiboPy v0.7.0
DOI: <Zenodo DOI for the v0.7.0 release>
Repository: https://github.com/Ahram-Ahn/MitoRiboPy/releases/tag/v0.7.0
Python: <Python version recorded in run_manifest.json>
External tools (cutadapt, bowtie2, umi_tools, pysam, pydeseq2):
  versions recorded in <run_root>/run_manifest.json under tool_versions.
```

A machine-readable [`CITATION.cff`](CITATION.cff) sits at the repository root for tools that consume the [Citation File Format](https://citation-file-format.github.io/) (GitHub's citation widget, Zenodo, Zotero, etc.).

---

## Known limitations

Reviewers tend to trust tools more when their boundaries are stated up front. MitoRiboPy is intentionally narrow:

- **Mitochondrial Ribo-seq, not nuclear.** The default references and the periodicity / coverage outputs assume mt-transcriptome-style alignment (one chromosome ≈ one transcript). Running it against a nuclear genome with cytoplasmic Ribo-seq is out of scope; `mt_index` and `mt_mrna_substring_patterns` are not generic Ribo-seq plumbing.
- **`rnaseq_mode: from_fastq` is exploratory.** The in-tree pyDESeq2 path runs on the **mt-mRNA subset only** (typically 13 transcripts) — not a full-transcriptome DE result. `mitoribopy all --strict` refuses this mode by default; the publication-safe route is `rnaseq_mode: de_table` with an external full-transcriptome DESeq2 / Xtail / Anota2Seq table.
- **Pseudo-replicates are not biological replicates.** `allow_pseudo_replicates_for_demo_not_publication: true` exists for demo / smoke runs only. Strict mode rejects it; do not report p-values from a pseudo-replicate run as biological evidence.
- **Offset selection reliability scales with depth and periodicity.** Per-sample offsets need both healthy 3-nt periodicity and enough reads per (sample, length). Below ~1 000 reads per length the per-sample pick is noisy — use the combined-sample offsets or pool replicates first.
- **Bicistronic mt-mRNAs need careful interpretation.** ATP8/ATP6 and ND4L/ND4 share nucleotides in different reading frames. The Fourier bundle ships dedicated `*_ATP86.png` / `*_ND4L4.png` panels for this; reviewers should inspect those panels rather than collapsing them into the combined metagene.
- **NUMT suppression depends on reference design and MAPQ.** Nuclear-mitochondrial transfers can recruit reads from real mt-RNA away from the mt reference. The default `mapq: 10` filter helps, but custom references should be built to suppress NUMTs explicitly (mask known NUMT regions, or use an mt-only transcriptome + contam index).
- **Metagene Fourier confidence depends on enough genes.** The score table carries `spectral_ratio_3nt(_local)_ci_{low,high}` (200-iteration bootstrap CI over genes, 90 % percentile by default) and `permutation_p` / `permutation_p_local` (200-iteration circular-shift null, Laplace-smoothed). The CI is skipped (NaN columns + `ci_method == "skipped_too_few_genes"`) when fewer than 3 qualifying per-gene tracks are available — a CI from < 3 genes would be misleadingly tight. The `snr_call` four-tier verdict is the human-readable headline; cite the CI bounds + permutation p in any publication-facing context.

For pipeline-level provenance and what survives a refactor, see [`docs/validation/`](docs/validation/) — in particular the TACO1-KO regression dataset, which is the biological signal that any periodicity refactor has to preserve.

---

## License

MIT. See [LICENSE](LICENSE).
