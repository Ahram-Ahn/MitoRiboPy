<p align="center">
  <img src="docs/banner.png" alt="MitoRiboPy - A Python package for mt-Ribo-seq" width="100%">
</p>

# MitoRiboPy

MitoRiboPy is a command-line pipeline for mitochondrial ribosome profiling
(mt-Ribo-seq). It takes raw FASTQ files through trimming, contaminant removal,
mitochondrial transcriptome alignment, RPF offset selection, P-site/A-site
profiles, codon usage, coverage plots, and optional RNA-seq integration for
translation efficiency (TE / delta-TE).

It is designed for focused mitochondrial analyses: human and yeast mt-mRNAs are
supported out of the box, and other mitochondrial transcriptomes can be supplied
with a custom annotation and codon table.

## When to use it

Use MitoRiboPy when you want to:

- process mt-Ribo-seq FASTQs into BAM/BED files and per-sample QC tables,
- estimate RPF offsets and generate P-site/A-site translation profiles,
- inspect 3-nt periodicity with metagene Fourier QC,
- compare mitochondrial RPF counts with paired RNA-seq,
- run mixed batches where adapters, UMIs, and deduplication choices may differ
  by sample,
- keep a reproducible run record with manifests, resolved configs, warnings,
  and figure QC.

MitoRiboPy is intentionally narrow. It is not a general nuclear Ribo-seq
pipeline, and its built-in RNA-seq-from-FASTQ mode runs pyDESeq2 on the
mt-mRNA subset only. For publication-grade differential expression, run DE on
the full transcriptome externally and pass that table back to MitoRiboPy.

## Pipeline At A Glance

![Pipeline overview](docs/diagrams/01_pipeline_overview.png)

| Command | Input | Output | Use it when |
|---|---|---|---|
| `mitoribopy align` | FASTQ | MAPQ-filtered BAM, BED6, read-count QC | You only need preprocessing and alignment. |
| `mitoribopy rpf` | BED or BAM | offsets, translation profiles, codon usage, coverage plots, periodicity QC | You already have aligned Ribo-seq reads. |
| `mitoribopy rnaseq` | RNA-seq + Ribo-seq inputs | `te.tsv`, `delta_te.tsv`, RNA/RPF comparison plots | You want TE / delta-TE. |
| `mitoribopy all` | one config file | full align + rpf + optional rnaseq run | You want the standard end-to-end workflow. |

Utility commands are also included:

- `mitoribopy validate-config` checks a YAML/JSON/TOML config before a long run.
- `mitoribopy validate-reference` checks a custom FASTA + annotation pair.
- `mitoribopy validate-figures` checks generated plots after a run.
- `mitoribopy periodicity` reruns metagene Fourier QC from a saved site table.
- `mitoribopy summarize` regenerates `SUMMARY.md` from a finished run.
- `mitoribopy benchmark` measures runtime, memory, and disk use.

For every flag, see the generated [CLI reference](docs/reference/cli.md) or run:

```bash
mitoribopy <command> --help
```

## Install

MitoRiboPy requires Python 3.10 or newer. The current documented interface is
v0.7.1.

### Option 1: install from PyPI

```bash
python -m pip install "mitoribopy>=0.7.1"
mitoribopy --version
```

Install the optional RNA-seq FASTQ workflow only if you need in-tree pyDESeq2:

```bash
python -m pip install "mitoribopy[fastq]>=0.7.1"
```

### Option 2: install from source

```bash
git clone https://github.com/Ahram-Ahn/MitoRiboPy.git
cd MitoRiboPy
python -m pip install -e .
mitoribopy --version
```

For development:

```bash
python -m pip install -e ".[dev]"
```

### External tools

Real runs also need standard bioinformatics tools on `PATH`:

| Tool | Needed for |
|---|---|
| `cutadapt` | adapter trimming, length/quality filtering |
| `bowtie2` and `bowtie2-build` | contaminant subtraction and mt-transcriptome alignment |
| `umi_tools` | UMI-aware deduplication when a sample has UMIs |
| `samtools` | optional output inspection |

The conda environment file installs the Python stack plus the external tools:

```bash
conda env create -f docs/environment/environment.yml
conda activate mitoribopy
```

Native Windows is not supported because the external bioinformatics tools are
not consistently available there. Use Linux, macOS, WSL2, or a container.

## Quick Check

After installation, confirm that the CLI is available:

```bash
mitoribopy --version
mitoribopy --help
```

The repository also includes a tiny synthetic smoke fixture:

```bash
cd examples/smoke
python generate_smoke_fastqs.py
mitoribopy all --config pipeline_config.smoke.yaml --output results/
```

The smoke fixture is a wiring check for install, FASTQ generation, alignment,
RPF execution, and output creation. It is not a biological validation dataset.
See [examples/smoke/README.md](examples/smoke/README.md) for current fixture
limitations and the optional `pytest -m smoke` workflow.

## First Real Run

For most users, start with `mitoribopy all` and one YAML config:

```bash
mitoribopy all --print-config-template --profile minimal > pipeline_config.yaml
$EDITOR pipeline_config.yaml
mitoribopy validate-config pipeline_config.yaml
mitoribopy all --config pipeline_config.yaml --output results/ --threads 8
```

Use `--dry-run` before launching a long job:

```bash
mitoribopy all --config pipeline_config.yaml --output results/ --threads 8 --dry-run
```

Use `--resume` after an interrupted run:

```bash
mitoribopy all --config pipeline_config.yaml --output results/ --threads 8 --resume
```

## Inputs You Need

A typical end-to-end run needs:

- Ribo-seq FASTQ files,
- a contaminant bowtie2 index, usually rRNA/tRNA or other sequences to remove,
- a mitochondrial transcriptome bowtie2 index,
- the FASTA used to build the mitochondrial transcriptome index,
- a sample sheet that names samples, assays, conditions, and FASTQ paths,
- optional RNA-seq FASTQs or an external full-transcriptome DE table.

The recommended layout is one sample sheet plus one pipeline config. The same
sample sheet can drive both Ribo-seq and RNA-seq stages.

```yaml
samples:
  table: samples.tsv

align:
  adapter_detection: auto
  library_strandedness: forward
  contam_index: input_data/indexes/rrna_contam
  mt_index: input_data/indexes/mt_tx
  mapq: 10
  min_length: 15
  max_length: 45
  dedup_strategy: auto

rpf:
  strain: h.sapiens
  fasta: input_data/human-mt-mRNA.fasta
  footprint_class: monosome
  align: stop
  offset_type: "5"
  offset_site: p
  offset_pick_reference: p_site
  offset_mode: per_sample
  analysis_sites: both
  min_5_offset: 10
  max_5_offset: 22
  min_3_offset: 10
  max_3_offset: 22
  plot_format: svg
```

Reference docs:

- [Inputs](docs/inputs.md) explains required files for each stage.
- [Sample sheet schema](docs/reference/sample_sheet_schema.md) lists every
  sample-sheet column.
- [Config schema](docs/reference/config_schema.md) documents the YAML shape.
- [Custom organisms](docs/custom_organisms.md) covers non-human/non-yeast
  mitochondrial transcriptomes.

## Common Workflows

### End-to-end from FASTQ

```bash
mitoribopy all --config pipeline_config.yaml --output results/full_run --threads 8
```

Set `align.max_parallel_samples` in the YAML when aligning many samples. The
orchestrator divides the thread budget across sample workers, while the joint
RPF stage remains serial because offset selection is pooled across samples.

### Alignment only

```bash
mitoribopy align \
  --fastq-dir fastqs/ \
  --contam-index indexes/rrna_contam \
  --mt-index indexes/mt_tx \
  --output results/align/ \
  --threads 8
```

Adapter detection is automatic by default. Pin `--adapter` only when detection
cannot identify your library, or pass `--pretrimmed` for already-trimmed FASTQs.

### RPF analysis from existing BEDs

```bash
mitoribopy rpf \
  --strain h.sapiens \
  --fasta references/human-mt-mRNA.fasta \
  --directory results/align/bed/ \
  --footprint-class monosome \
  --offset-mode per_sample \
  --analysis-sites both \
  --read-counts-file results/align/read_counts.tsv \
  --output results/rpf/ \
  --plot-format svg
```

For disome libraries, use `--footprint-class disome`. For custom organisms,
use `--strain custom` with `--annotation-file` and `--codon-table-name`.

### Translation efficiency with external DE results

This is the recommended publication route. Run RNA-seq differential expression
on the full transcriptome with DESeq2, Xtail, Anota2Seq, or another validated
workflow, then combine it with MitoRiboPy RPF counts:

```bash
mitoribopy rnaseq \
  --de-table external_full_transcriptome_de.tsv \
  --ribo-dir results/full_run/rpf/ \
  --reference-gtf references/gencode_or_refseq.gtf \
  --gene-id-convention hgnc \
  --base-sample WT \
  --compare-sample KO \
  --output results/full_run/rnaseq/
```

### Translation efficiency from raw FASTQs

This mode is useful for exploration and demos. It runs pyDESeq2 only on the
mitochondrial mRNA subset, so do not treat its p-values as full-transcriptome
publication statistics.

```bash
python -m pip install "mitoribopy[fastq]>=0.7.1"

mitoribopy rnaseq \
  --rna-fastq input_data/rna_seq/ \
  --ribo-fastq input_data/ribo_seq/ \
  --reference-fasta references/human-mt-mRNA.fasta \
  --gene-id-convention bare \
  --condition-map samples.tsv \
  --base-sample control \
  --compare-sample knockdown \
  --output results/rnaseq/ \
  --align-threads 8
```

## Publication Mode

For manuscript, preprint, or shared-dataset runs, use strict mode:

```bash
mitoribopy validate-config pipeline_config.yaml --strict
mitoribopy all \
  --config pipeline_config.yaml \
  --output results/full_run \
  --threads 8 \
  --strict \
  --progress jsonl \
  --progress-file results/full_run/progress.jsonl
```

Strict mode:

- validates the config before running,
- rejects ambiguous or non-publication-safe align choices,
- writes the resolved `canonical_config.yaml`,
- runs figure validation after the pipeline finishes,
- promotes warning-level figure QC issues to failures,
- writes `SUMMARY.md`, `summary_qc.tsv`, `warnings.tsv`, `outputs_index.tsv`,
  `run_manifest.json`, and optional progress JSONL.

For publication TE / delta-TE, use `rnaseq_mode: de_table` with an external
full-transcriptome DE table. The `from_fastq` mode is exploratory and is
rejected by strict mode unless explicitly overridden for a non-publication run.

More detail: [RNA-seq TE boundaries](docs/rnaseq_te.md),
[TE numerics](docs/te_numerics.md), and
[periodicity QC](docs/reference/periodicity.md).

## What To Inspect First

After a run, start with these files:

1. `SUMMARY.md` - human-readable summary of what ran and what was produced.
2. `warnings.tsv` - structured warnings with suggested actions.
3. `summary_qc.tsv` - per-sample QC roll-up across stages.
4. `align/kit_resolution.tsv` - adapter, UMI, and dedup decisions per sample.
5. `align/read_counts.tsv` - read-count funnel from input through alignment.
6. `rpf/offset_diagnostics/plots/offset_drift_<align>.svg` - offset drift by sample.
7. `rpf/qc/fourier_period3_score_combined.tsv` - 3-nt periodicity verdicts.
8. `rpf/coverage_profile_plots/` - frame-colored density plots.
9. `rpf/rpf_counts.tsv` - per-sample, per-gene RPF counts.
10. `rnaseq/te.tsv` and `rnaseq/delta_te.tsv` - TE outputs when RNA-seq ran.
11. `figure_qc.tsv` - plot validation results when figure QC ran.

Be cautious with downstream interpretation if you see:

- low-confidence adapter detection,
- high contaminant fraction,
- low mt-mRNA alignment fraction,
- fallback offset selection rather than a clear offset peak,
- `broken` or `no_signal` periodicity calls,
- pseudo-replicate mode,
- low gene-ID match rate between DE and RPF tables,
- any `fail` row in `figure_qc.tsv`.

For the full output tree and column definitions, see
[output schema](docs/reference/output_schema.md).

## Documentation Map

Start with the workflow guides when setting up a run:

- [End-to-end FASTQ tutorial](docs/tutorials/01_end_to_end_fastq.md) - a complete `mitoribopy all` walk-through from FASTQs to RPF outputs.
- [RNA-seq integration tutorial](docs/tutorials/02_rnaseq_integration.md) - TE / delta-TE workflows, including the external-DE publication route.
- [HPC / cluster tutorial](docs/tutorials/05_hpc_cluster_run.md) - SLURM/LSF execution, thread accounting, scratch storage, and resume.
- [Custom organism workflow](docs/custom_organisms.md) - non-human/non-yeast references, annotation CSVs, and codon tables.

Use these references when preparing inputs or debugging outputs:

- [Input reference](docs/inputs.md) - required files for each stage.
- [Sample sheet schema](docs/reference/sample_sheet_schema.md) - canonical `samples.tsv` columns and validation rules.
- [Config schema](docs/reference/config_schema.md) - orchestrator YAML sections and auto-wiring behavior.
- [CLI reference](docs/reference/cli.md) - generated flag reference for every subcommand.
- [Output schema](docs/reference/output_schema.md) - file inventory and column definitions.
- [Warning and error codes](docs/reference/warning_codes.md) - stable `warnings.tsv` codes and remediations.

Use these for interpretation, operations, and release context:

- [RNA-seq TE boundaries](docs/rnaseq_te.md) - when `from_fastq` is exploratory and when `de_table` is publication-grade.
- [TE numerics](docs/te_numerics.md) - equations behind `te.tsv` and `delta_te.tsv`.
- [Periodicity QC](docs/reference/periodicity.md) - metagene Fourier method and `snr_call` interpretation.
- [Benchmarking](docs/benchmarking.md) - `mitoribopy benchmark` outputs for runtime, memory, and disk sizing.
- [Validation notes](docs/validation/) - regression and biological validation materials.
- [Release notes](docs/release-notes/) - version-by-version upgrade notes.
- [Documentation index](docs/README.md) - full docs tree.

## Development

Run the default test suite:

```bash
PYTHONPATH=src pytest
```

Useful subsets:

```bash
PYTHONPATH=src pytest -k offset
PYTHONPATH=src pytest -x --tb=short
PYTHONPATH=src pytest tests/test_align_sample_resolve.py -v
```

The generated CLI reference should not be edited by hand:

```bash
PYTHONPATH=src python docs/generate_cli_reference.py --check
```

## Citation

Cite the released software version, not only the package name.

Recommended citation for MitoRiboPy v0.7.1:

```text
Ahn A. MitoRiboPy: Mitochondrial ribosome profiling analysis pipeline.
Version 0.7.1. 2026. https://github.com/Ahram-Ahn/MitoRiboPy
```

BibTeX:

```bibtex
@software{ahn_mitoribopy_2026,
  author = {Ahn, Ahram},
  title = {MitoRiboPy: Mitochondrial ribosome profiling analysis pipeline},
  version = {0.7.1},
  date = {2026-05-02},
  url = {https://github.com/Ahram-Ahn/MitoRiboPy}
}
```

A machine-readable [CITATION.cff](CITATION.cff) is included for tools that
consume the Citation File Format. Runtime provenance, including Python and
external-tool versions, is recorded in each run's `run_manifest.json`.

## Known Limitations

- MitoRiboPy targets mitochondrial Ribo-seq, not nuclear genome-wide Ribo-seq.
- The built-in `rnaseq_mode: from_fastq` path is exploratory because it runs
  DE only on the mt-mRNA subset.
- Pseudo-replicates are for demos and smoke runs only, not biological evidence.
- Offset selection needs enough depth and healthy 3-nt periodicity.
- Bicistronic mt-mRNAs such as ATP8/ATP6 and ND4L/ND4 need careful
  interpretation.
- NUMT suppression depends on reference design and MAPQ thresholds.
- Fourier confidence intervals are skipped when too few qualifying genes are
  available.

## License

MIT. See [LICENSE](LICENSE).
