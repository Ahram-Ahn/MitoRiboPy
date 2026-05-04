# Documentation

The top-level [`README.md`](../README.md) is the standalone overview for
installation, first runs, common workflows, release citation, and known
limitations. The files under `docs/` provide task-specific tutorials,
generated references, output schemas, environment recipes, validation notes,
and release notes.

## Workflow guides

- [`tutorials/01_end_to_end_fastq.md`](tutorials/01_end_to_end_fastq.md) - step-by-step `mitoribopy all` run from raw FASTQ to translation-profile and codon-usage outputs.
- [`tutorials/02_rnaseq_integration.md`](tutorials/02_rnaseq_integration.md) - RNA-seq integration, external DE-table mode, reference-consistency gate, and TE / delta-TE outputs.
- [`tutorials/05_hpc_cluster_run.md`](tutorials/05_hpc_cluster_run.md) - SLURM / LSF execution, thread accounting, scratch storage, disk budgeting, and hash-validated resume.
- [`custom_organisms.md`](custom_organisms.md) - non-human/non-yeast mitochondrial transcriptomes, annotation CSVs, codon tables, and footprint classes.

## Input and command references

- [`inputs.md`](inputs.md) - required files for each subcommand and the recommended sample-sheet workflow.
- [`reference/sample_sheet_schema.md`](reference/sample_sheet_schema.md) - canonical per-project sample sheet: required and optional columns, validation rules, and conflicts with per-stage inputs.
- [`reference/config_schema.md`](reference/config_schema.md) - orchestrator YAML schema: `samples:`, `execution:`, `periodicity:`, `align:`, `rpf:`, and `rnaseq:` sections.
- [`reference/cli.md`](reference/cli.md) - generated CLI reference for every subcommand and flag.

## Output and interpretation references

- [`reference/output_schema.md`](reference/output_schema.md) - column-by-column schema for every TSV / CSV / JSON written under a run root.
- [`reference/warning_codes.md`](reference/warning_codes.md) - stable warning/error codes written to `warnings.tsv`.
- [`reference/periodicity.md`](reference/periodicity.md) - metagene Fourier periodicity QC, `snr_call` tiers, and output interpretation.
- [`rnaseq_te.md`](rnaseq_te.md) - publication boundary between `rnaseq_mode: de_table` and exploratory `rnaseq_mode: from_fastq`.
- [`te_numerics.md`](te_numerics.md) - equations behind `te.tsv` and `delta_te.tsv`.

## Operations and release materials

- [`benchmarking.md`](benchmarking.md) - `mitoribopy benchmark` command, output schema, and runtime/disk sizing guidance.
- [`environment/environment.yml`](environment/environment.yml) - bioconda environment file with the Python stack and external tools.
- [`environment/Dockerfile`](environment/Dockerfile) - container build for reproducible runs.
- [`diagrams/`](diagrams/) - pipeline diagrams for the full workflow and each major stage.
- [`release-notes/`](release-notes/) - version-by-version upgrade notes. The full changelog is [`CHANGELOG.md`](../CHANGELOG.md).
- [`validation/`](validation/) - regression, biological validation, UMI dedup, offset selection, and RNA-seq TE validation notes.
