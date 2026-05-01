# Documentation

The repo's top-level [`README.md`](../README.md) is the canonical, standalone reference: pipeline overview, every CLI flag, input formats, custom-organism workflow, and example invocations. The files under `docs/` are the deeper material — tutorials, an at-a-glance CLI reference, environment recipes, and the diagram set.

## For users

- [`tutorials/01_end_to_end_fastq.md`](tutorials/01_end_to_end_fastq.md) — step-by-step walk through `mitoribopy all` from raw FASTQ to per-sample translation profile and codon usage.
- [`tutorials/02_rnaseq_integration.md`](tutorials/02_rnaseq_integration.md) — the optional `mitoribopy rnaseq` stage: DE-table integration, SHA256 reference-consistency gate, TE / ΔTE outputs.
- [`tutorials/05_hpc_cluster_run.md`](tutorials/05_hpc_cluster_run.md) — running on SLURM / LSF: thread accounting, scratch storage, disk budgeting, hash-validated resume.
- [`reference/cli.md`](reference/cli.md) — at-a-glance per-subcommand flag list. Run `mitoribopy <subcommand> --help` for the full per-flag help text.
- [`reference/sample_sheet_schema.md`](reference/sample_sheet_schema.md) — canonical per-project sample sheet: required / optional columns, validation rules, conflicts with per-stage inputs.
- [`reference/config_schema.md`](reference/config_schema.md) — orchestrator YAML schema: `samples:` / `align:` / `rpf:` / `rnaseq:` sections, the `rnaseq.mode` key, auto-wired cross-stage defaults.
- [`environment/environment.yml`](environment/environment.yml) — bioconda environment file with every external tool (cutadapt, bowtie2, umi_tools, samtools, …).
- [`environment/Dockerfile`](environment/Dockerfile) — container build for reproducible runs.
- [`diagrams/`](diagrams/) — four large-format PNG diagrams of the pipeline:
  - [`01_pipeline_overview.png`](diagrams/01_pipeline_overview.png) — horizontal LR view of all three stages.
  - [`02_align_stage.png`](diagrams/02_align_stage.png) — internals of `mitoribopy align`.
  - [`03_rpf_stage.png`](diagrams/03_rpf_stage.png) — internals of `mitoribopy rpf`.
  - [`04_rnaseq_stage.png`](diagrams/04_rnaseq_stage.png) — internals of the optional `mitoribopy rnaseq` stage.
  - Regenerate with `python diagrams/render_diagrams.py` (matplotlib only; no Node / mermaid-cli required).

## For maintainers

- [`release-notes/`](release-notes/) — version-by-version highlights aimed at users upgrading. The authoritative changelog is at [`CHANGELOG.md`](../CHANGELOG.md) at the repo root.
- [`validation/synthetic_mini.md`](validation/synthetic_mini.md) — known-answer integration test on a programmatically-generated fixture; runs in pytest.
- [`validation/taco1_ko_regression.md`](validation/taco1_ko_regression.md) — biological regression test plan (TACO1-KO polyproline stalling) the package must pass before tagging a release candidate.
- [`validation/public_dataset_reanalysis.md`](validation/public_dataset_reanalysis.md) — reproducible smoke recipe for re-running MitoRiboPy on a small public mt-Ribo-seq dataset.
- [`validation/offset_selection_validation.md`](validation/offset_selection_validation.md) — how the per-sample 5'/3' offset is picked, with confidence-label expectations.
- [`validation/umi_dedup_validation.md`](validation/umi_dedup_validation.md) — UMI source labels and the conservative inference policy.
- [`validation/rnaseq_te_validation.md`](validation/rnaseq_te_validation.md) — TE / ΔTE arithmetic, the `de_table` vs `from_fastq` mode split, the pseudo-replicate gate, and the reference-consistency gate.
