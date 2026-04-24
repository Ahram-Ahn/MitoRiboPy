# Documentation

Pointers to the package documentation. The repo's top-level [`README.md`](../README.md) is the canonical, standalone description of the current version (v0.4.1) including pipeline overview, every CLI flag, input formats, and example invocations.

## For users

- [`tutorials/01_end_to_end_fastq.md`](tutorials/01_end_to_end_fastq.md) — step-by-step walk through `mitoribopy all` from raw FASTQ to per-sample translation profile and codon usage.
- [`tutorials/02_rnaseq_integration.md`](tutorials/02_rnaseq_integration.md) — the optional `mitoribopy rnaseq` stage: DE-table integration, SHA256 reference-consistency gate, TE / ΔTE outputs.
- [`reference/cli.md`](reference/cli.md) — at-a-glance per-subcommand flag list. Run `mitoribopy <subcommand> --help` for the full per-flag help text.
- [`environment/environment.yml`](environment/environment.yml) — bioconda environment file with every external tool (cutadapt, bowtie2, umi_tools, samtools, …).
- [`environment/Dockerfile`](environment/Dockerfile) — container build for reproducible runs.

## For maintainers

- [`release-notes/`](release-notes/) — version-by-version highlights aimed at users upgrading. The authoritative changelog is at [`CHANGELOG.md`](../CHANGELOG.md) at the repo root.
- [`validation/taco1_ko_regression.md`](validation/taco1_ko_regression.md) — biological regression test plan (TACO1-KO polyproline stalling) the package must pass before tagging a release candidate.
- [`diagrams/`](diagrams/) — Mermaid sources for the pipeline architecture diagram.
