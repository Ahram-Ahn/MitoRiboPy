# RNA-seq + Ribo-seq translation efficiency: validation

## Two modes, two evidence levels

`mitoribopy rnaseq` runs in one of two modes:

| Mode | Inputs | When to use | Evidence quality |
|---|---|---|---|
| `de_table` | external DE table + prior `mitoribopy rpf` outputs | publication, when you've already run DESeq2 / Xtail / Anota2Seq externally on the **full transcriptome** | publication-grade |
| `from_fastq` | raw RNA-seq + Ribo-seq FASTQs + transcriptome FASTA | exploratory: tutorials, smoke tests, quick mt-only scoping | **exploratory only** |

The mode is set via `rnaseq.mode:` in the YAML or via
`--rnaseq-mode {de_table,from_fastq,none}` on the CLI. A run that
infers `from_fastq` from the inputs (no explicit mode) emits a
stderr banner reminding the user that the in-tree pyDESeq2 path is
exploratory.

## What the synthetic-mini test asserts

`tests/test_synthetic_mini.py` constructs WT and KO samples with
known counts and verifies the TE / ΔTE arithmetic against the
package's documented formula:

```
TE     = (RPF + δ) / (mRNA + δ)        # δ = 0.5 (small-count smoothing)
ΔTE    = log2(RPF_ratio) - log2(mRNA_ratio)
```

Concrete check: KO doubles `MT-ND2` reads from 100 to 200 per
sample; the synthetic DE table holds mRNA log2FC at 0; therefore
`delta_te_log2(MT-ND2) == log2(200.5/100.5) ≈ 0.997`, while the
held-equal transcripts must report `delta_te_log2 == 0`.

## Pseudo-replicate gate

The from-FASTQ flow detects conditions with `n == 1` up front and
**fails** unless the user opts in with
`--allow-pseudo-replicates-for-demo-not-publication`. With the opt-in,
each n=1 condition's FASTQ is split by record parity into rep1 / rep2;
the resulting padj / p-values come from mechanical halves of the same
library and are explicitly **not** biologically defensible. The run
writes an `EXPLORATORY.md` sidecar to the rnaseq output dir and a
`pseudo_replicate_mode: true` flag into `run_settings.json` so a
reviewer can never miss the caveat.

This gate is exercised by `tests/test_rnaseq_split_replicates.py`
and by `tests/test_rnaseq_cli_from_fastq.py::*pseudo*`.

## P0.2 reuse path

In `mitoribopy all`, the from-FASTQ rnaseq stage **reuses** the
rpf stage's `rpf_counts.tsv` instead of independently re-aligning
the Ribo FASTQs (the rpf stage already did the work). The orchestrator
threads `--upstream-rpf-counts <path>` to the rnaseq subcommand;
to force a second alignment pass, set
`rnaseq.recount_ribo_fastq: true` in the YAML.

The reuse is exercised end-to-end by
`tests/test_all_rnaseq_reuse.py::test_run_from_fastq_reuse_path_skips_alignment`
(needs the `[fastq]` extra installed; runs in the CI matrix).

## Reference-consistency gate (de_table mode)

The de_table flow does NOT re-align anything. To catch the case
where the external DE table was computed against a different
reference than the one used by the rpf stage, the rnaseq subcommand
hashes the `--reference-gtf` (or accepts `--reference-checksum`) and
verifies it matches the SHA256 the rpf run recorded in its
`run_settings.json`. On mismatch the run hard-fails with a
`MISMATCH` error pointing at both checksums.

Exercised by `tests/test_rnaseq_cli.py::test_rnaseq_reference_mismatch_hard_fails`.
