# RNA-seq integration & translation efficiency: publication boundaries

This page is the user-facing reference for the rnaseq stage of
MitoRiboPy. Its purpose is to draw a clear line between the two
supported flows and to document why one of them is **not**
publication-grade. For the regression-test details that pin the gates
described here see
[`docs/validation/rnaseq_te_validation.md`](validation/rnaseq_te_validation.md).

---

## TL;DR — pick the right mode for your claim

| Want to publish | Use | Why |
|---|---|---|
| Translation efficiency / ΔTE in a manuscript | **`rnaseq_mode: de_table`** with a full-transcriptome external DE table (DESeq2, Xtail, Anota2Seq) | Full-transcriptome dispersion estimation. The mt-mRNA subset is too small to fit dispersion well on its own. |
| Tutorial run, smoke test, mt-only scoping | `rnaseq_mode: from_fastq` | Convenience. Self-contained, no external pipeline needed — at the cost of dispersion that uses only the 13 mt-mRNAs. |
| One-off, quick sanity check on padj / log2FC | `rnaseq_mode: de_table` if you have a DE table; `from_fastq` otherwise — but **do not cite** `from_fastq` p-values | The point estimates of `from_fastq` log2FC are usually fine; the p-values are the danger. |

> **Publication warning.** The built-in `from_fastq` mode runs
> pyDESeq2 on the **mitochondrial mRNA subset only**. With 13 mt
> transcripts (or fewer, after filtering), dispersion estimation is
> unreliable, the trend curve frequently does not converge, and the
> resulting `padj` values can look real even when the underlying
> evidence does not justify them. Use it for tutorials, smoke runs,
> and exploratory scoping — not for a manuscript figure.
>
> For manuscript-grade DE / TE / ΔTE, run a full-transcriptome
> DESeq2 / Xtail / Anota2Seq externally, save the per-gene table, and
> feed it to MitoRiboPy via `rnaseq_mode: de_table` plus
> `rnaseq.de_table: <path>`. The mt-mRNA rows from that
> full-transcriptome fit are the ones to cite.

---

## How `--strict` enforces the boundary

`mitoribopy all --strict` and `mitoribopy validate-config --strict`
both refuse two configurations that lead to publication-unsafe output:

1. **Pseudo-replicates.** Any of these spellings being truthy is a
   strict-mode hard error before any stage runs:

   ```yaml
   rnaseq:
     allow_pseudo_replicates: true                                  # the parser dest
     # OR
     allow_pseudo_replicates_for_demo_not_publication: true         # the long alias
   ```

   Pseudo-replicates split a single library's FASTQ records into
   parity-defined halves to satisfy pyDESeq2's two-replicate
   minimum. The resulting `padj` is a function of how many records
   landed on the even side, not of biological variance.

2. **`rnaseq_mode: from_fastq` without an explicit override.** Under
   strict mode the orchestrator refuses `from_fastq` unless the user
   sets:

   ```yaml
   rnaseq:
     rnaseq_mode: from_fastq
     allow_exploratory_from_fastq_in_strict: true   # opt-in, not for publication
   ```

   The override exists so a tutorial config can still run under
   `--strict` for CI purposes; it is named to be unmistakably
   non-publication.

Both gates are validated by
[`tests/test_all_strict_rnaseq_gates.py`](../tests/test_all_strict_rnaseq_gates.py)
and (since v0.7.1) duplicated in `mitoribopy validate-config --strict`
so a CI script that only runs the validator catches the same
misconfigurations without needing to start the orchestrator.

---

## What `from_fastq` actually does

When you do want to use it (e.g. tutorial, mt-only scoping):

1. cutadapt + bowtie2 alignment of RNA and Ribo FASTQs against the
   mt-transcriptome FASTA.
2. Per-transcript read counting on each side (mt mRNAs only).
3. pyDESeq2 fit on the joint matrix — **on the mt-mRNA rows only**.
4. TE / ΔTE / plots.

The pyDESeq2 fit emits:

* `mean_disp = NaN` warnings when the mt-only sample is too small.
* `dispersion trend curve fitting did not converge` for low-replicate
  designs; the code falls back to a mean-based dispersion trend.
* `SmallSampleWarning` from pyDESeq2 when any subset is below its
  minimum-sample threshold; downstream contrast statistics for that
  subset are NaN.

These warnings are *not* spurious — they are the symptoms of the
mt-only design and are the reason `--strict` refuses the mode.

For the orchestrated path, `mitoribopy all` reuses the rpf stage's
`rpf_counts.tsv` instead of re-aligning the Ribo FASTQs (P0.2 reuse
path); set `rnaseq.recount_ribo_fastq: true` to force a second pass.

---

## What `de_table` mode expects

```yaml
rnaseq:
  rnaseq_mode: de_table
  de_table: /path/to/external_de.tsv     # full-transcriptome DESeq2 / Xtail / Anota2Seq
  ribo_dir: results/rpf/                  # the rpf stage's output, used for ribo counts
  reference_gtf: /path/to/reference.fa    # SHA256-checked against rpf-side checksum
  condition_map: /path/to/conditions.tsv
  condition_a: control
  condition_b: knockdown
```

The DE table needs the standard DESeq2 columns (`baseMean`,
`log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`) plus a gene-id
column matching the convention in `rpf_counts.tsv` (set
`gene_id_convention: hgnc | ensembl | refseq | bare` to match).
Aliases are accepted; see
[`docs/reference/sample_sheet_schema.md`](reference/sample_sheet_schema.md)
and the rnaseq subcommand `--help` for the full mapping.

A **reference-consistency gate** hashes the `--reference-gtf` (or
accepts a precomputed `--reference-checksum`) and refuses to run when
it disagrees with the SHA256 the rpf stage stamped into its
`run_settings.json`. This catches the silent-corruption case where
the external DE was computed against a different reference build.

---

## Output files

The two modes produce a partially overlapping output set:

| File | `de_table` | `from_fastq` | What it carries |
|---|---|---|---|
| `te.tsv` | ✓ | ✓ | per-gene TE = (RPF + δ) / (mRNA + δ), δ = 0.5 |
| `delta_te.tsv` | ✓ | ✓ | log2(RPF_ratio) − log2(mRNA_ratio) |
| `te_volcano.png` / `.svg` | ✓ | ✓ | volcano of the contrast |
| `de_table.tsv` | ✓ (passthrough) | ✓ (in-tree fit) | per-gene DE statistics |
| `rna_counts.tsv` | — | ✓ | per-transcript RNA counts (mt only) |
| `rpf_counts_matrix.tsv` | — | ✓ when no upstream reuse | per-transcript RPF counts (mt only) |
| `EXPLORATORY.md` | — | ✓ when pseudo-replicate mode | sidecar listing fields not to cite |
| `run_settings.json` | ✓ | ✓ | mode, padj policy, pseudo-replicate flag, reference checksum |

When `from_fastq` runs without enough replicates for the chosen
contrast, the in-tree fit returns `padj = NA` rather than synthesising
a publication-shaped value.

---

## See also

* [`docs/validation/rnaseq_te_validation.md`](validation/rnaseq_te_validation.md)
  — the validation evidence and the regression tests that pin the
  gates described above.
* [`docs/reference/sample_sheet_schema.md`](reference/sample_sheet_schema.md)
  — the unified sample-sheet schema (replaces the legacy
  `condition_map` for orchestrator runs).
* [`README.md`](../README.md) — the run-mode quickstart and the
  publication-readiness profile.
