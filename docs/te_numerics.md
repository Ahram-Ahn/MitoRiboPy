# What the numbers mean — RNA, RPF, TE, ΔTE

Translation efficiency is one of the most-asked-for and most-misread
quantities in Ribo-seq. This page pins each quantity to a specific
output file and a specific equation so reviewers don't have to guess
what `te.tsv` means versus a coverage-ratio plot.

For the column-by-column schema of every output mentioned here, see
[`docs/reference/output_schema.md`](reference/output_schema.md). For
the publication-grade boundary between the `de_table` and `from_fastq`
modes (and the `--strict` gates that protect it), see
[`docs/rnaseq_te.md`](rnaseq_te.md).

---

## The four quantities

| # | Quantity | Definition | What it is NOT |
|---|---|---|---|
| 1 | **RNA abundance** | Per-sample, per-transcript RNA-seq counts (or a normalised abundance estimate downstream of pyDESeq2). | Not a coverage profile, not a Ribo-seq quantity. |
| 2 | **RPF abundance** | Per-sample, per-transcript ribosome footprint counts after offset-aware filtering. | Not raw aligned-read counts (the read-length auto-filter prunes noise bins) and not a per-position density. |
| 3 | **TE (translation efficiency)** | Per-sample, per-gene RPF abundance normalised by RNA abundance. See equation below. | Not an inferential statistic — TE is a point estimate per (sample, gene). |
| 4 | **ΔTE (delta-TE)** | Per-gene log2 change in TE between two conditions. See equation below. | Not the ratio of two TEs computed sample-by-sample; it is a contrast on the log scale. |

---

## TE equation (`te.tsv`)

For sample *s* and gene *g*:

```
TE(s, g) = (RPF_count(s, g) + δ) / (mRNA_abundance(g) + δ)
```

where δ is the package's pseudocount (default 0.5; see
[`src/mitoribopy/rnaseq/te.py`](../src/mitoribopy/rnaseq/te.py)). The
pseudocount goes on **both** numerator and denominator to avoid
div-by-zero and log(0). Per-gene mRNA abundance is the DE table's
`baseMean` column (DESeq2 / Xtail / Anota2Seq) or a pyDESeq2-derived
estimate in the from-FASTQ flow. Genes missing from the DE table are
skipped — TE is undefined without an mRNA denominator.

Output: `<output>/rnaseq/te.tsv`. Schema details are in
[`docs/reference/output_schema.md`](reference/output_schema.md#rnaseqtetsv-schema-20).

---

## ΔTE equation (`delta_te.tsv`)

For gene *g*, comparing condition *B* to condition *A* (with *A* =
`--base-sample` / `--condition-a`):

```
ΔTE_log2(g) = log2(RPF_FC) - log2(mRNA_FC)

         where  RPF_FC  = mean(RPF_B) / mean(RPF_A)   (with pseudocount on each mean)
                mRNA_FC = exp2(mRNA_log2fc(g))         (from the DE table)
```

ΔTE is **not** computed as `TE(B) / TE(A)` per sample. It is a
contrast on log-fold-changes: each side's denominator (RNA
abundance) and numerator (RPF abundance) are aggregated separately
before the ratio. This is the standard formulation used by Xtail,
Anota2Seq, and downstream summary statistics, and it is what removes
the spurious sample-level pseudocount asymmetry you'd get from a
per-sample TE ratio.

Output: `<output>/rnaseq/delta_te.tsv`. The `note` column carries
any qualifier MitoRiboPy attached to the row:

| `note` value | Meaning |
|---|---|
| `""` (empty) | Both RPF and mRNA log2FC available; ΔTE is computed from replicate-based means. |
| `single_replicate_no_statistics` | No condition map / contrast was supplied (or only one replicate per condition); ΔTE row carries only the mRNA log2FC and `rpf_log2fc=None`. |
| `insufficient_ribo_replicates` | Condition map present but the gene has zero Ribo counts in one of the two conditions; `rpf_log2fc=None` and ΔTE is `None`. |
| `missing_from_de_table` | Gene is in the Ribo counts but absent from the DE table; everything except `gene` is `None`. |

---

## What is NOT in te.tsv / delta_te.tsv

These two files are **gene-level summary tables**. If you want any
of the following, look elsewhere:

| You want | Look at |
|---|---|
| Per-position P-site / A-site density | `<output>/translation_profile/<sample>/footprint_density/<transcript>_footprint_density.csv` |
| Per-codon usage by sample | `<output>/translation_profile/<sample>/codon_usage/` |
| Periodicity / 3-nt phasing diagnostic | The metagene Fourier bundle under `<output>/rpf/qc/` — full reference at [`docs/reference/periodicity.md`](reference/periodicity.md). |
| The descriptive coverage-normalised RPF metric (the legacy "rna-seq ratio" the v0.2.x module emitted) | **Removed in this refactor.** Use `te.tsv` for sample-level TE, or compute your own ratio from the per-position density CSVs above. |

---

## Pseudo-replicate runs

When a from-FASTQ flow run is launched with
`--allow-pseudo-replicates-for-demo-not-publication`, both `te.tsv`
and `delta_te.tsv` are still written, but `run_settings.json`
carries `pseudo_replicate_mode: true` and an `EXPLORATORY.md`
sidecar lists every column you should NOT cite (padj, "significant"
markers, dispersion estimates). The TE and ΔTE point estimates
remain readable as exploratory; only the inferential statistics are
unsafe. See [`docs/rnaseq_te.md`](rnaseq_te.md) for the full
publication-boundary reference and the `--strict` gates that refuse
this mode by default.

---

## See also

* [`docs/rnaseq_te.md`](rnaseq_te.md) — publication-boundary
  reference (mode comparison, strict-mode gates).
* [`docs/reference/output_schema.md`](reference/output_schema.md)
  — column-by-column schemas for `te.tsv` and `delta_te.tsv`.
* [`docs/reference/periodicity.md`](reference/periodicity.md) — the
  metagene Fourier QC bundle and its statistical hardening.
