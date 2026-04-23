# Tutorial 02 &mdash; RNA-seq integration for translation-efficiency analysis

`mitoribopy rnaseq` consumes a pre-computed differential-expression
(DE) table from DESeq2, Xtail, or Anota2Seq plus a prior
`mitoribopy rpf` run and produces translation-efficiency (TE) tables,
delta-TE tables, and diagnostic plots.

This tutorial continues from Tutorial 01. It assumes:

- You have run `mitoribopy rpf` (or `mitoribopy all` through the rpf
  stage) and `results/rpf/rpf_counts.tsv` and
  `results/rpf/run_settings.json` exist.
- Your RNA-seq data has been aligned to the **identical** transcript
  reference as the Ribo-seq side and passed through DESeq2 / Xtail /
  Anota2Seq externally. MitoRiboPy does not run DE itself because the
  13-mt-mRNA universe violates the shrinkage assumptions those methods
  need; run DE on the full transcriptome and subset.

## Step 1 &mdash; Export the DE table

Most users run this in R:

```r
# DESeq2 example
res <- results(dds)
write.table(
  as.data.frame(res) |> tibble::rownames_to_column("gene_id"),
  file = "de.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
```

Xtail and Anota2Seq outputs also work as-is; `mitoribopy rnaseq
--de-format auto` detects the format from the column headers.

## Step 2 &mdash; Choose a gene-ID convention

The `--gene-id-convention` flag is **required** (no default). It tells
rnaseq how to match your DE table's `gene_id` column against the
mt-mRNA registry:

| Convention      | Example human ID    |
|-----------------|---------------------|
| `hgnc`          | `MT-ND1`            |
| `mt_prefixed`   | `MT-ND1`            |
| `bare`          | `ND1`               |
| `ensembl`       | `ENSG00000198888`   |
| `refseq`        | `YP_003024026.1`    |

For yeast, only `bare` and `hgnc` are meaningful (they resolve to the
same set: `COX1`, `COX2`, `COX3`, `COB`, `ATP6`, `ATP8`, `ATP9`,
`VAR1`). Pass `--organism y`.

## Step 3 &mdash; Set up the reference-consistency gate

The rpf run you completed in Tutorial 01 wrote
`results/rpf/run_settings.json` with a `reference_checksum` &mdash; the
SHA-256 of the FASTA you passed to `-f / --fasta`. rnaseq will refuse
to run unless the rnaseq-side reference hashes to the same value.

Pass either:

- `--reference-gtf <path>` &mdash; rnaseq hashes the file locally and
  compares against the recorded hash, OR
- `--reference-checksum <sha256>` &mdash; useful when the reference
  file is not on the rnaseq host.

## Step 4 &mdash; (Optional) Condition map for replicate-based &Delta;TE

For replicate-based Ribo-seq log2FC, write a TSV that maps Ribo-seq
sample names to conditions:

```text
sample      condition
ctrl_1      control
ctrl_2      control
kd_1        knockdown
kd_2        knockdown
```

Pass it with `--condition-map <path> --condition-a control
--condition-b knockdown`. Without this mapping, rnaseq computes a
point-estimate &Delta;TE using only the DE table's mRNA log2FC and
emits rows with a `single_replicate_no_statistics` note.

## Step 5 &mdash; Run the integration

```bash
$ mitoribopy rnaseq \
    --de-table de.tsv \
    --gene-id-convention hgnc \
    --ribo-dir results/rpf \
    --reference-gtf references/human_mt_transcriptome.fa \
    --condition-map samples.tsv \
    --condition-a control \
    --condition-b knockdown \
    --output results/rnaseq
```

Outputs:

```text
results/rnaseq/
  te.tsv              # one row per (sample, gene) with rpf_count, mrna_abundance, te
  delta_te.tsv        # one row per gene: mrna_log2fc, rpf_log2fc, delta_te_log2, padj, note
  plots/
    mrna_vs_rpf.png        # four-quadrant log2FC scatter
    delta_te_volcano.png   # delta-TE (log2) vs -log10(padj)
  run_settings.json
```

## Step 6 &mdash; Interpret

- **Four-quadrant scatter.** Each dot is one mt-mRNA. Upper right =
  translation and transcription both up. Upper left = translation up
  despite mRNA decrease (buffered-up). Lower right = mRNA increase
  without matching translation (buffered-down). These are the classic
  mt-translation dysregulation phenotypes.
- **delta-TE volcano.** Genes far from `x=0` with low `padj` are
  candidates for translation-specific regulation. With single-replicate
  data the y axis is just zeros and `note` columns read
  `single_replicate_no_statistics`; use the table's
  `delta_te_log2` column directly and validate orthogonally.

## Common pitfalls

- **Missing mt-mRNAs.** If fewer than 13 human (or 8 yeast) mt-mRNAs
  match the DE table's gene IDs, rnaseq emits a WARNING naming every
  missing gene. Most common cause: wrong `--gene-id-convention`. Second
  most common: your DE table was subsetted to nuclear-only genes before
  export.
- **Reference MISMATCH.** If the rnaseq-side reference does not match
  the rpf-side hash, rnaseq exits 2 with a prominent error. Re-align
  one side so both use the same transcript set, or verify the correct
  `--reference-gtf` was passed.
- **Zero-count genes.** A 0.5 Laplace pseudocount is added to both
  RPF and mRNA values so genes with genuinely zero coverage don't
  break the math. Filter those genes out of the analysis downstream;
  their TE is uninformative.
- **DE statistics over mt only.** Do not run DESeq2 / Xtail / Anota2Seq
  on a 13-gene universe; the shrinkage / dispersion assumptions need a
  large background. Run DE on the full transcriptome and subset.
