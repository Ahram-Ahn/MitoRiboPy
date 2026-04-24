# Tutorial 02 — RNA-seq integration for translation-efficiency analysis

`mitoribopy rnaseq` consumes a pre-computed differential-expression (DE) table from DESeq2, Xtail, or Anota2Seq plus a prior `mitoribopy rpf` run and produces translation-efficiency (TE) tables, ΔTE tables, and diagnostic plots.

This tutorial continues from [Tutorial 01](01_end_to_end_fastq.md). It assumes:

- You have run `mitoribopy rpf` (or `mitoribopy all` through the `rpf` stage) and `results/rpf/rpf_counts.tsv` and `results/rpf/run_settings.json` exist.
- Your RNA-seq data has been aligned to the **identical** transcript reference as the Ribo-seq side and passed through DESeq2 / Xtail / Anota2Seq externally. MitoRiboPy does not run DE itself because the 13-mt-mRNA universe (8 in yeast) violates the shrinkage assumptions those methods need; run DE on the full transcriptome and subset.

---

## Step 1 — Export the DE table

Most users run this in R. DESeq2 example:

```r
res <- results(dds)
write.table(
  as.data.frame(res) |> tibble::rownames_to_column("gene_id"),
  file = "de.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
```

Xtail and Anota2Seq outputs also work as-is; `mitoribopy rnaseq --de-format auto` (default) detects the format from the column headers.

For arbitrary DE outputs, switch to `--de-format custom` and supply the column names with `--de-gene-col`, `--de-log2fc-col`, `--de-padj-col`, and (optionally) `--de-basemean-col`.

---

## Step 2 — Choose a gene-ID convention

The `--gene-id-convention` flag is **required** (no default). It tells `rnaseq` how to match your DE table's gene-ID column against the mt-mRNA registry:

| Convention | Example human ID |
|---|---|
| `hgnc` | `MT-ND1` |
| `mt_prefixed` | `MT-ND1` |
| `bare` | `ND1` |
| `ensembl` | `ENSG00000198888` |
| `refseq` | `YP_003024026.1` |

For yeast (`--organism y`), only `bare` and `hgnc` are meaningful; both resolve to `COX1, COX2, COX3, COB, ATP6, ATP8, ATP9, VAR1`.

If fewer than 13 (human) or 8 (yeast) mt-mRNAs match the DE table's IDs, `mitoribopy rnaseq` emits a WARNING naming every missing gene. The most common cause is the wrong `--gene-id-convention`; the second most common is that the DE table was subsetted to nuclear-only genes before export.

---

## Step 3 — Set up the reference-consistency gate

The rpf run from Tutorial 01 wrote `results/rpf/run_settings.json` with a `reference_checksum` — the SHA-256 of the FASTA you passed to `-f / --fasta`. `mitoribopy rnaseq` will refuse to run unless the rnaseq-side reference hashes to the same value.

Pass exactly one of:

- `--reference-gtf <path>` — `mitoribopy rnaseq` hashes the file locally and compares against the recorded hash, OR
- `--reference-checksum <sha256>` — useful when the reference file is not on the rnaseq host.

When the rnaseq-side reference does not match the rpf-side hash, `mitoribopy rnaseq` exits 2 with a `MISMATCH` error banner and writes no outputs. Re-align one side so both use the same transcript set, or verify the correct `--reference-gtf` was passed.

---

## Step 4 — (Optional) Condition map for replicate-based ΔTE

For replicate-based Ribo-seq log2FC, write a TSV that maps Ribo-seq sample names to conditions:

```text
sample      condition
ctrl_1      control
ctrl_2      control
kd_1        knockdown
kd_2        knockdown
```

Pass it with `--condition-map samples.tsv --condition-a control --condition-b knockdown`.

Without this mapping, `mitoribopy rnaseq` computes a point-estimate ΔTE using only the DE table's mRNA log2FC and emits rows with a `single_replicate_no_statistics` note.

---

## Step 5 — Run the integration

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
  te.tsv              # one row per (sample, gene): rpf_count, mrna_abundance, te
  delta_te.tsv        # one row per gene: mrna_log2fc, rpf_log2fc, delta_te_log2, padj, note
  plots/
    mrna_vs_rpf.png        # four-quadrant log2FC scatter
    delta_te_volcano.png   # ΔTE (log2) vs -log10(padj)
  run_settings.json
```

You can also drive `rnaseq` from the same YAML as the rest of the pipeline by adding an `rnaseq:` section and running through `mitoribopy all`:

```yaml
# Add to your pipeline_config.yaml
rnaseq:
  de_table: de.tsv
  gene_id_convention: hgnc
  reference_gtf: references/human_mt_transcriptome.fa
  condition_map: samples.tsv
  condition_a: control
  condition_b: knockdown
```

`mitoribopy all` auto-wires `rnaseq.ribo_dir` to `<run_root>/rpf/`, so you don't need to set it.

---

## Step 6 — Interpret

- **Four-quadrant scatter (`mrna_vs_rpf.png`).** Each dot is one mt-mRNA. Upper right = translation and transcription both up. Upper left = translation up despite mRNA decrease (buffered-up). Lower right = mRNA increase without matching translation (buffered-down). Lower left = both down. These are the classic mt-translation dysregulation phenotypes.
- **ΔTE volcano (`delta_te_volcano.png`).** Genes far from `x=0` with low `padj` are candidates for translation-specific regulation. With single-replicate data (no `--condition-map`) the y axis carries only zeros and the `note` column reads `single_replicate_no_statistics`; in that case use the table's `delta_te_log2` column directly and validate orthogonally.
- **Zero-count genes.** A 0.5 Laplace pseudocount is added to both RPF and mRNA values so genes with genuinely zero coverage don't break the math. Filter those genes out of the analysis downstream — their TE is uninformative.

---

## Common pitfalls

- **DE statistics over mt only.** Do not run DESeq2 / Xtail / Anota2Seq on a 13-gene universe; the shrinkage / dispersion assumptions need a large background. Run DE on the full transcriptome and subset for input to MitoRiboPy.
- **Mixed gene-ID conventions.** If your DE table uses Ensembl IDs but you pass `--gene-id-convention hgnc`, the per-row match will silently fail for every gene. Check the WARNING line.
- **Reference MISMATCH.** The SHA256 gate is intentionally strict; even a single trailing newline change in the FASTA flips the hash. Re-align one side with the exact file the other side uses.
