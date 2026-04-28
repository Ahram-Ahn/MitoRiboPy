# Tutorial 02 — RNA-seq integration for translation-efficiency analysis

`mitoribopy rnaseq` produces translation-efficiency (TE) tables, ΔTE tables, and diagnostic plots from paired RNA-seq + Ribo-seq data. There are two ways to run it:

| Mode | Trigger flag | Use when |
|---|---|---|
| **A — Pre-computed DE** *(original)* | `--de-table de.tsv` | You already ran DESeq2 / Xtail / Anota2Seq externally on the full transcriptome and have the results table. |
| **B — From-FASTQ** *(added in v0.5.0)* | `--rna-fastq …` | You have raw RNA-seq FASTQs and want the subcommand to handle alignment + counting + pyDESeq2 itself. |

The two modes are **mutually exclusive** — passing both flags exits with code 2. Both modes fall through into the **same** TE / ΔTE / plot path, so the downstream output (`te.tsv`, `delta_te.tsv`, plots, `run_settings.json`) is identical in shape; Mode B additionally writes the intermediate counts matrices and a generated `de_table.tsv`.

---

## Mode A — Pre-computed DE table

This mode is the original v0.3.0 design and the right choice when:

- You already have a DE table from the full-transcriptome analysis (which is the right place to run DE — the 13-mt-mRNA universe is too small for shrinkage / dispersion estimation).
- You want strict reference-consistency enforcement against a prior `mitoribopy rpf` run.

### Step A1 — Export the DE table

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

### Step A2 — Choose a gene-ID convention

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

### Step A3 — Set up the reference-consistency gate

The rpf run wrote `results/rpf/run_settings.json` with a `reference_checksum` — the SHA-256 of the FASTA you passed to `-f / --fasta`. `mitoribopy rnaseq` will refuse to run unless the rnaseq-side reference hashes to the same value.

Pass exactly one of:

- `--reference-gtf <path>` — `mitoribopy rnaseq` hashes the file locally and compares against the recorded hash, OR
- `--reference-checksum <sha256>` — useful when the reference file is not on the rnaseq host.

When the rnaseq-side reference does not match the rpf-side hash, `mitoribopy rnaseq` exits 2 with a `MISMATCH` error banner and writes no outputs.

### Step A4 — (Optional) Condition map for replicate-based ΔTE

For replicate-based Ribo-seq log2FC, write a TSV mapping Ribo-seq sample names to conditions:

```text
sample      condition
ctrl_1      control
ctrl_2      control
kd_1        knockdown
kd_2        knockdown
```

Pass it with `--condition-map samples.tsv --condition-a control --condition-b knockdown`. Without this mapping, `mitoribopy rnaseq` computes a point-estimate ΔTE using only the DE table's mRNA log2FC and emits rows with a `single_replicate_no_statistics` note.

### Step A5 — Run

```bash
mitoribopy rnaseq \
  --de-table de.tsv \
  --gene-id-convention hgnc \
  --ribo-dir results/rpf \
  --reference-gtf references/human_mt_transcriptome.fa \
  --condition-map samples.tsv \
  --condition-a control \
  --condition-b knockdown \
  --output results/rnaseq
```

---

## Mode B — From-FASTQ

This mode is right when:

- You have raw RNA-seq + Ribo-seq FASTQs and a transcriptome FASTA, and want the subcommand to handle the whole RNA-side pipeline itself.
- You are comfortable running pyDESeq2 on the small mt-mRNA gene set (acceptable for *exploratory* TE analysis when the goal is "let the math run end-to-end on the same library you aligned"; for publication-grade DE statistics still run DESeq2 / Xtail / Anota2Seq externally on the full transcriptome and use Mode A).

### Step B1 — Install the `[fastq]` extra

pyDESeq2 is a **soft optional dependency**. The pre-computed-DE flow does not need it. For Mode B:

```bash
pip install 'mitoribopy[fastq]'
```

### Step B2 — Lay out your FASTQs

`mitoribopy rnaseq` auto-detects SE vs PE from filename mate tokens, in this precedence order:

1. bcl2fastq: `_R1_001` / `_R2_001` (e.g. `A_S1_L001_R1_001.fastq.gz`)
2. read1/read2: `_read1` / `_read2`
3. R1/R2: `_R1` / `_R2`
4. dot-numbered: `.1.` / `.2.` (e.g. `A.1.fastq.gz`)
5. underscore-numbered: `_1` / `_2`

Both R1 + R2 present → paired sample. Lone R1 → SE sample with the stem name. Lone R2 → SE sample with `_R2` appended to the stem name (so the missing-mate is visible at a glance). `--rna-fastq` and `--ribo-fastq` accept files OR directories (directories are auto-globbed for `.fq`, `.fq.gz`, `.fastq`, `.fastq.gz`).

### Step B3 — Pick a contrast

Mode B uses `--condition-a` and `--condition-b` as the pyDESeq2 contrast. Any conditions in your `--condition-map` that are NOT one of these two are still fitted by DESeq2 (because they contribute to dispersion estimation) but do not appear in the contrast results.

### Step B4 — Auto pseudo-replicates for n=1 designs

Real-world experiments frequently have only one library per condition. pyDESeq2 needs **at least two samples per condition** to estimate dispersion, otherwise it errors out at the `Fitting dispersions` step. To prevent that, when a condition has exactly one sample, `mitoribopy rnaseq` automatically stream-splits its FASTQ by record parity (record N → `rep1` if even, `rep2` if odd) so pyDESeq2 sees n=2.

This is a **mechanical workaround**, not biological replication: the dispersion estimates will be artificially low because the two halves share the exact same library minus the read order. The subcommand prints one stderr WARNING per split so you are never surprised:

```text
[mitoribopy rnaseq] WARNING: RNA-seq condition 'KO' has only 1 sample
('KO_rnaseq'); auto-splitting reads by record parity into pseudo-
replicates 'KO_rnaseq_rep1' / 'KO_rnaseq_rep2'. These are mechanical
halves of the same library, NOT biological replicates — DESeq2
dispersion estimates will be artificially low. Pass
--no-auto-pseudo-replicate to disable.
```

The augmented condition map (original entries + rep1 / rep2 entries) is persisted to `<output>/condition_map.augmented.tsv`. Pass `--no-auto-pseudo-replicate` if you have biological replicates already named correctly in your condition map.

### Step B5 — Run

```bash
mitoribopy rnaseq \
  --rna-fastq input_data/rna_seq/        \
  --ribo-fastq input_data/ribo_seq/      \
  --reference-fasta references/human-mt-mRNA.fasta \
  --gene-id-convention bare              \
  --condition-map samples.tsv            \
  --condition-a control                  \
  --condition-b knockdown                \
  --output results/rnaseq                \
  --align-threads 8
```

The bowtie2 index is content-addressed (`<workdir>/bt2_cache/transcriptome_<sha12_of_fasta>`) so repeated runs against the same FASTA do not rebuild. Pass `--bowtie2-index <prefix>` to skip the build entirely; `--reference-fasta` is still required for the SHA256 we record under `from_fastq.reference_checksum` in `run_settings.json`.

### Mode B caveats

- **PE + UMI is `NotImplementedError`.** Preprocess UMIs into the read name before invoking from-FASTQ mode, or pass `--de-table` with your own pre-computed DE results.
- **DE statistics over 13 genes are noisy.** Mode B is for exploratory pipeline runs and validation. For publication-grade DE, run DESeq2 / Xtail / Anota2Seq externally on the full transcriptome and use Mode A.
- **Reference-consistency gate is skipped.** Mode B has no upstream rpf hash to compare against. Instead the FASTA SHA256 is recorded in `from_fastq.reference_checksum` so a future Mode A `mitoribopy rnaseq` against the same data can verify consistency.

---

## Common to both modes — Outputs

```text
results/rnaseq/
  te.tsv                          # one row per (sample, gene): rpf_count, mrna_abundance, te
  delta_te.tsv                    # one row per gene: mrna_log2fc, rpf_log2fc, delta_te_log2, padj, note
  plots/
    mrna_vs_rpf.png               # four-quadrant log2FC scatter
    delta_te_volcano.png          # ΔTE (log2) vs -log10(padj)
    ma.png                        # log10(baseMean) vs log2FoldChange, coloured by padj
    te_bar_by_condition.png       # log2(TE) per gene, bars grouped by condition + SE
    te_heatmap.png                # gene × sample log2(TE) heatmap (RdBu, centred at 0)
    sample_pca.png                # Mode B only — PC1 vs PC2 from log1p counts
  run_settings.json               # full provenance + (Mode B) per-sample alignment stats
  # --- Mode B only ---
  de_table.tsv                    # pyDESeq2 result, DESeq2 schema; reloads via load_de_table
  rna_counts.tsv                  # wide gene × sample RNA-seq counts
  rpf_counts.tsv                  # long-format Ribo-seq counts
  rpf_counts_matrix.tsv           # wide gene × sample Ribo-seq counts
  condition_map.augmented.tsv     # original + auto-generated rep1 / rep2 names
```

### How to read each plot

- **`mrna_vs_rpf.png`** — four-quadrant scatter. Upper right = translation and transcription both up. Upper left = translation up despite mRNA decrease (buffered-up). Lower right = mRNA increase without matching translation (buffered-down). Lower left = both down.
- **`delta_te_volcano.png`** — genes far from `x=0` with low `padj` are candidates for translation-specific regulation.
- **`ma.png`** — DESeq2 diagnostic. Look for direction (up / down) as a function of expression level; low-expression high-LFC outliers (e.g. `ND6` in mt-Ribo-seq) often dominate the volcano but should be sanity-checked here.
- **`te_bar_by_condition.png`** — the primary biological readout. Bars per (gene, condition) with SE error bars across replicates; `log2(TE)` so y=0 means "no change relative to TE=1". This is usually the figure to put in a presentation.
- **`te_heatmap.png`** — compact summary; columns are sample-ordered by condition so replicates cluster. Useful for spotting outlier samples or genes that move in the same direction across the entire experiment.
- **`sample_pca.png`** *(Mode B)* — PC1 vs PC2 from `log1p` counts. PC1 typically separates RNA from Ribo (different absolute count scales); PC2 typically separates conditions within an assay. Use it as the "did my libraries cluster as expected?" QC.

---

## YAML / `mitoribopy all` integration

Both modes are also drivable from the `rnaseq:` section of a `mitoribopy all` config. See [examples/templates/pipeline_config.example.yaml](../../examples/templates/pipeline_config.example.yaml) for the full schema. The relevant Mode B keys mirror the CLI flags:

```yaml
rnaseq:
  # Mode B trigger:
  rna_fastq:
    - input_data/rna_seq/
  ribo_fastq:
    - input_data/ribo_seq/
  reference_fasta: input_data/human-mt-mRNA.fasta
  gene_id_convention: bare
  condition_map: samples.tsv
  condition_a: control
  condition_b: knockdown
  align_threads: 8
  # bowtie2_index: cache/bt2_cache/transcriptome_5ca397907373   # optional
  # no_auto_pseudo_replicate: false                             # default
```

`mitoribopy all` auto-wires `rnaseq.ribo_dir` to `<run_root>/rpf/` for Mode A, but in Mode B `--ribo-dir` / `--ribo-counts` are not used — Ribo-seq counts are produced from `--ribo-fastq`.

---

## Common pitfalls

- **Mode A: DE statistics over mt only.** Do not run DESeq2 / Xtail / Anota2Seq on a 13-gene universe. Run DE on the full transcriptome and subset for input to MitoRiboPy.
- **Mode A: Mixed gene-ID conventions.** If your DE table uses Ensembl IDs but you pass `--gene-id-convention hgnc`, the per-row match silently fails for every gene. Check the WARNING line.
- **Mode A: Reference MISMATCH.** The SHA256 gate is intentionally strict; even a single trailing newline in the FASTA flips the hash. Re-align one side with the exact file the other side uses.
- **Mode B: pyDESeq2 not installed.** The subcommand emits a clear runtime error pointing at the `[fastq]` extra. Install with `pip install 'mitoribopy[fastq]'`.
- **Mode B: PE + UMI.** Currently `NotImplementedError`. Preprocess UMIs into the read name first, or use Mode A.
- **Mode B: pseudo-replicate inflation.** Auto-split halves are NOT biological replicates and DESeq2's padj will be artificially small. Read the WARNING lines.
