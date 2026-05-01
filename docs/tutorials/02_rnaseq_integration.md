# Tutorial 02 — RNA-seq integration for translation-efficiency analysis

`mitoribopy rnaseq` produces translation-efficiency (TE) tables, ΔTE tables, and diagnostic plots from paired RNA-seq + Ribo-seq data. There are two ways to drive it; the **default flow** is the right starting point for most users.

| Flow | Trigger flag | Use when |
|---|---|---|
| **Default** — from raw FASTQ | `--rna-fastq …` | You have raw RNA-seq FASTQs and want the subcommand to handle alignment + counting + pyDESeq2 itself. |
| **Alternative** — bring your own DE table | `--de-table de.tsv` | You already ran DESeq2 / Xtail / Anota2Seq externally on the full transcriptome and have the results table — the recommended path for publication-grade DE statistics. |

The two flows are **mutually exclusive** — passing both flags exits with code 2. Both fall through into the **same** TE / ΔTE / plot path, so the downstream output (`te.tsv`, `delta_te.tsv`, plots, `run_settings.json`) is identical in shape; the default flow additionally writes intermediate counts matrices, a generated `de_table.tsv`, and the sample PCA plot.

---

## Default flow — from raw FASTQ

### Step 1 — Install the `[fastq]` extra

pyDESeq2 is a **soft optional dependency**. The alternative `--de-table` flow does not need it. For the default flow, install once:

```bash
pip install 'mitoribopy[fastq]'
```

### Step 2 — Lay out your FASTQs

`mitoribopy rnaseq` auto-detects SE vs PE from filename mate tokens, in this precedence order:

1. bcl2fastq: `_R1_001` / `_R2_001` (e.g. `A_S1_L001_R1_001.fastq.gz`)
2. read1/read2: `_read1` / `_read2`
3. R1/R2: `_R1` / `_R2`
4. dot-numbered: `.1.` / `.2.` (e.g. `A.1.fastq.gz`)
5. underscore-numbered: `_1` / `_2`

Both R1 + R2 present → paired sample. Lone R1 → SE sample with the stem name. Lone R2 → SE sample with `_R2` appended to the stem name (so the missing-mate is visible at a glance). `--rna-fastq` and `--ribo-fastq` accept files OR directories (directories are auto-globbed for `.fq`, `.fq.gz`, `.fastq`, `.fastq.gz`).

### Step 3 — Pick a contrast

The default flow uses `--base-sample` (reference / denominator) and `--compare-sample` (comparison / numerator) as the pyDESeq2 contrast. The legacy `--condition-a` / `--condition-b` flags are still accepted as aliases. Any conditions in your `--condition-map` that are NOT one of these two are still fitted by DESeq2 (because they contribute to dispersion estimation) but do not appear in the contrast results.

### Step 4 — n=1 designs (publication-safe default + opt-in fallback)

pyDESeq2 needs **at least two samples per condition** to estimate dispersion, otherwise it errors out at the `Fitting dispersions` step. As of v0.5.2 the default behaviour for n=1 designs is to **fail fast** — the run exits with code 2 and tells you the offending condition(s). This is the publication-safe default: a tutorial run that silently inflates significance is worse than a clear error.

```text
[mitoribopy rnaseq] ERROR: the following condition(s) have only 1 sample
and pyDESeq2 cannot fit dispersion on n=1 designs: 'KO'.
  Resolve by ONE of:
    1. supply biological replicates (recommended);
    2. run external DE and use the --de-table flow;
    3. pass --allow-pseudo-replicates-for-demo-not-publication
       for a non-publication exploratory run (FASTQ-record
       parity halves; padj/p-values are NOT biologically
       defensible).
```

If you genuinely just want a tutorial / smoke-test run, opt in with `--allow-pseudo-replicates-for-demo-not-publication`. The subcommand then stream-splits the FASTQ by record parity (record N → `rep1` if even, `rep2` if odd) and continues, but it stamps `pseudo_replicate_mode: true` in `run_settings.json`, writes an `EXPLORATORY.md` sidecar listing the outputs that are not biologically defensible, and brackets the run with a loud stderr banner. The augmented condition map (original entries + rep1 / rep2 entries) lands at `<output>/condition_map.augmented.tsv`.

The pre-v0.5.2 `--no-auto-pseudo-replicate` flag is still accepted but is now a deprecated no-op (the safe default already does what it asked for). Drop it from your scripts.

### Step 5 — Run

```bash
mitoribopy rnaseq \
  --rna-fastq input_data/rna_seq/        \
  --ribo-fastq input_data/ribo_seq/      \
  --reference-fasta references/human-mt-mRNA.fasta \
  --gene-id-convention bare              \
  --condition-map samples.tsv            \
  --base-sample control                  \
  --compare-sample knockdown             \
  --output results/rnaseq                \
  --align-threads 8
```

The bowtie2 index is content-addressed (`<workdir>/bt2_cache/transcriptome_<sha12_of_fasta>`) so repeated runs against the same FASTA do not rebuild. Pass `--bowtie2-index <prefix>` to skip the build entirely; `--reference-fasta` is still required for the SHA256 we record under `from_fastq.reference_checksum` in `run_settings.json`.

### Default-flow caveats

- **PE + UMI is `NotImplementedError`.** Preprocess UMIs into the read name first (e.g. via `umi_tools extract`), or use the alternative `--de-table` flow with your own pre-computed DE results.
- **DE statistics over 13 mt-mRNAs are noisy.** The default flow is for exploratory pipeline runs and validation. For publication-grade DE, run DESeq2 / Xtail / Anota2Seq externally on the full transcriptome and use the alternative flow.
- **Reference-consistency gate is skipped.** The default flow has no upstream rpf hash to compare against. Instead the FASTA SHA256 is recorded in `from_fastq.reference_checksum` so a future alternative-flow run against the same data can verify consistency.

---

## Alternative flow — bring your own DE table

This flow is the right choice when:

- You already have a DE table from the full-transcriptome analysis (which is the right place to run DE — the 13-mt-mRNA universe is too small for shrinkage / dispersion estimation).
- You want strict reference-consistency enforcement against a prior `mitoribopy rpf` run.

### Step 1 — Export the DE table

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

### Step 2 — Choose a gene-ID convention

The `--gene-id-convention` flag is **required** in both flows. It tells `rnaseq` how to match your DE table's gene-ID column against the mt-mRNA registry:

| Convention | Example human ID |
|---|---|
| `hgnc` | `MT-ND1` |
| `mt_prefixed` | `MT-ND1` |
| `bare` | `ND1` |
| `ensembl` | `ENSG00000198888` |
| `refseq` | `YP_003024026.1` |

For yeast (`--organism y`), only `bare` and `hgnc` are meaningful; both resolve to `COX1, COX2, COX3, COB, ATP6, ATP8, ATP9, VAR1`.

If fewer than 13 (human) or 8 (yeast) mt-mRNAs match the DE table's IDs, `mitoribopy rnaseq` emits a WARNING naming every missing gene. The most common cause is the wrong `--gene-id-convention`; the second most common is that the DE table was subsetted to nuclear-only genes before export.

### Step 3 — Set up the reference-consistency gate

The rpf run wrote `results/rpf/run_settings.json` with a `reference_checksum` — the SHA-256 of the FASTA you passed to `-f / --fasta`. `mitoribopy rnaseq --de-table` will refuse to run unless the rnaseq-side reference hashes to the same value.

Pass exactly one of:

- `--reference-gtf <path>` — `mitoribopy rnaseq` hashes the file locally and compares against the recorded hash, OR
- `--reference-checksum <sha256>` — useful when the reference file is not on the rnaseq host.

When the rnaseq-side reference does not match the rpf-side hash, `mitoribopy rnaseq` exits 2 with a `MISMATCH` error banner and writes no outputs.

### Step 4 — (Optional) Condition map for replicate-based ΔTE

For replicate-based Ribo-seq log2FC, write a TSV mapping Ribo-seq sample names to conditions:

```text
sample      condition
ctrl_1      control
ctrl_2      control
kd_1        knockdown
kd_2        knockdown
```

Pass it with `--condition-map samples.tsv --base-sample control --compare-sample knockdown` (or the legacy `--condition-a` / `--condition-b`; both spellings work). Without this mapping, `mitoribopy rnaseq` computes a point-estimate ΔTE using only the DE table's mRNA log2FC and emits rows with a `single_replicate_no_statistics` note. (The condition map is **required** in the default flow; here it is optional.)

### Step 5 — Run

```bash
mitoribopy rnaseq \
  --de-table de.tsv \
  --gene-id-convention hgnc \
  --ribo-dir results/rpf \
  --reference-gtf references/human_mt_transcriptome.fa \
  --condition-map samples.tsv \
  --base-sample control \
  --compare-sample knockdown \
  --output results/rnaseq
```

---

## Common to both flows — Outputs

```text
results/rnaseq/
  te.tsv                          # one row per (sample, gene): rpf_count, mrna_abundance, te
  delta_te.tsv                    # one row per gene: mrna_log2fc, rpf_log2fc, delta_te_log2, padj, note
  plots/
    # WT-vs-X comparison plots — titles include `<base> vs <compare>` so
    # reviewers know the contrast direction at a glance
    mrna_vs_rpf.png               # four-quadrant log2FC scatter (mRNA vs RPF)
    delta_te_volcano.png          # ΔTE (log2) vs -log10(padj)
    ma.png                        # log10(baseMean) vs log2FoldChange, coloured by padj
    de_volcano_mrna.png           # mRNA DE volcano: red sig-up / blue sig-down /
                                  # grey n.s.; threshold guides at padj=0.05, |L2FC|=1
    de_volcano_rpf.png            # default flow only — Ribo DE volcano from a second
                                  # pyDESeq2 fit on the Ribo-seq subset
    te_bar_by_condition.png       # log2(TE) per gene, bars grouped by condition + SE
    te_heatmap.png                # gene × sample log2(TE) heatmap (RdBu, centred at 0)
    te_compare_scatter.png        # per-gene mean log2(TE) in <base> (x) vs <compare> (y)
                                  # with identity line; needs --condition-map + base/compare
    te_log2fc_bar.png             # sorted bar of log2(TE_compare / TE_base) per gene
    sample_pca.png                # default flow only — PC1 vs PC2 from log1p counts
  run_settings.json               # full provenance + (default flow) per-sample alignment stats
  # --- default flow only ---
  de_table.tsv                    # mRNA pyDESeq2 result, DESeq2 schema; reloads via load_de_table
  rpf_de_table.tsv                # Ribo-seq pyDESeq2 result, same schema; drives de_volcano_rpf;
                                  # absent when the Ribo subset has fewer than two condition levels
  rna_counts.tsv                  # wide gene × sample RNA-seq counts
  rpf_counts.tsv                  # long-format Ribo-seq counts
  rpf_counts_matrix.tsv           # wide gene × sample Ribo-seq counts
  condition_map.augmented.tsv     # original + auto-generated rep1 / rep2 names
```

### Publication style — applies to every plot

Every plot in this set renders under a shared publication style:

- **300 dpi PNG + editable-text SVG sidecar.** Every `*.png` ships with a sibling `*.svg` (same stem) so figures land in Illustrator without becoming path soup (`svg.fonttype = none`).
- **Okabe-Ito colour-blind-safe palette.** Vermillion = up, blue = down, light grey = n.s. — the same semantics across volcanos, MA, and TE plots so the colour reads consistently.
- **White-bbox gene labels with leader lines** + a small smart placer that tries 8 candidate positions per label and prefers above-the-point — handles the small mt-mRNA universe cleanly without dragging in `adjustText`.
- **Consistent contrast labelling.** Every plot title carries `<base> vs <compare>` (driven by `--base-sample` / `--compare-sample`) so a reviewer who lands on a single PNG knows the direction without reading the surrounding markdown.

### How to read each plot

- **`mrna_vs_rpf.png`** — four-quadrant scatter with faint corner captions ("co-regulated up / buffered up / co-regulated down / buffered down") and points coloured by quadrant; identity `y = x` line dashed. Upper right = translation and transcription both up. Upper left = translation up despite mRNA decrease (buffered-up). Lower right = mRNA increase without matching translation (buffered-down). Lower left = both down.
- **`delta_te_volcano.png`** — genes far from `x = 0` with low `padj` are candidates for translation-specific regulation. Threshold guides at `padj < 0.05` and `|ΔTE| ≥ 1`; stat box reports `n_up / n_down / n_total`. TE-up = vermillion, TE-down = blue, n.s. = grey.
- **`ma.png`** — DESeq2 diagnostic. log10(baseMean) (x) vs log2FoldChange (y), sig up = vermillion, sig down = blue, n.s. = grey. Look for direction (up / down) as a function of expression level; low-expression high-LFC outliers (e.g. `ND6` in mt-Ribo-seq) often dominate the volcano but should be sanity-checked here. Sig genes wear bbox-bg labels with leader lines.
- **`de_volcano_mrna.png`** — the WT-vs-X mRNA differential expression volcano. `log2FoldChange` (x) vs `-log10(padj)` (y); points coloured **vermillion** for sig up, **blue** for sig down, **grey** otherwise (Okabe-Ito). Threshold guides at `padj < 0.05` and `|log2FC| ≥ 1`. Stat box (lower-right) reports counts; legend (lower-left) only lists bins that actually have points. Read alongside `delta_te_volcano.png`: if a gene moves in `de_volcano_mrna` but not in `delta_te_volcano`, the change is transcription-driven; the reverse means translation-specific regulation; both moving the same way means co-regulation.
- **`de_volcano_rpf.png`** *(default flow)* — same shape as the mRNA volcano but driven by a second pyDESeq2 fit on the Ribo-seq subset (`rpf_de_table.tsv`). A reviewer can pair it with the mRNA volcano to localise the layer at which a gene's TE shift is happening: footprint-level change (RPF moves alone), mRNA-level (mRNA moves alone), or both. Skipped with a stderr WARNING when the Ribo subset has fewer than two condition levels.
- **`te_bar_by_condition.png`** — the primary biological readout. Bars per (gene, condition) with SE error bars across replicates AND **every replicate's individual TE value drawn as a black dot jittered along x within the bar's footprint**, so the within-group spread is visible by eye and a reviewer can never mistake a 2-rep mean for a many-rep one. `log2(TE)` so y=0 means "no change relative to TE=1". Usually the figure to put in a presentation.
- **`te_heatmap.png`** — compact summary; columns are sample-ordered by condition so replicates cluster, with a coloured **condition strip** above the columns spelling the assignment (KO / RESC / WT etc.) in white-bold so reviewers do not have to read sample names. Cell-value annotations get auto-contrast (white on saturated, black on faint). Useful for spotting outlier samples or genes that move in the same direction across the entire experiment.
- **`te_compare_scatter.png`** — per-gene mean `log2(TE)` in the base condition (x-axis) vs the compare condition (y-axis), with the identity `y = x` line dashed; gene-mean dots are coloured by direction (vermillion above identity, blue below, grey near-identity). Per-replicate values are drawn as a faint cloud behind the means so the within-group spread is visible. **Pearson r** + above/below identity counts in stat box. Emitted whenever `--condition-map` + `--base-sample` + `--compare-sample` are all set (in either flow).
- **`te_log2fc_bar.png`** — sorted bar of `log2(TE_compare / TE_base)` per gene with the **numeric log2FC printed at each bar's tip**. Bars above zero (vermillion) mean TE up in the compare condition; below zero (blue) means down. Direction-coded legend names both branches even when only one is realised in your data. The fastest-readable summary of TE direction across the gene set — great for picking the hits to chase. Same conditions as `te_compare_scatter`.
- **`sample_pca.png`** *(default flow)* — PC1 vs PC2 from `log1p` counts; condition = colour, assay = marker shape. PC1 typically separates RNA from Ribo (different absolute count scales); PC2 typically separates conditions within an assay. Use it as the "did my libraries cluster as expected?" QC.

---

## YAML / `mitoribopy all` integration

Both flows are also drivable from the `rnaseq:` section of a `mitoribopy all` config. See [examples/templates/pipeline_config.example.yaml](../../examples/templates/pipeline_config.example.yaml) for the full schema. The relevant default-flow keys mirror the CLI flags:

```yaml
rnaseq:
  # Default-flow trigger:
  rna_fastq:
    - input_data/rna_seq/
  ribo_fastq:
    - input_data/ribo_seq/
  reference_fasta: input_data/human-mt-mRNA.fasta
  gene_id_convention: bare
  condition_map: samples.tsv
  base_sample: control                  # alias for condition_a; pick one spelling
  compare_sample: knockdown             # alias for condition_b
  align_threads: 8
  # bowtie2_index: cache/bt2_cache/transcriptome_5ca397907373   # optional
  # allow_pseudo_replicates_for_demo_not_publication: false     # default;
  #   set true ONLY for tutorials / smoke tests on n=1 designs
```

`mitoribopy all` auto-wires `rnaseq.ribo_dir` to `<run_root>/rpf/` for the alternative flow; in the default flow `--ribo-dir` / `--ribo-counts` are not used — Ribo-seq counts are produced from `--ribo-fastq`.

---

## Common pitfalls

- **Default flow: pyDESeq2 not installed.** The subcommand emits a clear runtime error pointing at the `[fastq]` extra. Install with `pip install 'mitoribopy[fastq]'`.
- **Default flow: PE + UMI.** Currently `NotImplementedError`. Preprocess UMIs into the read name first, or use the alternative flow.
- **Default flow: pseudo-replicate inflation.** When the n=1 fallback is opted into via `--allow-pseudo-replicates-for-demo-not-publication`, the auto-split halves are NOT biological replicates and DESeq2's padj will be artificially small. Treat the run as exploratory; the `EXPLORATORY.md` sidecar in the output dir lists which fields not to cite.
- **Alternative flow: DE statistics over mt only.** Do not run DESeq2 / Xtail / Anota2Seq on a 13-gene universe. Run DE on the full transcriptome and subset for input to MitoRiboPy.
- **Alternative flow: Mixed gene-ID conventions.** If your DE table uses Ensembl IDs but you pass `--gene-id-convention hgnc`, the per-row match silently fails for every gene. Check the WARNING line.
- **Alternative flow: Reference MISMATCH.** The SHA256 gate is intentionally strict; even a single trailing newline in the FASTA flips the hash. Re-align one side with the exact file the other side uses.
