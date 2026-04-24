# Tutorial 01 — End-to-end mt-Ribo-seq from FASTQ to codon usage

This tutorial walks through a full mt-Ribo-seq analysis with **MitoRiboPy v0.4.1**, starting from raw FASTQ and ending with per-sample translation-profile and codon-usage outputs.

Every command below is a shell command (lines starting with `$`). The outputs are sketched so you know what to expect at each stage. The tutorial uses a human library mix as the running example; adjust `--strain`, the `-rpf` window, and the FASTA path to match your own data.

> **Prerequisites.** Install MitoRiboPy (see the repo [README](../../README.md)) and the external tools it shells out to. The quickest way is the bioconda environment shipped under `docs/environment/environment.yml`:
>
> ```bash
> $ conda env create -f docs/environment/environment.yml
> $ conda activate mitoribopy
> $ mitoribopy --version   # expect "MitoRiboPy 0.4.1"
> ```
>
> `cutadapt`, `bowtie2`, `bowtie2-build`, and `umi_tools` (when at least one sample has UMIs) must be on `$PATH`. `samtools`, `fastqc`, and `picard` are optional.

---

## Step 0 — Data layout

Create a workspace with the following shape:

```text
<project_root>/
  fastqs/
    ctrl_1.fq.gz
    ctrl_2.fq.gz
    kd_1.fq.gz
    kd_2.fq.gz
  references/
    human_mt_transcriptome.fa      # one FASTA record per mt-mRNA
    human_rrna_contam.fa           # 12S, 16S, cytoplasmic rRNA spike
  pipeline_config.yaml
```

FASTA headers in `human_mt_transcriptome.fa` must match the `sequence_name` column of the annotation CSV. MitoRiboPy ships annotation CSVs for human and yeast under `src/mitoribopy/data/`.

Mixed batches are fine — you can have raw FASTQs, pre-trimmed FASTQs (e.g. SRA-deposited), UMI samples, and non-UMI samples all in the same `fastqs/` directory. Each sample is resolved independently.

---

## Step 1 — Build the bowtie2 indexes

```bash
$ mkdir -p indexes
$ bowtie2-build --quiet references/human_rrna_contam.fa       indexes/rrna_contam
$ bowtie2-build --quiet references/human_mt_transcriptome.fa  indexes/mt_tx
```

Each command produces six `*.bt2` files sharing a common prefix. The prefix (e.g. `indexes/mt_tx`) is what you pass to MitoRiboPy, not an individual `.bt2` file.

---

## Step 2 — Write a shared pipeline config

Config files can be YAML, JSON, or TOML. The following YAML is the most common shape for `mitoribopy all`:

```yaml
align:
  # `auto` (default) detects the kit per sample by scanning each FASTQ.
  # An explicit preset becomes a per-sample fallback used only when
  # detection fails. Pre-trimmed FASTQs are auto-detected and routed
  # through cutadapt with no -a flag.
  kit_preset: auto                # auto | illumina_smallrna | illumina_truseq |
                                  # illumina_truseq_umi | qiaseq_mirna |
                                  # pretrimmed | custom
  adapter_detection: auto         # auto | off | strict
  library_strandedness: forward
  fastq: fastqs/                  # directory string (auto-globs *.fq, *.fq.gz, ...)
  contam_index: indexes/rrna_contam
  mt_index: indexes/mt_tx
  mapq: 10
  min_length: 15
  max_length: 45
  dedup_strategy: auto            # umi-tools per sample if UMI, else skip

rpf:
  strain: h
  fasta: references/human_mt_transcriptome.fa
  rpf: [29, 34]
  align: stop
  offset_type: "5"
  offset_site: p                  # selection coordinate space (P-site)
  offset_pick_reference: p_site
  offset_mode: per_sample         # each sample uses its own offsets
  analysis_sites: both            # write BOTH P-site and A-site outputs
  min_5_offset: 10
  max_5_offset: 22
  min_3_offset: 10
  max_3_offset: 22
  offset_mask_nt: 5
  plot_format: svg
  merge_density: true

# Optional rnaseq section — see Tutorial 02.
# rnaseq:
#   de_table: de.tsv
#   gene_id_convention: hgnc
#   reference_gtf: references/human_mt_transcriptome.fa
```

> Every key under a section maps to the corresponding subcommand's CLI flag with hyphens turned to underscores (so `kit_preset` → `--kit-preset`, `library_strandedness` → `--library-strandedness`). Booleans emit the bare flag (`true`) or are omitted entirely (`false`); `null` values are dropped.

You can also generate a fresh, fully-commented template:

```bash
$ mitoribopy all --print-config-template > pipeline_config.yaml
```

---

## Step 3 — Dry-run the full pipeline

```bash
$ mitoribopy all --config pipeline_config.yaml --output results/ --dry-run
[all] dry-run: planned actions
  1. align: --kit-preset auto --library-strandedness forward --fastq-dir fastqs/ ...
  2. rpf: --strain h --fasta references/human_mt_transcriptome.fa -rpf 29 34 ...
  3. write manifest to results/run_manifest.json
```

The dry-run prints the exact argv each subcommand will receive, and auto-wires `rpf`'s `--directory` to `results/align/bed/` and `--read_counts_file` to `results/align/read_counts.tsv` so `rpf` ingests the BED6 and counts produced by `align`.

---

## Step 4 — Run end-to-end

```bash
$ mitoribopy all --config pipeline_config.yaml --output results/ --threads 8
```

On completion, inspect:

```text
results/
  align/
    mitoribopy.log
    read_counts.tsv             # per-stage counts (input → trimmed → contam → mt → MAPQ → dedup)
    kit_resolution.tsv          # per-sample kit + dedup decisions; the v0.4 provenance spine
    run_settings.json           # includes per_sample[] block
    trimmed/                    # *.trimmed.fq.gz, *.cutadapt.json
    contam_filtered/            # *.nocontam.fq.gz
    aligned/                    # *.bam, *.mapq.bam
    deduped/                    # *.dedup.bam
    bed/                        # *.bed (strand-aware BED6)
  rpf/
    mitoribopy.log
    rpf_counts.tsv              # per-sample per-gene RPF counts → feeds rnaseq
    run_settings.json           # includes reference_checksum (SHA256 of --fasta)
    plots_and_csv/
      offset_stop.csv           # COMBINED enrichment summary (diagnostic)
      p_site_offsets_stop.csv   # COMBINED selected offsets (diagnostic)
      offset_drift_stop.svg     # per-sample drift comparison; READ THIS FIRST
      offset_*.svg              # heatmaps + line plots
      per_sample/
        ctrl_1/                 # per-sample enrichment + selected offsets
        ctrl_2/
        kd_1/
        kd_2/
    translation_profile_p/      # P-site outputs (analysis_sites=both)
      ctrl_1/
        footprint_density/      # *_footprint_density.csv (P/A/E columns)
        translating_frame/      # frame_usage_total.csv, frame_usage_by_transcript.csv
        codon_usage/            # codon_usage_<gene>.csv, codon_usage_total.csv, ...
        debug_csv/
      ...
    translation_profile_a/      # A-site outputs (analysis_sites=both); same shape
      ...
    coverage_profile_plots_p/   # per-site coverage plots
      read_coverage_rpm/, p_site_coverage_rpm/, *_codon/, *_frame/
    coverage_profile_plots_a/
      ...
  run_manifest.json             # composed provenance
```

---

## Step 5 — Interpret the outputs

The shortest path through the outputs:

1. **Did each sample get the right kit?**
   Open `results/align/kit_resolution.tsv`:

   ```text
   sample    fastq           applied_kit         adapter                                umi_length  dedup_strategy  detected_kit       match_rate  source
   ctrl_1    fastqs/ctrl_1.fq.gz   illumina_truseq    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA      0           skip            illumina_truseq    0.9912      detected
   kd_1      fastqs/kd_1.fq.gz     pretrimmed         (none)                                 0           skip            (none)             0.0001      inferred_pretrimmed
   ```

   `source=detected` means the scanner identified the kit. `source=inferred_pretrimmed` means no adapter signal was found at all — typical for SRA-deposited inputs. `source=user_fallback` means detection failed and the user-supplied `--kit-preset` was used.

2. **Did each pipeline stage make biological sense?**
   Open `results/align/read_counts.tsv`. The invariants must hold:

   - `rrna_aligned + post_rrna_filter == post_trim`
   - `mt_aligned + unaligned_to_mt == post_rrna_filter`

   If they don't, the alignment step has a bug or the inputs were mis-formatted.

3. **Are the per-sample offsets consistent?**
   Open `results/rpf/plots_and_csv/offset_drift_stop.svg`. Each sample has its own per-read-length 5' and 3' offset bar; the combined diagnostic is overlaid as a dashed line. Outliers are visible by eye in seconds. If a single sample drifts by more than ~2 nt at most read lengths, inspect its enrichment heatmap under `plots_and_csv/per_sample/<sample>/` — usually low coverage is the cause, not a real biological signal.

4. **Does the canonical 12–15 nt P-site peak show up?**
   Open `results/rpf/plots_and_csv/offset_stop.svg` (heatmap + line plot). A sharp peak at the canonical 12–15 nt 5' offset confirms library quality. Diffuse / flat enrichment usually means the wrong strand, the wrong RPF window, or low coverage.

5. **What's the codon-usage profile?**
   For each sample, the totals are at:

   - `results/rpf/translation_profile_p/<sample>/codon_usage/codon_usage_total.csv` — overall P-site codon occupancy.
   - `results/rpf/translation_profile_a/<sample>/codon_usage/codon_usage_total.csv` — overall A-site codon occupancy.

   Compare the two side by side to inspect both sites independently.

---

## Step 6 — Resume a failed run

If `align` finished but `rpf` bailed on a config typo, fix the config and re-run with `--resume`:

```bash
$ mitoribopy all --config pipeline_config.yaml --output results/ --resume
```

`--resume` skips a stage when its sentinel output file already exists (`align/read_counts.tsv`, `rpf/rpf_counts.tsv`, `rnaseq/delta_te.tsv`). Individual stages can also be force-skipped via `--skip-align`, `--skip-rpf`, `--skip-rnaseq`.

---

## Step 7 — Iterating on rpf parameters

When you only want to tweak rpf parameters (offset windows, RPF range, plot format, …) without re-running alignment, point `mitoribopy rpf` directly at the existing BED dir:

```bash
$ mitoribopy rpf \
    -s h \
    -f references/human_mt_transcriptome.fa \
    --directory results/align/bed/ \
    -rpf 29 34 \
    --align stop \
    --offset_type 5 \
    --offset_site p \
    --analysis_sites both \
    --min_5_offset 11 --max_5_offset 18 \
    --min_3_offset 11 --max_3_offset 18 \
    --plot_format svg \
    --read_counts_file results/align/read_counts.tsv \
    --rpm_norm_mode mt_mrna \
    --output results/rpf_v2/
```

This skips the expensive trim + align + dedup steps entirely.

---

## Strandedness + ND5/ND6 note

Human mt-ND5 and mt-ND6 are transcribed from opposite strands and their 3' ends overlap at the genome level. On Path A (transcriptome reference, the MitoRiboPy default) every mt-mRNA is its own FASTA record, so the overlap is not a concern at alignment time. `--library-strandedness forward` (the default) additionally passes `--norc` to bowtie2 so reverse-complement hits are rejected. The resulting `bed/<sample>.bed` should show only `+` strand rows; any `-` strand rows are worth investigating because they suggest either a library-prep mismatch or a reverse-stranded kit that should have been run with `--library-strandedness reverse`.

---

## What's next

- For the optional translation-efficiency integration with paired RNA-seq, continue to [Tutorial 02 — RNA-seq integration](02_rnaseq_integration.md).
- For the full per-flag reference, see [docs/reference/cli.md](../reference/cli.md).
- For changes per release, see [docs/release-notes/](../release-notes/).
