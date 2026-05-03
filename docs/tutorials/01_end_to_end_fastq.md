# Tutorial 01 — End-to-end mt-Ribo-seq from FASTQ to codon usage

This tutorial walks through a full mt-Ribo-seq analysis with **MitoRiboPy**, starting from raw FASTQ and ending with per-sample translation-profile and codon-usage outputs.

![Pipeline overview](../diagrams/01_pipeline_overview.png)

Every command below is a shell command (lines starting with `$`). Outputs are sketched so you know what to expect at each stage. The tutorial uses a human library mix as the running example; adjust `--strain`, the `-rpf` window, and the FASTA path to match your own data.

> **Prerequisites.** Install MitoRiboPy with `pip install mitoribopy` or from the bioconda environment shipped under `docs/environment/environment.yml`:
>
> ```bash
> $ conda env create -f docs/environment/environment.yml
> $ conda activate mitoribopy
> $ mitoribopy --version
> ```
>
> `cutadapt`, `bowtie2`, `bowtie2-build`, and `umi_tools` (when at least one sample has UMIs) must be on `$PATH`. `samtools` and `fastqc` are optional.

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

FASTA headers in `human_mt_transcriptome.fa` must match the `sequence_name` column of the annotation CSV. MitoRiboPy ships annotation CSVs for `h.sapiens` and `s.cerevisiae` under `src/mitoribopy/data/`.

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
  # Adapter auto-detection (since v0.7.1) is the default and the only
  # kit-name input vector — the user-facing `kit_preset:` knob was
  # removed in v0.7.1. Pin the 3' sequence with `adapter:` when
  # detection cannot identify the library; declare already-trimmed
  # FASTQs with `pretrimmed: true`. The two are mutually exclusive.
  # Detection still names the matched adapter family in
  # kit_resolution.tsv (`detected_kit` / `applied_kit`) for provenance.
  # adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
  # pretrimmed: true
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
  strain: h.sapiens               # h.sapiens | s.cerevisiae | custom
  fasta: references/human_mt_transcriptome.fa
  footprint_class: monosome       # short | monosome | disome | custom
  rpf: [29, 34]                   # h.sapiens monosome window
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
  codon_density_window: true      # smooth codon-density with +/-1 nt window

# Optional rnaseq section -- see Tutorial 02.
# rnaseq:
#   de_table: de.tsv
#   gene_id_convention: hgnc
#   reference_gtf: references/human_mt_transcriptome.fa
```

> Every key under a section maps to the corresponding subcommand's CLI flag with hyphens turned to underscores (so `library_strandedness` → `--library-strandedness`, `adapter_detection` → `--adapter-detection`). Booleans emit the bare flag (`true`) or are omitted entirely (`false`); `null` values are dropped.

You can also start from the exhaustive copy-and-edit template under `examples/templates/`:

```bash
$ cp examples/templates/pipeline_config.example.yaml pipeline_config.yaml
```

or get the curated minimal template from the CLI:

```bash
$ mitoribopy all --print-config-template > pipeline_config.yaml
```

---

## Step 3 — Dry-run the full pipeline

```bash
$ mitoribopy all --config pipeline_config.yaml --output results/ --dry-run
[all] dry-run: planned actions
  1. align: --adapter-detection auto --library-strandedness forward --fastq-dir fastqs/ ...
  2. rpf:   --strain h.sapiens --fasta references/human_mt_transcriptome.fa -rpf 29 34 ...
  3. write manifest to results/run_manifest.json
```

The dry-run prints the exact argv each subcommand will receive, and auto-wires `rpf`'s `--directory` to `results/align/bed/` and `--read_counts_file` to `results/align/read_counts.tsv` so `rpf` ingests the BED6 and counts produced by `align`.

---

## Step 4 — Run end-to-end

```bash
$ mitoribopy all --config pipeline_config.yaml --output results/ --threads 8
```

For multi-sample batches you can align several samples concurrently. `mitoribopy all` does not take `--max-parallel-samples` directly — set it under the `align:` section of your YAML (it is forwarded to the align stage automatically):

```yaml
align:
  # ... existing keys ...
  max_parallel_samples: 4
```

Then re-run with the same `--threads 8`:

```bash
$ mitoribopy all --config pipeline_config.yaml --output results/ --threads 8
```

`--threads` is auto-divided across workers (here, 4 workers × 2 threads/tool ≈ 8 cores total), and the joint `rpf` stage stays serial because offset selection is pooled across all samples. If you only want to run the standalone align stage, the same flag is available on the CLI:

```bash
$ mitoribopy align --max-parallel-samples 4 --threads 8 ...
```

On completion, inspect:

```text
results/
  align/
    mitoribopy.log
    read_counts.tsv             # per-stage counts (input -> trimmed -> contam -> mt -> MAPQ -> dedup)
    kit_resolution.tsv          # per-sample kit + dedup decisions
    run_settings.json           # includes per_sample[] block
    sample_overrides.tsv        # only when align.samples: was set in YAML
    .sample_done/               # per-sample resume markers; --resume skips finished samples
      ctrl_1.json
      kd_1.json ...
    trimmed/                    # *.cutadapt.json (per-sample)
                                # NOTE: *.trimmed.fq.gz is deleted as soon as
                                # contam-filter consumes it; pass
                                # --keep-intermediates to retain it.
    contam_filtered/            # empty unless --keep-intermediates
    aligned/                    # *.mapq.bam (kept)
                                # NOTE: pre-MAPQ *.bam is deleted unless
                                # --keep-intermediates.
    deduped/                    # *.dedup.bam ONLY when at least one sample uses
                                # umi-tools. For dedup=skip samples the orchestrator
                                # wires mapq.bam straight into BED conversion --
                                # no duplicate dedup.bam is written.
    bed/                        # *.bed (strand-aware BED6) -- input to rpf
  rpf/
    mitoribopy.log
    rpf_counts.tsv              # per-sample per-gene RPF counts -> feeds rnaseq
    run_settings.json           # includes reference_checksum (SHA256 of --fasta)
    offset_diagnostics/         # renamed from plots_and_csv/ in v0.4.4
      csv/
        offset_stop.csv         # COMBINED enrichment summary
        p_site_offsets_stop.csv # COMBINED selected offsets
        per_sample_offset/      # renamed from per_sample/
          ctrl_1/
            offset_stop.csv
            p_site_offsets_stop.csv
            offset_applied.csv  # exact offsets row applied to ctrl_1 downstream
          ctrl_2/ ...
      plots/
        offset_drift_stop.svg   # per-sample drift comparison; READ THIS FIRST
        offset_*.svg            # heatmaps + line plots
    translation_profile/        # FLAT (v0.4.4): site is filename prefix
      <sample>/
        footprint_density/      # <gene>_footprint_density.csv (cols: Position,
                                #   Nucleotide, A_site, P_site)
                                # + <gene>_p_site_depth.png and/or
                                #   <gene>_a_site_depth.png
        translating_frame/      # p_site_frame_usage_total.csv,
                                # a_site_frame_usage_total.csv,
                                # *_by_transcript.csv + matching plots
        codon_usage/            # p_site_codon_usage_<gene>.csv,
                                # p_site_codon_usage_total.csv,
                                # a_site_codon_usage_<gene>.csv,
                                # a_site_codon_usage_total.csv + plots
    coverage_profile_plots/     # FLAT (v0.4.4): site is filename prefix
      read_coverage_rpm/        # gated by --read_coverage_rpm
      read_coverage_raw/        # gated by --read_coverage_raw
      read_coverage_rpm_codon/
      read_coverage_raw_codon/
      p_site_density_rpm/
      p_site_density_raw/
      p_site_density_rpm_codon/
      p_site_density_raw_codon/
      p_site_density_rpm_frame/       # frame-coloured CDS-only nt plots (overlay);
      p_site_density_raw_frame/       # frame-0 dominance is the canonical mt-Ribo-seq QC
      p_site_density_rpm_frame_split/ # v0.6.2: per-frame split companion (3 sub-rows
      p_site_density_raw_frame_split/ #         per sample, shared y-axis)
      a_site_density_*/               # A-site mirror folders (analysis_sites in {a,both})
    qc/                               # metagene Fourier QC bundle (always written)
      fourier_spectrum_combined.tsv         # per-(sample, length, gene_set, region) amplitude curve
      fourier_period3_score_combined.tsv    # spectral_ratio_3nt + snr_call tier
      periodicity.metadata.json             # knobs that produced the tables
      fourier_spectrum/<sample>/<sample>_<length>nt_{combined,ATP86,ND4L4}.{png,svg}
      metagene_{start,stop}.tsv             # P-site metagene tables
      metagene_{start,stop}_p_site.svg      # 3-nt periodicity plots
    codon_correlation/          # if --cor_plot
      p_site/<base>_vs_<sample>_*.{csv,svg,png}
      a_site/<base>_vs_<sample>_*.{csv,svg,png}
    igv_tracks/                 # if --igv_export
      <sample>/<sample>_{p_site,a_site}.bedgraph
  run_manifest.json             # composed provenance
```

---

## Step 5 — Interpret the outputs

The shortest path through the outputs:

1. **Did each sample get the right kit?**
   Open `results/align/kit_resolution.tsv`:

   ```text
   sample    fastq                applied_kit         adapter                                umi_length  dedup_strategy  detected_kit       match_rate  source
   ctrl_1    fastqs/ctrl_1.fq.gz  illumina_truseq     AGATCGGAAGAGCACACGTCTGAACTCCAGTCA      0           skip            illumina_truseq    0.9912      detected
   kd_1      fastqs/kd_1.fq.gz    pretrimmed          (none)                                 0           skip            (none)             0.0001      inferred_pretrimmed
   ```

   `source=detected` means the scanner identified the kit. `source=inferred_pretrimmed` means no adapter signal was found at all — typical for SRA-deposited inputs. `source=user_fallback` means detection failed and a user-supplied `--adapter <SEQ>` (or `--pretrimmed`) was used. `source=per_sample_override:*` means an `align.samples:` YAML override pinned the adapter / pretrimmed flag explicitly for this sample.

2. **Did each pipeline stage make biological sense?**
   Open `results/align/read_counts.tsv`. The invariants must hold:

   - `rrna_aligned + post_rrna_filter == post_trim`
   - `mt_aligned + unaligned_to_mt == post_rrna_filter`

   If they don't, the alignment step has a bug or the inputs were mis-formatted.

3. **Are the per-sample offsets consistent?**
   Open `results/rpf/offset_diagnostics/plots/offset_drift_stop.svg`. Each sample has its own per-read-length 5' and 3' offset bar; the combined diagnostic is overlaid as a dashed line. Outliers are visible by eye in seconds. If a single sample drifts by more than ~2 nt at most read lengths, inspect its enrichment heatmap under `offset_diagnostics/csv/per_sample_offset/<sample>/` — usually low coverage is the cause, not a real biological signal. The `offset_applied.csv` in that same folder confirms which offset row downstream actually used for the sample.

4. **Does the canonical 12–15 nt P-site peak show up?**
   Open `results/rpf/offset_diagnostics/plots/offset_stop.svg` (heatmap + line plot). A sharp peak at the canonical 12–15 nt 5' offset confirms library quality. Diffuse / flat enrichment usually means the wrong strand, the wrong RPF window, or low coverage.

5. **Does the CDS show frame-0 dominance?**
   Open `results/rpf/coverage_profile_plots/p_site_density_rpm_frame/<gene>_p-site_density_(rpm,_cds_frame_coloring).svg`. Bars are coloured by reading frame (0 / 1 / 2 relative to CDS start). A healthy library shows ~70–90% of the CDS density on frame 0 with much smaller frames 1 and 2; flat or jittery frames suggest contamination, poor offset selection, or a low-complexity region.

   When the overlay's tallest frame is ambiguous (typical for fused overlapping ORFs like `ATP86`, where the ATP6 ORF is +2 nt offset from ATP8), open the per-frame split companion at `results/rpf/coverage_profile_plots/p_site_density_rpm_frame_split/<gene>_*.svg` instead. Three sub-rows per sample (frame 0, +1, +2) sharing the y-axis make low-frame signal visible even when another frame stacks at the same CDS position.

6. **What's the periodicity QC verdict?**
   Open `results/rpf/qc/fourier_period3_score_combined.tsv` and filter on `gene_set=combined`. The `snr_call` column is the headline tier per `(sample, read_length, region)`: `excellent ≥ 10×`, `healthy ≥ 5×`, `modest ≥ 2×`, `broken < 2×`. A `broken` call means the offsets did not produce codon phasing for that read length; codon-occupancy interpretation is unsafe in that case. The three figures under `qc/fourier_spectrum/<sample>/` (`*_combined.png`, `*_ATP86.png`, `*_ND4L4.png`) report the same metric in-figure on each panel.

   To re-score periodicity with a different window or codon-skip without re-running offset selection:

   ```bash
   $ mitoribopy periodicity \
       --site-table results/rpf/qc/site_table.tsv \
       --output     results/rpf/qc/standalone_periodicity \
       --site p \
       --fourier-window-nt 99
   ```

7. **What's the codon-usage profile?**
   For each sample, the totals are in the same flat folder:

   - `results/rpf/translation_profile/<sample>/codon_usage/p_site_codon_usage_total.csv` — overall P-site codon occupancy.
   - `results/rpf/translation_profile/<sample>/codon_usage/a_site_codon_usage_total.csv` — overall A-site codon occupancy.

   Compare the two side by side to inspect both sites independently.

8. **Want a publication-ready codon scatter?**
   Set `cor_plot: true` and `base_sample: <reference_sample>` in your YAML (or pass `--cor_plot --base_sample <ref>`). With `analysis_sites: both`, you get parallel outputs at `results/rpf/codon_correlation/p_site/` and `.../a_site/`.

9. **Want to view footprint density in IGV?**
   Set `igv_export: true` (or pass `--igv_export`). BedGraph tracks land at `results/rpf/igv_tracks/<sample>/<sample>_{p_site,a_site}.bedgraph` — drop them into IGV alongside your reference FASTA.

---

## Step 6 — Resume a failed run

If `align` finished but `rpf` bailed on a config typo, fix the config and re-run with `--resume`:

```bash
$ mitoribopy all --config pipeline_config.yaml --output results/ --resume
```

`--resume` works at two levels:

1. **Stage level** — each stage is skipped when its final sentinel exists: `align/read_counts.tsv`, `rpf/rpf_counts.tsv`, `rnaseq/delta_te.tsv`.
2. **Per-sample inside align** — every sample that finishes successfully writes `<output>/align/.sample_done/<sample>.json` containing its read-count row. On a subsequent `--resume` (when `read_counts.tsv` is missing because the previous run crashed mid-batch), samples whose marker is present are reloaded from JSON instead of re-run; the markers from completed samples are merged with freshly-processed ones to produce the aggregated `read_counts.tsv`. A 50-sample run that died at sample 30 picks up at sample 31.

Individual stages can also be force-skipped via `--skip-align`, `--skip-rpf`, `--skip-rnaseq`.

---

## Step 7 — Iterating on rpf parameters

When you only want to tweak rpf parameters (offset windows, RPF range, plot format, …) without re-running alignment, point `mitoribopy rpf` directly at the existing BED dir:

```bash
$ mitoribopy rpf \
    -s h.sapiens \
    -f references/human_mt_transcriptome.fa \
    --directory results/align/bed/ \
    -rpf 29 34 \
    --align stop \
    --offset-type 5 \
    --offset-site p \
    --analysis-sites both \
    --min-5-offset 11 --max-5-offset 18 \
    --min-3-offset 11 --max-3-offset 18 \
    --plot-format svg \
    --read-counts-file results/align/read_counts.tsv \
    --rpm-norm-mode mt_mrna \
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
- For non-human / non-yeast organisms, see the [Custom organisms](../../README.md#custom-organisms) section of the README.
