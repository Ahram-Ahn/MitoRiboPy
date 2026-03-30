# MitoRiboPy: Mitochondrial Ribo-seq Analysis Pipeline

MitoRiboPy processes mitochondrial ribosome profiling (Ribo-seq) BED files and generates:
- Read-length QC summaries
- Start/stop codon offset enrichment tables and plots
- Selected offsets by read length
- Footprint density tracks (P-site, A-site, E-site)
- Frame usage summaries
- Codon usage summaries
- IGV-style RPM/raw coverage plots
- Optional VARNA-ready exports
- Optional codon-correlation plots between samples
- Optional RNA-seq integration

This repository is currently in transition from script-style analysis code to the formal Python package **MitoRiboPy**. The full packaging scheme is documented in:
- [`docs/PACKAGE_REFACTOR_SCHEME.md`](docs/PACKAGE_REFACTOR_SCHEME.md)

---

## 1) Current Status

What is stable now:
- Human and yeast pipelines run from CLI with explicit arguments.
- P-site vs A-site behavior is implemented and validated.
- mt-mRNA-only RPM normalization is available.
- Offset selection strategy now supports canonical P-site-linked picking for robust A/P coupling.
- Phase II package migration is complete for active pipeline modules.
- Phase III readability cleanup is complete for the package path used by the CLI.

What is in progress:
- Expanding tests and CI.
- Formal public API docs.

---

## 2) Requirements

- Python 3.10+
- pandas
- numpy
- matplotlib
- seaborn
- biopython
- scipy

If you are using a clean environment:

```bash
python -m venv .venv
source .venv/bin/activate
pip install pandas numpy matplotlib seaborn biopython scipy
```

---

## 3) Quick Start

### Human example (stop-aligned, P-site, mt-mRNA normalization)

```bash
python main.py \
  -s h \
  -f human_riboseq_mt_entire_genome_withND6.fasta \
  -rpf 29 34 \
  --directory human_test_input \
  --align stop \
  --offset_type 5 \
  --offset_site p \
  --offset_pick_reference p_site \
  --min_offset 10 \
  --max_offset 22 \
  --read_counts_file human_read_counts.txt \
  --read_counts_sample_col sample \
  --read_counts_reads_col read_count \
  --read_counts_reference_col reference \
  --rpm_norm_mode mt_mrna \
  --mrna_ref_patterns mt_genome \
  --codon_overlap_mode full \
  --order_samples NT3 NT5 \
  --cap_percentile 0.999 \
  --plot_format svg \
  -m \
  --output analysis_results_human
```

### Enable optional VARNA and codon correlation

```bash
python main.py ... -v --cor_plot --base_sample NT3
```

### Phase I package entrypoint

During Phase I, package CLI calls are forwarded to the current legacy parser:

```bash
PYTHONPATH=src python -m mitoribopy --help
PYTHONPATH=src python -m mitoribopy -s h -f human.fa --directory input_dir
```

After editable install, use:

```bash
mitoribopy --help
```

---

## 4) Core CLI Options

### Input and organism
- `--fasta`: FASTA file path (required)
- `--directory`: BED input directory
- `--strain y|h`: yeast or human
- `-rpf <min> <max>`: read-length filter range

### Offset analysis
- `--align start|stop`: codon anchor for offset enrichment
- `--offset_type 5|3`: downstream offset type
- `--offset_site p|a`: interpret selected offsets as P-site or A-site
- `--offset_pick_reference p_site|selected_site`:
  - `p_site` (recommended): pick offsets in canonical P-site space, then transform for A-site if needed
  - `selected_site` (legacy): pick directly in selected site space
- `--min_offset`, `--max_offset`: absolute offset filter range for selection
- `--codon_overlap_mode full|any`: full codon overlap recommended

### Normalization and read-count table
- `--read_counts_file`: count table path
- `--read_counts_sample_col`: sample column name
- `--read_counts_reads_col`: read-count column name
- `--rpm_norm_mode total|mt_mrna`:
  - `total`: denominator = total reads in count file
  - `mt_mrna`: denominator = rows matching mt-mRNA reference patterns
- `--read_counts_reference_col`: reference column name for `mt_mrna` mode
- `--mrna_ref_patterns`: list of substrings for mt-mRNA row matching

### Outputs and plotting
- `--output`: output directory
- `--plot_format png|pdf|svg`
- `--plot_dir`: subdirectory for offset CSV/plots
- `--line_plot_style combined|separate`
- `--cap_percentile`: upper clipping percentile for coverage plots
- `-m/--merge_density`: merge adjacent nt into codon-center density

### Optional modules
- `-v/--varna`: generate VARNA exports
- `--cor_plot --base_sample <sample>`: codon-correlation plots
- `--use_rna_seq`: run RNA-seq module (requires RNA args)

---

## 5) Input File Formats

### BED files
Expected columns (minimum 3):
1. `chrom` (transcript name)
2. `start` (0-based)
3. `end` (end-exclusive)
Additional columns are tolerated.

### FASTA
- FASTA headers must match transcript names in BED/annotation logic.

### Read-count file
Delimited text (CSV/TSV autodetected). Typical columns:
- sample column (`sample`, `Sample`, etc.)
- read-count column (`read_count`, `Reads`, etc.)
- optional reference column (`reference`) for `mt_mrna` mode

Example:

```text
sample,reference,read_count
NT3,mt_genome.aligned,1685015
NT3,sapiens_mito_rRNA.aligned,3895562
...
```

---

## 6) Output Structure

Typical run directory:

```text
<output>/
  plots_and_csv/
    offset_<align>.csv
    p_site_offsets_<align>.csv
    offset_enrichment_heatmap_<align>.<fmt>
    offset_enrichment_combined_lineplots_<align>.<fmt>
    read_length_summary.csv
    unfiltered_read_length_summary_15_50.csv
    unfiltered_heatmap_15_50_count.svg
    unfiltered_heatmap_15_50_RPM.svg

  <sample>/
    footprint_density/
      <transcript>_footprint_density.csv
      <transcript>_P_site_depth.png or <transcript>_A_site_depth.png
    translating_frame/
      frame_usage_total.csv
      frame_usage_by_transcript.csv
      frame_usage_total_plot.png
      frame_usage_by_transcript_plot.png
    codon_usage/
      codon_usage_total.csv
      a_site_codon_usage_total.csv
      codon_usage_total_plot.png

  igv_style_plots/
    read_coverage_rpm/
    p_site_coverage_rpm/
    read_coverage_raw/
    p_site_coverage_raw/

  varna/                      # if --varna
  codon_correlation/          # if --cor_plot
  rna_seq_results/            # if --use_rna_seq
```

---

## 7) Coordinate and Biological Assumptions

- BED is treated as standard 0-based, end-exclusive.
- Read length is `end - start`.
- `start_codon` is derived from transcript annotation UTR lengths.
- `stop_codon` is represented as the first base of stop codon in 0-based indexing.
- P/A/E site assignment uses selected offsets and `offset_type` rules.

---

## 8) P-site vs A-site Behavior

- `--offset_site p`: selected offset corresponds to P-site position.
- `--offset_site a`: selected offset corresponds to A-site position.
- In downstream analysis, primary output site follows `offset_site`.

### Recommended robust mode
Use:
- `--offset_pick_reference p_site`

Why:
- Offset picking happens in canonical P-site geometry first.
- A-site selected offsets are transformed from P-site picks.
- This preserves codon-consistent +3/-3 A/P relationship.

Legacy behavior can be restored with:
- `--offset_pick_reference selected_site`

Detailed refinement notes:
- [`OFFSET_PICKING_REFINEMENT_2026-03-29.txt`](OFFSET_PICKING_REFINEMENT_2026-03-29.txt)

---

## 9) Reproducible Debug Runs

For fast, deterministic testing, create subsampled BED files:

```bash
python subsample_bed.py \
  --input NT3.bed \
  --output debug_subsample/NT3_subsample.bed \
  --n 200000 \
  --seed 42
```

---

## 10) Configuration Files

You can run with a JSON config file and override from CLI.

Template:
- `pipeline_config.example.json`

Use:

```bash
python main.py --config pipeline_config.example.json --fasta <path_to_fasta>
```

---

## 11) Troubleshooting

### No reads survive filtering
- Check `-rpf` range and input BED coordinate validity.
- Inspect `<output>/plots_and_csv/read_length_summary.csv`.

### RPM looks too small/large
- Confirm `--rpm_norm_mode` choice.
- Confirm sample names in BED and read-count table match.
- For mt-only denominator, set `--rpm_norm_mode mt_mrna` and correct `--mrna_ref_patterns`.

### A-site/P-site offsets look inconsistent
- Use `--offset_pick_reference p_site`.
- Verify `--min_offset/--max_offset` are biologically appropriate for your data.

### Missing plots for some transcripts
- Confirm transcript names match between BED and FASTA.

---

## 12) Development and Contribution

Current best practice for edits:
- Keep behavior-preserving refactors separate from scientific-logic changes.
- Add docstrings and type hints for public functions.
- Prefer explicit names over short abbreviations.
- Keep plotting style centralized in `plot_style.py` until package split completes.

Planned packaging roadmap:
- [`docs/PACKAGE_REFACTOR_SCHEME.md`](docs/PACKAGE_REFACTOR_SCHEME.md)

---

## 13) Citation and Versioning (Planned)

For public package release, this repository will include:
- Semantic versioning
- Changelog
- Citation metadata
- DOI-friendly release notes
