# Periodicity QC reference

Two complementary signals defend a mt-Ribo-seq run's 3-nt periodicity
claim:

1. **Frame-fraction QC** — pooled and per-read-length tables of how
   often the assigned P-site lands in frame 0, 1, 2 across the CDS.
   Outputs: `frame_counts_by_sample_length.tsv`,
   `frame_counts_by_gene.tsv`, `gene_periodicity.tsv`, `qc_summary.tsv`,
   `qc_summary.md`.
2. **Fourier amplitude spectrum** — the Wakigawa-style discrete
   Fourier transform on the 100-nt window upstream of (ORF) and
   downstream of (3' UTR) every canonical stop codon, computed
   per `(sample, read_length, gene, region)`. The defensible claim is
   "ORF shows a sharp peak at period = 3 nt AND the matched 3' UTR
   does not." Outputs: `fourier_spectrum.tsv`,
   `fourier_period3_score.tsv`, `fourier_period3_summary.tsv`, plus
   per-`(sample, read_length)` two-panel overlay plots under
   `fourier_spectrum/<sample>/`.

This document covers (2). The frame-fraction QC is documented in the
v0.6.2 release notes.

## What replaced what

| v0.6.x (removed) | v0.7.0+ (current) |
|---|---|
| `fft_period3_power.tsv` (single scalar per `(sample, gene)`) | `fourier_spectrum.tsv` (full amplitude curve, period 2-10 nt) + `fourier_period3_score.tsv` (per-`(sample, read_length, gene, region)` scalar) + `fourier_period3_summary.tsv` (per-`(sample, length, region)` combined summary, overlap-upstream genes split out) |
| `--fft-period3` / `--periodicity-fft-period3` | `--fourier-spectrum` / `--periodicity-fourier-spectrum` (default ON), `--fourier-window-nt`, `--fourier-period-min`, `--fourier-period-max`, `--fourier-no-plots` |
| `compute_fft_period3=True` (Python kwarg) | `compute_fourier_spectrum=True` |
| `calculate_fft_period3_ratio()`, `build_fft_period3_table()` | `mitoribopy.analysis.fourier_spectrum.compute_amplitude_spectrum()`, `build_fourier_spectrum_table()`, `build_period3_score_table()`, `build_period3_summary_table()` |

Migration recipe: scripts that used to grab the single-scalar
`fft_period3_power` column should read either
`fourier_period3_score.tsv` (new column: `period3_score`) or, when a
single defensible per-sample number is needed,
`fourier_period3_summary.tsv` → `median_period3_score_combined`.

## Window convention (Wakigawa et al., Mol Cell 2025)

For a transcript with annotated `stop_codon` (the 0-based first nt of
the stop codon trinucleotide):

| Region | nt range (transcript-local) | Width |
|---|---|---|
| `orf`  | `[stop_codon - W, stop_codon)` | W nt strictly upstream of the stop codon |
| `utr3` | `[stop_codon + 3, stop_codon + 3 + W)` | W nt strictly downstream of the stop trinucleotide |

Default W = 100 nt (the published Wakigawa value). The stop codon
trinucleotide itself is never included in either panel. Genes whose
window falls off the transcript end are dropped silently.

## Read assignment

The Wakigawa figure displays the **A-site** position of each read
(`A-site = P-site + 3` nt). The default is preserved here: `--site a`.
For libraries where the P-site is the more diagnostic anchor, pass
`--site p`.

## Spectrum normalization

`amplitude = |FFT[k]| * 2 / N` for `k > 0` — the standard one-sided
amplitude scale that recovers the peak amplitude of a pure sinusoid
from a finite signal. The DC component (`k = 0`) is dropped via
mean-subtraction before the FFT. This is depth-sensitive on raw
counts; downstream consumers that want depth-normalised scores can
divide by `n_sites`.

## Output schema

### `fourier_spectrum.tsv`

One row per `(sample, read_length, gene, region, period_nt)` cell —
the long-format spectrum a plotting tool can re-render directly.

| Column | Type | Description |
|---|---|---|
| `sample` | string | Sample id |
| `read_length` | int | Footprint length (nt) — Wakigawa publishes one figure per length |
| `gene` | string | Transcript id (matches the FASTA / annotation `transcript` field) |
| `transcript_id` | string | Same value as `gene` (kept for downstream tooling that expects both) |
| `region` | `orf` \| `utr3` | Window relative to the canonical stop codon |
| `n_sites` | int | Total assigned sites in the window (depth proxy) |
| `n_nt` | int | Window width in nt (= `--fourier-window-nt`) |
| `period_nt` | float | FFT bin period in nt — non-integer because the bin spacing is `N / k` (e.g. for N=100 the period closest to 3 nt is `100/33 ≈ 3.03`) |
| `amplitude` | float | One-sided amplitude (see normalization above) |
| `is_overlap_upstream_orf` | bool | `True` for `MT-ATP8` / `MT-ND4L` (genes whose stop codon sits inside another ORF) — see Overlap-upstream policy below |

### `fourier_period3_score.tsv`

One row per `(sample, read_length, gene, region)`. Replaces the
legacy single-scalar `fft_period3_power` column.

| Column | Type | Description |
|---|---|---|
| `sample`, `read_length`, `gene`, `transcript_id`, `region`, `n_sites`, `n_nt`, `is_overlap_upstream_orf` | (as above) | |
| `period3_amplitude` | float | Amplitude at the FFT bin closest to period 3 nt |
| `period3_score` | float | `period3_amplitude / mean(amplitude in reference bins 4, 5, 6, 7, 8, 9, 10 nt)`. `NaN` when the reference mean is 0 |

### `fourier_period3_summary.tsv`

One row per `(sample, read_length, region)`. Combines genes EXCLUDING
overlap-upstream and reports the excluded genes separately.

| Column | Type | Description |
|---|---|---|
| `sample`, `read_length`, `region` | | |
| `n_genes_combined` | int | Number of non-overlap-upstream genes scored |
| `median_period3_score_combined` | float | Median `period3_score` over those genes — the **headline number** for "is this sample's translation periodicity healthy?" |
| `max_period3_score_combined` | float | Max `period3_score` — useful for spotting one heroic gene that dominates the median |
| `n_genes_overlap_upstream` | int | Number of overlap-upstream genes (typically 2 for human: ATP8 + ND4L) |
| `median_period3_score_overlap_upstream` | float | Median `period3_score` over the overlap-upstream genes — these typically also score high in the 3' UTR panel because the "3' UTR" window is actually inside another ORF |
| `max_period3_score_overlap_upstream` | float | Max over the overlap-upstream genes |

## Plot bundle

Per `(sample, read_length)`, two figures land under
`<output>/rpf/qc/fourier_spectrum/<sample>/`:

* `<sample>_<read_length>nt_combined.{png,svg}` — Wakigawa-style
  stacked two-panel overlay (ORF top, 3' UTR bottom), one trace per
  gene EXCEPT the overlap-upstream pair. This is the publication-grade
  figure.
* `<sample>_<read_length>nt_overlap_upstream.{png,svg}` — same layout
  but only `MT-ATP8` and `MT-ND4L` (or whatever overlap-upstream genes
  are present). Skipped when none are in the input. Its purpose is
  diagnostic: a reviewer can confirm the overlap-upstream 3' UTR
  panel really does show period-3 signal (because the "3' UTR" window
  is inside the next ORF), which validates the rationale for excluding
  these genes from the combined panel.

Each PNG carries the canonical per-plot `.metadata.json` sidecar so
`mitoribopy validate-figures` can score it without re-running
matplotlib.

## Overlap-upstream policy

The annotated stop codon of `MT-ATP8` sits inside the `MT-ATP6` ORF;
the annotated stop codon of `MT-ND4L` sits inside the `MT-ND4` ORF.
For these two genes the 100-nt window labelled "3' UTR" is not really
3' UTR — it's another ORF being translated, so the 3' UTR panel will
show period-3 signal that has nothing to do with stop-codon
termination behaviour. Wakigawa et al. handle this by excluding
both genes from the combined figure entirely.

This package preserves Wakigawa's exclusion in the **combined data**
(combined panel of the plot, `*_combined` columns of the summary
table) but **always scores the genes individually** so a reviewer can
look at them deliberately:

* per-gene rows in `fourier_spectrum.tsv` and
  `fourier_period3_score.tsv` include `MT-ATP8` and `MT-ND4L` with
  `is_overlap_upstream_orf=true`;
* the combined-summary TSV reports them in `*_overlap_upstream`
  columns separate from the `*_combined` columns;
* a separate `*_overlap_upstream.png` figure renders just these
  genes.

The recognised gene names are case-insensitive and treat `_` as `-`:
`MT-ATP8`, `MTATP8`, `ATP8`, `MT-ND4L`, `MTND4L`, `ND4L`. Fused-FASTA
spellings (`ATP86`, `ND4L4`) are deliberately NOT excluded — those
entries cover the whole combined region and have a legitimate 3' UTR
downstream of the natural transcript-end stop codon.

To override the exclusion (analyse the overlap-upstream genes as part
of the combined data anyway), bypass the predicate by editing the
score / summary tables directly with pandas — the per-gene rows have
the boolean flag set so the filtering is one line of code.

## Reading the figure

A healthy mt-Ribo-seq library produces:

1. A **sharp peak at period 3 nt in the ORF panel** for every gene
   trace (the published amplitude reaches ~5-12 in the Wakigawa
   figures, but the absolute scale is depth-dependent — compare across
   genes within a sample, not across samples).
2. A **flat 3' UTR panel** with no period-3 peak (any signal there
   means non-translating mitoribosome footprints landed in the 3'
   UTR — which is the diagnostic for the published `mtRF1L KO`
   stop-codon read-through story).
3. The **overlap-upstream panel** shows period-3 in BOTH panels (the
   "3' UTR" window is really inside the next ORF, so it carries
   translating-ribosome signal).

If the ORF panel is flat or the 3' UTR panel matches the ORF panel
amplitude, the offset assignment for that read length is suspect.
Cross-check `frame_counts_by_sample_length.tsv` → `qc_call` for the
same length.
