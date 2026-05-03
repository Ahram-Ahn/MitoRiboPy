# Periodicity QC reference

The mt-Ribo-seq 3-nt periodicity QC bundle is a metagene Fourier
analysis: per-gene tracks are normalised, aggregated, and a single
direct DFT is evaluated on the metagene.

The bundle consists of:

* **Metagene start / stop profiles** — `metagene_start.tsv`,
  `metagene_stop.tsv`, plus the start- and stop-anchored `.svg` plots.
  These are the community-standard "show me 3-nt phasing around the
  start / stop codon" plots and remain useful as a sanity check for
  offset assignment.
* **Metagene Fourier spectrum** — three figures per
  `(sample, read_length)` and two TSVs per run: per-(sample, length,
  gene_set, region) amplitude curve over period 2-10 nt, plus the
  derived 3-nt spectral ratio with an `snr_call` tier. This is the
  headline QC artefact.

The legacy frame-fraction QC bundle (`qc_summary.tsv`,
`frame_counts_*.tsv`, `gene_periodicity.tsv`,
`frame_fraction_heatmap.svg`, `read_length_periodicity_barplot.svg`,
`gene_phase_score_dotplot.svg`) was retired in v0.8.0. The
`spectral_ratio_3nt` + `snr_call` columns now carry the headline QC
story; the per-frame breakdown is no longer reported.

## Why aggregate-then-DFT instead of per-gene overlay

The previous v0.7.x implementation overlaid one FFT trace per gene per
panel, which produced a noisy, hard-to-read plot. Per-gene DFT is
structurally underpowered for several reasons:

* **Initiation pile-ups create flat broadband spectra.** The first
  ~5 codons after every AUG carry massive ribosome occupancy
  (initiation stalling). The Fourier transform of an impulse is a
  flat spectrum across all frequencies — which smears period-3 energy
  uniformly across periods 2-10. Skipping the first 5 codons (15 nt)
  of every transcript window removes that contaminant.
* **Stop-codon stalling does the same on the 3' end.** Especially
  severe for noncanonical-stop ORFs (MT-CO1: AGA; MT-ND6: AGG).
  Including the stop trinucleotide adds another impulse, so the
  last codon (3 nt) before the stop is dropped.
* **Per-gene 100 nt is structurally underpowered.** 33 codons of
  signal per gene is too short for clean DFT on noisy data.
  Aggregating ~9 ORFs into one metagene before the DFT averages out
  gene-specific noise while reinforcing codon-locked 3-nt phasing.
* **Coverage magnitude dominates raw overlays.** The highest-
  expression gene (MT-ND6 in many libraries) buries every other
  trace. Per-gene unit-mean normalisation equalises contributions
  before aggregation.
* **Hann window + mean-centering kill spectral leakage.** Without
  them, the DC component leaks into nearby frequencies and creates
  an upward amplitude ramp toward periods 8-10 (the "rising tail"
  that the legacy figure showed).

## Window convention

For a transcript with annotated `start_codon` and `stop_codon`
(each the 0-based first nt of its trinucleotide) and window width `W`:

| Region | nt range | Notes |
|---|---|---|
| `orf_start` | `[start_codon + 15, start_codon + 15 + W)` | Skip 5 codons of initiation peak |
| `orf_stop`  | `[stop_codon - 3 - W, stop_codon - 3)` | Skip 1 codon before stop |

Default `W = 99 nt = 33 codons`. Multiple of 3 — no period-3 leakage.

A 3' UTR negative-control window is intentionally NOT computed: human
mt-mRNA 3' UTRs are typically too short for a meaningful 99-nt window,
and the user-facing artefact is cleaner without a sparse third panel.
Override the window via `--periodicity-fourier-window-nt`.

## Per-gene processing pipeline

For each (sample, read_length, transcript, region):

```
1. Build P-/A-site coverage vector of length W in window-local coords.
2. Filter: skip if mean(coverage) < 0.1 or sum(coverage) < 30.
3. x = coverage / mean(coverage)        # unit-mean normalisation
4. x = x - mean(x)                      # mean-centre (kills DC leakage)
5. x = x * np.hanning(W)                # Hann taper
```

## Aggregation: three gene_sets per (sample, length)

Three figures get rendered per (sample, read_length):

### `combined`

Element-wise mean of the per-gene normalised tracks across every
transcript that is **not** part of the ATP8/ATP6 or ND4L/ND4 overlap
pair. Typically 9 mt-mRNAs in human (`MT-CO1`, `MT-CO2`, `MT-CO3`,
`MT-CYB`, `MT-ND1`, `MT-ND2`, `MT-ND3`, `MT-ND5`, `MT-ND6`).

### `ATP86` (junction-bracketed)

The ATP8 / ATP6 bicistronic transcript carries two ORFs in different
reading frames that overlap at nt 162-205 of the ATP86 transcript.
The figure has TWO panels, each interrogating a different frame at
the junction:

* **`orf_stop` panel — ATP8 frame.** Window =
  `[ATP8.stop - 3 - 99, ATP8.stop - 3)` = `[103, 202)`. Captures the
  end of ATP8's reading frame, including the start of the bicistronic
  overlap. A period-3 peak here means ribosomes are reading in ATP8's
  frame.
* **`orf_start` panel — ATP6 frame.** Window =
  `[ATP6.start + 15, ATP6.start + 15 + 99)` = `[177, 276)`. Captures
  the start of ATP6's reading frame, just past the bicistronic
  overlap. A period-3 peak here means ribosomes are reading in
  ATP6's frame.

The two windows OVERLAP at `[177, 202)` — exactly the bicistronic
overlap region. The ATP86 figure is the only place in the bundle
where this region is interrogated separately in each frame.

### `ND4L4` (junction-bracketed)

Same idea for the ND4L / ND4 bicistronic pair. The 4-nt overlap is
too short to fall inside both windows; they flank the junction:

* **`orf_stop` panel — ND4L frame.** Window =
  `[ND4L.stop - 3 - 99, ND4L.stop - 3)` = `[192, 291)`. Captures the
  end of ND4L's reading frame, ending just before the 4-nt overlap.
* **`orf_start` panel — ND4 frame.** Window =
  `[ND4.start + 15, ND4.start + 15 + 99)` = `[305, 404)`. Captures
  the start of ND4's reading frame, beginning just after the overlap.

## Direct DFT at exact period 3.0

Standard `np.fft.rfft` evaluates at frequency bins `k / N` for
integer `k`, which means the period-3 bin lands at
`N / round(N / 3)` — never exactly 3 except when `N` is a multiple
of 3. We sidestep that by computing the DFT directly at arbitrary
periods:

```
z         = sum_n x[n] * exp(-2j * pi * n / period)
amplitude = |z| / sqrt(sum(x ** 2) + 1e-12)
```

Amplitude is in `[0, 1]`: 1.0 means `x` is a pure sinusoid at that
period; 0 means no projection onto that frequency. This is comparable
across windows of different length (unlike rfft amplitude).

The headline scalar:

```
spectral_ratio_3nt = amp(3.0) / median( amp(p) for p in 2..10
                                        excluding 2.8 <= p <= 3.2 )
```

`snr_call` tier:

| ratio | snr_call | meaning |
|---|---|---|
| ≥ 10× | `excellent` | Pristine library — every key length is well-phased. |
| ≥ 5×  | `healthy` | Publication-defensible periodicity. |
| ≥ 2×  | `modest` | Real but weak — suggests poor offset assignment or low-quality library. |
| < 2×  | `broken` | Likely frame-shifted or wrong offset. |
| NaN   | `no_signal` | Empty input — usually a depth issue. |

## Output schema

### `fourier_spectrum_combined.tsv`

One row per `(sample, read_length, gene_set, region, period_nt)`.

| Column | Type | Description |
|---|---|---|
| `sample` | string | Sample id |
| `read_length` | int | Footprint length (nt) |
| `gene_set` | `combined` \| `ATP86` \| `ND4L4` | See aggregation above |
| `region` | `orf_start` \| `orf_stop` | See window convention above |
| `n_genes` | int | Number of qualifying transcripts in the metagene |
| `n_sites_total` | int | Sum of site counts across qualifying transcripts |
| `n_nt` | int | Window width in nt (= `--periodicity-fourier-window-nt`) |
| `period_nt` | float | Period evaluated by direct DFT (grid: 2.0..10.0 step 0.05) |
| `amplitude` | float | Normalised amplitude in [0, 1] (see Direct DFT) |

### `fourier_period3_score_combined.tsv` (schema 1.1)

One row per `(sample, read_length, gene_set, region)`.

| Column | Type | Description |
|---|---|---|
| `sample`, `read_length`, `gene_set`, `region`, `n_genes`, `n_sites_total`, `n_nt` | (as above) | |
| `amp_at_3nt` | float | Direct DFT amplitude at exactly period 3.0 |
| `background_amp_median` | float | Median amplitude in [2..10] nt excluding (2.8, 3.2) |
| `spectral_ratio_3nt` | float | `amp_at_3nt / background_amp_median` |
| `local_background_amp_median` | float | Median amplitude in the codon-scale neighbourhood [4..6] nt |
| `spectral_ratio_3nt_local` | float | `amp_at_3nt / local_background_amp_median` — robust to long-period structure |
| `snr_call` | string | Tier (see table above) for the global ratio |
| `snr_call_local` | string | Tier for the local-background ratio |
| `amp_3nt_ci_low`, `amp_3nt_ci_high` | float | Bootstrap-over-genes percentile CI on `amp_at_3nt` (v0.9.0+) |
| `spectral_ratio_3nt_ci_low`, `spectral_ratio_3nt_ci_high` | float | Bootstrap CI on `spectral_ratio_3nt` (v0.9.0+) |
| `spectral_ratio_3nt_local_ci_low`, `spectral_ratio_3nt_local_ci_high` | float | Bootstrap CI on the local ratio (v0.9.0+) |
| `permutation_p` | float | Laplace-smoothed empirical p-value from a per-gene circular-shift null on the global ratio (v0.9.0+) |
| `permutation_p_local` | float | Same for the local ratio (v0.9.0+) |
| `n_bootstrap`, `n_permutations`, `ci_alpha` | int / float | Audit fields recording the exact procedure (defaults: 200 / 200 / 0.10) |
| `ci_method` | string | `percentile_over_genes` (live), `skipped_too_few_genes` (< 3 qualifying tracks), `disabled` (stats turned off via `--no-stats`) |
| `null_method` | string | `circular_shift_per_gene` (live) or `disabled` |
| `transcripts` | string | Semicolon-joined list of transcripts that contributed |

When fewer than `MIN_GENES_FOR_BOOTSTRAP_CI = 3` qualifying per-gene
tracks are available, the bootstrap is skipped (NaN CI columns +
`ci_method == "skipped_too_few_genes"`) rather than emitted at
misleadingly tight precision. The point estimate `spectral_ratio_3nt`
is still emitted; the user retains the ability to look at the
metagene figure and judge it visually.

## Plot bundle

Per `(sample, read_length)`, three figures land under
`<output>/rpf/qc/fourier_spectrum/<sample>/`:

* `<sample>_<length>nt_combined.{png,svg}` — single trace per panel
  drawn from the combined gene_set metagene.
* `<sample>_<length>nt_ATP86.{png,svg}` — junction-bracketed ATP86
  analysis (ATP8 frame top, ATP6 frame bottom).
* `<sample>_<length>nt_ND4L4.{png,svg}` — junction-bracketed ND4L4
  analysis (ND4L frame top, ND4 frame bottom).

Each PNG carries the canonical per-plot `.metadata.json` sidecar so
`mitoribopy validate-figures` can score it without re-running
matplotlib. The sidecar's `panel_layout` field records which
transcript drives each panel (auditability).

In-figure annotations on each panel show the `spectral_ratio_3nt`
value and the `snr_call` tier so a reviewer can read the headline
number off the plot.

## Statistical hardening (shipped in v0.9.0)

Two companion estimates are computed alongside `spectral_ratio_3nt`
and live in the score table (see the column list above):

* **Bootstrap CI over genes** — resample the per-gene normalised
  tracks with replacement (default 200 iterations), recompute the
  metagene + `amp_at_3nt` per resample, report the two-sided
  percentile interval (default 90 % CI). Surfaced as
  `amp_3nt_ci_{low,high}`, `spectral_ratio_3nt_ci_{low,high}`, and
  `spectral_ratio_3nt_local_ci_{low,high}`. Tells you "how sensitive
  is this peak to which genes happen to carry it."
* **Circular-shift permutation null** — independently shift each
  per-gene RAW coverage track by a random offset in `[0, W-1]`,
  re-apply the unit-mean → linear-detrend → Hann-window pipeline,
  re-aggregate, recompute the spectral ratio (default 200 iterations).
  Report the Laplace-smoothed empirical p-value
  `(k+1)/(n_perm+1)` as `permutation_p` and `permutation_p_local`.
  Tells you "is this peak unlikely under a same-magnitude signal
  with randomised codon-phase coherence across genes."

Both are deterministic given the seed (default 42, recorded in
`periodicity.metadata.json`). Neither changes the value of
`spectral_ratio_3nt`; they tell you how much to trust it.

Tune them per run via:

```
mitoribopy periodicity --bootstrap-n 1000 --permutations-n 1000 \
                       --ci-alpha 0.05 --random-seed 7
mitoribopy rpf --periodicity-fourier-bootstrap-n 1000 \
               --periodicity-fourier-permutations-n 1000 \
               --periodicity-fourier-ci-alpha 0.05 \
               --periodicity-fourier-random-seed 7
```

In `mitoribopy all`, the same knobs live under the YAML
`periodicity:` section as `fourier_bootstrap_n`,
`fourier_permutations_n`, `fourier_ci_alpha`, `fourier_random_seed`,
and `no_fourier_stats` (set the last to `true` to skip both for a
fast smoke run).

## Reading the figure

A healthy mt-Ribo-seq library produces:

1. A **sharp peak at period 3 nt in BOTH ORF panels** of the
   `combined` figure — `spectral_ratio_3nt >= 5`. Agreement between
   `orf_start` and `orf_stop` is the strongest evidence that
   elongation engages frame-3 reading promptly after initiation AND
   keeps phasing through to termination.
2. The **ATP86 panels** show period-3 peaks in EACH of the two
   frames (ATP8 frame and ATP6 frame), confirming the bicistronic
   transcript is translated in both reading frames.
3. The **ND4L4 panels** show similar two-frame translation evidence.

If the `orf_start` panel is flat but `orf_stop` is healthy, the
offset assignment for that read length is suspect — most often the
5'-end-to-P-site offset is correct for elongation footprints but
wrong for initiation footprints. Cross-check `metagene_start_p_site.svg`
for the same length.

## Method summary

In one paragraph: per-(sample, read_length, transcript) coverage is
sliced into a `W`-nt window anchored downstream of the start codon
(`orf_start`) or upstream of the stop codon (`orf_stop`), with the
first 5 codons after AUG and the last codon before stop dropped.
Each per-gene window is divided by its own mean (unit-mean
normalisation), mean-centred, and Hann-windowed. Qualifying tracks
are averaged element-wise into a metagene; a direct DFT is then
evaluated on the metagene at every period in `[2.0, 10.0]` (step
0.05) and at exactly period 3.0. The headline scalar is
`spectral_ratio_3nt = amp(3.0) / median(background)`, reported with
an `snr_call` tier. `W = 99 nt` (33 codons; multiple of 3 so the
period-3 bin lands on the integer) is the default; pass
`--periodicity-fourier-window-nt N` to override.
