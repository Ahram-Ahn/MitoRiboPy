# Synthetic-mini known-answer integration test

A purely synthetic fixture used by [`tests/test_synthetic_mini.py`](../../tests/test_synthetic_mini.py) to verify that the four core analytical stages of MitoRiboPy agree on a known answer. The fixture is generated programmatically inside the test (no binary blobs in the repo); this document describes how it is constructed and what it is designed to assert.

## Fixture design

| Property | Value |
|---|---|
| Transcripts | `ND1`, `ND2`, `ND3` from the bundled `human_annotation_df` |
| Conditions | `WT` (samples `WT_R1`, `WT_R2`) and `KO` (samples `KO_R1`, `KO_R2`) |
| Read length | 30 nt |
| Reads per WT sample, per transcript | 100 |
| Reads per KO sample, per transcript | 100 (`ND1`, `ND3`) and **200** on `ND2` |
| Read placement | 60% anchored at the stop codon; 40% tiled across the CDS in 3-nt steps |
| Strand | All reads on `+` (matches the sense-oriented mt-transcriptome FASTA) |

The 60/40 split is deliberate. The anchor reads provide a sharp peak at one offset value so the offset-selection routine has an unambiguous winner; the tile reads provide enough density across the CDS that the frame-summary computation is meaningful instead of degenerate.

## Read-start arithmetic

The package's stop-anchored 5'-offset arithmetic is

```
5' offset = stop_codon - read_start - 2
```

so to make the package report a 5' offset of **12**, anchor reads are placed at `read_start = stop_codon - 14`. Tile reads are placed at `read_start = codon_start - 12` so that, when the selected offset (12) is applied, every tile read's P-site lands on a codon start (frame 0 relative to that transcript).

## What the test asserts

| Stage | Assertion | Why this is the right invariant |
|---|---|---|
| Offset enrichment + selection | `Most Enriched 5' Offset == 12` for the combined run AND every per-sample run; `confidence_5 == "high"` for the combined fit | Confirms the offset selector picks the correct value out of the dominant peak and assigns the documented `high` label when the peak is clean. |
| Frame summary | `frame_1 == 0` for every sample; `frame_0 + frame_2 == 1.0` | All reads land at codon-aligned positions or at `stop - 2`, so the off-by-one phase (frame 1) must be empty. Frame 0 vs frame 2 mass depends on each transcript's CDS-length parity — the off-by-one being empty is the strongest negative control. |
| Frame summary, n_reads | WT samples report 300; KO samples report 400 | KO doubles `ND2` from 100 to 200 reads per sample. |
| Start-aligned metagene | Frame-1 mass == 0 across the whole 60-nt window | Same invariant as above, evaluated on the metagene profile. |
| TE math | `TE(WT, ND1) == (100 + 0.5) / (100 + 0.5) == 1.0`; `TE(KO, ND2) == (200 + 0.5) / (100 + 0.5)` | Direct check of the documented `TE = (RPF + δ) / (mRNA + δ)` equation with δ = 0.5. |
| ΔTE math | `delta_te_log2(ND2) == log2(200.5/100.5)`; `delta_te_log2(ND1) == 0`; `delta_te_log2(ND3) == 0` | KO doubles `ND2` RPF; the synthetic DE table holds mRNA log2FC at 0; ΔTE for ND2 is therefore the log2 of the RPF ratio; the held-equal transcripts must be zero. |
| Output artefacts | `frame_summary.tsv` exists and lists every sample | Catches accidental output-path regressions. |

## What the test does NOT cover

- The real `align` stage (cutadapt + bowtie2). The fixture starts from BED-equivalent in-memory dataframes; trim/align is verified separately in [`tests/test_align_*.py`](../../tests/).
- The from-FASTQ rnaseq flow (pyDESeq2). The DE table is hand-built; the from-FASTQ flow is exercised by [`tests/test_rnaseq_cli_from_fastq.py`](../../tests/test_rnaseq_cli_from_fastq.py).
- Plotting. The metagene Fourier path is invoked with `plot=False` so the synthetic test stays headless and fast; the plotting code path is exercised by the dedicated suites in [`tests/test_periodicity_cli.py`](../../tests/test_periodicity_cli.py) and the figure-validation tests in [`tests/test_validate_figures_cli.py`](../../tests/test_validate_figures_cli.py).
- Periodicity-bundle statistics (bootstrap CI over genes, circular-shift permutation null, `spectral_ratio_3nt_local`, `snr_call` tiering, dominant-fraction bug fix). Those have dedicated unit tests in [`tests/test_fourier_spectrum.py`](../../tests/test_fourier_spectrum.py), [`tests/test_fourier_stats.py`](../../tests/test_fourier_stats.py), and [`tests/test_fourier_edge_cases.py`](../../tests/test_fourier_edge_cases.py). The legacy `test_periodicity_qc.py` / `test_periodicity_refinements.py` modules were retired in v0.8.0 along with the frame-fraction QC bundle.

## Why this fixture is in the repo

The architectural review asked for a tiny known-answer integration test that lets a future reviewer (or CI in a year) prove that the four numerical claims in the README's "What the numbers mean" section still hold. This fixture is the smallest dataset that exercises offset selection + frame QC + metagene + TE + ΔTE in one pass on real (bundled) annotation, with hand-computable answers at every stage.
