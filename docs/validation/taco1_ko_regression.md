# TACO1-KO polyproline stalling regression — validation plan

`TACO1` encodes a mt-translational activator of `MT-CO1`. TACO1-KO libraries from prior publications (e.g. Weraarpachai et al. 2009; Richman et al. 2016; Soto et al. 2022) display a characteristic polyproline-stalling phenotype: ribosomes pile up at proline codons in `MT-CO1` due to impaired release-factor recycling downstream of the TACO1 activator lesion.

This file records the regression test MitoRiboPy must pass before tagging a release candidate.

---

## Acceptance criteria

1. **Primary signal.** In the TACO1-KO condition, P-site codon occupancy at proline codons (`CCU`, `CCC`, `CCA`, `CCG`) in `MT-CO1` exceeds the control condition by at least the ratio reported in the original publication (typically ≥ 1.5×).

2. **Specificity.** The polyproline signal is `MT-CO1`-specific; other mt-mRNAs do not show a matching proline-codon occupancy increase.

3. **Signal preservation under safe dedup defaults.** Running with `--dedup-strategy auto` (the default; resolves to `umi-tools` for UMI samples and `skip` for non-UMI samples) preserves the polyproline signal. Running the same library with `--dedup-strategy mark-duplicates` (gated behind the long confirmation flag) collapses the polyproline signal to near-baseline. This is the empirical proof of the warning that coordinate-only dedup destroys codon-occupancy signal on low-complexity mt-Ribo-seq libraries.

4. **Per-sample resolution stability** (v0.4+). When the TACO1-KO and control libraries use different kits (e.g. one TruSeq and one NEBNext UMI), per-sample kit detection (`--kit-preset auto`, default) and per-sample dedup resolution still preserve the per-sample polyproline signal — the criterion above must hold for both samples regardless of kit chemistry.

5. **Per-sample offset stability** (v0.4+). With `--offset_mode per_sample` (default), per-sample offset selection must produce 5'/3' offsets within ±1 nt of the combined-across-samples offsets at every read length where coverage is sufficient. The drift plot at `offset_diagnostics/plots/offset_drift_<align>.svg` (renamed from `plots_and_csv/` in v0.4.4) is the primary check; the per-sample audit at `offset_diagnostics/csv/per_sample_offset/<sample>/offset_applied.csv` confirms which row was actually applied downstream.

6. **Reference-consistency gate.** When the downstream `rnaseq` run uses a different reference FASTA than the rpf-side recorded checksum, `mitoribopy rnaseq` exits 2 with the `MISMATCH` error banner (no outputs written).

---

## Command sequence

```bash
# 1. End-to-end for control and TACO1-KO samples together (safe defaults).
$ mitoribopy all \
    --config validation/taco1_ko_pipeline.yaml \
    --output validation/taco1_ko_results

# 2. Alternate run with mark-duplicates to demonstrate signal loss.
$ mitoribopy all \
    --config validation/taco1_ko_mark_duplicates.yaml \
    --output validation/taco1_ko_mark_duplicates_results
```

The second config differs from the first only in `align.dedup_strategy: mark-duplicates` plus `align.confirm_mark_duplicates: true` — every other parameter is held constant so any signal difference is attributable to the dedup step.

---

## Inputs required

- Control + TACO1-KO FASTQs (typically ≥ 3 biological replicates each).
- Published FASTQ identifiers (e.g. from GEO / SRA) plus kit chemistry metadata so per-sample auto-detection resolves correctly. With `kit_preset: auto` (the v0.4+ default) the only metadata you actually need is "is this raw or pre-trimmed?" — pre-trimmed SRA inputs are auto-resolved to the `pretrimmed` kit.
- Human mt-transcriptome FASTA + bowtie2 index.
- DE table comparing control vs TACO1-KO mRNA levels (from paired RNA-seq or estimated from the RPF pre-offset counts; the rnaseq stage's SHA256 gate enforces reference consistency regardless).

---

## Status

**Pending dataset delivery.** The validation cannot execute in the MitoRiboPy development sandbox because the tutorial-scale FASTQs and the TACO1-KO dataset are not part of the repository. This document encodes the exact acceptance criteria so the test is reproducible once the data is in place on the maintainer's workstation.

When the run completes, append the observed proline-codon occupancy ratios (KO vs control; safe-dedup vs mark-duplicates) to this file and link any resulting figures from `docs/validation/figures/`. A release candidate should not be tagged until the primary and specificity criteria pass.
