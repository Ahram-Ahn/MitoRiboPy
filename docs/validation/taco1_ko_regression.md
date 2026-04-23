# TACO1-KO polyproline stalling regression &mdash; validation plan

`TACO1` encodes a mt-translational activator of `MT-CO1`. TACO1-KO
libraries from prior publications (e.g. Weraarpachai et al. 2009;
Richman et al. 2016; Soto et al. 2022) display a characteristic
polyproline-stalling phenotype: ribosomes pile up at proline codons
in `MT-CO1` due to impaired release-factor recycling downstream of the
TACO1 activator lesion.

This file records the Phase 7.7 regression test MitoRiboPy v0.3.0
must pass before tagging a release candidate.

## Acceptance criteria

1. **Primary signal.** In the TACO1-KO condition, P-site codon
   occupancy at proline codons (`CCU`, `CCC`, `CCA`, `CCG`) in
   `MT-CO1` exceeds the control condition by at least the ratio
   reported in the original publication (typically >= 1.5x).
2. **Specificity.** The polyproline signal is MT-CO1-specific; other
   mt-mRNAs do not show a matching proline-codon occupancy increase.
3. **Signal preservation under dedup.** Running the same library
   with `--dedup-strategy mark-duplicates` collapses the polyproline
   signal to near-baseline levels (this is the empirical proof of the
   Phase-3.5 warning that coordinate-only dedup destroys
   codon-occupancy signal). The safe defaults (`skip` for no-UMI,
   `umi-tools` for UMI) preserve the signal.
4. **Reference-consistency gate.** When the downstream rnaseq run
   uses a different reference FASTA than the rpf-side recorded
   checksum, `mitoribopy rnaseq` exits 2 with the `MISMATCH` error
   banner (no outputs written).

## Command sequence

```bash
# 1. End-to-end for control and TACO1-KO samples together.
$ mitoribopy all \
    --config validation/taco1_ko_pipeline.yaml \
    --output validation/taco1_ko_results

# 2. Alternate run with mark-duplicates to demonstrate signal loss.
$ mitoribopy all \
    --config validation/taco1_ko_mark_duplicates.yaml \
    --output validation/taco1_ko_mark_duplicates_results
```

The second config differs from the first only in its
`align.dedup_strategy` setting (`mark-duplicates` plus the
confirmation flag) &mdash; every other parameter is held constant so
any signal difference is attributable to the dedup step.

## Inputs required

- Control + TACO1-KO FASTQs (typically >= 3 biological replicates each).
- Published FASTQ identifiers (e.g. from GEO / SRA) along with kit
  chemistry metadata so `--kit-preset` resolves to the correct
  adapter + UMI settings.
- Human mt-transcriptome FASTA + bowtie2 index.
- DE table comparing control vs TACO1-KO mRNA levels (from paired
  RNA-seq or estimated from the RPF pre-offset counts; the rnaseq
  stage's SHA256 gate enforces reference consistency regardless).

## Status

**Pending dataset delivery.** The validation cannot execute in the
MitoRiboPy development sandbox because the tutorial-scale FASTQs and
the TACO1-KO dataset are not part of the repository. This document
encodes the exact acceptance criteria so the test is reproducible once
the data is in place on the maintainer's workstation.

When the run completes, append the observed proline-codon occupancy
ratios (KO vs control; safe-dedup vs mark-duplicates) to this file and
link any resulting figures from `docs/validation/figures/`. The v0.3.0
release candidate should not be tagged until the primary and
specificity criteria pass.
