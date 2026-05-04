# MitoRiboPy — architecture history

This document collects design decisions that previously lived inside
module docstrings as "Phase 3 / Phase 6 / v0.3.0 refactor" annotations.
The implementation in `src/mitoribopy/` describes current behaviour;
this file is the audit trail of how we got here.

For a full version-by-version log of user-visible changes see
[CHANGELOG.md](../../CHANGELOG.md). For per-release notes see
[release-notes/](../release-notes/).

---

## v0.3.0 refactor — phase outline

The v0.3.0 refactor split the package into a small set of focused
sub-packages. Internal phase numbers are no longer surfaced in
module docstrings; this section preserves the mapping.

| Internal phase    | Scope                                                                |
| ----------------- | -------------------------------------------------------------------- |
| Phase 0           | CLI / API surface audit                                              |
| Phase 1           | Package layout (`align/`, `analysis/`, `cli/`, `io/`, `pipeline/`)   |
| Phase 2           | `mitoribopy.config` + canonical YAML loader                          |
| Phase 3 (Path A)  | `mitoribopy align`: cutadapt → contam → mt-transcriptome bowtie2     |
| Phase 3.1         | Fail-loud tool checks (`align/tool_check.py`)                        |
| Phase 3.3         | Path A vs Path B decision (transcriptome reference, `--norc`/`--nofw`) |
| Phase 3.6         | pysam-only BAM I/O (no samtools/bedtools PATH dependency)            |
| Phase 4           | rpf BAM input support (`io/bam_reader.py`, BAM → BED6)               |
| Phase 5           | rnaseq integration (DE table + reference SHA gate, TE / ΔTE)         |
| Phase 6           | `mitoribopy all` orchestrator + manifest                             |
| Phase 7           | Diagnostics (offsets, coverage-frame, periodicity)                   |

## Path A vs Path B (mt-transcriptome vs full mtDNA)

We chose **Path A** — align to one FASTA record per mt-mRNA, not to the
full mtDNA contig — for two reasons:

1. *ND5 / ND6 antisense overlap* on the L-strand of human mtDNA makes
   genome-level strand assignment ambiguous. With one FASTA record per
   mt-mRNA, every read is unambiguously on the sense strand of *that*
   transcript and `--norc`/`--nofw` enforces it at alignment time.
2. *Annotation symmetry.* Yeast mt-introns drop out cleanly when the
   reference is the spliced mt-mRNA set. A genomic reference would need
   per-strain custom annotation maps.

The cost is one extra reference-prep step (build the transcriptome
FASTA) which `mitoribopy validate-reference` audits.

## P-site / A-site offset selection (rpf)

`mitoribopy.analysis.offset_selection` picks per-length, per-sample
offsets in canonical P-site space (not "reported-site" space) and
records the picked vector under `offset_applied.csv`. The legacy
"selected_site" pick reference is rewritten to "reported_site" by
`migrate-config`.

## Refactor-4 (publication-readiness, v0.6.0)

Major user-facing changes that landed in v0.6.0:

* `mitoribopy run` / `mitoribopy <bare flags>` removed.
* `mitoribopy rpf --help` reports as `mitoribopy rpf`.
* `mitoribopy all --strict` is the single publication-mode switch.
* Unified top-level `samples:` sample sheet.
* `from_fastq` / `de_table` / `rna_only` / `none` modes for rnaseq.
* Pseudo-replicate fallback is opt-in only via the
  `--allow-pseudo-replicates-for-demo-not-publication` flag.

## Fourth-edit cleanup (v0.6.2)

* Top-level `execution:` block; resource plan written to
  `<run_root>/resource_plan.json`.
* Canonical `dedup_strategy: umi_coordinate`; legacy `umi-tools` is
  rewritten by `migrate-config`.
* Historical periodicity QC bundle extended with `qc_summary.tsv`,
  `qc_summary.md`, gene-level frame fractions, optional phase score
  and FFT period-3 power. The current package has since replaced this
  with the metagene Fourier bundle documented in
  [`docs/reference/periodicity.md`](../reference/periodicity.md).
* Stable warning / error code registry under
  `docs/reference/warning_codes.md`.
* Release checklist under `docs/developer/release_checklist.md` is
  now a hard gate for every PyPI cut.

### Refinements layered on the v0.6.2 periodicity bundle

The bundle introduced in v0.6.2 received four targeted refinements
during the cycle (no rewrites; no architecture changes):

1. `build_qc_summary` `best_read_length_dominant_fraction` is now
   a direct `max(f0, f1, f2)` from the chosen row; previously a
   tortured `A and B or C` ternary collapsed to the wrong value when
   the dominant frame differed from the expected frame.
2. `build_qc_summary` best-length pick is depth-aware: the candidate
   pool is restricted to rows clearing `min_reads_per_length` whenever
   any qualify, so a 5-read length with lucky frame-0 dominance can no
   longer crown a 5,000-read length at 0.75.
3. `exclude_start_codons` / `exclude_stop_codons` became first-class
   parameters of `compute_frame_summary`,
   `compute_frame_summary_by_length`, `build_gene_periodicity`,
   `run_periodicity_qc`, and `run_periodicity_qc_bundle`. Defaults
   stay at 0 in the pipeline path to preserve historical numbers; the
   standalone CLI applies the spec defaults of 6 / 3 uniformly.
4. `gene_periodicity.tsv` shipped an `is_overlap_pair` column in that
   historical bundle, flagging known human mt-mRNA overlap regions.
   The current Fourier path handles ATP8/ATP6 and ND4L/ND4 through the
   dedicated `ATP86` and `ND4L4` gene sets instead.

### Per-frame split coverage plot

The `_plot_frame_colored` overlay was complemented with a
`_plot_frame_split` companion that stacks 3 sub-rows per sample
(frame 0, +1, +2) sharing the y-axis. The sibling
`{p_site,a_site}_density_{rpm,raw}_frame_split/` directories ship the
same `coverage_plot.metadata.json` sidecar as the overlay, so a
reviewer reads one definition of "frame 0" applied to both.
