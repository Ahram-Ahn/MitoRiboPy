# Release Notes

Versioned release summaries for MitoRiboPy. The authoritative line-by-line changelog lives in [`CHANGELOG.md`](../../CHANGELOG.md) at the repo root; the entries here are higher-level highlights aimed at users upgrading from one version to the next.

- [v0.5.0](v0.5.0.md): `mitoribopy rnaseq` is now from-FASTQ-by-default — cutadapt + bowtie2 + pyDESeq2 run inside the subcommand with auto SE/PE / adapter / UMI detection and auto-pseudo-replicate fallback for n=1 conditions; the original "bring your own DE table" path stays available as the alternative flow (`--de-table`); new `[fastq]` optional-dependency extra; four new TE / DE / QC plots (MA, TE bar by condition, TE heatmap, sample PCA)
- [v0.4.5](v0.4.5.md): `mark-duplicates` (picard) dedup option removed; `--rpf_min_count_frac` auto-filter for offset selection; compact `[COVERAGE]` log; templates moved to `examples/templates/`
- [v0.4.4](v0.4.4.md): unified per-site output layout (P-site / A-site filename prefixes replace `p/` / `a/` subdirs)
- [v0.4.2](v0.4.2.md): per-sample UMI / kit overrides, unified per-site output layout, strain rename (`h.sapiens` / `s.cerevisiae`), `short` footprint class, per-sample resume markers, intermediate-file cleanup, publication-quality codon-correlation plot, PyPI install, diagram redesign, working `mitoribopy align --config`
- [v0.4.1](v0.4.1.md): pre-trimmed FASTQ handling, kit-preset consolidation by adapter family, detection tuning flags
- [v0.4.0](v0.4.0.md): per-sample adapter detection + dedup, per-sample offsets, `--analysis_sites both`, polymorphic `align.fastq` YAML key
- [v0.3.0](v0.3.0.md): full FASTQ → translation-efficiency pipeline (`align`, `rnaseq`, `all` subcommands); strain `vm` / `ym` + `--footprint_class`; SHA256 reference-consistency gate
- [v0.2.0](v0.2.0.md): standalone package pipeline, packaged reference data, custom-organism support, and naming cleanup
- [v0.1.0](v0.1.0.md): first package-oriented public release with CI, tests, and release metadata
