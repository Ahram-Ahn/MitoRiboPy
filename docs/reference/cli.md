# CLI reference

> **GENERATED FILE — do not hand-edit.** Regenerate with
> `python docs/generate_cli_reference.py`. CI runs the same script in
> `--check` mode to fail when this file drifts from the live argparse
> parsers in `src/mitoribopy/`.

This document is the canonical machine reference for every flag in
every `mitoribopy` subcommand. For prose, examples, and decision
trees see [the README](../../README.md) and the tutorials under
[`docs/tutorials/`](../tutorials/). For the column-by-column output
catalogue (every TSV, CSV, and JSON the package writes, with units
and coordinate spaces) see
[`docs/reference/output_schema.md`](output_schema.md).

Generated against MitoRiboPy v0.7.1.

## Subcommand summary

| Subcommand | What it does |
|---|---|
| [`mitoribopy align`](#mitoribopy-align) | Preprocess FASTQ inputs: cutadapt trim + bowtie2 contaminant subtraction + bowtie2 mt-transcriptome alignment + MAPQ filter + dedup + BAM->BED6. Produces drop-in inputs for 'mitoribopy rpf'. |
| [`mitoribopy rpf`](#mitoribopy-rpf) | Ribo-seq analysis from BED/BAM inputs. |
| [`mitoribopy rnaseq`](#mitoribopy-rnaseq) | Translation efficiency (TE / delta-TE) from paired RNA-seq + Ribo-seq. Default flow: pass --rna-fastq + --ribo-fastq + --reference-fasta and the subcommand runs trimming, bowtie2 alignment, per-transcript counting, and pyDESeq2 itself before emitting te.tsv, delta_te.tsv, and plots. Alternative: pass --de-table from a prior external DESeq2 / Xtail / Anota2Seq run together with --ribo-dir; this path is mutually exclusive with --rna-fastq and enforces a SHA256 reference-consistency gate. |
| [`mitoribopy all`](#mitoribopy-all) | End-to-end orchestrator: align + rpf, plus rnaseq when the config carries an 'rnaseq' section configured for either flow (from-FASTQ via 'rna_fastq' + 'reference_fasta', or external-DE via 'de_table'). Writes a composed run_manifest.json with tool versions, parameters, and input/output hashes across all three stages. |
| [`mitoribopy periodicity`](#mitoribopy-periodicity) | Quantify 3-nt periodicity by running the metagene Fourier analysis on a pre-assigned site table. |
| [`mitoribopy migrate-config`](#mitoribopy-migrateconfig) | Rewrite legacy MitoRiboPy YAML keys to their canonical names. Input is read from a path; output is written to stdout (the change log goes to stderr). Use to upgrade old pipeline configs without manually hunting down every renamed key. |
| [`mitoribopy validate-config`](#mitoribopy-validateconfig) | Pre-flight a MitoRiboPy YAML / JSON / TOML config: parse, canonicalise legacy keys, check file paths and mutually-exclusive sections, and resolve rnaseq.mode against supplied inputs. Exit code is 0 on success, 2 when at least one error was found. |
| [`mitoribopy validate-reference`](#mitoribopy-validatereference) | Pre-flight a custom mitochondrial reference: check that the FASTA and annotation CSV are consistent (matching transcript IDs, matching lengths, CDS divisible by 3, valid start / stop codons under the selected codon table). |
| [`mitoribopy validate-figures`](#mitoribopy-validatefigures) | Mechanically validate every plot under a finished MitoRiboPy run root: check label / legend / stat-box overlap, label clipping, point counts vs source TSV, SVG text editability, PNG dpi, and metadata sidecar coverage. Writes <RUN_DIR>/figure_qc.tsv. Exit 0 / 1 / 2 (all pass / warn-only / fail; --strict upgrades warn → fail). |
| [`mitoribopy summarize`](#mitoribopy-summarize) | Regenerate SUMMARY.md and summary_qc.tsv from a finished MitoRiboPy run by reading the run_manifest.json and per-stage TSV outputs. Useful for re-rendering summaries on archival runs without re-executing any pipeline stage. |
| [`mitoribopy benchmark`](#mitoribopy-benchmark) | Time and disk-measure a full `mitoribopy all` run, optionally after pre-subsampling each FASTQ to N reads. Produces benchmark.tsv and benchmark_summary.md at the run root for tuning thread counts, disk budgets, and per-stage wall time. |

---

## `mitoribopy align`

```text
usage: mitoribopy align [-h] [--config CONFIG] [--dry-run] [--threads N]
                        [--log-level {DEBUG,INFO,WARNING,ERROR}]
                        [--fastq-dir DIR] [--fastq PATH]
                        [--contam-index BT2_PREFIX] [--mt-index BT2_PREFIX]
                        [--output DIR] [--adapter SEQ] [--pretrimmed]
                        [--umi-length N] [--umi-position {5p,3p,both}]
                        [--umi-length-5p N] [--umi-length-3p N]
                        [--sample-overrides TSV] [--keep-intermediates]
                        [--tmpdir TMPDIR] [--allow-count-invariant-warning]
                        [--strict-publication-mode] [--resume]
                        [--adapter-detection MODE] [--adapter-detect-reads N]
                        [--adapter-detect-min-rate FRAC]
                        [--adapter-detect-min-len N]
                        [--adapter-detect-pretrimmed-threshold FRAC]
                        [--no-pretrimmed-inference]
                        [--library-strandedness {forward,reverse,unstranded}]
                        [--min-length NT] [--max-length NT] [--quality Q]
                        [--mapq Q] [--seed N]
                        [--dedup-strategy {auto,umi_coordinate,umi-tools,umi_tools,skip}]
                        [--umi-dedup-method {unique,percentile,cluster,adjacency,directional}]
                        [--max-parallel-samples N|auto] [--single-sample-mode]
                        [--memory-gb GB|auto]

Preprocess FASTQ inputs: cutadapt trim + bowtie2 contaminant subtraction + bowtie2 mt-transcriptome alignment + MAPQ filter + dedup + BAM->BED6. Produces drop-in inputs for 'mitoribopy rpf'.

options:
  -h, --help                          show this help message and exit

Shared options:
  --config CONFIG                     Configuration file (.json, .yaml, .yml, or .toml). CLI arguments override values read from the file.
  --dry-run                           Print planned actions and exit without executing.
  --threads N                         Preferred thread count for external tools and BLAS libraries.
  --log-level {DEBUG,INFO,WARNING,ERROR}
                                      Python logging level for MitoRiboPy console output. [default: INFO]

Inputs:
  --fastq-dir DIR                     Directory containing input FASTQ files (*.fq, *.fq.gz, *.fastq, *.fastq.gz). Pass --fastq-dir OR --fastq (or both).
  --fastq PATH                        Individual FASTQ input; repeatable.
  --contam-index BT2_PREFIX           bowtie2 index prefix of the contaminant panel (rRNA, tRNA, any nuclear spike-ins to subtract). Required for non-dry-run invocations. Build with 'bowtie2-build contaminants.fa <prefix>'.
  --mt-index BT2_PREFIX               bowtie2 index prefix of the mt-transcriptome (one FASTA record per mt-mRNA, header matching annotation sequence_name). Required for non-dry-run invocations.
  --output DIR                        Output directory for BAM/BED/read_counts (required).

Library prep:
  --adapter SEQ                       3' adapter sequence. By default the pipeline auto-detects the adapter from the head of each FASTQ; pass --adapter <SEQ> when detection cannot identify your library or when you want to pin a specific sequence. Mutually exclusive with --pretrimmed.
  --pretrimmed                        Declare that the input FASTQ has already been adapter-trimmed (e.g. SRA-deposited data). cutadapt skips the -a flag and only enforces length and quality filtering. Auto-detection also infers this when no known adapter signature is present; pass this flag to assert it explicitly. Mutually exclusive with --adapter.
  --umi-length N                      UMI length in nt. Overrides the kit preset's default. For --umi-position=both this MUST equal --umi-length-5p + --umi-length-3p (it is the canonical concatenated QNAME UMI length umi_tools dedups on).
  --umi-position {5p,3p,both}         UMI position within the insert (overrides kit preset). '5p' / '3p' are single-end UMIs; 'both' is a dual-end UMI library (e.g. xGen Duplex, Twist) — supply --umi-length-5p and --umi-length-3p in that mode.
  --umi-length-5p N                   Per-end 5' UMI length in nt. Used only when --umi-position=both; ignored otherwise.
  --umi-length-3p N                   Per-end 3' UMI length in nt. Used only when --umi-position=both; ignored otherwise.
  --sample-overrides TSV              Path to a TSV with per-sample overrides for adapter / pretrimmed / umi_length / umi_position / dedup_strategy. Required header columns: 'sample' plus at least one of the override columns. The 'sample' value must match the FASTQ basename with the .fq[.gz] / .fastq[.gz] suffix removed. Empty cells (or NA / None / null) fall through to the global CLI default for that field, so a single sample can override only its UMI without restating the rest. Useful for mixed-UMI / mixed-adapter batches.
  --keep-intermediates                Keep the per-step intermediate files (trimmed FASTQ, contam-filtered FASTQ, pre-MAPQ BAM). By default these are deleted as soon as the next step has consumed them, since they are large, regenerable, and not needed by any downstream stage. Pass this flag when debugging a sample or comparing per-step intermediate counts.
  --tmpdir TMPDIR                     Optional override for the directory used for per-step scratch files (trimmed FASTQ, contam-filtered FASTQ, intermediate BAMs). Defaults to a subdirectory of --output. Set this to a fast local SSD when running on a cluster with slow shared storage, or to a pre-mounted tmpfs to avoid hitting disk altogether for short runs.
  --allow-count-invariant-warning     DEVELOPER / DEBUG ONLY. Demote read_counts.tsv invariant violations from errors to warnings. The default is to fail the run on any violation, since a real violation indicates a bug somewhere upstream. NEVER use for a publication run (`mitoribopy all --strict` will reject this flag too).
  --strict-publication-mode           Reject runs that rely on inferred-rather-than-declared metadata: inferred-pretrimmed kits and ambiguous adapter detection (confidence margin < 0.10). Use when preparing a publication run so the strict checks fail loud rather than a single sample silently picking the wrong defaults.
  --resume                            Skip samples that have already completed in a previous invocation against this --output directory. Each completed sample writes a small JSON file under <output>/.sample_done/; on resume, samples whose JSON is present and parses are reloaded instead of re-run. Use this after a crash or kill mid-batch to avoid redoing the samples that already finished. The orchestrator ('mitoribopy all --resume') sets this automatically when the align stage's read_counts.tsv is missing.
  --adapter-detection MODE            Per-sample adapter detection policy. 'auto' (default) scans every input FASTQ and picks the matching adapter family per sample; samples whose scan fails fall back to --adapter when supplied, or to 'pretrimmed' when no fallback is set and the data looks already-trimmed. 'strict' scans and HARD-FAILS on any sample whose detected adapter conflicts with an explicit --adapter or where no adapter can be identified. 'off' skips the scan and trusts --adapter / --pretrimmed for every sample (one of those is required). [default: auto]
  --adapter-detect-reads N            Number of FASTQ reads to scan per sample during adapter auto-detection. Increase for noisy libraries where the first 5000 reads have unusual adapter distributions; decrease for a faster pre-flight pass on cleaner data. [default: 5000]
  --adapter-detect-min-rate FRAC      Minimum fraction of scanned reads that must contain an adapter prefix for the kit to be considered detected. Lower for sparsely-adapted libraries (e.g. 0.10); raise for stricter calls. [default: 0.3]
  --adapter-detect-min-len N          Adapter prefix length used as the search needle (nt). Default 12. Lower (e.g. 8) tolerates noisy adapter regions; raise (e.g. 16) for stricter matches. [default: 12]
  --adapter-detect-pretrimmed-threshold FRAC
                                      When EVERY kit's match rate is at or below this value, the FASTQ is classified as already adapter-trimmed and resolved to the 'pretrimmed' kit (cutadapt skips the -a flag). Default 0.05 (5%). [default: 0.05]
  --no-pretrimmed-inference           Disable the auto-fallback to 'pretrimmed' when adapter detection finds no known kit. With this flag, detection failure with no --adapter / --pretrimmed fallback raises an error instead. [default: True]
  --library-strandedness {forward,reverse,unstranded}
                                      Library strandedness. 'forward' (default, NEBNext/TruSeq small-RNA kits) enforces --norc at alignment time; 'reverse' enforces --nofw; 'unstranded' leaves bowtie2 permissive. [default: forward]
  --min-length NT                     Minimum read length kept after trimming (mt-RPF default 15). [default: 15]
  --max-length NT                     Maximum read length kept after trimming (mt-RPF default 45). [default: 45]
  --quality Q                         cutadapt -q Phred+33 3' quality trim threshold. [default: 20]

Alignment:
  --mapq Q                            MAPQ threshold for the post-alignment filter (default 10). The main reason for this filter is NUMT cross-talk suppression. [default: 10]
  --seed N                            bowtie2 --seed value (deterministic output). [default: 42]

Deduplication:
  --dedup-strategy {auto,umi_coordinate,umi-tools,umi_tools,skip}
                                      Canonical: auto | umi_coordinate | skip. 'auto' (default) -> umi_coordinate when UMIs are present, else skip. umi_coordinate collapses reads on (coordinate, UMI); the implementation calls into umi_tools but the statistical operation is coordinate+UMI dedup. The legacy aliases 'umi-tools' / 'umi_tools' are accepted and rewritten to 'umi_coordinate' (canonical_config.yaml records the canonical name). Coordinate-only mark-duplicates is not supported: it destroys codon-occupancy signal on low-complexity mt-Ribo-seq libraries. [default: auto]
  --umi-dedup-method {unique,percentile,cluster,adjacency,directional}
                                      umi_tools --method. 'unique' (default) collapses only on exact coord+UMI match; other methods may over-collapse in low-complexity mt regions. [default: unique]

Execution:
  --max-parallel-samples N|auto       Number of samples to process concurrently. Accepts an integer or 'auto' (default). 'auto' picks min(n_samples, threads/min_per_sample, memory_gb/est_per_sample) so a modern CPU is not left idle on a multi-sample run. With --threads T, each worker uses max(1, T // N) tool threads so total CPU use stays around T. Pass --max-parallel-samples 1 (or --single-sample-mode) for legacy serial behaviour. Per-sample work (cutadapt + bowtie2 + dedup + BAM->BED) is embarrassingly parallel; the joint 'mitoribopy rpf' stage is unaffected.
  --single-sample-mode                Force serial execution (alias for --max-parallel-samples 1). Use when you want one sample's logs interleaved cleanly or when memory pressure rules out concurrency.
  --memory-gb GB|auto                 Total memory budget (in GiB) the auto scheduler may use. Accepts a float or 'auto' (no memory cap). Used only when --max-parallel-samples auto.
```

## `mitoribopy rpf`

```text
usage: mitoribopy rpf [-h] [--config CONFIG] -f REF_FASTA [-s STRAIN]
                      [-d BED_DIR] [--bam-mapq Q] [-rpf MIN_LEN MAX_LEN]
                      [--footprint-class {short,monosome,disome,custom}]
                      [--annotation-file ANNOTATION.csv]
                      [--codon-tables-file CODON_TABLES.json]
                      [--codon-table-name TABLE_NAME]
                      [--start-codons CODON [CODON ...]]
                      [--atp8-atp6-baseline {ATP6,ATP8}]
                      [--nd4l-nd4-baseline {ND4,ND4L}] [-a {start,stop}]
                      [-r NT] [--min-offset NT] [--max-offset NT]
                      [--rpf-min-count-frac FRAC] [--min-5-offset NT]
                      [--max-5-offset NT] [--min-3-offset NT]
                      [--max-3-offset NT] [--offset-mask-nt NT]
                      [--offset-pick-reference {reported_site,p_site,selected_site}]
                      [--offset-type {5,3}] [--offset-site {p,a}]
                      [--analysis-sites {p,a,both}]
                      [--codon-overlap-mode {full,any}] [-p NT]
                      [--offset-mode {per_sample,combined}] [-o OUTPUT_DIR]
                      [--downstream-dir DOWNSTREAM_DIR] [--plot-dir PLOT_DIR]
                      [-fmt {png,pdf,svg}] [--x-breaks NT [NT ...]]
                      [--line-plot-style {combined,separate}]
                      [--cap-percentile CAP_PERCENTILE] [-m]
                      [--order-samples ORDER_SAMPLES [ORDER_SAMPLES ...]]
                      [--read-counts-file COUNTS_TABLE]
                      [--read-counts-sample-col READ_COUNTS_SAMPLE_COL]
                      [--read-counts-reads-col READ_COUNTS_READS_COL]
                      [--unfiltered-read-length-range MIN_LEN MAX_LEN]
                      [--rpm-norm-mode {total,mt_mrna}]
                      [--read-counts-reference-col READ_COUNTS_REFERENCE_COL]
                      [--mt-mrna-substring-patterns PATTERN [PATTERN ...]]
                      [--structure-density]
                      [--structure-density-norm-perc STRUCTURE_DENSITY_NORM_PERC]
                      [--cor-plot] [--base-sample BASE_SAMPLE]
                      [--cor-mask-method {percentile,fixed,none}]
                      [--cor-mask-percentile COR_MASK_PERCENTILE]
                      [--cor-mask-threshold COR_MASK_THRESHOLD]
                      [--cor-metric {log2_density_rpm,log2_rpm,linear,raw_count}]
                      [--cor-regression {theil_sen,ols,none}]
                      [--cor-support-min-raw COR_SUPPORT_MIN_RAW]
                      [--cor-label-top-n COR_LABEL_TOP_N]
                      [--cor-pseudocount COR_PSEUDOCOUNT]
                      [--cor-raw-panel {qc_only,off}]
                      [--read-coverage-raw | --no-read-coverage-raw | --read_coverage_raw | --no-read_coverage_raw]
                      [--read-coverage-rpm | --no-read-coverage-rpm | --read_coverage_rpm | --no-read_coverage_rpm]
                      [--igv-export | --no-igv-export | --igv_export | --no-igv_export]
                      [--periodicity-enabled | --no-periodicity-enabled | --periodicity_enabled | --no-periodicity_enabled]
                      [--periodicity-fourier-window-nt PERIODICITY_FOURIER_WINDOW_NT]
                      [--periodicity-metagene-nt PERIODICITY_METAGENE_NT]
                      [--periodicity-metagene-normalize {per_gene_unit_mean,none}]
                      [--periodicity-fourier-bootstrap-n PERIODICITY_FOURIER_BOOTSTRAP_N]
                      [--periodicity-fourier-permutations-n PERIODICITY_FOURIER_PERMUTATIONS_N]
                      [--periodicity-fourier-ci-alpha PERIODICITY_FOURIER_CI_ALPHA]
                      [--periodicity-fourier-random-seed PERIODICITY_FOURIER_RANDOM_SEED]
                      [--periodicity-no-fourier-stats]

Run the MitoRiboPy Ribo-seq analysis stage on BED / BAM inputs.
This subcommand filters reads, estimates P-site / A-site offsets,
and then generates translation-profile, codon, coverage-profile,
and optional RNA-seq / structure-density outputs.

options:
  -h, --help                          show this help message and exit

Core Inputs:
  --config CONFIG                     Optional config file (.yaml, .yml, .json, or .toml; format auto-detected by extension). CLI arguments override values from this file.
  -f REF_FASTA, --fasta REF_FASTA     Reference FASTA file used to build the transcript/annotation context.
  -s STRAIN, --strain STRAIN          Reference preset (organism). Built-in strains ship a complete
                                      annotation table and codon table, so the user only needs to
                                      supply the FASTA reference. Other organisms must use
                                      'custom' and supply their own annotation + codon-table
                                      choice (see --annotation_file / --codon_table_name).
                                        h.sapiens     human mt (default; vertebrate_mitochondrial
                                                      codon table; ATG/ATA start codons)
                                        s.cerevisiae  budding yeast mt (yeast_mitochondrial codon
                                                      table; ATG start codon)
                                        custom        any other organism (mouse, fly, A. thaliana,
                                                      fungi, ...); requires --annotation_file plus
                                                      a --codon_table_name picked from the
                                                      built-in NCBI Genetic Codes list (see
                                                      --codon_table_name help) or a
                                                      --codon_tables_file you supply.
                                        h, y          deprecated short aliases for h.sapiens and
                                                      s.cerevisiae; still accepted for one cycle. [default: h.sapiens]
  -d BED_DIR, --directory BED_DIR     Directory containing Ribo-seq input files.
                                      Both .bed and .bam are accepted; BAM files are auto-converted
                                      to BED6 under <output>/bam_converted/ via pysam. [default: current working directory]
  --bam-mapq Q, --bam_mapq Q          MAPQ threshold applied to BAM inputs before BAM->BED6 conversion.
                                      Set to 0 to disable. Default 10 is the same NUMT-suppression
                                      default used by 'mitoribopy align'. [default: 10]
  -rpf MIN_LEN MAX_LEN                Inclusive read-length filter range, for example: -rpf 29 34 [default: short    h.sapiens / s.cerevisiae: 16-24; monosome h.sapiens: 28-34, s.cerevisiae: 37-41; disome   h.sapiens: 50-70, s.cerevisiae: 60-90]
  --footprint-class {short,monosome,disome,custom}, --footprint_class {short,monosome,disome,custom}
                                      Expected ribosome-protected fragment class. Selects the
                                      default --rpf and --unfiltered_read_length_range windows for
                                      the chosen --strain when the user does not pass them
                                      explicitly. Built-in defaults exist for h.sapiens and
                                      s.cerevisiae; for --strain custom you must also pass -rpf.
                                        short     truncated RNase products (~16-24 nt). Sit just
                                                  below the canonical monosome window; useful for
                                                  mapping context-dependent pausing and as a QC
                                                  indicator of digest aggressiveness.
                                        monosome  single-ribosome footprint (default).
                                                  h.sapiens 28-34 nt, s.cerevisiae 37-41 nt.
                                        disome    collided-ribosome footprint. h.sapiens 50-70 nt,
                                                  s.cerevisiae 60-90 nt. Used for ribosome-stalling /
                                                  queueing studies (e.g. eIF5A-depletion).
                                        custom    no biological default; supply -rpf and
                                                  --unfiltered_read_length_range yourself. [default: monosome]
  --annotation-file ANNOTATION.csv, --annotation_file ANNOTATION.csv
                                      Per-transcript annotation CSV. Required when --strain
                                      custom; the built-in human (h.sapiens) and yeast
                                      (s.cerevisiae) tables are used otherwise. The CSV must have
                                      one row per transcript with the columns: transcript,
                                      sequence_name, display_name, sequence_aliases, l_utr5,
                                      l_cds, l_utr3, start_codon, stop_codon. See the 'Custom
                                      organisms' section of the README for a worked example.
  --codon-tables-file CODON_TABLES.json, --codon_tables_file CODON_TABLES.json
                                      Optional path to a codon-table JSON file. Supports either
                                      (a) ONE flat 64-codon mapping {codon: amino_acid}, or
                                      (b) multiple named tables {table_name: {codon: amino_acid}}.
                                      Combine with --codon_table_name when (b). Use this only when
                                      your organism's genetic code is not already in the built-in
                                      list — see --codon_table_name help.
  --codon-table-name TABLE_NAME, --codon_table_name TABLE_NAME
                                      Codon-table name to load from built-ins or from
                                      --codon_tables_file. The built-in tables are the NCBI
                                      Genetic Codes (ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
                                      with NCBI numbers converted to descriptive names. Pick the
                                      one matching your organism's MITOCHONDRIAL code (or NUCLEAR
                                      code for ciliate / candida / etc.):
                                        vertebrate_mitochondrial   (NCBI #2)  mouse, fish, frog,
                                                                             any non-human vertebrate
                                        yeast_mitochondrial        (NCBI #3)  S. cerevisiae and
                                                                             close relatives
                                        mold_mitochondrial         (NCBI #4)  Neurospora, Aspergillus,
                                                                             Trichoderma, mycoplasmas
                                        invertebrate_mitochondrial (NCBI #5)  fruit fly, mosquito,
                                                                             nematodes
                                        echinoderm_mitochondrial   (NCBI #9)  sea urchin, starfish
                                        ascidian_mitochondrial     (NCBI #13) sea squirts
                                        alternative_flatworm_mitochondrial (NCBI #14) flatworms
                                        trematode_mitochondrial    (NCBI #21) trematodes
                                        standard                   (NCBI #1)  use this when the
                                                                             organism uses the
                                                                             standard genetic code
                                                                             (most plants, including
                                                                             A. thaliana mt + plastid)
                                      Full list of 27 built-in tables: alternative_flatworm_mitochondrial, alternative_yeast_nuclear, ascidian_mitochondrial, bacterial, balanophoraceae_plastid, blastocrithidia_nuclear, blepharisma_macronuclear, candidate_division_sr1, cephalodiscidae_mitochondrial, chlorophycean_mitochondrial, ciliate_nuclear, condylostoma_nuclear, echinoderm_mitochondrial, euplotid_nuclear, invertebrate_mitochondrial, karyorelict_nuclear, mesodinium_nuclear, mold_mitochondrial, pachysolen_tannophilus_nuclear, peritrich_nuclear, pterobranchia_mitochondrial, scenedesmus_obliquus_mitochondrial, standard, thraustochytrium_mitochondrial, trematode_mitochondrial, vertebrate_mitochondrial, yeast_mitochondrial. [default: h.sapiens: vertebrate_mitochondrial, s.cerevisiae: yeast_mitochondrial]
  --start-codons CODON [CODON ...], --start_codons CODON [CODON ...]
                                      Allowed translation start codons used for codon classification.
                                      Defaults to strain presets; custom organisms can override them here. [default: y: ATG, h: ATG ATA, custom: ATG]
  --atp8-atp6-baseline {ATP6,ATP8}, --atp8_atp6_baseline {ATP6,ATP8}
                                      Preferred baseline name when resolving the shared ATP8/ATP6 bicistronic transcript.
                                      Titles remain ATP8/ATP6. [default: ATP6]
  --nd4l-nd4-baseline {ND4,ND4L}, --nd4l_nd4_baseline {ND4,ND4L}
                                      Preferred baseline name when resolving the shared ND4L/ND4 bicistronic transcript.
                                      Titles remain ND4L/ND4. [default: ND4]

Offset Enrichment and Selection:
  -a {start,stop}, --align {start,stop}
                                      Anchor offset enrichment around the start or stop codon. [default: start]
  -r NT, --range NT                   Plot offsets from -range to +range around the chosen anchor codon. [default: 20]
  --min-offset NT, --min_offset NT    Backward-compatible shared minimum absolute offset used only when the 5'/3' end-specific bounds are not provided. [default: 11]
  --max-offset NT, --max_offset NT    Backward-compatible shared maximum absolute offset used only when the 5'/3' end-specific bounds are not provided. [default: 20]
  --rpf-min-count-frac FRAC, --rpf_min_count_frac FRAC
                                      Drop read-length bins from the RPF window whose total count
                                      across all samples is below FRAC x the most-enriched length.
                                      FRAC=0 disables the filter. Default 0.20 keeps only read
                                      lengths with >=20% of the dominant-length count, so noisy
                                      low-count bins do not pollute offset selection. Pair with a
                                      wide --rpf range (e.g. 27 36 for human) to let the data
                                      decide which lengths actually carry signal. [default: 0.2]
  --min-5-offset NT, --min_5_offset NT
                                      Recommended minimum absolute 5' offset considered during offset selection. [default: same as --min_offset]
  --max-5-offset NT, --max_5_offset NT
                                      Recommended maximum absolute 5' offset considered during offset selection. [default: same as --max_offset]
  --min-3-offset NT, --min_3_offset NT
                                      Recommended minimum absolute 3' offset considered during offset selection. [default: same as --min_offset]
  --max-3-offset NT, --max_3_offset NT
                                      Recommended maximum absolute 3' offset considered during offset selection. [default: same as --max_offset]
  --offset-mask-nt NT, --offset_mask_nt NT
                                      Mask near-anchor bins from -N..-1 and +1..+N in offset summaries,
                                      line plots, and heatmaps. The first visible bins are -(N+1) and +(N+1). [default: 5]
  --offset-pick-reference {reported_site,p_site,selected_site}, --offset_pick_reference {reported_site,p_site,selected_site}
                                      Which coordinate space the best offset is chosen in.
                                        p_site         (default) pick the offset in canonical P-site
                                                       space first, then convert into the space
                                                       named by --offset_site if that is 'a'.
                                                       Use this when comparing across samples with
                                                       different offset_site choices.
                                        reported_site  pick the offset directly in the same space
                                                       named by --offset_site (no P<->A conversion).
                                                       Use this when you want what you see in the
                                                       enrichment table for --offset_site to be the
                                                       exact value picked.
                                        selected_site  DEPRECATED alias for 'reported_site'. [default: p_site]
  --offset-type {5,3}, --offset_type {5,3}
                                      Which read end defines downstream site placement after offsets are selected:
                                        5  measure offsets from the read 5' end
                                        3  measure offsets from the read 3' end [default: 5]
  --offset-site {p,a}, --offset_site {p,a}
                                      Coordinate space for the SELECTED OFFSETS table. p = P-site,
                                      a = A-site. This controls the column values in
                                      p_site_offsets_<align>.csv only. To control which downstream
                                      outputs (codon usage, coverage plots) are generated, use
                                      --analysis_sites. [default: p]
  --analysis-sites {p,a,both}, --analysis_sites {p,a,both}
                                      Which downstream outputs to generate per sample.
                                        both  (default) write both P-site and A-site codon usage and
                                              coverage plots, side by side under per-site subdirs.
                                        p     P-site outputs only.
                                        a     A-site outputs only.
                                      Independent of --offset_site, which only controls the offset
                                      selection coordinate space. [default: both]
  --codon-overlap-mode {full,any}, --codon_overlap_mode {full,any}
                                      How reads count toward the codon-level offset-enrichment
                                      table at each anchor codon (start or stop, per --align).
                                        full  (default) the read must span ALL 3 nt of the anchor
                                              codon. Example: anchor codon at positions 101-103;
                                              a 30-nt read at 100-129 counts (read covers 101, 102,
                                              103); a read at 102-131 does NOT count (does not
                                              cover position 101). Strict and the right default for
                                              ribo-seq footprints (>= ~26 nt).
                                        any   any 1+ nt overlap with the anchor codon counts. The
                                              102-131 read above WOULD count. Use only for very
                                              short reads relative to a codon (rare for ribo-seq). [default: full]
  -p NT, --psite-offset NT, --psite_offset NT
                                      Use one fixed offset for every read length and sample.
                                      This bypasses enrichment-based offset selection.
  --offset-mode {per_sample,combined}, --offset_mode {per_sample,combined}
                                      How offsets drive the downstream translation-profile and coverage plots.
                                        per_sample  (default) run enrichment + offset selection per
                                                    sample; each sample uses its own offsets for
                                                    downstream outputs. Offset drift across samples
                                                    is surfaced in offset_drift_<align>.svg.
                                        combined    select one offset table from every sample pooled
                                                    together and apply it uniformly; matches the
                                                    v0.3.x behavior. Use when you have very low
                                                    coverage per sample and need the pooled signal. [default: per_sample]

Read-Count Normalization:
  --read-counts-file COUNTS_TABLE, --read_counts_file COUNTS_TABLE
                                      Read-count table for RPM normalization (.csv, .tsv, or .txt; delimiter auto-detected). [default: read_counts_summary.txt]
  --read-counts-sample-col READ_COUNTS_SAMPLE_COL, --read_counts_sample_col READ_COUNTS_SAMPLE_COL
                                      Optional sample column override.
                                      If omitted, matching is case-insensitive and then falls back to column 1.
  --read-counts-reads-col READ_COUNTS_READS_COL, --read_counts_reads_col READ_COUNTS_READS_COL
                                      Optional read-count column override.
                                      If omitted, names like reads/read_count/counts are matched first, then column 3.
  --unfiltered-read-length-range MIN_LEN MAX_LEN, --unfiltered_read_length_range MIN_LEN MAX_LEN
                                      Inclusive read-length range used for the unfiltered QC summary and heatmaps.
                                      Broaden this when you want to inspect longer footprints such as disomes. [default: [15, 50]]
  --rpm-norm-mode {total,mt_mrna}, --rpm_norm_mode {total,mt_mrna}
                                      RPM denominator used for unfiltered and coverage-profile normalization:
                                        total    sum all rows in the read-count table
                                        mt_mrna  sum only rows whose reference matches
                                                 --mt_mrna_substring_patterns [default: total]
  --read-counts-reference-col READ_COUNTS_REFERENCE_COL, --read_counts_reference_col READ_COUNTS_REFERENCE_COL
                                      Optional reference column override.
                                      Needed for --rpm_norm_mode mt_mrna if auto-detection fails; otherwise falls back to column 2.
  --mt-mrna-substring-patterns PATTERN [PATTERN ...], --mt_mrna_substring_patterns PATTERN [PATTERN ...], --mrna_ref_patterns PATTERN [PATTERN ...]
                                      When --rpm_norm_mode mt_mrna is selected, the RPM denominator
                                      uses only rows from the read-count table whose value in the
                                      reference column (set via --read_counts_reference_col, or
                                      auto-detected) contains ANY of these substrings. Default
                                      matches reference names like 'mt_genome', 'mt-mrna', and
                                      'mt_mrna'. Example: in a read-count table with rows for
                                      'rrna', 'trna', 'mt_mrna_ND1', and 'mt_mrna_COX1', the
                                      default pattern keeps the latter two and drops the
                                      (r/t)RNA rows. The legacy flag name '--mrna_ref_patterns'
                                      is kept as a deprecated alias. [default: ['mt_genome', 'mt-mrna', 'mt_mrna']]

Outputs and Plotting:
  -o OUTPUT_DIR, --output OUTPUT_DIR  Base output directory for all pipeline results, including mitoribopy.log. [default: analysis_results]
  --downstream-dir DOWNSTREAM_DIR, --downstream_dir DOWNSTREAM_DIR
                                      Per-sample subdirectory name for frame and codon analyses. [default: footprint_density]
  --plot-dir PLOT_DIR, --plot_dir PLOT_DIR
                                      Shared subdirectory name for offset CSV files and plots. [default: offset_diagnostics]
  -fmt {png,pdf,svg}, --plot-format {png,pdf,svg}, --plot_format {png,pdf,svg}
                                      File format used for saved plots. [default: png]
  --x-breaks NT [NT ...], --x_breaks NT [NT ...]
                                      Optional custom x-axis tick marks for offset line plots.
  --line-plot-style {combined,separate}, --line_plot_style {combined,separate}
                                      Draw offset line plots in one combined panel or separate 5'/3' panels. [default: combined]
  --cap-percentile CAP_PERCENTILE, --cap_percentile CAP_PERCENTILE
                                      Upper percentile cap for coverage-style plots (0 < value <= 1). [default: 0.999]
  -m, --codon-density-window, --codon_density_window, --merge_density
                                      When set, the codon-level density at each codon centre is summed with its +/-1 nt neighbours (a 3-nt sliding window) before being written into the codon coverage / usage tables and plots. Smooths short-window noise around the codon centre; does NOT collapse reading frames despite the historical flag name '--merge_density' (which is kept as a deprecated alias).
  --order-samples ORDER_SAMPLES [ORDER_SAMPLES ...], --order_samples ORDER_SAMPLES [ORDER_SAMPLES ...]
                                      Optional sample order for Ribo-seq plots and aggregated outputs.

Optional Modules:
  --structure-density, --structure_density
                                      Generate structure-density exports from footprint-density tables.
  --structure-density-norm-perc STRUCTURE_DENSITY_NORM_PERC, --structure_density_norm_perc STRUCTURE_DENSITY_NORM_PERC
                                      Upper percentile used to cap and scale structure-density values. [default: 0.99]
  --cor-plot, --cor_plot              Generate codon-correlation plots.
  --base-sample BASE_SAMPLE, --base_sample BASE_SAMPLE
                                      Reference sample for codon-correlation comparisons.
  --cor-mask-method {percentile,fixed,none}, --cor_mask_method {percentile,fixed,none}
                                      Masking rule for extreme codon-correlation outliers. [default: percentile]
  --cor-mask-percentile COR_MASK_PERCENTILE, --cor_mask_percentile COR_MASK_PERCENTILE
                                      Percentile cutoff used when --cor_mask_method percentile. [default: 0.99]
  --cor-mask-threshold COR_MASK_THRESHOLD, --cor_mask_threshold COR_MASK_THRESHOLD
                                      Fixed absolute cutoff used when --cor_mask_method fixed.
  --cor-metric {log2_density_rpm,log2_rpm,linear,raw_count}, --cor_metric {log2_density_rpm,log2_rpm,linear,raw_count}
                                      Primary codon-correlation metric. log2_density_rpm (default) log-transforms RPM-normalised codon density; raw_count reproduces the legacy depth-dominated scatter (and emits W_CODON_RAW_COUNT_PRIMARY). [default: log2_density_rpm]
  --cor-regression {theil_sen,ols,none}, --cor_regression {theil_sen,ols,none}
                                      Robust regression method drawn through the codon-correlation scatter. theil_sen (default) is median-based and tolerant of outliers; ols reproduces the legacy linear fit. [default: theil_sen]
  --cor-support-min-raw COR_SUPPORT_MIN_RAW, --cor_support_min_raw COR_SUPPORT_MIN_RAW
                                      Minimum raw value required in BOTH samples for a codon to be marked include_primary. Low-support codons remain in the metrics TSV but are dimmed in the figure and excluded from the residual ranking. [default: 10]
  --cor-label-top-n COR_LABEL_TOP_N, --cor_label_top_n COR_LABEL_TOP_N
                                      Number of codons to label on the scatter (by support-aware label score). [default: 10]
  --cor-pseudocount COR_PSEUDOCOUNT, --cor_pseudocount COR_PSEUDOCOUNT
                                      Additive pseudocount before log2. 'auto' or float. [default: auto]
  --cor-raw-panel {qc_only,off}, --cor_raw_panel {qc_only,off}
                                      Where to write the raw-count panel when metric=raw_count. 'qc_only' (default) routes it to raw_count_qc/; 'off' skips it. [default: qc_only]
  --read-coverage-raw, --no-read-coverage-raw, --read_coverage_raw, --no-read_coverage_raw
                                      Write read-coverage plots in raw counts under coverage_profile_plots/read_coverage_raw[_codon]/. Use --no-read_coverage_raw to skip. [default: True]
  --read-coverage-rpm, --no-read-coverage-rpm, --read_coverage_rpm, --no-read_coverage_rpm
                                      Write read-coverage plots in RPM under coverage_profile_plots/read_coverage_rpm[_codon]/. Use --no-read_coverage_rpm to skip. [default: True]
  --igv-export, --no-igv-export, --igv_export, --no-igv_export
                                      Export per-sample BedGraph tracks (P-site / A-site) under <output>/igv_tracks/<sample>/, suitable for opening in IGV.
  --periodicity-enabled, --no-periodicity-enabled, --periodicity_enabled, --no-periodicity_enabled
                                      Skip the periodicity QC step entirely when set to false. [default: True]
  --periodicity-fourier-window-nt PERIODICITY_FOURIER_WINDOW_NT, --periodicity_fourier_window_nt PERIODICITY_FOURIER_WINDOW_NT
                                      Window (nt) for the Fourier metagene per region (orf_start, orf_stop). Default: 99 (33 codons; multiple of 3 for clean period-3 bin alignment).
  --periodicity-metagene-nt PERIODICITY_METAGENE_NT, --periodicity_metagene_nt PERIODICITY_METAGENE_NT
                                      Window (nt) up/downstream of start/stop codons for the metagene_start.tsv / metagene_stop.tsv plots. Default: 300 (100 codons).
  --periodicity-metagene-normalize {per_gene_unit_mean,none}, --periodicity_metagene_normalize {per_gene_unit_mean,none}
                                      How metagene_{start,stop}.tsv aggregates per-transcript signals. 'per_gene_unit_mean' (default) divides each transcript's per-position density by its own mean before averaging, which removes the depth-weighting bias where one high-expression transcript dominates the metagene shape. 'none' preserves the older raw-position-count sum for users that need depth-weighted metagene values. [default: per_gene_unit_mean]
  --periodicity-fourier-bootstrap-n PERIODICITY_FOURIER_BOOTSTRAP_N, --periodicity_fourier_bootstrap_n PERIODICITY_FOURIER_BOOTSTRAP_N
                                      Bootstrap iterations for the percentile CI on the metagene Fourier ratio. Default: 200. Set to 0 to disable the CI without disabling the permutation null.
  --periodicity-fourier-permutations-n PERIODICITY_FOURIER_PERMUTATIONS_N, --periodicity_fourier_permutations_n PERIODICITY_FOURIER_PERMUTATIONS_N
                                      Circular-shift permutations for the empirical null on the metagene Fourier ratio. Default: 200. Set to 0 to disable the null without disabling the CI.
  --periodicity-fourier-ci-alpha PERIODICITY_FOURIER_CI_ALPHA, --periodicity_fourier_ci_alpha PERIODICITY_FOURIER_CI_ALPHA
                                      Two-sided alpha for the Fourier percentile bootstrap CI. Default: 0.10 (90% CI).
  --periodicity-fourier-random-seed PERIODICITY_FOURIER_RANDOM_SEED, --periodicity_fourier_random_seed PERIODICITY_FOURIER_RANDOM_SEED
                                      RNG seed for the Fourier bootstrap + permutation draws. Default: 42. Recorded in periodicity.metadata.json so a reviewer can reproduce the exact CI / p-value bounds.
  --periodicity-no-fourier-stats, --periodicity_no_fourier_stats
                                      Skip the bootstrap CI + circular-shift permutation null on the metagene Fourier ratio. Faster, but the score table loses CI / permutation_p columns.

Examples:
  mitoribopy rpf --strain h.sapiens --fasta ref.fa --directory beds \
                 -rpf 29 34 --align stop --offset-type 5 \
                 --offset-site p --offset-mask-nt 5 --output out \
                 --codon-density-window
  mitoribopy rpf ... --min-5-offset 10 --max-5-offset 22 \
                     --min-3-offset 12 --max-3-offset 30
  mitoribopy rpf --strain custom --fasta ref.fa --directory beds \
                 --annotation-file ann.csv \
                 --codon-tables-file codon_tables.json \
                 --codon-table-name standard

Tip: prefer the end-specific 5'/3' offset bounds.
Use --min-offset / --max-offset only as shared fallback bounds.
Underscore-style flags (--offset_type, --min_5_offset, ...) are
still accepted for one transition cycle but no longer shown here.
```

## `mitoribopy rnaseq`

```text
usage: mitoribopy rnaseq [-h] [--config CONFIG] [--dry-run] [--threads N]
                         [--log-level {DEBUG,INFO,WARNING,ERROR}]
                         [--rnaseq-mode MODE] [--sample-sheet PATH]
                         [--rna-fastq PATH [PATH ...]]
                         [--ribo-fastq PATH [PATH ...]]
                         [--reference-fasta PATH] [--bowtie2-index PREFIX]
                         [--workdir DIR] [--align-threads N]
                         [--recount-ribo-fastq] [--align-only]
                         [--allow-pseudo-replicates-for-demo-not-publication]
                         [--gene-id-convention {ensembl,refseq,hgnc,mt_prefixed,bare}]
                         [--organism {h.sapiens,s.cerevisiae,h,y,human,yeast}]
                         [--condition-map PATH] [--condition-a NAME]
                         [--condition-b NAME] [--base-sample NAME]
                         [--compare-sample NAME] [--output DIR]
                         [--de-table PATH]
                         [--de-format {auto,deseq2,xtail,anota2seq,custom}]
                         [--de-gene-col NAME] [--de-log2fc-col NAME]
                         [--de-padj-col NAME] [--de-basemean-col NAME]
                         [--ribo-dir DIR] [--ribo-counts PATH]
                         [--reference-gtf PATH] [--reference-checksum SHA256]

Translation efficiency (TE / delta-TE) from paired RNA-seq + Ribo-seq. Default flow: pass --rna-fastq + --ribo-fastq + --reference-fasta and the subcommand runs trimming, bowtie2 alignment, per-transcript counting, and pyDESeq2 itself before emitting te.tsv, delta_te.tsv, and plots. Alternative: pass --de-table from a prior external DESeq2 / Xtail / Anota2Seq run together with --ribo-dir; this path is mutually exclusive with --rna-fastq and enforces a SHA256 reference-consistency gate.

options:
  -h, --help                          show this help message and exit
  --rnaseq-mode MODE                  Explicit rnaseq stage mode. Recommended over relying on input-presence inference. Modes:
                                        de_table   PUBLICATION: consume external DE table + prior rpf run (--de-table + --ribo-dir).
                                        from_fastq EXPLORATORY: run cutadapt + bowtie2 + pyDESeq2 on the mt-mRNA subset (--rna-fastq or --sample-sheet).
                                        none       stage configured but inert (no inputs).
                                      When omitted, the mode is inferred from supplied inputs.

Shared options:
  --config CONFIG                     Configuration file (.json, .yaml, .yml, or .toml). CLI arguments override values read from the file.
  --dry-run                           Print planned actions and exit without executing.
  --threads N                         Preferred thread count for external tools and BLAS libraries.
  --log-level {DEBUG,INFO,WARNING,ERROR}
                                      Python logging level for MitoRiboPy console output. [default: INFO]

Inputs (from raw FASTQ — default flow):
  --sample-sheet PATH                 Unified sample sheet TSV (one row per FASTQ; columns documented in mitoribopy.sample_sheet). When provided, --rna-fastq, --ribo-fastq, and --condition-map are derived from the sheet and passing any of those flags alongside --sample-sheet is an error. Pairing between RNA-seq and Ribo-seq is by sample_id; index-based pairing is no longer supported.
  --rna-fastq PATH [PATH ...]         RNA-seq FASTQ files or directories. The default driver of the rnaseq subcommand: pass these and the rest of the pipeline (trim → bowtie2 → count → pyDESeq2 → TE / ΔTE / plots) runs end-to-end. Mutually exclusive with --de-table and --sample-sheet.
  --ribo-fastq PATH [PATH ...]        Ribo-seq FASTQ files or directories. When omitted the run short-circuits after writing de_table.tsv (mode='from-fastq-rna-only' in run_settings.json).
  --reference-fasta PATH              Transcriptome FASTA used to build (or reuse) a content-addressed bowtie2 index. Required for the default flow. The SHA256 is recorded under from_fastq.reference_checksum in run_settings.json.
  --bowtie2-index PREFIX              Pre-built bowtie2 index prefix. When set, bowtie2-build is skipped; --reference-fasta is still required for hashing.
  --workdir DIR                       Scratch directory for trim / index / BAM artefacts. Defaults to <output>/work.
  --align-threads N                   Threads passed to cutadapt and bowtie2. Separate from --threads (which caps BLAS / pyDESeq2 thread pools). [default: 4]
  --recount-ribo-fastq                When 'mitoribopy all' has just produced rpf_counts.tsv from the rpf stage, the rnaseq from-FASTQ flow defaults to REUSING those counts instead of re-aligning the Ribo FASTQs (rpf already did the work; re-running cutadapt + bowtie2 + counting is duplicate effort and a source of subtle inconsistencies between the two stages' counts). Pass --recount-ribo-fastq to force a second pass over the Ribo FASTQs from inside the rnaseq stage.
  --align-only                        Run only the FASTQ-trim + bowtie2 + counts phase and exit before DESeq2 / TE / plots. Used by `mitoribopy all` to kick off RNA-seq alignment in parallel with the Ribo-seq align stage so total wall time is dominated by the slower of the two assays instead of the sum. Writes `rna_counts.tsv` (and `rpf_counts.tsv` if Ribo FASTQs are supplied) to the output dir; subsequent rnaseq invocations with --upstream-rna-counts pick up the prebuilt matrix.
  --allow-pseudo-replicates-for-demo-not-publication
                                      Opt INTO the auto-pseudo-replicate fallback for conditions with only 1 sample. Without this flag, an n=1 design is a hard error (exit 2) — the safe default for publication-grade DE. With this flag, each n=1 condition is split into rep1 / rep2 by FASTQ record parity so pyDESeq2 has n>=2 to fit dispersion on, BUT the resulting p-values and padj are NOT biologically defensible (the two halves are mechanical subsamples of one library). Use ONLY for demos / tutorials. The run_settings.json records pseudo_replicate_mode=true and an EXPLORATORY.md sidecar is written to the output dir.

Gene identifiers:
  --gene-id-convention {ensembl,refseq,hgnc,mt_prefixed,bare}
                                      Gene identifier convention used in the FASTA / DE table. REQUIRED (no default): mismatched conventions silently produce zero-match runs.
  --organism {h.sapiens,s.cerevisiae,h,y,human,yeast}
                                      Organism for the mt-mRNA registry. Use 'h.sapiens' for human, 's.cerevisiae' for budding yeast. Short / spelled-out forms ('h', 'y', 'human', 'yeast') are accepted as synonyms. [default: h.sapiens]

Conditions (required in default flow; optional in --de-table flow):
  --condition-map PATH                TSV with columns 'sample' and 'condition'. Required by the default flow (drives the pyDESeq2 contrast). In the --de-table flow it is optional and enables a replicate-based Ribo log2FC for delta-TE; without it the delta-TE rows carry only the mRNA log2FC and a 'single_replicate_no_statistics' note.
  --condition-a NAME                  Reference condition (denominator of the WT-vs-X contrast). Used as the baseline in DE / TE comparison plots. ``--base-sample`` is the preferred spelling (matches the rpf config's ``base_sample`` key); both flags resolve to the same internal field.
  --condition-b NAME                  Comparison condition (numerator of the contrast). ``--compare-sample`` is the preferred alias.
  --base-sample NAME                  Reference condition for the contrast and for all per-gene comparison plots (mRNA / RPF DE volcanoes, TE compare scatter, TE log2FC bar). Alias for ``--condition-a``; named to mirror the rpf config's ``base_sample`` key. When both ``--base-sample`` and ``--condition-a`` are provided they must agree.
  --compare-sample NAME               Comparison condition (alias for ``--condition-b``). When both ``--compare-sample`` and ``--condition-b`` are provided they must agree.

Output:
  --output DIR                        Output directory for te.tsv, delta_te.tsv, and plots.

Alternative inputs (bring your own DE table; mutually exclusive with --rna-fastq):
  --de-table PATH                     Pre-computed DESeq2 / Xtail / Anota2Seq results table (CSV or TSV). Use this when you already ran DE externally on the full transcriptome — the recommended path for publication-grade DE statistics. Requires --ribo-dir / --ribo-counts and one of --reference-gtf / --reference-checksum. Mutually exclusive with --rna-fastq.
  --de-format {auto,deseq2,xtail,anota2seq,custom}
                                      DE table format. 'auto' detects from column headers. [default: auto]
  --de-gene-col NAME
  --de-log2fc-col NAME
  --de-padj-col NAME
  --de-basemean-col NAME

Ribo-seq inputs (--de-table flow):
  --ribo-dir DIR                      Directory produced by a prior 'mitoribopy rpf' run; expected to contain rpf_counts.tsv and run_settings.json / run_manifest.json with a recorded reference_checksum.
  --ribo-counts PATH                  Explicit path to rpf_counts.tsv. Defaults to <ribo-dir>/rpf_counts.tsv when --ribo-dir is set.

Reference-consistency gate (--de-table flow; exactly one):
  --reference-gtf PATH                Reference GTF / FASTA used by RNA-seq; we hash this and verify it matches the hash recorded in the rpf run's manifest. EXACTLY ONE of --reference-gtf / --reference-checksum must be provided in the --de-table flow.
  --reference-checksum SHA256         Precomputed SHA-256 of the shared reference (use instead of --reference-gtf when the file is not on this host).
```

## `mitoribopy all`

```text
usage: mitoribopy all [-h] [--config CONFIG] [--dry-run] [--threads N]
                      [--log-level {DEBUG,INFO,WARNING,ERROR}] [--output DIR]
                      [--resume] [--force-resume] [--skip-align] [--skip-rpf]
                      [--skip-rnaseq] [--manifest PATH]
                      [--show-stage-help STAGE] [--print-config-template]
                      [--profile {minimal,publication,exhaustive}]
                      [--print-canonical-config] [--strict] [--progress MODE]
                      [--progress-file PATH] [--no-progress]

End-to-end orchestrator: align + rpf, plus rnaseq when the config carries an 'rnaseq' section configured for either flow (from-FASTQ via 'rna_fastq' + 'reference_fasta', or external-DE via 'de_table'). Writes a composed run_manifest.json with tool versions, parameters, and input/output hashes across all three stages.

options:
  -h, --help                          show this help message and exit
  --output DIR                        Run root directory. Each stage writes under <output>/align/, <output>/rpf/, and <output>/rnaseq/ when the stage is active.
  --resume                            Skip stages whose expected output already exists: align is skipped when <output>/align/read_counts.tsv is present; rpf when <output>/rpf/rpf_counts.tsv is present; rnaseq when <output>/rnaseq/delta_te.tsv is present. The skip decision is gated by a hash check against the prior run_manifest.json (config_source_sha256, sample_sheet_sha256, reference_checksum, mitoribopy_version, schema_version); edits to any of those fields force the affected stage(s) to re-run unless --force-resume is also set.
  --force-resume                      Like --resume, but bypass the hash guard. Use only when you know the stage outputs are still valid for the new config (e.g. you edited a comment-only line). Also honoured via the MITORIBOPY_FORCE_RESUME=1 environment variable for CI scripts that cannot easily change argv.
  --skip-align                        Skip the align stage even when an [align] section exists.
  --skip-rpf                          Skip the rpf stage.
  --skip-rnaseq                       Skip the rnaseq stage even when an [rnaseq] section exists.
  --manifest PATH                     Manifest filename (relative to --output). [default: run_manifest.json]
  --show-stage-help STAGE             Print the full help for one stage and exit. Useful because 'mitoribopy all --help' only shows orchestrator-level flags.
  --print-config-template             Print a commented YAML config template covering every stage (align / rpf / rnaseq) with sensible defaults, then exit. Pipe this into a file to start a new project: 'mitoribopy all --print-config-template > pipeline_config.yaml'. Pair with --profile to pick which template to emit (default: minimal).
  --profile {minimal,publication,exhaustive}
                                      Template profile for --print-config-template. 'minimal' (default): the curated 80-line template with the most-edited keys. 'publication': adds publication-readiness defaults (--strict, rnaseq_mode=de_table, fourier_bootstrap_n=200) so a methods-paper run uses the right gates. 'exhaustive': prints the full annotated example from examples/templates/pipeline_config.example.yaml — every single flag with its default and a one-line comment. [default: minimal]
  --print-canonical-config            Load --config, apply every auto-wiring + sample-sheet expansion that 'mitoribopy all' would normally apply, then print the resulting canonical config to stdout (YAML if PyYAML is available, JSON otherwise) and exit. Useful for diffing your input config against what was actually executed: 'mitoribopy all --print-canonical-config --config pipeline_config.yaml --output results/'. The same blob is embedded in run_manifest.json under 'config_canonical' on real runs.
  --strict                            Publication-safe mode. A single switch that forwards strictness to every stage and post-run validation:
                                        * align: --strict-publication-mode (fail on non-default policies that would invalidate a publication run);
                                        * config: --strict on the up-front validate-config pass (treat any deprecated-key rewrite or unknown stage key as a hard error);
                                        * rnaseq: refuse 'allow_pseudo_replicates_for_demo_not_publication: true' and refuse 'rnaseq_mode: from_fastq' (the mt-mRNA-only DE path) unless the user explicitly opts back in via 'allow_exploratory_from_fastq_in_strict: true';
                                        * figures: --strict on the post-run validate-figures pass (promote warn-only QC findings to fail);
                                        * summary: warning rows in warnings.tsv are mirrored as WARN bullets in SUMMARY.md.
                                      Recommended for any run that backs a manuscript figure.
  --progress MODE                     Progress display mode. 'auto' (default) picks a tqdm bar on an interactive TTY and one-line plain logs on HPC / non-TTY streams. 'plain' forces stable [PROGRESS] log lines. 'bar' forces tqdm bars (falls back to plain if tqdm is missing). 'rich' degrades to bar in this version. 'jsonl' streams machine-readable JSON events to stderr (or to --progress-file). 'off' silences all progress output, but --progress-file still attaches if set. [default: auto]
  --progress-file PATH                Append every typed progress event to this file as JSONL (one JSON object per line). Independent of --progress: an explicit --progress-file always attaches, even when --progress=off. Default location is <output>/progress.jsonl when --output is set; pass explicitly to override.
  --no-progress                       Shortcut for --progress=off.

Shared options:
  --config CONFIG                     Configuration file (.json, .yaml, .yml, or .toml). CLI arguments override values read from the file.
  --dry-run                           Print planned actions and exit without executing.
  --threads N                         Preferred thread count for external tools and BLAS libraries.
  --log-level {DEBUG,INFO,WARNING,ERROR}
                                      Python logging level for MitoRiboPy console output. [default: INFO]

This subcommand owns only the orchestrator flags. Stage-specific options live
inside the config file sections whose keys match the subcommand flags
(with dashes replaced by underscores, e.g. '--adapter' -> 'adapter').
  align:  keys for 'mitoribopy align --help'
  rpf:    keys for 'mitoribopy rpf --help'
  rnaseq: keys for 'mitoribopy rnaseq --help'

Start a new project:
  mitoribopy all --print-config-template > pipeline_config.yaml
  # edit the file to point at your FASTQs / indexes, then:
  mitoribopy all --config pipeline_config.yaml --output results/

Inspect a stage's full flag list:
  mitoribopy all --show-stage-help align
  mitoribopy all --show-stage-help rpf
  mitoribopy all --show-stage-help rnaseq
```

## `mitoribopy periodicity`

```text
usage: mitoribopy periodicity [-h] --site-table PATH --output DIR
                              [--site {p,a}] [--fourier-window-nt N]
                              [--drop-codons-after-start N]
                              [--drop-codons-before-stop N]
                              [--min-mean-coverage X] [--min-total-counts N]
                              [--no-plots] [--no-stats] [--bootstrap-n B]
                              [--permutations-n P] [--ci-alpha A]
                              [--random-seed S]

Quantify 3-nt periodicity by running the metagene Fourier analysis on a pre-assigned site table.

options:
  -h, --help                   show this help message and exit
  --site-table PATH            Per-read site table; required columns: sample, gene, transcript_id, read_length, site_type, site_pos, cds_start, cds_end. `count` is optional (defaults to 1).
  --output DIR                 Directory for the Fourier QC bundle.
  --site {p,a}                 Ribosomal site to score. Default: p. [default: p]
  --fourier-window-nt N        Window size (nt) per region. Default: 99 (33 codons). Must be a multiple of 3 for clean period-3 bin alignment. [default: 99]
  --drop-codons-after-start N  Codons after the AUG to skip in the orf_start window. Default: 5 (skip the initiation peak). [default: 5]
  --drop-codons-before-stop N  Codons before the stop codon to skip in the orf_stop window. Default: 1 (skip the termination peak). [default: 1]
  --min-mean-coverage X        Skip per-gene windows whose mean coverage is below X. Default: 0.1. [default: 0.1]
  --min-total-counts N         Skip per-gene windows whose total site count is below N. Default: 30. [default: 30]
  --no-plots                   Skip the per-(sample, read_length) figures; TSVs still written. [default: True]

statistical hardening (bootstrap CI + permutation null):
  --no-stats                   Skip the bootstrap CI + circular-shift permutation null. Faster, but the score table loses amp_3nt_ci_*, spectral_ratio_3nt(_local)_ci_*, and permutation_p columns. Use only for smoke runs — the publication-facing path keeps stats on. [default: True]
  --bootstrap-n B              Bootstrap iterations for the percentile CI over genes. Default: 200. [default: 200]
  --permutations-n P           Circular-shift permutations for the empirical null on the spectral ratios. Default: 200. [default: 200]
  --ci-alpha A                 Two-sided alpha for the percentile bootstrap CI. Default: 0.1 (90% CI). [default: 0.1]
  --random-seed S              RNG seed for the bootstrap and permutation draws. Default: 42. The seed is recorded in periodicity.metadata.json so a reviewer can reproduce the exact CI / p-value bounds. [default: 42]

Frame definition: (site_pos - cds_start) mod 3. The Fourier analysis is anchored at the start codon (orf_start) and the stop codon (orf_stop). site_pos must be transcript-oriented (forward-strand-relative) and 0-based.
```

## `mitoribopy migrate-config`

```text
usage: mitoribopy migrate-config [-h] [--in-place] PATH

Rewrite legacy MitoRiboPy YAML keys to their canonical names. Input is read from a path; output is written to stdout (the change log goes to stderr). Use to upgrade old pipeline configs without manually hunting down every renamed key.

positional arguments:
  PATH        Path to the legacy YAML / JSON / TOML config file.

options:
  -h, --help  show this help message and exit
  --in-place  Overwrite the input file in place (with a `.bak` backup) instead of writing the canonical config to stdout.
```

## `mitoribopy validate-config`

```text
usage: mitoribopy validate-config [-h] [--no-path-checks] [--strict]
                                  [--stage {align,rpf,rnaseq}]
                                  PATH

Pre-flight a MitoRiboPy YAML / JSON / TOML config: parse, canonicalise legacy keys, check file paths and mutually-exclusive sections, and resolve rnaseq.mode against supplied inputs. Exit code is 0 on success, 2 when at least one error was found.

positional arguments:
  PATH                        Path to the YAML / JSON / TOML config to validate.

options:
  -h, --help                  show this help message and exit
  --no-path-checks            Skip the on-disk existence checks for path-shaped values. Useful when validating a config on a different host than where the run will execute (e.g. CI-side validation).
  --strict                    Treat any legacy-key rewrite as a warning AND make the validator exit 2 when at least one warning fired. Use for publication-grade configs that should already be canonical.
  --stage {align,rpf,rnaseq}  Pin the per-stage parser to validate against when the file is a flat per-stage config (no align: / rpf: / rnaseq: wrapper — the layout the standalone subcommands accept via `mitoribopy <stage> --config <file>`). Without this flag the validator infers the stage by matching keys against each per-stage parser's argparse dests; the flag exists so you can disambiguate when inference is ambiguous and so CI scripts can declare intent. Ignored for orchestrator-shape configs.
```

## `mitoribopy validate-reference`

```text
usage: mitoribopy validate-reference [-h] --fasta FASTA --annotation
                                     ANNOTATION
                                     [--codon-table {standard,vertebrate_mitochondrial,yeast_mitochondrial}]

Pre-flight a custom mitochondrial reference: check that the FASTA and annotation CSV are consistent (matching transcript IDs, matching lengths, CDS divisible by 3, valid start / stop codons under the selected codon table).

options:
  -h, --help                          show this help message and exit
  --fasta FASTA                       Path to the custom mt-transcriptome FASTA.
  --annotation ANNOTATION             Path to the per-transcript annotation CSV (the same schema accepted by `mitoribopy rpf --annotation_file`).
  --codon-table {standard,vertebrate_mitochondrial,yeast_mitochondrial}
                                      Codon table for start / stop codon validation. Default is vertebrate_mitochondrial (NCBI translation table 2). [default: vertebrate_mitochondrial]
```

## `mitoribopy validate-figures`

```text
usage: mitoribopy validate-figures [-h] [--strict] [--require-png-dpi DPI]
                                   [--out PATH]
                                   RUN_DIR

Mechanically validate every plot under a finished MitoRiboPy run root: check label / legend / stat-box overlap, label clipping, point counts vs source TSV, SVG text editability, PNG dpi, and metadata sidecar coverage. Writes <RUN_DIR>/figure_qc.tsv. Exit 0 / 1 / 2 (all pass / warn-only / fail; --strict upgrades warn → fail).

positional arguments:
  RUN_DIR                The output directory of a finished `mitoribopy all` run.

options:
  -h, --help             show this help message and exit
  --strict               Treat warn-level conditions (low PNG dpi, non-editable SVG text) as fail. Hard contract violations (point-count mismatch, label-overlap > 0) are always fail regardless.
  --require-png-dpi DPI  Minimum acceptable PNG dpi (default 300). [default: 300]
  --out PATH             Override the figure_qc.tsv output path. Defaults to <RUN_DIR>/figure_qc.tsv.
```

## `mitoribopy summarize`

```text
usage: mitoribopy summarize [-h] [--manifest NAME] RUN_DIR

Regenerate SUMMARY.md and summary_qc.tsv from a finished MitoRiboPy run by reading the run_manifest.json and per-stage TSV outputs. Useful for re-rendering summaries on archival runs without re-executing any pipeline stage.

positional arguments:
  RUN_DIR          The output directory of a finished `mitoribopy all` run.

options:
  -h, --help       show this help message and exit
  --manifest NAME  Manifest filename within RUN_DIR. [default: run_manifest.json]
```

## `mitoribopy benchmark`

```text
usage: mitoribopy benchmark [-h] --config PATH --output DIR [--subsample N]
                            [--threads N] [--seed SEED]

Time and disk-measure a full `mitoribopy all` run, optionally after pre-subsampling each FASTQ to N reads. Produces benchmark.tsv and benchmark_summary.md at the run root for tuning thread counts, disk budgets, and per-stage wall time.

options:
  -h, --help     show this help message and exit
  --config PATH  Pipeline YAML / JSON / TOML, exactly as accepted by `mitoribopy all --config`.
  --output DIR   Run root for the benchmarked execution. Will be created if missing; existing contents are not deleted.
  --subsample N  If set, reservoir-sample each Ribo / RNA FASTQ to N reads before running. Subsampled copies are written under <output>/.benchmark_subsamples/ and the canonical config is rewritten to point at them. Requires reservoir-friendly FASTQ inputs (no streaming sources).
  --threads N    Forwarded to `mitoribopy all --threads`. Recorded in the benchmark.tsv `threads` column.
  --seed SEED    Random seed for the FASTQ subsampler (no effect when --subsample is omitted). [default: 42]
```

