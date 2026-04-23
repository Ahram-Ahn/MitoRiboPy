# MitoRiboPy v0.3.0-rc1 &mdash; release readiness report

## Status

**Ready for PR review on `refactor/v0.3.0-preprocessing` -> `main`.**
Push access and PR creation are blocked in this sandbox (HTTP 403 on
both the local git proxy and the GitHub MCP integration). The 19
commits below sit on the local branch; Ahram can push them from a
workstation where his GitHub credentials are configured, then tag
`v0.3.0-rc1` against the tip and open the draft PR.

## Final test matrix

- `python -m pytest` &mdash; **189 passed, 1 skipped** on Python 3.11.15.
  The single skipped test is `tests/test_align_integration.py` (marked
  `@pytest.mark.requires_tools`); it auto-skips when `cutadapt`,
  `bowtie2`, `bowtie2-build`, or `samtools` is not on `$PATH`. It must
  run successfully in a bioconda environment before tagging v0.3.0
  final.
- `python -m build` &mdash; `mitoribopy-0.3.0.tar.gz` +
  `mitoribopy-0.3.0-py3-none-any.whl` build cleanly.

### Artifact hashes (SHA-256)

```
5be12522cd6407784d8b845a271f7192db3c09bf97bc25f0117dc9a01018854c  mitoribopy-0.3.0-py3-none-any.whl
ab679aafa8209740b06e6c844176816035d8de135c7c5b41aca7f2f441a72a57  mitoribopy-0.3.0.tar.gz
```

### Wheel contents spot-check

- All Phase-3 modules present: `mitoribopy/align/{__init__,_types,align,
  bam_utils,contam,dedup,read_counts,tool_check,trim}.py`.
- All Phase-5 modules present: `mitoribopy/rnaseq/{__init__,_types,
  counts,de_loader,gene_ids,plots,reference_gate,te}.py`.
- Subcommand CLI package: `mitoribopy/cli/{__init__,align,all_,common,
  rnaseq,rpf}.py`.
- Legacy files correctly absent: no `main.py`, no
  `mitoribopy/config/{defaults,loader,models}.py`, no `_legacy/`.
- Data files shipped: `data/codon_tables.json`,
  `data/{human,yeast}_annotation.csv`.
- Entry point: `mitoribopy = mitoribopy.cli:main`.

## Commit history on `refactor/v0.3.0-preprocessing`

19 commits layered on top of `7babcf5` (v0.2.0 tip of main):

```
de6f0d4 docs: v0.3.0 tutorials, CLI reference, and validation plan
a07da61 feat: mitoribopy all end-to-end orchestrator with provenance manifest
4faf922 feat: mitoribopy rnaseq for DE-based translation-efficiency integration
4ff7bce feat: rpf accepts BAM input in addition to BED
172d1e4 docs: v0.3.0 Phase 3 CHANGELOG entries and retire align stub tests
797194d test: mitoribopy align end-to-end integration test (requires_tools)
3ec44e1 docs: bioconda environment.yml and Dockerfile for mitoribopy align deps
b5d3b32 feat: wire mitoribopy align CLI end-to-end
7c2fb4b feat: add mitoribopy align - per-sample read-counts provenance table
fcd12f1 feat: add mitoribopy align - BAM utilities (MAPQ filter, BAM->BED6) via pysam
e52c264 feat: add mitoribopy align - dedup (umi_tools / skip / mark-duplicates)
ca4a3eb feat: add mitoribopy align - bowtie2 mt-transcriptome alignment (Path A)
b74a2f0 feat: add mitoribopy align - bowtie2 contaminant subtraction
86eb2ad feat: add mitoribopy align - trimming via cutadapt with kit presets
74e5eff feat: add mitoribopy align - tool availability checks and shared types
1d62362 feat: introduce subcommand CLI (align/rpf/rnaseq/all) with rpf as default fallback
3ea81b7 breaking: remove main.py shim, vestigial config scaffolding, and --varna* aliases
ea23103 docs: v0.3.0 Phase 0 baseline audit (inventory, CLI surface, dep graph)
```

(Plus the Phase 8 readiness commit this file lands in.)

## Biological invariants preserved

Verified through targeted unit tests and the Phase-3 strand-aware design:

1. Vertebrate + yeast mitochondrial codon tables unchanged
   (`tests/test_codon_tables.py`; `data/codon_tables.json` untouched
   since pre-v0.3.0).
2. Bicistronic ATP8/ATP6 and ND4L/ND4 handling preserved
   (`tests/test_reference_data.py`; `--atp8_atp6_baseline` and
   `--nd4l_nd4_baseline` flags untouched).
3. Frame-0 reference + P/A-site shift behavior preserved
   (`tests/test_integration_pipeline.py::test_offset_selection_
   psite_vs_asite_shifts_by_three_nt` and related).
4. mt-RPF length range 15-45 nt enforced at trim, align, and the
   rpf `-rpf MIN MAX` filter.
5. ND5 / ND6 antisense overlap resolved by construction on Path A
   (transcriptome reference, one FASTA record per mt-mRNA) plus
   `--library-strandedness forward` -> bowtie2 `--norc` at alignment
   time.
6. Codon-occupancy signal preservation: safe-by-default dedup
   (Approach D; `mark-duplicates` gated behind
   `--i-understand-mark-duplicates-destroys-mt-ribo-seq-signal`).
7. NUMT cross-talk suppression: `--mapq 10` default at align and on
   BAM inputs to rpf.
8. Ribo-seq / RNA-seq reference consistency: SHA-256 gate in rnaseq
   against rpf-recorded `reference_checksum`. Hard fail on mismatch.

## Outstanding items before v0.3.0 final tag

1. **Push the branch and open the draft PR.** Blocked on sandbox
   credentials; see status section above.
2. **Run the `requires_tools` integration test on a bioconda env.**
   Expected to pass on the maintainer's workstation; documented in
   `docs/environment/environment.yml`.
3. **Execute the TACO1-KO polyproline regression.** Plan +
   acceptance criteria in `docs/validation/taco1_ko_regression.md`;
   awaiting dataset delivery.
4. **Optional:** run the CI matrix across Python 3.10 / 3.11 / 3.12.
   Local run only exercised 3.11; CI config updates (if needed) are
   outside the v0.3.0 scope and can be deferred to v0.3.1.

## Recommended tag + PR commands

When push access is restored:

```bash
git push -u origin refactor/v0.3.0-preprocessing
git tag -a v0.3.0-rc1 -m "MitoRiboPy v0.3.0 release candidate 1" \
        refactor/v0.3.0-preprocessing
git push origin v0.3.0-rc1
gh pr create \
  --title "MitoRiboPy v0.3.0 refactor" \
  --body "See CHANGELOG.md [0.3.0] section" \
  --draft \
  --base main \
  --head refactor/v0.3.0-preprocessing
```

## Revert plan

Every commit on the branch is independent and revertable:
`git revert <sha>` cleanly undoes the change. Phase-3 commits are
particularly isolated because each one adds a single `src/mitoribopy/
align/*.py` file plus its tests. Phase-5 rnaseq and Phase-6 `all` are
single-commit stages and can be reverted wholesale without affecting
the rest of the refactor.
