# UMI deduplication: validation and conservative defaults

## What gets called a UMI

For each sample, the resolved UMI length comes from one of four
sources, recorded in `align/kit_resolution.tsv`'s `umi_source` column:

| `umi_source` value | Meaning |
|---|---|
| `declared` | The sample sheet supplied `umi_length` and `umi_position` for this sample. **Recommended for publication.** |
| `preset_default` | The chosen kit preset's own `umi_length` was used (e.g. `illumina_truseq_umi` declares an 8 nt 5' UMI). |
| `inferred` | Per-sample R1 entropy detector landed the value (rnaseq from-FASTQ path only; align does not infer). Triggers a `UMI_INFERRED_NO_DECLARATION` warning into `warnings.tsv`. |
| `none` | `umi_length == 0`; dedup is skipped. |

## Why we snap inferred UMIs DOWN

`mitoribopy.rnaseq.umi_detect` measures Shannon entropy across the
first / last 16 positions of R1 and identifies the longest run of
high-entropy positions at either end. The result is then **snapped
down** to the nearest of `{0, 4, 6, 8, 10, 12}`:

* a wrong-too-long UMI silently corrupts dedup (real bases get
  trimmed, distinct fragments collapse);
* a wrong-too-short UMI just means some duplicate reads survive,
  which is recoverable downstream.

The conservative snap is therefore the safe direction, and is the
basis for the package's claim that an inferred UMI is *advisory*.

## Publication-mode discipline

For publication-grade runs, declare every sample's UMI explicitly in
the sample sheet:

```
sample_id    assay    condition    fastq_1            adapter                                    umi_length    umi_position    dedup_strategy
WT_Ribo_1    ribo     WT           ribo/WT.fq.gz      AGATCGGAAGAGCACACGTCTGAACTCCAGTCA          8             5p              umi-tools
KO_Ribo_1    ribo     KO           ribo/KO.fq.gz      AGATCGGAAGAGCACACGTCTGAACTCCAGTCA          8             5p              umi-tools
```

> The `kit_preset` column was removed in v0.7.1; declare the 3'
> adapter sequence explicitly via `adapter` (or set `pretrimmed: true`).
> The detector still names the matched adapter family in
> `kit_resolution.tsv` (`detected_kit` / `applied_kit`) for provenance.

When a row lacks `umi_length` and the rnaseq detector subsequently
infers one, the run still completes (the inferred value is applied)
but a structured warning lands in `warnings.tsv`:

```
component    severity    sample        code                              message
RNASEQ       warn        WT_RNA_1      UMI_INFERRED_NO_DECLARATION       sample 'WT_RNA_1': UMI length 8 was inferred from R1 entropy; the sample sheet did not declare umi_length / umi_position. ...
```

A run with any `UMI_INFERRED_NO_DECLARATION` warnings should not be
labelled "publication-grade" — fix the sheet and re-run with
`--resume` (the hash guard will allow the resume because the sheet
edit is the intended change, but `--strict` mode in
`validate-config` will catch the legacy case).

## Tests

* `tests/test_align_sample_resolve.py::*umi*` — exercises every
  `umi_source` outcome of `_resolve_one`.
* `tests/test_rnaseq_umi_detect.py` — checks the entropy detector's
  snap-DOWN behaviour and `length=0` fallback.
* `tests/test_warnings_log.py::test_mitoribopy_all_writes_warnings_tsv_and_manifest`
  — end-to-end check that the inference warning lands in
  `warnings.tsv` and in the manifest.
