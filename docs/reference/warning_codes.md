# Warning / error codes

Every record written to `<run_root>/warnings.tsv` carries a stable
`code` value. The codes are owned by
[`src/mitoribopy/io/warning_codes.py`](../../src/mitoribopy/io/warning_codes.py);
this page is a human-readable mirror.

The naming convention is:

| Prefix | Severity | Behaviour |
| ------ | -------- | --------- |
| `E_`   | error    | Promoted to a non-zero exit under `mitoribopy all --strict`; otherwise the run finishes and the manifest carries the error. |
| `W_`   | warn     | Surfaced in `SUMMARY.md` and `warnings.tsv`. The publication recipe expects `--strict` to gate on these. |

Use the codes in CI / paper-figure scripts that grep the manifest:

```bash
jq -r '.warnings[] | "\(.severity)\t\(.code)\t\(.message)"' \
  results/run_manifest.json \
  | sort -u
```

## Config / sample sheet

| Code | Severity | Stage | Summary | Remediation |
| ---- | -------- | ----- | ------- | ----------- |
| `E_CONFIG_UNKNOWN_KEY` | error | all | Top-level YAML key not in the canonical schema. | Remove the key, or run `mitoribopy migrate-config` to rewrite legacy spellings into the canonical form. |
| `E_CONFIG_DEPRECATED_KEY_STRICT` | error | all | A deprecated key was rewritten by the canonicaliser, but `--strict` was set. | Update the YAML to use the canonical key name; re-run with `--strict`. |
| `E_SAMPLE_SHEET_DUPLICATE_ID` | error | all | Two rows in the unified sample sheet share a sample_id. | Make `sample_id` unique across every row. |
| `E_SAMPLE_SHEET_UNKNOWN_COLUMN` | error | all | Sample sheet has a column not in the documented schema. | Drop the column, or move it to a non-colliding name. |

## Align preprocessing

| Code | Severity | Stage | Summary | Remediation |
| ---- | -------- | ----- | ------- | ----------- |
| `E_FASTQ_MATE_PAIRING_FAILED` | error | align | Could not pair R1/R2 FASTQs from filename. | Rename to `<sample>_R1.fq.gz` / `<sample>_R2.fq.gz`, or list pairs explicitly. |
| `E_ADAPTER_AMBIGUOUS` | error | align | Adapter auto-detection saw two candidate kits above the threshold. | Pin the 3' adapter explicitly with `--adapter <SEQ>` (CLI) or `adapter: <SEQ>` (YAML / sample sheet). |
| `E_UMI_DECLARED_BUT_NOT_FOUND` | error | align | Config declared a UMI but cutadapt extracted none. | Verify `umi_length` / `umi_position`; pre-trimmed reads need `--pretrimmed` (or `pretrimmed: true` in YAML). |
| `W_UMI_SHORT_COLLISION_RISK` | warn | align | UMI shorter than 8 nt has a high collision rate. | Treat duplicate-fraction estimates as conservative; prefer ≥ 8 nt UMIs. |
| `W_UMI_HIGH_DUPLICATE_FRACTION` | warn | align | Coordinate+UMI dedup removed > 80 % of reads. | Re-check library complexity before biological inference. |
| `W_UMI_DEDUP_SKIPPED_WITH_UMI_PRESENT` | warn | align | UMI present but `dedup_strategy: skip` set. | Switch to `dedup_strategy: umi_coordinate` (or `auto`). |
| `W_UMI_METHOD_DIRECTIONAL_ON_SHORT_UMI` | warn | align | `directional` clustering on a < 8 nt UMI over-collapses. | Use `umi_dedup_method: unique` for short UMIs. |

## Rpf / periodicity / codon correlation

| Code | Severity | Stage | Summary | Remediation |
| ---- | -------- | ----- | ------- | ----------- |
| `W_CODON_RAW_COUNT_PRIMARY` | warn | rpf | Codon-correlation rendered with `--cor-metric raw_count`. | Re-run with `--cor-metric log2_density_rpm` for publication. |
| `W_PERIODICITY_LOW_DEPTH` | warn | rpf | Sample has < `min_reads_per_length` assigned P-sites. | Treat as `low_depth`; do not use it to validate the offset choice. |
| `W_PERIODICITY_WEAK` | warn | rpf | 3-nt spectral ratio below the warn threshold. | Inspect `fourier_period3_score_combined.tsv` (`snr_call` column) and the per-(sample, read_length) Fourier overlay PNGs under `rpf/qc/fourier_spectrum/`; consider tightening the read-length window. |

## Rnaseq / TE

| Code | Severity | Stage | Summary | Remediation |
| ---- | -------- | ----- | ------- | ----------- |
| `W_RNASEQ_FROM_FASTQ_EXPLORATORY` | warn | rnaseq | Mode = `from_fastq` (mt-mRNA-only pyDESeq2). | For publication, run DESeq2 / Xtail / Anota2Seq externally; re-run with `--rnaseq-mode de_table`. |
| `W_RNASEQ_PSEUDO_REPLICATE` | warn | rnaseq | n=1 condition split into FASTQ-record-parity pseudo-replicates. | Use only for demo / smoke tests. |
| `E_REFERENCE_CHECKSUM_MISMATCH` | error | rnaseq | rpf and rnaseq references disagree (different SHA256). | Rebuild the reference once and re-run both stages. |
| `E_TEMPLATE_INVALID` | error | all | A shipped template failed `validate-config --strict`. | File a MitoRiboPy issue. |

## Adding new codes

1. Add a `WarningCode(...)` entry to `src/mitoribopy/io/warning_codes.py`.
2. Add a row above (this file is hand-curated; the schema is intentionally short).
3. Call `mitoribopy.io.warnings_log.record(code="<NEW_CODE>", ...)`
   somewhere in the codebase.
4. The unit tests under `tests/test_warnings_log.py` will refuse a
   record whose `code` is not in the registry.
