# MitoRiboPy smoke fixture

Tiny end-to-end test that runs the full `mitoribopy all` pipeline in
under 30 seconds on a fresh install. Three synthetic mt-mRNAs, two
"samples" (one WT, one KO), no UMIs, no contaminants — designed to
exercise every stage (trim → align → BED → offsets → translation
profile → coverage → metagene Fourier QC) without needing real-world
input data.

This is **not** a biological-correctness test. The TACO1-KO regression
in [`docs/validation/taco1_ko_regression.md`](../../docs/validation/taco1_ko_regression.md)
is the biology gate. The smoke fixture exists so a user installing
MitoRiboPy can confirm in one command that the package is wired up
end-to-end on their machine.

---

## Prerequisites

The smoke test shells out to the same external tools the real `align`
stage uses. They must be on `$PATH`:

* `cutadapt` (>= 4.0)
* `bowtie2` and `bowtie2-build` (>= 2.4)
* `samtools` (>= 1.15) — optional, only used by inspection commands

The `[fastq]` extra is **not** required: the smoke fixture exercises the
align + rpf path only (no rnaseq stage).

---

## Run it

```bash
cd examples/smoke
python generate_smoke_fastqs.py    # writes *.fastq.gz + builds bowtie2 index
mitoribopy all --config pipeline_config.smoke.yaml --output results/
```

Expected wall-clock: ~10-30 seconds on a 2024 laptop.

---

## What's checked

After `mitoribopy all` exits 0, every file listed in
[`expected_outputs.txt`](expected_outputs.txt) must exist under
`results/` and be non-empty. The pytest harness
([`tests/test_smoke_fixture.py`](../../tests/test_smoke_fixture.py))
re-runs the same flow under `pytest -m smoke` when the external tools
are available.

---

## Files

| File | Purpose |
|---|---|
| [`generate_smoke_fastqs.py`](generate_smoke_fastqs.py) | Synthesises FASTQs + reference + bowtie2 index. Deterministic (seed 42). |
| [`human_mt_tiny.fa`](human_mt_tiny.fa) | Three synthetic mt-mRNA-like contigs (`MT-CO1`, `MT-ND1`, `MT-ATP6`), each ~600 nt. |
| [`samples.tsv`](samples.tsv) | Two-row sample sheet (`WT_smoke_1`, `KO_smoke_1`), both `assay=ribo`. |
| [`samples.condition_map.tsv`](samples.condition_map.tsv) | Maps both samples to their condition labels. |
| [`pipeline_config.smoke.yaml`](pipeline_config.smoke.yaml) | Config drives `mitoribopy all`. No rnaseq section. |
| [`expected_outputs.txt`](expected_outputs.txt) | List of files that must exist after a successful run. |
| `results/` | Run output (gitignored). |

---

## Why no UMIs?

The smoke fixture sets `umi_position: none` so it does not depend on
`umi_tools`. The UMI / dedup pipeline has its own dedicated tests
under `tests/test_align_dedup.py` and `tests/test_align_dual_umi.py`.
