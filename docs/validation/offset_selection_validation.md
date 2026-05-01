# Offset selection: validation evidence

How MitoRiboPy chooses the per-sample 5'/3' P-site offset, and why the
choice is defensible.

## Algorithm overview

* For each read length in the configured RPF window (default 29–34 nt
  for human monosomes; configurable via `rpf.rpf` and
  `rpf.footprint_class`), the rpf stage builds a histogram of the
  distance from the 5' end (or 3' end, per `offset_type`) of each
  read to the nearest stop codon (or start codon, per `align`).
* The selected offset is the value of the dominant peak inside the
  configured `[min_5_offset, max_5_offset]` (or `_3_`) window.
* When `offset_mode: per_sample` (default), every sample's offset is
  picked independently; with `offset_mode: combined` the offsets are
  picked from the pooled histogram across all samples.

## Confidence reporting

For every (sample, read length) pair the rpf stage emits an offset
table whose schema includes:

* `selected_offset` — the chosen value;
* `peak_count` — read count at `selected_offset`;
* `second_peak_count` — read count at the second-best offset;
* `peak_margin` — `peak_count - second_peak_count`;
* `confidence` — `high` / `medium` / `low`, derived from the margin
  and the entropy of the histogram.

A `confidence: low` row is the explicit signal that the per-sample
fallback to `offset_mode: combined` should be considered.

## Known-answer validation

The synthetic-mini fixture
([`synthetic_mini.md`](synthetic_mini.md)) constructs reads whose
5' end sits exactly 14 nt 5' of the stop codon. The rpf stage's
arithmetic is:

```
5' offset = stop_codon - read_start - 2
```

so the correct selected offset is `12`, and the test asserts:

* `Most Enriched 5' Offset == 12` for every per-sample selection AND
  for the combined fit;
* `confidence_5 == "high"` for the combined fit (the constructed peak
  is sharp by design, with 60% of mass at one offset).

These invariants are enforced by `tests/test_offset_per_sample.py`
and `tests/test_synthetic_mini.py`; failures there indicate that the
selector's arithmetic or its confidence labelling has drifted.

## Public-dataset spot-check

The full reanalysis recipe
([`public_dataset_reanalysis.md`](public_dataset_reanalysis.md))
includes the offset table in the artefacts that should be inspected
before publication. Typical defensible values for human mt-RPF
monosome libraries:

| read length | typical 5' offset | confidence |
|---|---|---|
| 29 nt | 12 | high |
| 30 nt | 12 | high |
| 31 nt | 12 | high |
| 32 nt | 13 | high–medium |
| 33 nt | 13 | medium |

Values outside this range, or `confidence: low` on multiple lengths,
are signals that the library was probably not size-selected to a
canonical monosome window — re-check the kit / size-selection step
before reporting.
