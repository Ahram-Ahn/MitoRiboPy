# Custom Reference Templates

Two small templates for `--strain custom` runs. The full schema and a
worked end-to-end example (mouse mitochondria) live in the
[Custom organisms](../../README.md#custom-organisms) section of the
top-level README — start there. This directory exists so you have a
syntactically valid file to copy and edit.

## Files

- [`annotation_template.csv`](annotation_template.csv) — example
  per-transcript annotation table. Required columns: `transcript`,
  `l_tr`, `l_utr5`, `l_utr3`. Optional columns: `l_cds` (computed
  if absent), `sequence_name`, `sequence_aliases`, `display_name`.
  The CDS bounds (`start_codon`, `stop_codon`) are computed by the
  loader and do NOT need to be supplied.
- [`codon_tables_template.json`](codon_tables_template.json) — example
  named codon-table JSON. Use this **only** when your organism's
  genetic code is not in the bundled NCBI Genetic Codes list. Most
  users can skip the file entirely and pass an NCBI name through
  `--codon_table_name` (e.g. `vertebrate_mitochondrial`,
  `mold_mitochondrial`, `invertebrate_mitochondrial`); see the
  Custom organisms section for the full picker.

## Quick recipe

```bash
$ cp examples/custom_reference/annotation_template.csv my_annotation.csv
$ # edit my_annotation.csv to match your transcripts

$ mitoribopy rpf \
    --strain custom \
    --fasta my_reference.fa \
    --annotation_file my_annotation.csv \
    --codon_table_name vertebrate_mitochondrial \
    --rpf 28 34 \
    --directory bed/ \
    --output results/
```

If you really do need a custom genetic code, drop the JSON in too:

```bash
$ cp examples/custom_reference/codon_tables_template.json my_codon_tables.json
$ # edit my_codon_tables.json
$ mitoribopy rpf ... \
    --codon_tables_file my_codon_tables.json \
    --codon_table_name your_table_name
```
