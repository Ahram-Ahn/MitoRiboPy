# Custom organisms, strain presets, and footprint classes

To analyse an organism other than human or yeast, run with
`--strain custom` and supply three things: a **reference FASTA** of
mt-mRNA transcripts, a **per-transcript annotation CSV**, and a
**codon-table choice**. Start codons can be left at the `[ATG]`
default or overridden with `--start_codons`.

For the input-side overview (what each subcommand needs on disk),
see [`docs/inputs.md`](inputs.md). For the output-side schema (what
the package emits and where), see
[`docs/reference/output_schema.md`](reference/output_schema.md).

---

## Strain presets

The strain preset selects the organism's mitochondrial annotation
and codon table. Two organisms ship complete reference data;
everything else uses `custom` and supplies its own files.

| Value | Organism | Codon table | Ships annotation? | Ships `-rpf` default? |
|---|---|---|:-:|:-:|
| `h.sapiens` (default) | *Homo sapiens* mt | `vertebrate_mitochondrial` (NCBI #2) | ✓ | ✓ |
| `s.cerevisiae` | *Saccharomyces cerevisiae* mt | `yeast_mitochondrial` (NCBI #3) | ✓ | ✓ |
| `custom` | Any other organism | user-supplied via `--codon_table_name` (built-in NCBI list) or `--codon_tables_file` | ✗ — pass `--annotation_file` | ✗ — pass `-rpf MIN MAX` |

`h` and `y` are also accepted as short synonyms for `h.sapiens` and
`s.cerevisiae`. Publication-grade configs should always use the
canonical long names; `validate-config --strict` rewrites the
shorthand under [`mitoribopy.config.migrate`](../src/mitoribopy/config/migrate.py).

---

## Footprint class (`--footprint_class`)

Pair `-s` with `--footprint_class` to pick sensible RPF and
unfiltered-length defaults. An explicit `-rpf MIN MAX` or
`--unfiltered_read_length_range MIN MAX` always wins over the
footprint-class default. Built-in defaults exist for `h.sapiens`
and `s.cerevisiae`; for `--strain custom` you must also pass
`-rpf`.

| Value | RPF window default | `--unfiltered_read_length_range` default | Use for |
|---|---|---|---|
| `short` | h.sapiens / s.cerevisiae: 16–24 | 10–30 | Truncated RNase products. Sit just below the canonical monosome window; useful for context-dependent pausing and as a QC indicator of digest aggressiveness. |
| `monosome` (default) | h.sapiens: 28–34, s.cerevisiae: 37–41 | 15–50 | Single-ribosome footprints. The standard mt-Ribo-seq class. |
| `disome` | h.sapiens: 50–70, s.cerevisiae: 60–90 | 40–100 | Collided-ribosome footprints. eIF5A-depletion, queueing, ribosome-stalling studies. |
| `custom` | user must pass `-rpf` | unchanged | Any non-standard footprint class. |

---

## Picking a codon table

MitoRiboPy bundles every NCBI Genetic Code as a named table. The
names match NCBI's organism-group labels; pick the one for your
organism's mitochondrial code (or nuclear code for ciliates /
*Candida* / similar).

| `--codon_table_name` | NCBI # | Use for |
|---|:---:|---|
| `standard` | 1 | Most plant mitochondria and plastids; many fungal nuclear genomes |
| `vertebrate_mitochondrial` | 2 | Mouse, rat, zebrafish, Xenopus, chicken — any non-human vertebrate mt |
| `yeast_mitochondrial` | 3 | *S. cerevisiae* and close relatives (matches `s.cerevisiae` preset) |
| `mold_mitochondrial` | 4 | *Neurospora*, *Aspergillus*, *Trichoderma*, mycoplasmas |
| `invertebrate_mitochondrial` | 5 | *Drosophila*, mosquito, *C. elegans* and other invertebrates |
| `echinoderm_mitochondrial` | 9 | Sea urchin, starfish |
| `ascidian_mitochondrial` | 13 | Sea squirts (tunicates) |
| `alternative_flatworm_mitochondrial` | 14 | Flatworms |
| `trematode_mitochondrial` | 21 | Trematodes |

Run `mitoribopy rpf --help` to see the full list of all 27 bundled
tables (NCBI #1–34, with a few historical numbers omitted by NCBI).
Sources: [NCBI Taxonomy Genetic Codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
— the JSON tables under
[`src/mitoribopy/data/codon_tables.json`](../src/mitoribopy/data/codon_tables.json)
are a faithful transcription of those NCBI tables.

If the table you need isn't in the built-in list (e.g. an unusual
reassignment in a non-model organism), pass
`--codon_tables_file your_tables.json --codon_table_name your_name`
to load your own JSON. Each entry maps the 64 codons to single-letter
amino-acid codes (use `*` for stop):

```json
{
  "your_organism_mitochondrial": {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TAA": "*", "TAG": "*", "TGA": "W",
    "...": "..."
  }
}
```

---

## Annotation CSV format

`--annotation_file` takes a CSV with one row per mt-mRNA transcript.
The loader requires four columns and accepts four optional columns;
CDS bounds (`start_codon`, `stop_codon`, `l_cds`) are computed from
the UTR lengths and do NOT need to be supplied.

| Column | Required? | Type | Meaning |
|---|:---:|---|---|
| `transcript` | yes | string | Display name used in plots and CSVs (e.g. `ND1`, `COX1`). Must be unique per row. |
| `l_tr` | yes | int | Total transcript length in nt (`l_utr5 + l_cds + l_utr3`). |
| `l_utr5` | yes | int | Length of the 5' UTR in nt (0 if none). |
| `l_utr3` | yes | int | Length of the 3' UTR in nt (0 if none). |
| `l_cds` | optional | int | Length of the CDS in nt. Computed as `l_tr − l_utr5 − l_utr3` when omitted. |
| `sequence_name` | optional | string | The exact FASTA header / BED chrom field for this transcript. Defaults to the `transcript` column when blank, but **must match your FASTA header** if those names differ. |
| `display_name` | optional | string | Human-readable label used in plot titles. Defaults to `transcript` when blank. |
| `sequence_aliases` | optional | string | Semicolon-separated list of legacy IDs that should also map to this transcript (e.g. `ATP86;ATP8_6` if older BEDs use those names). Leave blank when there are none. |

Minimal example for a hypothetical 3-transcript mouse mt-mRNA reference:

```csv
transcript,sequence_name,l_tr,l_utr5,l_utr3
ND1,mouse_ND1,957,0,0
COX1,mouse_COX1,1545,0,0
ATP6,mouse_ATP6,681,0,0
```

The `sequence_name` field MUST match the FASTA header you pass to
`-f` AND the BED `chrom` field produced by `mitoribopy align` (so
the bowtie2 indexes have to be built from the same FASTA).
[`examples/custom_reference/annotation_template.csv`](../examples/custom_reference/annotation_template.csv)
is a complete worked example.

---

## Putting it together: a mouse run

```bash
# Build the bowtie2 indexes from your mouse mt-mRNA FASTA + an rRNA decoy.
bowtie2-build mouse-mt-mRNA.fasta indexes/mt
bowtie2-build mouse_rrna.fa       indexes/rrna

# Run the full pipeline.
mitoribopy all --config mouse_pipeline.yaml --output mouse_results/
```

```yaml
# mouse_pipeline.yaml
align:
  # adapter auto-detection runs by default; add `adapter: <SEQ>` only
  # when detection cannot identify the library, or `pretrimmed: true`
  # for already-trimmed FASTQs.
  fastq: input_data/seq
  contam_index: indexes/rrna
  mt_index: indexes/mt

rpf:
  strain: custom
  fasta: mouse-mt-mRNA.fasta
  annotation_file: mouse_annotation.csv      # CSV described above
  codon_table_name: vertebrate_mitochondrial # NCBI #2
  rpf: [28, 34]                              # mouse mt-monosome window
  align: stop
  offset_type: "5"
  offset_site: p
```

For other organisms the recipe is identical — only the codon table
name and the annotation CSV change.

---

## Bicistronic transcript pairs

The two overlapping mt-mRNA pairs (`ATP8/ATP6` and `ND4L/ND4`) are
kept as paired display names in plot titles. Choose which member's
coordinates seed the merged sequence with
`--atp8_atp6_baseline ATP8|ATP6` and
`--nd4l_nd4_baseline ND4L|ND4` (defaults: `ATP6` and `ND4`). Legacy
FASTA / BED identifiers `ATP86` and `ND4L4` are also recognised
through the built-in alias map. The metagene Fourier QC bundle
treats these pairs as their own `gene_set` ("ATP86" and "ND4L4"
panels) — see
[`docs/reference/periodicity.md`](reference/periodicity.md).

---

## Built-in references

MitoRiboPy ships packaged reference data for two organisms:

- *Homo sapiens* mt-translation using the
  `vertebrate_mitochondrial` codon table (`-s h.sapiens`)
- *Saccharomyces cerevisiae* mt-translation using the
  `yeast_mitochondrial` codon table (`-s s.cerevisiae`)

Built-in annotation tables are stored as CSV and the bundled NCBI
Genetic Codes as JSON under
[`src/mitoribopy/data/`](../src/mitoribopy/data/). 27 codon tables
are available out of the box; the picker for non-human, non-yeast
organisms lives above.
