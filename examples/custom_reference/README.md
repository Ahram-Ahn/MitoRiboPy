# Custom Reference Templates

This directory contains small, valid templates for running `mitoribopy` with `--strain custom`.

Files:
- `annotation_template.csv`: example CDS annotation table
- `codon_tables_template.json`: example named codon-table JSON

How to use:
1. Copy the files and edit them for your organism.
2. Point `--annotation_file` to the edited CSV.
3. Point `--codon_tables_file` to the edited JSON.
4. Choose the table with `--codon_table_name`.

Notes:
- `transcript` is the logical CDS name used in codon and frame summaries.
- `sequence_name` is the FASTA/BED sequence ID that the CDS maps onto.
- `sequence_aliases` can list alternate FASTA/BED names separated by semicolons.
- `display_name` controls plot titles and can be shared by bicistronic CDS rows.
