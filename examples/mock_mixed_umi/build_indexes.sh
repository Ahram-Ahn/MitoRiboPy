#!/usr/bin/env bash
# Build the bowtie2 indexes the mock_mixed_umi smoke test needs.
# Run this once after `python make_mock.py`.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INDEXES="${HERE}/indexes"
mkdir -p "${INDEXES}"

bowtie2-build --quiet "${HERE}/human-mt-mRNA.mini.fasta" "${INDEXES}/mt"
bowtie2-build --quiet "${HERE}/human_rrna.mini.fa"      "${INDEXES}/rrna"

echo "Built bowtie2 indexes under ${INDEXES}/"
ls "${INDEXES}"
