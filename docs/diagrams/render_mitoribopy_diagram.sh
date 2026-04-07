#!/usr/bin/env bash
set -euo pipefail

# Render Mermaid diagram source into SVG using the public Mermaid CLI.
# Requires Node.js + npx.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
INPUT_FILE="${ROOT_DIR}/docs/diagrams/mitoribopy_package_flow.mmd"
OUTPUT_FILE="${1:-${ROOT_DIR}/docs/diagrams/mitoribopy_package_flow.svg}"

if [[ ! -f "${INPUT_FILE}" ]]; then
  echo "[diagram] Missing input: ${INPUT_FILE}" >&2
  echo "[diagram] Run: python docs/diagrams/generate_mitoribopy_diagram.py" >&2
  exit 1
fi

if ! command -v npx >/dev/null 2>&1; then
  echo "[diagram] npx was not found." >&2
  echo "[diagram] Install Node.js, then run this script again." >&2
  exit 1
fi

npx -y @mermaid-js/mermaid-cli -i "${INPUT_FILE}" -o "${OUTPUT_FILE}"
echo "[diagram] Rendered SVG: ${OUTPUT_FILE}"
