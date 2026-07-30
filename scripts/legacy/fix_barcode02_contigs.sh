#!/usr/bin/env bash
# Renomme les contigs dupliqués dans polished.fasta de barcode02
# puis relance pharokka (step 04) sur ce barcode uniquement
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

POLISHED="${RESULTS_DIR}/03_polished/barcode02_PAK_P4/polished.fasta"
FIXED="${RESULTS_DIR}/03_polished/barcode02_PAK_P4/polished_fixed.fasta"

[[ -f "${POLISHED}" ]] || { echo "ERROR: ${POLISHED} not found"; exit 1; }

echo "Input contigs:"
grep '^>' "${POLISHED}" | head -20

# Rename: add _1, _2 etc to duplicates, and just number all headers sequentially
python3 - "${POLISHED}" "${FIXED}" << 'PYEOF'
import sys
from collections import defaultdict

infile, outfile = sys.argv[1], sys.argv[2]
seen = defaultdict(int)
with open(infile) as fi, open(outfile, 'w') as fo:
    for line in fi:
        if line.startswith('>'):
            name = line.strip().lstrip('>')
            seen[name] += 1
            if seen[name] > 1:
                new_name = f"{name}_dup{seen[name]}"
            else:
                new_name = name
            fo.write(f">{new_name}\n")
        else:
            fo.write(line)
PYEOF

echo ""
echo "Fixed contigs:"
grep '^>' "${FIXED}"

# Backup original and replace
cp "${POLISHED}" "${POLISHED}.bak"
cp "${FIXED}"    "${POLISHED}"
echo ""
echo "Replaced polished.fasta (backup at polished.fasta.bak)"
echo ""
echo "Now delete the failed annotation and rerun step 04:"
echo "  rm -rf ${RESULTS_DIR}/04_annotated/barcode02_PAK_P4"
echo "  bash submit_all.sh --step 04"
