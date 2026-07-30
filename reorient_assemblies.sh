#!/usr/bin/env bash
# =============================================================================
# reorient_assemblies.sh
# Relance dnaapler sur les assembly.fasta existants (sans re-faire Flye/Unicycler).
# A exécuter depuis le dossier phage_nanopore/ sur le cluster.
# Supprime ensuite les polished.fasta pour forcer la re-polissage en aval.
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
activate_env
require_tool dnaapler

THREADS="${THREADS:-8}"

echo "========================================================"
echo "  Reorientation dnaapler — tous les barcodes"
echo "========================================================"

while IFS=$'\t' read -r barcode sample_name _rest; do
    [[ "$barcode" == "barcode" ]] && continue   # skip header
    [[ -z "$barcode" ]] && continue

    ASM_DIR="${RESULTS_DIR}/02_assembly/${barcode}_${sample_name}"
    ASM="${ASM_DIR}/assembly.fasta"
    DNAAP_OUT="${ASM_DIR}/dnaapler"

    if [[ ! -f "$ASM" ]]; then
        echo "[SKIP] $barcode — assembly.fasta introuvable"
        continue
    fi

    echo ""
    echo "── $barcode ($sample_name) ──"
    echo "   Input : $ASM"

    rm -rf "$DNAAP_OUT"
    dnaapler all \
        -i "$ASM" \
        -o "$DNAAP_OUT" \
        -p assembly \
        -t "$THREADS" \
        2>&1

    REORIENTED="${DNAAP_OUT}/assembly_reoriented.fasta"
    if [[ -f "$REORIENTED" ]]; then
        cp "$REORIENTED" "$ASM"
        echo "   ✓ Réorienté → $ASM"
    else
        echo "   ⚠ dnaapler n'a pas produit de fichier réorienté — orientation inchangée"
    fi

    # Supprimer polished.fasta pour forcer le re-polissage en aval
    POLISHED="${RESULTS_DIR}/03_polished/${barcode}_${sample_name}/polished.fasta"
    if [[ -f "$POLISHED" ]]; then
        rm -f "$POLISHED"
        echo "   → polished.fasta supprimé (re-polissage requis)"
    fi

done < <(grep -v '^#' "${RESULTS_DIR}/sample_sheet.tsv" | tail -n +2 | cut -f1,2)

echo ""
echo "========================================================"
echo "  Terminé. Lance ensuite :"
echo "  bash submit_all.sh --resume --steps 03,04,06,07,08,09,09b,09c,10"
echo "========================================================"
