#!/usr/bin/env bash
# =============================================================================
# clean_for_reassembly.sh
# -----------------------------------------------------------------------------
# Removes the DOWNSTREAM results (steps 03 -> 10) that were produced from the
# old contaminated assembly, so the pipeline regenerates them from the new
# host-depleted assembly instead of skipping (every step self-skips when its
# output already exists unless FORCE=true).
#
# What it DELETES (per selected barcode + global aggregates):
#   03_polished/<label>           04_annotated/<label>
#   07_modbase/<label>            09_breseq/<label>
#   09c_snippy/<label>_vs_*       09d_breseq_internal/<label>_vs_*
#   06_summary/  08_methyl/  09_breseq/gdtools/  10_final_report/   (global)
#
# What it KEEPS (inputs / raw data — never touched):
#   00_merged/   01_qc/   02_assembly/   sample_sheet.tsv   refs/
#   (02_assembly is left alone: your new assembly.fasta lives there; the
#    assemble script already backed up the contaminated one.)
#
# SAFE BY DEFAULT: prints what it would remove and does nothing. Add --yes to
# actually delete.
#
# USAGE
#   bash clean_for_reassembly.sh                 # dry-run, all barcodes
#   bash clean_for_reassembly.sh --yes           # really delete, all barcodes
#   bash clean_for_reassembly.sh --barcodes barcode03,barcode04,barcode05 --yes
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

DO_IT=false
BARCODES_ARG=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --yes)       DO_IT=true; shift ;;
        --barcodes)  BARCODES_ARG="$2"; shift 2 ;;
        -h|--help)   sed -n '2,40p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *)           die "Unknown argument: $1 (use --help)" ;;
    esac
done

if [[ -n "$BARCODES_ARG" ]]; then
    IFS=',' read -r -a BARCODES <<< "$BARCODES_ARG"
else
    mapfile -t BARCODES < <(list_barcodes)
fi

echo "RESULTS_DIR: ${RESULTS_DIR}"
echo "Mode       : $([[ "$DO_IT" == true ]] && echo 'DELETE (--yes)' || echo 'dry-run (add --yes to delete)')"
echo ""

TARGETS=()
for BARCODE in "${BARCODES[@]}"; do
    SAMPLE_NAME=$(get_sample_field "$BARCODE" "sample_name")
    LABEL="${BARCODE}_${SAMPLE_NAME}"
    TARGETS+=( "${RESULTS_DIR}/03_polished/${LABEL}" )
    TARGETS+=( "${RESULTS_DIR}/04_annotated/${LABEL}" )
    TARGETS+=( "${RESULTS_DIR}/07_modbase/${LABEL}" )
    TARGETS+=( "${RESULTS_DIR}/09_breseq/${LABEL}" )
    # comparison dirs are named <label>_vs_<ref> — glob them
    while IFS= read -r d; do TARGETS+=( "$d" ); done < <(
        find "${RESULTS_DIR}/09c_snippy" -maxdepth 1 -type d -name "${LABEL}_vs_*" 2>/dev/null
        find "${RESULTS_DIR}/09d_breseq_internal" -maxdepth 1 -type d -name "${LABEL}_vs_*" 2>/dev/null
    )
done
# global aggregates (mix all samples → always regenerate)
TARGETS+=( "${RESULTS_DIR}/06_summary" )
TARGETS+=( "${RESULTS_DIR}/08_methyl" )
TARGETS+=( "${RESULTS_DIR}/09_breseq/gdtools" )
TARGETS+=( "${RESULTS_DIR}/10_final_report" )

removed=0
for t in "${TARGETS[@]}"; do
    if [[ -e "$t" ]]; then
        sz=$(du -sh "$t" 2>/dev/null | cut -f1)
        if [[ "$DO_IT" == true ]]; then
            rm -rf "$t" && echo "  removed  ${t}  (${sz})"
        else
            echo "  would remove  ${t}  (${sz})"
        fi
        removed=$((removed+1))
    fi
done

echo ""
if [[ "$removed" -eq 0 ]]; then
    echo "Nothing to remove (already clean)."
elif [[ "$DO_IT" == true ]]; then
    echo "Done — ${removed} path(s) removed. Now relaunch:  bash submit_all.sh"
else
    echo "${removed} path(s) would be removed. Re-run with --yes to delete."
fi
