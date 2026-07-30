#!/usr/bin/env bash
# =============================================================================
# refresh_sample_sheet.sh
# =============================================================================
# Scans $ILLUMINA_DIR for R1/R2 FASTQ files and updates (or generates) the
# illumina_r1 / illumina_r2 columns in the sample sheet.
#
# Run once when Illumina reads are placed in a new location (new run dir,
# data transfer, etc.).  All other columns (barcode, sample_name, mode) are
# preserved exactly.
#
# Usage:
#   bash refresh_sample_sheet.sh              # update in place (backup kept)
#   bash refresh_sample_sheet.sh --dry-run    # print what would be written
#
# Reads are found by searching $ILLUMINA_DIR/barcodeXX_*/ for files matching
# *R1*.fastq.gz, *R1*.fq.gz (and R2 equivalents).  Works regardless of the
# per-sample naming convention (OM104_S7_L001, merged_R1, CS049_nxq, etc.).
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

SAMPLE_SHEET="${RESULTS_DIR}/sample_sheet.tsv"
[[ -f "$SAMPLE_SHEET" ]] || die "Sample sheet not found: ${SAMPLE_SHEET}"

log "Refreshing Illumina paths in: ${SAMPLE_SHEET}"
log "Searching under: ${ILLUMINA_DIR}"

TMPFILE=$(mktemp)
trap 'rm -f "$TMPFILE"' EXIT

# ── Preserve comment lines verbatim ──────────────────────────────────────────
grep '^#' "$SAMPLE_SHEET" >> "$TMPFILE" || true

# ── Process header and data rows ─────────────────────────────────────────────
# Read 5 columns (barcode sample_name mode illumina_r1 illumina_r2).
while IFS=$'\t' read -r barcode sample mode r1 r2; do
    # Header row
    if [[ "$barcode" == "barcode" ]]; then
        printf "barcode\tsample_name\tmode\tillumina_r1\tillumina_r2\n" >> "$TMPFILE"
        continue
    fi

    # Find the Illumina directory for this barcode (barcodeXX_*)
    ill_dir=$(find "$ILLUMINA_DIR" -maxdepth 1 -type d -name "${barcode}_*" 2>/dev/null \
              | sort | head -1)

    if [[ -z "$ill_dir" ]]; then
        warn "  ${barcode}: no directory found under ${ILLUMINA_DIR} — keeping existing path"
        printf "%s\t%s\t%s\t%s\t%s\n" \
            "$barcode" "$sample" "$mode" "$r1" "$r2" >> "$TMPFILE"
        continue
    fi

    # Find R1 and R2 — handles mixed naming conventions
    new_r1=$(find "$ill_dir" -maxdepth 1 -type f \
             \( -name "*.fastq.gz" -o -name "*.fq.gz" \) 2>/dev/null \
             | grep -i "R1" | sort | head -1 || true)
    new_r2=$(find "$ill_dir" -maxdepth 1 -type f \
             \( -name "*.fastq.gz" -o -name "*.fq.gz" \) 2>/dev/null \
             | grep -i "R2" | sort | head -1 || true)

    if [[ -z "$new_r1" || -z "$new_r2" ]]; then
        warn "  ${barcode}: could not find R1/R2 in ${ill_dir} — keeping existing path"
        new_r1="$r1"; new_r2="$r2"
    else
        log "  ${barcode}: R1=$(basename "$new_r1")  R2=$(basename "$new_r2")"
    fi

    printf "%s\t%s\t%s\t%s\t%s\n" \
        "$barcode" "$sample" "$mode" "$new_r1" "$new_r2" >> "$TMPFILE"

done < <(grep -v '^#' "$SAMPLE_SHEET")

# ── Write or preview ─────────────────────────────────────────────────────────
if [[ "$DRY_RUN" == "true" ]]; then
    echo ""
    echo "── Dry-run: would write the following sample sheet ──"
    cat "$TMPFILE"
else
    cp "$SAMPLE_SHEET" "${SAMPLE_SHEET}.bak"
    cp "$TMPFILE" "$SAMPLE_SHEET"
    log "Done. Backup saved to ${SAMPLE_SHEET}.bak"
fi
