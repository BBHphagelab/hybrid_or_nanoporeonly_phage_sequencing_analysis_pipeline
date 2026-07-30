#!/usr/bin/env bash
# =============================================================================
# Step 00 — Setup: detect samples, merge ONT FASTQs, generate sample sheet
# =============================================================================
# Phage pipeline: samples are discovered from FASTQ_PASS_DIR (ONT reads).
# No Illumina directory is required; hybrid mode is optional.
#
# Sample naming (in order of priority):
#   1. sample_names.tsv  in RUN_DIR (two-column TSV: barcode<TAB>name)
#      → create this file manually before running step 00 to assign proper names
#   2. Auto-generated name from barcode number (e.g. barcode03 → phage03)
#
# Illumina reads (optional hybrid polishing):
#   If ILLUMINA_DIR exists and contains barcodeXX_<name>/ subfolders with
#   R1/R2 FASTQ files, mode is set to "hybrid" for those samples.
#   Otherwise mode is "ont".
#
# Outputs:
#   $RESULTS_DIR/sample_sheet.tsv
#   $RESULTS_DIR/00_merged/<barcode>_<name>.fastq.gz
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/config.sh"

step_banner "00 — Setup & FASTQ Merge"
activate_env

mkdir -p "${RESULTS_DIR}/00_merged" "${RESULTS_DIR}/logs"

SAMPLE_SHEET="${RESULTS_DIR}/sample_sheet.tsv"

# ---------------------------------------------------------------------------
# Helper: find R1/R2 Illumina files in a directory (multiple naming styles)
# ---------------------------------------------------------------------------
find_illumina_reads() {
    local dir="$1" strand="$2"
    local files=()
    while IFS= read -r -d '' f; do
        files+=("$f")
    done < <(find "$dir" -maxdepth 1 \
        \( -name "*_${strand}*.fastq.gz" \
           -o -name "*_${strand}*.fq.gz" \
           -o -name "*_${strand}_*.fastq" \
           -o -name "${strand}_*.fastq" \
           -o -name "${strand}_*.fastq.gz" \
        \) \
        -not -name "merged_*" \
        -print0 | sort -z)
    echo "${files[@]:-}"
}

# ---------------------------------------------------------------------------
# 1. Build sample name map (priority: SAMPLE_NAMES_MAP in config > sample_names.tsv)
# ---------------------------------------------------------------------------
declare -A CUSTOM_NAMES=()

# Source 1 — SAMPLE_NAMES_MAP variable from config.sh
# Format: "barcodeXX:name barcodeXX:name ..."
if [[ -n "${SAMPLE_NAMES_MAP:-}" ]]; then
    log "Reading sample names from SAMPLE_NAMES_MAP (config.sh)..."
    for pair in ${SAMPLE_NAMES_MAP}; do
        bc="${pair%%:*}"
        name="${pair#*:}"
        CUSTOM_NAMES["$bc"]="$name"
        log "  ${bc} → ${name}"
    done
fi

# Source 2 — fallback: sample_names.tsv file in RUN_DIR (two-column TSV: barcode<TAB>name)
NAMES_FILE="${RUN_DIR}/sample_names.tsv"
if [[ -f "${NAMES_FILE}" ]]; then
    log "Reading sample names from ${NAMES_FILE}..."
    while IFS=$'\t' read -r bc name; do
        [[ "$bc" == "#"* || -z "$bc" ]] && continue
        # config takes priority — only add if not already set
        [[ -z "${CUSTOM_NAMES[$bc]:-}" ]] && CUSTOM_NAMES["$bc"]="$name"
        log "  ${bc} → ${name}  (from file)"
    done < "${NAMES_FILE}"
fi

# ---------------------------------------------------------------------------
# 2. Discover samples from FASTQ_PASS_DIR
# ---------------------------------------------------------------------------
log "Scanning ${FASTQ_PASS_DIR} for barcode directories..."

[[ -d "${FASTQ_PASS_DIR}" ]] \
    || die "FASTQ_PASS_DIR not found: ${FASTQ_PASS_DIR}"

declare -a BARCODES=()
declare -A SAMPLE_NAMES=()
declare -A ILLUMINA_R1=()
declare -A ILLUMINA_R2=()
declare -A MODES=()

while IFS= read -r -d '' bc_dir; do
    folder_name=$(basename "$bc_dir")
    # Accept barcodeXX directories
    [[ "$folder_name" =~ ^barcode[0-9]{2,}$ ]] || {
        warn "Skipping unexpected directory: ${folder_name}"
        continue
    }
    barcode="$folder_name"

    # Assign sample name
    if [[ -n "${CUSTOM_NAMES[$barcode]:-}" ]]; then
        sample_name="${CUSTOM_NAMES[$barcode]}"
    else
        # Auto-generate: phage01, phage02, …
        num="${barcode#barcode}"
        num="${num#0}"   # strip leading zero for arithmetic
        [[ -z "$num" ]] && num=0
        sample_name="$(printf 'phage%02d' "$num")"
    fi

    # Count FASTQ chunks — skip empty barcodes (no-template control, etc.)
    n_chunks=$(find "$bc_dir" -name "*.fastq.gz" 2>/dev/null | wc -l)
    if [[ "$n_chunks" -eq 0 ]]; then
        warn "  [${barcode}] No FASTQ chunks found — skipping"
        continue
    fi

    # Quick size check: MinKNOW creates barcode dirs even for empty wells.
    # A real sample should have at least ~1 MB of compressed reads.
    total_bytes=$(find "$bc_dir" -name "*.fastq.gz" -exec du -sb {} + 2>/dev/null \
                  | awk '{s+=$1} END{print s+0}')
    if [[ "$total_bytes" -lt 1048576 ]]; then
        warn "  [${barcode}] Only ${total_bytes} bytes of FASTQ data — likely empty barcode (MinKNOW artifact), skipping"
        continue
    fi

    BARCODES+=("$barcode")
    SAMPLE_NAMES["$barcode"]="$sample_name"
    ILLUMINA_R1["$barcode"]="-"
    ILLUMINA_R2["$barcode"]="-"
    MODES["$barcode"]="ont"

    log "  Found: ${barcode} → ${sample_name}  (${n_chunks} FASTQ chunks)"

done < <(find "${FASTQ_PASS_DIR}" -maxdepth 1 -mindepth 1 -type d -print0 | sort -z)

[[ ${#BARCODES[@]} -eq 0 ]] && die "No barcode directories found in ${FASTQ_PASS_DIR}"
log "Total samples detected: ${#BARCODES[@]}"

# ---------------------------------------------------------------------------
# 3. Optional: detect Illumina reads for hybrid mode
# ---------------------------------------------------------------------------
if [[ -d "${ILLUMINA_DIR:-}" ]]; then
    log "Illumina directory found — checking for hybrid-mode reads..."
    for bc in "${BARCODES[@]}"; do
        bc_dir=$(find "${ILLUMINA_DIR}" -maxdepth 1 -type d -name "${bc}_*" 2>/dev/null | head -1 || true)
        [[ -z "$bc_dir" ]] && continue

        r1_files=( $(find_illumina_reads "$bc_dir" "R1") )
        r2_files=( $(find_illumina_reads "$bc_dir" "R2") )

        if [[ ${#r1_files[@]} -gt 0 && -n "${r1_files[0]:-}" ]]; then
            MODES["$bc"]="hybrid"
            ILLUMINA_R1["$bc"]="${r1_files[*]}"
            ILLUMINA_R2["$bc"]="${r2_files[*]}"
            log "  [${bc}] Illumina reads found → hybrid mode"
        fi
    done
else
    log "ILLUMINA_DIR not present — all samples will run in ONT-only mode."
fi

# ---------------------------------------------------------------------------
# 4. Write sample_sheet.tsv
# ---------------------------------------------------------------------------
if [[ -f "${SAMPLE_SHEET}" && "${FORCE:-false}" != "true" ]]; then
    log "sample_sheet.tsv already exists — preserving (use FORCE=true to regenerate)"
    log "  Path: ${SAMPLE_SHEET}"
else
    log "Writing ${SAMPLE_SHEET}..."
    {
        echo "# Nanopore phage pipeline — sample sheet"
        echo "# mode        : 'ont' (long reads only) or 'hybrid' (set automatically if Illumina reads are found)"
        echo "# illumina_r1 / illumina_r2: set automatically if ILLUMINA_DIR is present (used for short-read polishing in step 03)"
        echo "# Each barcode is analysed independently."
        printf "barcode\tsample_name\tmode\tillumina_r1\tillumina_r2\n"
        for bc in "${BARCODES[@]}"; do
            printf "%s\t%s\t%s\t%s\t%s\n" \
                "$bc" \
                "${SAMPLE_NAMES[$bc]}" \
                "${MODES[$bc]}" \
                "${ILLUMINA_R1[$bc]}" \
                "${ILLUMINA_R2[$bc]}"
        done
    } > "${SAMPLE_SHEET}"
    log "  sample_sheet.tsv written."
fi

# ---------------------------------------------------------------------------
# 5. Merge Illumina reads if multiple files per sample
# ---------------------------------------------------------------------------
ILLUMINA_MERGED=false
for bc in "${BARCODES[@]}"; do
    [[ "${MODES[$bc]}" != "hybrid" ]] && continue
    r1_files=( ${ILLUMINA_R1[$bc]} )
    r2_files=( ${ILLUMINA_R2[$bc]} )
    [[ ${#r1_files[@]} -le 1 ]] && continue

    bc_dir=$(find "${ILLUMINA_DIR}" -maxdepth 1 -type d -name "${bc}_*" | head -1)
    merged_r1="${bc_dir}/merged_R1.fastq.gz"
    merged_r2="${bc_dir}/merged_R2.fastq.gz"

    if [[ -f "$merged_r1" && "${FORCE:-false}" != "true" ]]; then
        log "  [${bc}] Merged Illumina already exists — skipping"
    else
        log "  [${bc}] Merging ${#r1_files[@]} R1 files → $(basename $merged_r1)"
        cat "${r1_files[@]}" > "$merged_r1"
        cat "${r2_files[@]}" > "$merged_r2"
        ILLUMINA_MERGED=true
    fi

    # Update sample sheet to point to merged files (cols: barcode, sample_name, mode, illumina_r1, illumina_r2)
    sed -i "s|^${bc}\t\([^\t]*\)\t\([^\t]*\)\t[^\t]*\t[^\t]*|${bc}\t\1\t\2\t${merged_r1}\t${merged_r2}|" \
        "${SAMPLE_SHEET}"
done

# ---------------------------------------------------------------------------
# 6. Concatenate ONT FASTQ chunks per barcode
# ---------------------------------------------------------------------------
log "Merging ONT FASTQ chunks per barcode..."

for bc in "${BARCODES[@]}"; do
    ont_source="${FASTQ_PASS_DIR}/${bc}"
    merged_ont="${RESULTS_DIR}/00_merged/${bc}_${SAMPLE_NAMES[$bc]}.fastq.gz"

    if [[ -f "$merged_ont" && "${FORCE:-false}" != "true" ]]; then
        log "  [${bc}] Merged ONT FASTQ exists — skipping"
        continue
    fi

    n_chunks=$(find "$ont_source" -name "*.fastq.gz" | wc -l)
    [[ "$n_chunks" -eq 0 ]] && die "[${bc}] No FASTQ files found in ${ont_source}"

    log "  [${bc}] Concatenating ${n_chunks} chunks → $(basename $merged_ont)"
    cat "${ont_source}"/*.fastq.gz > "$merged_ont"

    n_reads=$(zcat "$merged_ont" | awk 'NR%4==1' | wc -l)
    log "  [${bc}] Total reads after merge: ${n_reads}"

    [[ "$n_reads" -lt 500 ]] && \
        warn "[${bc}] Very few reads (${n_reads}) — this barcode may be a contaminant or empty well."
done

log ""
log "Step 00 complete."
log "  Samples     : ${#BARCODES[@]}"
log "  Sample sheet: ${SAMPLE_SHEET}"
