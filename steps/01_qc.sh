#!/usr/bin/env bash
# =============================================================================
# Step 01 — Read QC: chopper filtering + NanoStat report
# =============================================================================
# Per sample (SLURM array: ARRAY_IDX maps to row in sample_sheet.tsv):
#   chopper → length + quality filter on ONT reads
#   NanoStat → QC report on raw and filtered reads
#   Coverage estimate → warn if below threshold
#   Illumina stats → read count + mean length (if illumina_r1/r2 paths are set,
#                    regardless of assembly mode — needed for hybrid short-read polishing, step 03)
#
# Input:  $RESULTS_DIR/00_merged/<barcode>_<sample>.fastq.gz
# Output: $RESULTS_DIR/01_qc/<barcode>_<sample>/
#           filtered.fastq.gz      ← cleaned reads for downstream steps
#           nanostat_raw.txt       ← stats before filtering
#           nanostat_filtered.txt  ← stats after filtering
#           coverage_estimate.txt  ← coverage estimate + ONT mean read length
#           illumina_stats.txt     ← Illumina read count + mean length (if applicable)
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/config.sh"

activate_env
require_tool chopper
require_tool NanoStat

SAMPLE_SHEET="${RESULTS_DIR}/sample_sheet.tsv"

# Resolve this array task's barcode
mapfile -t BARCODES < <(list_barcodes)
if [[ "${ARRAY_IDX:-${SLURM_ARRAY_TASK_ID:-0}}" -ge "${#BARCODES[@]}" ]]; then
    echo "[SKIP] ARRAY_IDX=${ARRAY_IDX:-${SLURM_ARRAY_TASK_ID:-0}} >= barcodes in sample_sheet (${#BARCODES[@]}) — nothing to do."
    exit 0
fi
ARRAY_IDX="${ARRAY_IDX:-${SLURM_ARRAY_TASK_ID:-0}}"
BARCODE="${BARCODES[${ARRAY_IDX}]}"
SAMPLE_NAME=$(get_sample_field "$BARCODE" "sample_name")
MODE=$(get_sample_field        "$BARCODE" "mode")
ILLUMINA_R1=$(get_sample_field "$BARCODE" "illumina_r1")
ILLUMINA_R2=$(get_sample_field "$BARCODE" "illumina_r2")

step_banner "01 — QC  [${BARCODE} / ${SAMPLE_NAME}]"

MERGED_ONT="${RESULTS_DIR}/00_merged/${BARCODE}_${SAMPLE_NAME}.fastq.gz"
OUT="${RESULTS_DIR}/01_qc/${BARCODE}_${SAMPLE_NAME}"
mkdir -p "$OUT"

check_file "$MERGED_ONT"

# Skip if already done
if [[ -f "${OUT}/filtered.fastq.gz" && "${FORCE:-false}" != "true" ]]; then
    log "  Filtered reads exist — skipping (delete ${OUT} or set FORCE=true to rerun)"
    exit 0
fi

# ---------------------------------------------------------------------------
# NanoStat on raw reads (before filtering)
# ---------------------------------------------------------------------------
log "  Running NanoStat on raw reads..."
NanoStat \
    --fastq "$MERGED_ONT" \
    --threads "${THREADS}" \
    > "${OUT}/nanostat_raw.txt" 2>&1

# Extract total bases and read count for coverage estimate (before filter)
raw_bases=$(grep -oP 'Total bases.*?:\s+\K[\d,]+' "${OUT}/nanostat_raw.txt" \
    | tr -d ',' | head -1 || echo "0")
raw_reads=$(grep -oP 'Number of reads.*?:\s+\K[\d,]+' "${OUT}/nanostat_raw.txt" \
    | tr -d ',' | head -1 || echo "0")
ont_raw_mean_len=$(grep -oP 'Mean read length.*?:\s+\K[\d,.]+' "${OUT}/nanostat_raw.txt" \
    | tr -d ',' | head -1 || echo "?")

log "  Raw: ${raw_reads} reads, ${raw_bases} bases, mean length ${ont_raw_mean_len} bp"

# ---------------------------------------------------------------------------
# Rough coverage estimate (using configured genome size)
# ---------------------------------------------------------------------------
# Convert genome size string (e.g. "100k", "6.4m", "2g") to bp
genome_bp=$(echo "${GENOME_SIZE_ESTIMATE}" \
    | awk '{
        val=$1
        if      (val ~ /[Kk]$/) { sub(/[Kk]$/, "", val); val=val*1000 }
        else if (val ~ /[Mm]$/) { sub(/[Mm]$/, "", val); val=val*1000000 }
        else if (val ~ /[Gg]$/) { sub(/[Gg]$/, "", val); val=val*1000000000 }
        printf "%.0f", val
    }')

coverage_raw=$(awk "BEGIN{printf \"%.1f\", ${raw_bases}/${genome_bp}}")
log "  Estimated raw coverage: ~${coverage_raw}×"

if awk "BEGIN{exit (${coverage_raw} >= ${MIN_COVERAGE_WARN})}"; then
    warn "  Coverage ~${coverage_raw}× is below warning threshold (${MIN_COVERAGE_WARN}×)"
    warn "  Assembly quality may be reduced. Check that this barcode is not empty."
fi
if awk "BEGIN{exit (${coverage_raw} >= ${MIN_COVERAGE_FAIL})}"; then
    die "  Coverage ~${coverage_raw}× is below minimum threshold (${MIN_COVERAGE_FAIL}×)"
fi

# ---------------------------------------------------------------------------
# chopper — quality and length filtering
# ---------------------------------------------------------------------------
log "  Running chopper (Q≥${MIN_READ_QUALITY}, len≥${MIN_READ_LENGTH})..."
zcat "$MERGED_ONT" \
    | chopper \
        --quality "${MIN_READ_QUALITY}" \
        --minlength "${MIN_READ_LENGTH}" \
        --threads "${THREADS}" \
    2>"${OUT}/chopper.log" \
    | gzip > "${OUT}/filtered.fastq.gz"

[[ -s "${OUT}/filtered.fastq.gz" ]] \
    || die "chopper produced empty output — check ${OUT}/chopper.log"

# ---------------------------------------------------------------------------
# NanoStat on filtered reads
# ---------------------------------------------------------------------------
log "  Running NanoStat on filtered reads..."
NanoStat \
    --fastq "${OUT}/filtered.fastq.gz" \
    --threads "${THREADS}" \
    > "${OUT}/nanostat_filtered.txt" 2>&1 || true

filt_bases=$(grep -oP 'Total bases.*?:\s+\K[\d,]+' "${OUT}/nanostat_filtered.txt" \
    | tr -d ',' | head -1 || echo "0")
filt_reads=$(grep -oP 'Number of reads.*?:\s+\K[\d,]+' "${OUT}/nanostat_filtered.txt" \
    | tr -d ',' | head -1 || echo "0")
coverage_filt=$(awk "BEGIN{printf \"%.1f\", ${filt_bases}/${genome_bp}}")

log "  Filtered: ${filt_reads} reads, ${filt_bases} bases, ~${coverage_filt}× coverage"

# Write coverage summary
{
    echo "barcode=${BARCODE}"
    echo "sample=${SAMPLE_NAME}"
    echo "genome_size_estimate=${GENOME_SIZE_ESTIMATE}"
    echo "raw_reads=${raw_reads}"
    echo "raw_bases=${raw_bases}"
    echo "raw_coverage_est=${coverage_raw}"
    echo "ont_raw_mean_len=${ont_raw_mean_len}"
    echo "filtered_reads=${filt_reads}"
    echo "filtered_bases=${filt_bases}"
    echo "filtered_coverage_est=${coverage_filt}"
} > "${OUT}/coverage_estimate.txt"

# ---------------------------------------------------------------------------
# Illumina read stats — if illumina_r1/r2 are set in sample sheet
# (checked regardless of assembly mode: Illumina reads may be used for short-read polishing (step 03)
#  in step 09 even when assembly was done in ONT-only mode)
# ---------------------------------------------------------------------------
if [[ "${ILLUMINA_R1}" != "-" && -n "${ILLUMINA_R1}" && -f "${ILLUMINA_R1}" ]]; then
    log "  Computing Illumina read stats (sampling first 50 000 reads from R1)..."

    # zcat -f handles both gzipped (.fastq.gz) and plain (.fastq) files
    # pipefail disabled: zcat receives SIGPIPE (exit 141) when head closes early — non-fatal
    set +o pipefail
    illum_stats=$(zcat -f "${ILLUMINA_R1}" 2>/dev/null \
        | head -200000 \
        | awk 'NR%4==2 { sum += length($0); n++ }
               END     { printf "%d\t%.1f\n", n, (n>0 ? sum/n : 0) }')
    set -o pipefail

    illum_r1_sampled=$(echo "${illum_stats}" | cut -f1)
    illum_mean_len=$(echo "${illum_stats}"    | cut -f2)
    illum_total_est=$(( illum_r1_sampled * 2 ))

    {
        echo "illumina_r1=${ILLUMINA_R1}"
        echo "illumina_r1_sampled=${illum_r1_sampled}"
        echo "illumina_total_est=${illum_total_est}"
        echo "illumina_mean_len=${illum_mean_len}"
    } > "${OUT}/illumina_stats.txt"

    log "  Illumina: ~${illum_total_est} reads (R1+R2 est.), mean length ${illum_mean_len} bp"
elif [[ "${ILLUMINA_R1}" != "-" && -n "${ILLUMINA_R1}" ]]; then
        warn "  Illumina R1 path set but file not found: ${ILLUMINA_R1} — illumina_stats.txt not generated"
fi

log "  Step 01 done."
