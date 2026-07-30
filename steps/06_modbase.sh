#!/usr/bin/env bash
# =============================================================================
# Step 06 — Modified base detection
# =============================================================================
# Detects DNA modifications (6mA, 5mCG, 5hmCG) from mod-tagged BAM files
# produced directly by MinKNOW (modification calling enabled at sequencing time).
#
# INPUT BAM FILES — MinKNOW layout
# ─────────────────────────────────
# MinKNOW with modification calling produces per-barcode BAMs in MOD_DEMUX_DIR.
# Two supported layouts (auto-detected):
#
#   Layout A — flat directory (one BAM per barcode):
#     MOD_DEMUX_DIR/
#       barcode01.bam
#       barcode02.bam
#
#   Layout B — subdirectory per barcode (multiple chunk BAMs):
#     MOD_DEMUX_DIR/
#       barcode01/
#         PAW12345_pass_barcode01_*.bam
#         PAW12345_pass_barcode01_*.bam   ← merged by this step
#       barcode02/
#         ...
#
# In layout B the chunks are merged with samtools merge before processing.
# The merged BAM is written to MOD_DEMUX_DIR/barcode01_merged.bam and reused
# on subsequent runs (not regenerated if it already exists).
#
# This step (per barcode):
#   1. Locate / merge mod-tagged BAM(s)
#   2. Align mod-tagged reads to the polished assembly (minimap2 -ax map-ont -y)
#   3. Sort and index the aligned BAM
#   4. Run modkit pileup → per-base methylation BED file
#
# Output: $RESULTS_DIR/06_modbase/<barcode>_<sample>/
#   aligned_mods.bam      ← aligned, sorted, indexed mod-tagged BAM
#   methylation.bed       ← per-base modification frequencies
#   modkit.log
#
# Modification codes in the BED (column 4):
#   a  = 6mA  (N6-methyladenine)
#   m  = 5mC  (5-methylcytosine)
#   h  = 5hmC (5-hydroxymethylcytosine)
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/config.sh"

SORT_TMPDIR="${TMPDIR:-/tmp}/samtools_modbase_$$"
mkdir -p "${SORT_TMPDIR}"
trap 'rm -rf "${SORT_TMPDIR}"' EXIT

# ---------------------------------------------------------------------------
# Skip entirely if MOD_DEMUX_DIR is not configured
# ---------------------------------------------------------------------------
if [[ -z "${MOD_DEMUX_DIR:-}" ]]; then
    echo "[SKIP] Step 06 — MOD_DEMUX_DIR not set in config.sh. Skipping."
    echo "       Set MOD_DEMUX_DIR to your MinKNOW bam_pass directory to enable."
    exit 0
fi

if [[ ! -d "${MOD_DEMUX_DIR}" ]]; then
    echo "[SKIP] Step 06 — MOD_DEMUX_DIR does not exist: ${MOD_DEMUX_DIR}"
    exit 0
fi

activate_env
export PATH="${MAMBA_ROOT_PREFIX}/envs/${MAMBA_ENV}/bin:${PATH}"
load_samtools_cluster   # samtools/1.21 required for bam2fq -T, sort, view

require_tool minimap2
require_tool samtools
require_tool modkit

mapfile -t BARCODES < <(list_barcodes)
if [[ "${ARRAY_IDX}" -ge "${#BARCODES[@]}" ]]; then
    echo "[SKIP] ARRAY_IDX=${ARRAY_IDX} >= barcodes in sample_sheet (${#BARCODES[@]}) — rien a faire."
    exit 0
fi
BARCODE="${BARCODES[${ARRAY_IDX}]}"
SAMPLE_NAME=$(get_sample_field "$BARCODE" "sample_name")

step_banner "06 — Modified bases  [${BARCODE} / ${SAMPLE_NAME}]"

POLISHED="${RESULTS_DIR}/03_polished/${BARCODE}_${SAMPLE_NAME}/polished.fasta"
OUT="${RESULTS_DIR}/06_modbase/${BARCODE}_${SAMPLE_NAME}"
mkdir -p "$OUT"

check_file "$POLISHED"

# ---------------------------------------------------------------------------
# Locate the mod-base BAM for this barcode
#
# MinKNOW naming conventions tested (in order):
#   A) flat:   MOD_DEMUX_DIR/barcode01.bam
#   B) subdir: MOD_DEMUX_DIR/barcode01/*.bam  (one or more chunks)
#   C) flat wildcard: MOD_DEMUX_DIR/*barcode01*.bam (any prefix/suffix)
# ---------------------------------------------------------------------------
MOD_BAM=""

# Layout A: flat file named exactly barcode01.bam
if [[ -f "${MOD_DEMUX_DIR}/${BARCODE}.bam" ]]; then
    MOD_BAM="${MOD_DEMUX_DIR}/${BARCODE}.bam"
    log "  Layout A — single BAM: $(basename "$MOD_BAM")"

# Layout B: subdirectory barcode01/ with one or more chunk BAMs
elif [[ -d "${MOD_DEMUX_DIR}/${BARCODE}" ]]; then
    CHUNK_DIR="${MOD_DEMUX_DIR}/${BARCODE}"
    mapfile -t CHUNKS < <(find "${CHUNK_DIR}" -maxdepth 1 -name "*.bam" | sort)

    if [[ ${#CHUNKS[@]} -eq 0 ]]; then
        warn "  Subdirectory ${CHUNK_DIR} found but contains no BAM files."
    elif [[ ${#CHUNKS[@]} -eq 1 ]]; then
        MOD_BAM="${CHUNKS[0]}"
        log "  Layout B — single chunk: $(basename "$MOD_BAM")"
    else
        # Multiple chunks — merge into a single BAM (cached for reuse)
        MERGED_BAM="${MOD_DEMUX_DIR}/${BARCODE}_merged.bam"
        if [[ -f "${MERGED_BAM}" ]]; then
            log "  Layout B — merged BAM already exists: $(basename "$MERGED_BAM")"
        else
            log "  Layout B — merging ${#CHUNKS[@]} chunks → $(basename "$MERGED_BAM")..."
            samtools merge -@ "${THREADS}" -f "${MERGED_BAM}" "${CHUNKS[@]}" \
                || die "samtools merge failed for ${BARCODE}"
            samtools index "${MERGED_BAM}"
            log "  Merge done."
        fi
        MOD_BAM="${MERGED_BAM}"
    fi

# Layout C: flat directory, any file matching *barcode01*.bam (fallback)
else
    MOD_BAM=$(find "${MOD_DEMUX_DIR}" -maxdepth 1 -name "*${BARCODE}*.bam" | head -1 || true)
    [[ -n "$MOD_BAM" ]] && log "  Layout C — wildcard match: $(basename "$MOD_BAM")"
fi

# Give up if still not found
if [[ -z "$MOD_BAM" || ! -f "$MOD_BAM" ]]; then
    warn "  Mod-base BAM not found for ${BARCODE} in ${MOD_DEMUX_DIR}"
    warn "  Searched:"
    warn "    A) ${MOD_DEMUX_DIR}/${BARCODE}.bam"
    warn "    B) ${MOD_DEMUX_DIR}/${BARCODE}/*.bam"
    warn "    C) ${MOD_DEMUX_DIR}/*${BARCODE}*.bam"
    warn "  Files present: $(ls "${MOD_DEMUX_DIR}"/*.bam 2>/dev/null | xargs -n1 basename | tr '\n' ' ' || echo '(none)')"
    warn "  Subdirs present: $(find "${MOD_DEMUX_DIR}" -maxdepth 1 -mindepth 1 -type d | xargs -n1 basename | tr '\n' ' ' || echo '(none)')"
    warn "  Skipping ${BARCODE}."
    exit 0
fi
log "  Mod-base BAM: $(basename "$MOD_BAM")"

ALIGNED_BAM="${OUT}/aligned_mods.bam"
METHYL_BED="${OUT}/methylation.bed"

if [[ -f "$METHYL_BED" && "${FORCE:-false}" != "true" ]]; then
    log "  Modbase results already exist — skipping (set FORCE=true to rerun)"
    exit 0
fi

# ---------------------------------------------------------------------------
# 1. Align mod-tagged reads to polished assembly
#    bam2fq -T MM,ML : preserve modification tags through the FASTQ conversion
#    minimap2 -y     : copy per-read tags from FASTQ header to aligned SAM
#    --secondary=no  : skip secondary alignments (modkit requires primary only)
# ---------------------------------------------------------------------------
log "  Aligning mod-tagged reads to polished assembly..."
log "    Input BAM:  ${MOD_BAM}"
log "    Reference:  $(basename "$POLISHED")"

samtools bam2fq -@ "${THREADS}" -T MM,ML "$MOD_BAM" \
| minimap2 \
    -ax map-ont \
    -y \
    --secondary=no \
    -t "${THREADS}" \
    "$POLISHED" \
    - \
| samtools sort \
    -@ "${THREADS}" \
    -T "${SORT_TMPDIR}/sort" \
    -o "$ALIGNED_BAM"

samtools index "$ALIGNED_BAM"

TOTAL_READS=$(samtools view -c "$ALIGNED_BAM")
MAPPED_READS=$(samtools view -c -F 4 "$ALIGNED_BAM")
log "  Alignment: ${MAPPED_READS} / ${TOTAL_READS} reads mapped"

# ---------------------------------------------------------------------------
# 2. modkit pileup — per-base modification frequencies
# ---------------------------------------------------------------------------
log "  Running modkit pileup..."

modkit pileup \
    --ref      "$POLISHED" \
    --threads  "${THREADS}" \
    "$ALIGNED_BAM" \
    "$METHYL_BED" \
    2>"${OUT}/modkit.log"

[[ -s "$METHYL_BED" ]] \
    || die "modkit produced empty output — check ${OUT}/modkit.log"

# ---------------------------------------------------------------------------
# 3. Summary
# ---------------------------------------------------------------------------
log "  Methylation summary (all coverage levels):"
awk '$4=="a"     {n++} END {printf "    6mA  positions: %d\n", n+0}' "$METHYL_BED"
awk '$4=="m"     {n++} END {printf "    5mC  positions: %d\n", n+0}' "$METHYL_BED"
# 4mC — the modification the all-context re-basecalling adds (modkit emits ChEBI
# code 21839 for 4mC; 'c' as a fallback). Counting it here so it is never masked.
awk '$4=="21839" || $4=="c" {n++} END {printf "    4mC  positions: %d\n", n+0}' "$METHYL_BED"
awk '$4=="h"     {n++} END {printf "    5hmC positions: %d\n", n+0}' "$METHYL_BED"
# Non-silent: dump EVERY mod code actually present so a new modification can never
# be hidden by a hard-coded list again (col4 = raw mod code).
log "  Mod codes present (code:count):"
cut -f4 "$METHYL_BED" | sort | uniq -c | awk '{printf "    %s : %s\n", $2, $1}'

log "  High-confidence sites (>= 8x coverage, >= 50% modification frequency):"
awk '$4=="a"     && $10>=8 && $11>=50 {n++} END {printf "    6mA  hi-conf: %d\n", n+0}' "$METHYL_BED"
awk '$4=="m"     && $10>=8 && $11>=50 {n++} END {printf "    5mC  hi-conf: %d\n", n+0}' "$METHYL_BED"
awk '($4=="21839"||$4=="c") && $10>=8 && $11>=50 {n++} END {printf "    4mC  hi-conf: %d\n", n+0}' "$METHYL_BED"

log "  Output: ${OUT}/"
log "  Step 06 done."
