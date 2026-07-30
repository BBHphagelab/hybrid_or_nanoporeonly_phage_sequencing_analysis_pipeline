#!/usr/bin/env bash
# =============================================================================
# Step 05 — Annotation: Pharokka
# =============================================================================
# Functional annotation of the polished phage genome with Pharokka.
#
# The assembly is already reoriented to terL and its contigs renamed to
# <SAMPLE>_ctgNNN in step 03, so polished.fasta is used directly as input
# (a single consistent reference across annotation and modbase/methylation).
#
# Pharokka database must be downloaded once:
#   pharokka_db_install.py --db_dir /path/to/pharokka_db
# Set PHAROKKA_DB in config.sh to point to the downloaded directory.
#
# Output: $RESULTS_DIR/05_annotated/<barcode>_<sample>/
#   <sample>.gff                         ← GFF annotation
#   <sample>.gbk                         ← GenBank flat file
#   <sample>.faa                         ← protein sequences
#   <sample>_cds_final_merged_output.tsv ← detailed gene table
#   <sample>_summary_output.tsv          ← genome-level summary
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/config.sh"

activate_env
require_tool pharokka.py

SAMPLE_SHEET="${RESULTS_DIR}/sample_sheet.tsv"
mapfile -t BARCODES < <(list_barcodes)
if [[ "${ARRAY_IDX}" -ge "${#BARCODES[@]}" ]]; then
    echo "[SKIP] ARRAY_IDX=${ARRAY_IDX} >= barcodes in sample_sheet (${#BARCODES[@]}) — rien a faire."
    exit 0
fi
BARCODE="${BARCODES[${ARRAY_IDX}]}"
SAMPLE_NAME=$(get_sample_field "$BARCODE" "sample_name")

step_banner "05 — Annotate (Pharokka)  [${BARCODE} / ${SAMPLE_NAME}]"

POLISHED="${RESULTS_DIR}/03_polished/${BARCODE}_${SAMPLE_NAME}/polished.fasta"
OUT="${RESULTS_DIR}/05_annotated/${BARCODE}_${SAMPLE_NAME}"

check_file "$POLISHED"
check_dir  "$PHAROKKA_DB"

GFF="${OUT}/${SAMPLE_NAME}.gff"
if [[ -f "${GFF}" && "${FORCE:-false}" != "true" ]]; then
    log "  Pharokka annotation exists — skipping (set FORCE=true to rerun)"
    exit 0
fi

# ── NFS workaround ────────────────────────────────────────────────────────────
PHAROKKA_LOG_TMP="${TMPDIR:-/tmp}/pharokka_${BARCODE}_${SAMPLE_NAME}_$$.log"
mkdir -p "$(dirname "$OUT")"
# Remove any previous output. On the CIFS/NFS results filesystem a plain
# `rm -rf "$OUT"` can transiently fail with "Directory not empty" (async-unlink
# race), which aborted the job under `set -e`. Move it aside atomically first,
# then best-effort delete; never abort here (pharokka --force recreates it).
if [[ -e "$OUT" ]]; then
    _old="${OUT}.old_$$"
    if mv "$OUT" "$_old" 2>/dev/null; then
        rm -rf "$_old" 2>/dev/null || true
    else
        rm -rf "$OUT" 2>/dev/null || true
    fi
fi

# ---------------------------------------------------------------------------
# Reorientation + contig renaming are done in step 03 (on polished.fasta), so
# polished.fasta is the single, consistently-oriented, collision-proof reference
# shared by annotation (here) and modbase/methylation (06/07).
# ---------------------------------------------------------------------------
PHAROKKA_INPUT="${POLISHED}"

# ---------------------------------------------------------------------------
# Stage 2 — Pharokka: functional annotation
# ---------------------------------------------------------------------------
log "  Stage 2 — Pharokka (functional annotation)..."
log "    Input:  $(basename ${PHAROKKA_INPUT})"
log "    DB:     ${PHAROKKA_DB}"
log "    Prefix: ${SAMPLE_NAME}"

# Recount contigs in the actual input (might differ from POLISHED if reoriented)
N_CONTIGS_INPUT=$(grep -c '^>' "${PHAROKKA_INPUT}" || true)
META_FLAG=""
if [[ "${N_CONTIGS_INPUT}" -gt 1 ]]; then
    META_FLAG="--meta"
    log "    Mode:   meta (${N_CONTIGS_INPUT} contigs)"
else
    log "    Mode:   single genome (${N_CONTIGS_INPUT} contig)"
fi

# Gene predictor: prodigal-gv for ALL samples (robust and consistent; Phanotate
# is known to crash with "return code 1" on some sequences — it killed barcode01).
# Configurable via PHAROKKA_GENE_PREDICTOR in config.sh.
GENE_PRED="${PHAROKKA_GENE_PREDICTOR:-prodigal-gv}"
log "    Gene predictor: ${GENE_PRED}"

pharokka.py \
    -i           "${PHAROKKA_INPUT}" \
    --outdir     "$OUT" \
    --database   "$PHAROKKA_DB" \
    --prefix     "${SAMPLE_NAME}" \
    --threads    "${THREADS}" \
    --force \
    --gene_predictor "${GENE_PRED}" \
    ${META_FLAG} \
    2>&1 | tee "$PHAROKKA_LOG_TMP"

# Copy log to output dir (pharokka has now created $OUT itself)
cp "$PHAROKKA_LOG_TMP" "${OUT}/pharokka.log" 2>/dev/null || true
rm -f "$PHAROKKA_LOG_TMP"

[[ -f "${GFF}" ]] \
    || die "Pharokka failed (${GENE_PRED}) — ${SAMPLE_NAME}.gff not found. See ${OUT}/pharokka.log"

# ---------------------------------------------------------------------------
# Report annotation statistics
# ---------------------------------------------------------------------------
n_cds=$(grep -c $'\tCDS\t' "${GFF}" 2>/dev/null || echo 0)

# Count annotated vs hypothetical from the CDS TSV
CDS_TSV="${OUT}/${SAMPLE_NAME}_cds_final_merged_output.tsv"
n_hypothetical=0
n_annotated=0
if [[ -f "${CDS_TSV}" ]]; then
    # Column "annot" is "hypothetical protein" or a real function
    n_hypothetical=$(awk -F'\t' 'NR>1 && tolower($0) ~ /hypothetical/ {n++} END {print n+0}' "${CDS_TSV}")
    n_annotated=$(( n_cds - n_hypothetical ))
fi

log "  Pharokka done:"
log "    Input used      : $(basename ${PHAROKKA_INPUT})"
log "    CDSs annotated  : ${n_cds}"
log "    Functionally annotated: ${n_annotated}"
log "    Hypothetical    : ${n_hypothetical}"

# Summary TSV (genome-level stats)
SUMMARY_TSV="${OUT}/${SAMPLE_NAME}_summary_output.tsv"
if [[ -f "${SUMMARY_TSV}" ]]; then
    log "  Genome summary:"
    tail -n +2 "${SUMMARY_TSV}" | while IFS=$'\t' read -r key val; do
        log "    ${key}: ${val}"
    done
fi

log "  Outputs in: ${OUT}"
log "  Step 05 done."
