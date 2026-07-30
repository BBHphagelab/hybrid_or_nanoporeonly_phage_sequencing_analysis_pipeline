#!/usr/bin/env bash
# =============================================================================
# Step 01b — OPTIONAL host-read depletion (pre-assembly)  [SLURM array]
# =============================================================================
# Runs ONLY if HOST_DEPLETION=true in config.sh. Removes host-derived ONT reads
# BEFORE assembly by competitive mapping of the QC'd reads against a combined
# reference = HOST_REF + optional phage references (HOST_DEPLETE_PHAGE_REFS). A read whose
# PRIMARY alignment lands on a host sequence is dropped; reads whose best hit is
# a phage sequence — or that don't map at all — are kept. This is orthogonal to
# the Filtlong/Flye coverage handling (depth) and to the post-assembly
# <MIN_CONTIG_LEN parasite filter (it prevents host reads from entering the
# assembly graph in the first place, avoiding host-phage chimeras).
#
# Input:  01_qc/<bc>_<sample>/filtered.fastq.gz
# Output: 01b_hostdepleted/<bc>_<sample>/
#           filtered.host_depleted.fastq.gz   ← used by step 02 if present
#           depletion_stats.txt               ← total / host / kept counts
# Step 02 automatically prefers the depleted reads when HOST_DEPLETION=true.
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/config.sh"

# Skip entirely unless enabled
if [[ "${HOST_DEPLETION:-false}" != "true" ]]; then
    echo "[SKIP] Step 01b — HOST_DEPLETION not 'true' in config.sh. Nothing to do."
    exit 0
fi

activate_env
load_samtools_cluster
require_tool minimap2
require_tool seqkit

SAMPLE_SHEET="${RESULTS_DIR}/sample_sheet.tsv"
mapfile -t BARCODES < <(list_barcodes)
IDX="${ARRAY_IDX:-${SLURM_ARRAY_TASK_ID:-0}}"
if [[ "${IDX}" -ge "${#BARCODES[@]}" ]]; then
    echo "[SKIP] ARRAY_IDX=${IDX} >= barcodes (${#BARCODES[@]}) — nothing to do."
    exit 0
fi
BARCODE="${BARCODES[${IDX}]}"
SAMPLE_NAME=$(get_sample_field "$BARCODE" "sample_name")
LABEL="${BARCODE}_${SAMPLE_NAME}"

step_banner "01b — Host depletion  [${BARCODE} / ${SAMPLE_NAME}]"

READS="${RESULTS_DIR}/01_qc/${LABEL}/filtered.fastq.gz"
OUT="${RESULTS_DIR}/01b_hostdepleted/${LABEL}"
DEPLETED="${OUT}/filtered.host_depleted.fastq.gz"
mkdir -p "$OUT"

# ── Input guard ──────────────────────────────────────────────────────────────
if [[ ! -s "$READS" ]]; then
    echo "[SKIP] ${LABEL}: filtered reads missing (${READS}) — step 01 did not complete."
    exit 0
fi
if [[ -s "$DEPLETED" && "${FORCE:-false}" != "true" ]]; then
    log "  Depleted reads exist — skipping (set FORCE=true to rerun)"
    exit 0
fi
if [[ -z "${HOST_REF:-}" || ! -f "${HOST_REF}" ]]; then
    warn "  HOST_REF not set or not found ('${HOST_REF:-}') — cannot deplete. Step 02 will use full reads."
    exit 0
fi

# ---------------------------------------------------------------------------
# Build combined reference: host + phage(s); record host contig names.
# ---------------------------------------------------------------------------
to_fasta() {  # $1 in  $2 out(append)
    local in="$1" out="$2"
    case "$in" in
        *.gbk|*.gbff|*.gb) any2fasta "$in" >> "$out" 2>/dev/null ;;
        *)                 cat "$in"        >> "$out" ;;
    esac
}

HOST_FA="${OUT}/host.fasta"; : > "$HOST_FA"
to_fasta "$HOST_REF" "$HOST_FA"
[[ -s "$HOST_FA" ]] || die "Could not read HOST_REF into FASTA: ${HOST_REF}"
HOST_NAMES="${OUT}/host_contig_names.txt"
grep '^>' "$HOST_FA" | sed 's/^>//; s/[[:space:]].*//' | sort -u > "$HOST_NAMES"
log "  Host contigs: $(wc -l < "$HOST_NAMES")"

PHAGE_FA="${OUT}/phage_union.fasta"; : > "$PHAGE_FA"
for v in ${HOST_DEPLETE_PHAGE_REFS:-}; do
    [[ -n "$v" && -f "$v" ]] && to_fasta "$v" "$PHAGE_FA"
done
[[ -s "$PHAGE_FA" ]] || warn "  No phage references (HOST_DEPLETE_PHAGE_REFS empty) — competitive depletion reduces to host-only mapping."

COMBINED="${OUT}/combined_ref.fasta"
cat "$HOST_FA" "$PHAGE_FA" > "$COMBINED"

# ---------------------------------------------------------------------------
# Competitive mapping → host-primary read names → keep the rest
# ---------------------------------------------------------------------------
LOCAL="${TMPDIR:-/tmp}/hostdep_${LABEL}_$$"
rm -rf "$LOCAL"; mkdir -p "$LOCAL"
trap 'rm -rf "$LOCAL"' EXIT

log "  Competitive mapping (minimap2 map-ont, host + phage)..."
minimap2 -ax map-ont -t "${THREADS}" "$COMBINED" "$READS" 2>"${OUT}/minimap2.log" \
    | samtools view -b -o "${LOCAL}/all.bam" - 2>>"${OUT}/minimap2.log"

# Primary, mapped alignments only (-F 0x904 = exclude unmapped + secondary + supplementary)
HOST_READS="${LOCAL}/host_reads.txt"
samtools view -F 0x904 "${LOCAL}/all.bam" \
    | awk -v OFS='\t' '{print $1, $3}' \
    | awk 'NR==FNR{h[$1]=1; next} ($2 in h){print $1}' "$HOST_NAMES" - \
    | sort -u > "$HOST_READS"

total=$(zcat "$READS" | awk 'NR%4==1' | wc -l)
n_host=$(wc -l < "$HOST_READS")

log "  Removing ${n_host} host reads / ${total} total..."
if [[ "${n_host}" -gt 0 ]]; then
    seqkit grep -v -f "$HOST_READS" "$READS" -o "$DEPLETED" 2>>"${OUT}/minimap2.log"
else
    cp "$READS" "$DEPLETED"
fi
kept=$(zcat "$DEPLETED" | awk 'NR%4==1' | wc -l)

{
    echo "barcode=${BARCODE}"
    echo "sample=${SAMPLE_NAME}"
    echo "host_ref=${HOST_REF}"
    echo "total_reads=${total}"
    echo "host_reads_removed=${n_host}"
    echo "kept_reads=${kept}"
    echo "pct_host=$(awk -v h="${n_host}" -v t="${total}" 'BEGIN{printf "%.2f", (t>0)?100*h/t:0}')"
} > "${OUT}/depletion_stats.txt"

# Keep small refs for record; drop the big combined to save space
rm -f "$COMBINED"
log "  Depletion done: kept ${kept} / ${total} reads ($(grep pct_host "${OUT}/depletion_stats.txt" | cut -d= -f2)% host removed)."
log "  Output: ${DEPLETED}"
log "  Step 01b done."
