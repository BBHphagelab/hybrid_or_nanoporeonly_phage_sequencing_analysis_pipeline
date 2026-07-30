#!/usr/bin/env bash
# =============================================================================
# Step 04 — Assembly QC: CheckV + coverage uniformity  [SLURM array]
# =============================================================================
# Validates that each polished assembly is not just present but CORRECT.
# Two independent checks per sample:
#
#   1. CheckV (if CHECKV_DB is set in config.sh)
#        completeness %, contamination %, CheckV quality tier, and any
#        host-derived regions / provirus boundaries. Flags residual host
#        contamination and over/under-assembled genomes.
#
#   2. Coverage uniformity (always)
#        re-maps the filtered ONT reads to the polished assembly (minimap2),
#        then computes per-base depth statistics. A uniform genome shows a
#        flat profile; localised drops or spikes betray a misassembly
#        (chimeric junction, collapsed repeat). Reported: mean depth, depth CV
#        (sd/mean), % of positions below 0.5× mean (low-cov windows), breadth.
#
# Input:  03_polished/<bc>_<sample>/polished.fasta
#         01_qc/<bc>_<sample>/filtered.fastq.gz
# Output: 04_assembly_qc/<bc>_<sample>/
#           checkv/quality_summary.tsv   (full CheckV output, if run)
#           coverage_per_base.tsv.gz     (samtools depth -a)
#           asm_qc_summary.tsv           (one-line parsed summary for step 09/10)
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/config.sh"

activate_env
require_tool minimap2
load_samtools_cluster

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

step_banner "04 — Assembly QC  [${BARCODE} / ${SAMPLE_NAME}]"

POLISHED="${RESULTS_DIR}/03_polished/${LABEL}/polished.fasta"
ONT_READS="${RESULTS_DIR}/01_qc/${LABEL}/filtered.fastq.gz"
OUT="${RESULTS_DIR}/04_assembly_qc/${LABEL}"
mkdir -p "$OUT"

# ── Input guard: skip cleanly if the polished assembly is missing ────────────
if [[ ! -s "$POLISHED" ]]; then
    echo "[SKIP] ${LABEL}: polished.fasta missing — step 03 did not complete. Nothing to do."
    exit 0
fi

SUMMARY="${OUT}/asm_qc_summary.tsv"
if [[ -f "$SUMMARY" && "${FORCE:-false}" != "true" ]]; then
    log "  Assembly QC summary exists — skipping (set FORCE=true to rerun)"
    exit 0
fi

n_contigs=$(grep -c '^>' "$POLISHED" || echo 0)
total_bp=$(grep -v '^>' "$POLISHED" | tr -d '\n' | wc -c)

# =============================================================================
# 1. CheckV (optional — needs CHECKV_DB)
# =============================================================================
checkv_completeness="NA"; checkv_contamination="NA"; checkv_quality="NA"; checkv_warnings="NA"
if [[ -n "${CHECKV_DB:-}" ]] && command -v checkv &>/dev/null; then
    # CheckV expects the versioned db dir (checkv-db-vX.Y); accept either the
    # parent or the versioned dir itself.
    CKDB="${CHECKV_DB}"
    if [[ ! -f "${CKDB}/genome_db/checkv_reps.dmnd" ]]; then
        _v=$(compgen -G "${CHECKV_DB}/checkv-db-*" 2>/dev/null | head -1 || true)
        [[ -n "$_v" ]] && CKDB="$_v"
    fi
    if [[ -f "${CKDB}/genome_db/checkv_reps.dmnd" ]]; then
        log "  Running CheckV (db: ${CKDB})..."
        CK_LOCAL="${TMPDIR:-/tmp}/checkv_${LABEL}_$$"
        rm -rf "$CK_LOCAL"
        set +e
        checkv end_to_end "$POLISHED" "$CK_LOCAL" -d "$CKDB" -t "${THREADS}" \
            2>&1 | tee "${OUT}/checkv.log"
        _ckexit=${PIPESTATUS[0]}
        set -e
        if [[ ${_ckexit} -eq 0 && -f "${CK_LOCAL}/quality_summary.tsv" ]]; then
            mkdir -p "${OUT}/checkv"
            cp "${CK_LOCAL}"/*.tsv "${OUT}/checkv/" 2>/dev/null || true
            # Aggregate across contigs: mean completeness, summed contamination-ish,
            # worst quality tier. quality_summary.tsv columns include:
            # contig_id, contig_length, gene_count, viral_genes, host_genes,
            # checkv_quality, completeness, contamination, warnings
            read -r checkv_completeness checkv_contamination checkv_quality checkv_warnings < <(
                python3 - "${CK_LOCAL}/quality_summary.tsv" <<'PYEOF'
import sys, csv
rows=list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
def fnum(x):
    try: return float(x)
    except: return None
comps=[fnum(r.get('completeness','')) for r in rows]; comps=[c for c in comps if c is not None]
conts=[fnum(r.get('contamination','')) for r in rows]; conts=[c for c in conts if c is not None]
tiers=[r.get('checkv_quality','') for r in rows if r.get('checkv_quality','')]
order=['Complete','High-quality','Medium-quality','Low-quality','Not-determined']
worst='NA'
for t in reversed(order):
    if t in tiers: worst=t; break
host=sum(int(r.get('host_genes','0') or 0) for r in rows)
comp = f"{max(comps):.1f}" if comps else "NA"
cont = f"{max(conts):.1f}" if conts else "NA"
warn = f"host_genes={host}"
print(comp, cont, worst, warn)
PYEOF
            )
            log "  CheckV: completeness=${checkv_completeness}%  contamination=${checkv_contamination}%  quality=${checkv_quality}  (${checkv_warnings})"
        else
            warn "  CheckV failed (exit=${_ckexit}) — see ${OUT}/checkv.log. Continuing with coverage QC only."
        fi
        rm -rf "$CK_LOCAL"
    else
        warn "  CHECKV_DB set but no usable db found under ${CHECKV_DB} — skipping CheckV."
    fi
else
    log "  CheckV skipped (CHECKV_DB not set or checkv not installed)."
fi

# =============================================================================
# 2. Coverage uniformity (always)
# =============================================================================
mean_depth="NA"; depth_cv="NA"; pct_low_cov="NA"; pct_covered="NA"
if [[ -s "$ONT_READS" ]]; then
    log "  Remapping ONT reads for coverage uniformity..."
    QC_LOCAL="${TMPDIR:-/tmp}/asmqc_${LABEL}_$$"
    rm -rf "$QC_LOCAL"; mkdir -p "$QC_LOCAL"
    set +e
    minimap2 -ax map-ont -t "${THREADS}" "$POLISHED" "$ONT_READS" 2>"${QC_LOCAL}/minimap2.log" \
        | samtools sort -@ "${THREADS}" -o "${QC_LOCAL}/aln.bam" 2>>"${QC_LOCAL}/minimap2.log"
    _mmexit=$?
    samtools index "${QC_LOCAL}/aln.bam" 2>>"${QC_LOCAL}/minimap2.log"
    set -e
    if [[ ${_mmexit} -eq 0 && -s "${QC_LOCAL}/aln.bam" ]]; then
        samtools depth -a "${QC_LOCAL}/aln.bam" > "${QC_LOCAL}/depth.tsv" 2>/dev/null || true
        if [[ -s "${QC_LOCAL}/depth.tsv" ]]; then
            read -r mean_depth depth_cv pct_low_cov pct_covered < <(
                python3 - "${QC_LOCAL}/depth.tsv" <<'PYEOF'
import sys, statistics
d=[]
for line in open(sys.argv[1]):
    p=line.rstrip('\n').split('\t')
    if len(p)>=3:
        try: d.append(int(p[2]))
        except: pass
if not d:
    print("NA NA NA NA"); sys.exit()
n=len(d); mean=sum(d)/n
sd=statistics.pstdev(d) if n>1 else 0.0
cv=(sd/mean) if mean>0 else 0.0
low=sum(1 for x in d if x < 0.5*mean)/n*100.0
covered=sum(1 for x in d if x>0)/n*100.0
print(f"{mean:.1f} {cv:.3f} {low:.2f} {covered:.2f}")
PYEOF
            )
            gzip -c "${QC_LOCAL}/depth.tsv" > "${OUT}/coverage_per_base.tsv.gz" 2>/dev/null || true
            log "  Coverage: mean=${mean_depth}×  CV=${depth_cv}  low-cov(<0.5×mean)=${pct_low_cov}%  breadth=${pct_covered}%"
            # Heuristic misassembly warning
            awk -v cv="${depth_cv}" -v low="${pct_low_cov}" 'BEGIN{
                if (cv+0 > 0.5)  print "[WARN]  Depth CV > 0.5 — non-uniform coverage, inspect for misassembly.";
                if (low+0 > 5)   print "[WARN]  >5% of positions below half-mean depth — possible chimeric junction.";
            }'
        fi
    else
        warn "  Remapping failed (exit=${_mmexit}) — see ${QC_LOCAL}/minimap2.log. Coverage QC unavailable."
    fi
    rm -rf "$QC_LOCAL"
else
    warn "  Filtered ONT reads missing — coverage QC skipped."
fi

# =============================================================================
# 3. One-line summary for step 09 (summary) and step 10 (final report)
# =============================================================================
{
    printf "label\tn_contigs\ttotal_bp\tcheckv_completeness\tcheckv_contamination\tcheckv_quality\tcheckv_warnings\tmean_depth\tdepth_cv\tpct_low_cov\tpct_covered\n"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$LABEL" "$n_contigs" "$total_bp" \
        "$checkv_completeness" "$checkv_contamination" "$checkv_quality" "$checkv_warnings" \
        "$mean_depth" "$depth_cv" "$pct_low_cov" "$pct_covered"
} > "$SUMMARY"

log "  Assembly QC summary: ${SUMMARY}"
log "  Step 04 done."
