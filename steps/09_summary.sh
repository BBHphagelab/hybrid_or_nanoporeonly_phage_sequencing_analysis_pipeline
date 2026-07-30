#!/usr/bin/env bash
# =============================================================================
# Step 09 — Run summary across all samples  (TSV + TXT + HTML)
# =============================================================================
# Gathers metrics from every prior step into a single overview. Modelled on the
# richer bacterial-pipeline summary (cards, methods section, colour-coded table,
# quality flags), adapted for the phage pipeline:
#
#   QC          raw / filtered reads, ONT coverage estimate, ONT mean length,
#               Illumina reads + mean length (if present)
#   Assembly    n contigs, total bp, N50
#   Assembly QC contig completeness/contamination/quality (CheckV, step 04) and
#               coverage uniformity (mean depth, depth CV, % low-cov, breadth)
#   Annotation  Pharokka CDS + hypothetical (step 05)
#   Methylation 6mA / 5mC positions (step 06, if run)
#
# Output: $RESULTS_DIR/09_summary/
#   run_summary.tsv | run_summary.txt | run_summary.html
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/config.sh"

activate_env
step_banner "09 — Summary"

SAMPLE_SHEET="${RESULTS_DIR}/sample_sheet.tsv"
[[ -f "${SAMPLE_SHEET}" ]] || die "Sample sheet not found: ${SAMPLE_SHEET}"
OUT="${RESULTS_DIR}/09_summary"
mkdir -p "$OUT"
TSV="${OUT}/run_summary.tsv"
TXT="${OUT}/run_summary.txt"
HTML="${OUT}/run_summary.html"

# ---------------------------------------------------------------------------
compute_n50() {
    [[ -s "$1" ]] || { echo "?"; return; }
    python3 - "$1" <<'PYEOF'
import sys
lens=[]; cur=0
for line in open(sys.argv[1]):
    if line.startswith('>'):
        if cur: lens.append(cur)
        cur=0
    else: cur+=len(line.strip())
if cur: lens.append(cur)
if not lens: print("?"); sys.exit()
lens.sort(reverse=True); tot=sum(lens); run=0
for l in lens:
    run+=l
    if run>=tot/2: print(l); break
PYEOF
}
getkey() { grep -m1 "^$2=" "$1" 2>/dev/null | cut -d= -f2- || true; }
field()  { local v; v=$(grep -iE "$2" "$1" 2>/dev/null | grep -oE '[0-9]+' | head -1); echo "${v:-0}"; }

# ---------------------------------------------------------------------------
# TSV header
# ---------------------------------------------------------------------------
{
printf "barcode\tsample_name\tmode\t"
printf "raw_reads\tfiltered_reads\tcoverage_est\tont_mean_len\t"
printf "illumina_reads\tillumina_mean_len\t"
printf "n_contigs\ttotal_bp\tN50\t"
printf "checkv_completeness\tcheckv_contamination\tcheckv_quality\t"
printf "mean_depth\tdepth_cv\tpct_low_cov\tpct_covered\t"
printf "pharokka_cds\tpharokka_hypothetical\t"
printf "6mA_sites\t5mC_sites\n"
} > "$TSV"

mapfile -t BARCODES < <(list_barcodes)
N_SAMPLES=${#BARCODES[@]}; N_ONT=0; N_HYBRID=0

for BARCODE in "${BARCODES[@]}"; do
    SAMPLE_NAME=$(get_sample_field "$BARCODE" "sample_name")
    MODE=$(get_sample_field "$BARCODE" "mode")
    LABEL="${BARCODE}_${SAMPLE_NAME}"
    [[ "$MODE" == "hybrid" ]] && N_HYBRID=$((N_HYBRID+1)) || N_ONT=$((N_ONT+1))

    # ── QC ────────────────────────────────────────────────────────────────
    QCF="${RESULTS_DIR}/01_qc/${LABEL}/coverage_estimate.txt"
    raw_reads=$(getkey "$QCF" raw_reads);            raw_reads="${raw_reads:-?}"
    filtered_reads=$(getkey "$QCF" filtered_reads);  filtered_reads="${filtered_reads:-?}"
    coverage_est=$(getkey "$QCF" filtered_coverage_est); coverage_est="${coverage_est:-?}"
    ont_mean_len=$(getkey "$QCF" ont_raw_mean_len);  ont_mean_len="${ont_mean_len:-?}"
    ILF="${RESULTS_DIR}/01_qc/${LABEL}/illumina_stats.txt"
    illumina_reads=$(getkey "$ILF" illumina_total_est); illumina_reads="${illumina_reads:--}"
    illumina_mean_len=$(getkey "$ILF" illumina_mean_len); illumina_mean_len="${illumina_mean_len:--}"

    # ── Assembly ──────────────────────────────────────────────────────────
    POL="${RESULTS_DIR}/03_polished/${LABEL}/polished.fasta"
    if [[ -s "$POL" ]]; then
        n_contigs=$(grep -c '^>' "$POL"); total_bp=$(grep -v '^>' "$POL" | tr -d '\n' | wc -c)
        N50=$(compute_n50 "$POL")
    else
        n_contigs="?"; total_bp="?"; N50="?"
    fi

    # ── Assembly QC (CheckV + coverage) ───────────────────────────────────
    AQC="${RESULTS_DIR}/04_assembly_qc/${LABEL}/asm_qc_summary.tsv"
    ck_comp="?"; ck_cont="?"; ck_qual="?"; mean_depth="?"; depth_cv="?"; pct_low="?"; pct_cov="?"
    if [[ -f "$AQC" ]]; then
        IFS=$'\t' read -r _l _nc _tb ck_comp ck_cont ck_qual _ckw mean_depth depth_cv pct_low pct_cov \
            < <(tail -n +2 "$AQC" | head -1)
    fi

    # ── Annotation ────────────────────────────────────────────────────────
    GFF="${RESULTS_DIR}/05_annotated/${LABEL}/${SAMPLE_NAME}.gff"
    CDS_TSV="${RESULTS_DIR}/05_annotated/${LABEL}/${SAMPLE_NAME}_cds_final_merged_output.tsv"
    pharokka_cds="?"; pharokka_hyp="?"
    [[ -f "$GFF" ]] && pharokka_cds=$(grep -c $'\tCDS\t' "$GFF" 2>/dev/null || echo 0)
    if [[ -f "$CDS_TSV" ]]; then
        pharokka_hyp=$(grep -ciE 'hypothetical|unknown function' "$CDS_TSV" 2>/dev/null || echo 0)
    fi

    # ── Methylation ───────────────────────────────────────────────────────
    BED="${RESULTS_DIR}/06_modbase/${LABEL}/methylation.bed"
    m6A="-"; m5C="-"
    if [[ -f "$BED" ]]; then
        m6A=$(awk '$4=="a"{n++} END{print n+0}' "$BED")
        m5C=$(awk '$4=="m"{n++} END{print n+0}' "$BED")
    fi

    _s() { printf '%s' "${1:-?}"; }
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$(_s "$BARCODE")" "$(_s "$SAMPLE_NAME")" "$(_s "$MODE")" \
        "$(_s "$raw_reads")" "$(_s "$filtered_reads")" "$(_s "$coverage_est")" "$(_s "$ont_mean_len")" \
        "$(_s "$illumina_reads")" "$(_s "$illumina_mean_len")" \
        "$(_s "$n_contigs")" "$(_s "$total_bp")" "$(_s "$N50")" \
        "$(_s "$ck_comp")" "$(_s "$ck_cont")" "$(_s "$ck_qual")" \
        "$(_s "$mean_depth")" "$(_s "$depth_cv")" "$(_s "$pct_low")" "$(_s "$pct_cov")" \
        "$(_s "$pharokka_cds")" "$(_s "$pharokka_hyp")" \
        "$(_s "$m6A")" "$(_s "$m5C")" \
        >> "$TSV"
done

# ---------------------------------------------------------------------------
# Plain-text report
# ---------------------------------------------------------------------------
{
    echo "Nanopore phage pipeline — run summary"
    echo "Run: ${RUN_NAME}    Generated: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Samples: ${N_SAMPLES}  (ont: ${N_ONT}, hybrid: ${N_HYBRID})"
    echo ""
    column -t -s $'\t' "$TSV"
} > "$TXT"

# ---------------------------------------------------------------------------
# HTML report (cards + methods + colour-coded table with QC flags)
# ---------------------------------------------------------------------------
{
cat <<HEAD
<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">
<title>Phage pipeline summary — ${RUN_NAME}</title>
<style>
 body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:24px;color:#222;background:#f6f7f9}
 h1{margin:0 0 4px} h2{color:#2c3e50;border-bottom:2px solid #e2e6ea;padding-bottom:4px;margin-top:1.6em}
 .meta{color:#666;font-size:13px;margin-bottom:1em}
 .cards{display:flex;gap:14px;flex-wrap:wrap;margin:1em 0}
 .card{background:#fff;border:1px solid #dde;border-radius:8px;padding:14px 20px;min-width:120px;text-align:center}
 .card .v{font-size:26px;font-weight:700;color:#2c3e50} .card .l{font-size:12px;color:#777}
 .methods{background:#fff;border:1px solid #dde;border-radius:8px;padding:14px 18px;font-size:12.5px;line-height:1.7}
 table{border-collapse:collapse;width:100%;background:#fff;font-size:12.5px}
 th,td{border:1px solid #e3e7ea;padding:5px 8px;text-align:right;white-space:nowrap}
 th{background:#2c3e50;color:#fff;position:sticky;top:0} td:first-child,td:nth-child(2),th:first-child,th:nth-child(2){text-align:left}
 tr:nth-child(even){background:#fafbfc}
 .good{color:#27ae60;font-weight:600}.warn{color:#e67e22;font-weight:600}.bad{color:#c0392b;font-weight:600}
 .flag{color:#c0392b;font-size:11px}
 .tbl-wrap{overflow-x:auto;border-radius:8px;border:1px solid #dde}
 code{background:#eef;padding:1px 4px;border-radius:3px}
</style></head><body>
<h1>Phage pipeline — run summary</h1>
<div class="meta"><b>Run:</b> ${RUN_NAME} &nbsp;|&nbsp; <b>Generated:</b> $(date '+%Y-%m-%d %H:%M:%S') &nbsp;|&nbsp; <b>Results:</b> <code>${RESULTS_DIR}</code></div>
<div class="cards">
 <div class="card"><div class="v">${N_SAMPLES}</div><div class="l">Samples</div></div>
 <div class="card"><div class="v">${N_ONT}</div><div class="l">ONT</div></div>
 <div class="card"><div class="v">${N_HYBRID}</div><div class="l">Hybrid polish</div></div>
</div>

<h2>Methods</h2>
<div class="methods">
<b>QC:</b> <code>chopper</code> (Q≥${MIN_READ_QUALITY}, len≥${MIN_READ_LENGTH} bp) + NanoStat.<br>
<b>Assembly:</b> <code>Autocycler</code> multi-assembler consensus (default); fallback <code>Flye --nano-hq</code> on Filtlong-capped reads (~${FILTLONG_TARGET_COV}×). Genome ~${GENOME_SIZE_ESTIMATE}; structure is ONT-defined.<br>
<b>Polishing:</b> <code>Medaka</code> (${MEDAKA_MODEL}); hybrid samples then <code>Polypolish --careful</code> + <code>Pypolca --careful</code> (Illumina, unique regions only).<br>
<b>Reorientation:</b> <code>PhageTerm</code> (fallback <code>dnaapler</code> terL), contigs renamed <code>&lt;sample&gt;_ctgNNN</code>.<br>
<b>Assembly QC:</b> <code>CheckV</code> (completeness/contamination/host) + read-remap coverage uniformity (mean depth, CV, % low-coverage, breadth).<br>
<b>Annotation:</b> <code>Pharokka</code> (${PHAROKKA_GENE_PREDICTOR}).<br>
<b>Methylation:</b> MinKNOW mod-calls → <code>modkit</code> + <code>MicrobeMod</code>.<br>
<b>QC flags:</b> [!contigs&gt;1] fragmented &nbsp; [!CV&gt;0.5] non-uniform coverage &nbsp; [!low&gt;5%] coverage dropouts &nbsp; [!compl&lt;90] CheckV incomplete &nbsp; [!cont&gt;5] CheckV contamination.
</div>

<h2>Per-sample metrics</h2>
<div class="tbl-wrap"><table>
<tr><th>barcode</th><th>sample</th><th>mode</th><th>filt.reads</th><th>cov×</th><th>contigs</th><th>total bp</th><th>N50</th>
<th>CheckV %compl</th><th>%cont</th><th>quality</th><th>mean dp</th><th>CV</th><th>%low</th><th>%breadth</th><th>CDS</th><th>6mA</th><th>5mC</th>
<th>flags</th></tr>
HEAD

# Body rows from TSV
tail -n +2 "$TSV" | while IFS=$'\t' read -r bc sample mode raw filt cov ontlen ilr ilm nctg tot n50 ckc ckct ckq md cv low cov2 cds hyp m6a m5c; do
    flags=""
    [[ "$nctg" =~ ^[0-9]+$ && "$nctg" -gt 1 ]] && flags="${flags}[!contigs>1] "
    awk "BEGIN{exit !(\"$cv\"+0>0.5 && \"$cv\"!=\"?\")}" 2>/dev/null && [[ "$cv" != "?" && "$cv" != "NA" ]] && flags="${flags}[!CV>0.5] "
    awk "BEGIN{exit !(\"$low\"+0>5 && \"$low\"!=\"?\")}" 2>/dev/null && [[ "$low" != "?" && "$low" != "NA" ]] && flags="${flags}[!low>5%] "
    awk "BEGIN{exit !(\"$ckc\"+0<90 && \"$ckc\"!=\"?\" && \"$ckc\"!=\"NA\")}" 2>/dev/null && [[ "$ckc" =~ ^[0-9.]+$ ]] && flags="${flags}[!compl<90] "
    awk "BEGIN{exit !(\"$ckct\"+0>5 && \"$ckct\"!=\"?\" && \"$ckct\"!=\"NA\")}" 2>/dev/null && [[ "$ckct" =~ ^[0-9.]+$ ]] && flags="${flags}[!cont>5] "
    flag_html=""; [[ -n "$flags" ]] && flag_html="<span class='flag'>${flags}</span>"
    printf '<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\n' \
        "$bc" "$sample" "$mode" "$filt" "$cov" "$nctg" "$tot" "$n50" \
        "$ckc" "$ckct" "$ckq" "$md" "$cv" "$low" "$cov2" "$cds" "$m6a" "$m5c" \
        "$flag_html"
done

echo "</table></div>"
echo "<p class='meta'>Full machine-readable table: <code>run_summary.tsv</code></p>"
echo "</body></html>"
} > "$HTML"

log "  Summary written:"
log "    TSV : ${TSV}"
log "    TXT : ${TXT}"
log "    HTML: ${HTML}"
log "  Step 09 done."
