#!/usr/bin/env bash
# =============================================================================
# Step 10 — Final HTML report  [single global job]
# =============================================================================
# Produces:
#   A) Per-sample report  10_final_report/<barcode>_<sample>_report.html
#        1. Assembly & QC      reads, coverage, contigs, N50
#        2. Assembly validation CheckV (completeness/contamination/quality) +
#                              coverage uniformity (mean depth, CV, %low, breadth)
#        3. Annotation         Pharokka CDS / hypothetical, links to gbk/gff
#        4. Methylation        6mA / 5mC positions (step 06), MicrobeMod report link
#   B) Global report      10_final_report/run_final_summary.html
#        one row per sample + links to per-sample reports and the step-09 summary
#
# Data sources:
#   01_qc/<l>/coverage_estimate.txt | illumina_stats.txt
#   03_polished/<l>/polished.fasta
#   04_assembly_qc/<l>/asm_qc_summary.tsv
#   05_annotated/<l>/<s>.gbk | <s>.gff | <s>_cds_final_merged_output.tsv
#   06_modbase/<l>/methylation.bed
#   07_methyl/methylation_report.html
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/config.sh"

activate_env
step_banner "10 — Final HTML report"

SAMPLE_SHEET="${RESULTS_DIR}/sample_sheet.tsv"
[[ -f "${SAMPLE_SHEET}" ]] || die "Sample sheet not found: ${SAMPLE_SHEET}"
OUT_DIR="${RESULTS_DIR}/10_final_report"
mkdir -p "${OUT_DIR}"
GEN=$(date '+%Y-%m-%d %H:%M:%S')

getkey() { grep -m1 "^$2=" "$1" 2>/dev/null | cut -d= -f2- || true; }
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

CSS='<style>
 body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:24px;color:#222;background:#f6f7f9}
 h1{margin:0 0 4px} h2{color:#2c3e50;border-bottom:2px solid #e2e6ea;padding-bottom:4px;margin-top:1.5em}
 .meta{color:#666;font-size:13px} .box{background:#fff;border:1px solid #dde;border-radius:8px;padding:12px 16px;margin:10px 0}
 table{border-collapse:collapse;background:#fff;font-size:13px;margin:6px 0} th,td{border:1px solid #e3e7ea;padding:5px 9px}
 th{background:#2c3e50;color:#fff;text-align:left} .good{color:#27ae60;font-weight:600}.warn{color:#e67e22;font-weight:600}.bad{color:#c0392b;font-weight:600}
 code{background:#eef;padding:1px 4px;border-radius:3px} a{color:#2563c9} .kv td:first-child{color:#555;font-weight:600}
 img{max-width:680px;border:1px solid #dde;border-radius:6px;margin-top:8px}
</style>'

mapfile -t BARCODES < <(list_barcodes)

# =============================================================================
# A) Per-sample reports
# =============================================================================
for BARCODE in "${BARCODES[@]}"; do
    SAMPLE_NAME=$(get_sample_field "$BARCODE" "sample_name")
    MODE=$(get_sample_field "$BARCODE" "mode")
    LABEL="${BARCODE}_${SAMPLE_NAME}"
    RPT="${OUT_DIR}/${LABEL}_report.html"

    QCF="${RESULTS_DIR}/01_qc/${LABEL}/coverage_estimate.txt"
    ILF="${RESULTS_DIR}/01_qc/${LABEL}/illumina_stats.txt"
    POL="${RESULTS_DIR}/03_polished/${LABEL}/polished.fasta"
    AQC="${RESULTS_DIR}/04_assembly_qc/${LABEL}/asm_qc_summary.tsv"
    GFF="${RESULTS_DIR}/05_annotated/${LABEL}/${SAMPLE_NAME}.gff"
    GBK="${RESULTS_DIR}/05_annotated/${LABEL}/${SAMPLE_NAME}.gbk"
    CDS_TSV="${RESULTS_DIR}/05_annotated/${LABEL}/${SAMPLE_NAME}_cds_final_merged_output.tsv"
    BED="${RESULTS_DIR}/06_modbase/${LABEL}/methylation.bed"

    raw=$(getkey "$QCF" raw_reads); filt=$(getkey "$QCF" filtered_reads)
    cov=$(getkey "$QCF" filtered_coverage_est); ontlen=$(getkey "$QCF" ont_raw_mean_len)
    ilr=$(getkey "$ILF" illumina_total_est); ilm=$(getkey "$ILF" illumina_mean_len)
    if [[ -s "$POL" ]]; then nctg=$(grep -c '^>' "$POL"); tot=$(grep -v '^>' "$POL" | tr -d '\n' | wc -c); n50=$(compute_n50 "$POL"); else nctg="?"; tot="?"; n50="?"; fi
    ck_comp="?"; ck_cont="?"; ck_qual="?"; md="?"; cv="?"; low="?"; breadth="?"
    [[ -f "$AQC" ]] && IFS=$'\t' read -r _l _nc _tb ck_comp ck_cont ck_qual _w md cv low breadth < <(tail -n +2 "$AQC" | head -1)
    cds="?"; hyp="?"
    [[ -f "$GFF" ]] && cds=$(grep -c $'\tCDS\t' "$GFF" 2>/dev/null || echo 0)
    [[ -f "$CDS_TSV" ]] && hyp=$(grep -ciE 'hypothetical|unknown function' "$CDS_TSV" 2>/dev/null || echo 0)
    m6a="-"; m5c="-"
    if [[ -f "$BED" ]]; then m6a=$(awk '$4=="a"{n++}END{print n+0}' "$BED"); m5c=$(awk '$4=="m"{n++}END{print n+0}' "$BED"); fi

    {
    echo "<!DOCTYPE html><html lang='en'><head><meta charset='UTF-8'><title>${LABEL} — report</title>${CSS}</head><body>"
    echo "<h1>${SAMPLE_NAME}</h1><div class='meta'>${BARCODE} · mode ${MODE} · run ${RUN_NAME} · ${GEN}</div>"

    echo "<h2>1. Assembly &amp; QC</h2><div class='box'><table class='kv'>"
    echo "<tr><td>Raw reads</td><td>${raw:-?}</td></tr><tr><td>Filtered reads</td><td>${filt:-?}</td></tr>"
    echo "<tr><td>ONT coverage est.</td><td>${cov:-?}×</td></tr><tr><td>ONT mean read len</td><td>${ontlen:-?} bp</td></tr>"
    [[ "${ilr:-}" != "" ]] && echo "<tr><td>Illumina reads (est.)</td><td>${ilr} (mean ${ilm} bp)</td></tr>"
    echo "<tr><td>Contigs</td><td>${nctg}</td></tr><tr><td>Total length</td><td>${tot} bp</td></tr><tr><td>N50</td><td>${n50}</td></tr>"
    echo "</table>"
    [[ -s "$POL" ]] && echo "<p>Assembly: <code>03_polished/${LABEL}/polished.fasta</code></p>"
    echo "</div>"

    echo "<h2>2. Assembly validation</h2><div class='box'>"
    if [[ -f "$AQC" ]]; then
        cvcls="good"; awk "BEGIN{exit !(\"$cv\"+0>0.5)}" 2>/dev/null && cvcls="warn"
        lowcls="good"; awk "BEGIN{exit !(\"$low\"+0>5)}" 2>/dev/null && lowcls="warn"
        compcls="good"; { awk "BEGIN{exit !(\"$ck_comp\"+0<90)}" 2>/dev/null && [[ "$ck_comp" =~ ^[0-9.]+$ ]]; } && compcls="warn"
        echo "<table class='kv'>"
        echo "<tr><td>CheckV completeness</td><td class='${compcls}'>${ck_comp}%</td></tr>"
        echo "<tr><td>CheckV contamination</td><td>${ck_cont}%</td></tr>"
        echo "<tr><td>CheckV quality</td><td>${ck_qual}</td></tr>"
        echo "<tr><td>Mean depth (remap)</td><td>${md}×</td></tr>"
        echo "<tr><td>Depth CV</td><td class='${cvcls}'>${cv}</td></tr>"
        echo "<tr><td>% positions &lt; 0.5× mean</td><td class='${lowcls}'>${low}%</td></tr>"
        echo "<tr><td>Breadth covered</td><td>${breadth}%</td></tr>"
        echo "</table><p class='meta'>High depth-CV or % low-coverage flags a possible misassembly (chimeric junction / collapsed repeat).</p>"
    else
        echo "<p class='meta'>No assembly-QC summary (step 04 not run or skipped).</p>"
    fi
    echo "</div>"

    echo "<h2>3. Annotation (Pharokka)</h2><div class='box'><table class='kv'>"
    echo "<tr><td>CDS</td><td>${cds}</td></tr><tr><td>Hypothetical</td><td>${hyp}</td></tr></table>"
    [[ -f "$GBK" ]] && echo "<p>GenBank: <code>05_annotated/${LABEL}/${SAMPLE_NAME}.gbk</code> · GFF: <code>${SAMPLE_NAME}.gff</code></p>"
    echo "</div>"

    echo "<h2>4. Methylation</h2><div class='box'>"
    if [[ -f "$BED" ]]; then
        echo "<table class='kv'><tr><td>6mA positions</td><td>${m6a}</td></tr><tr><td>5mC positions</td><td>${m5c}</td></tr></table>"
        [[ -f "${RESULTS_DIR}/07_methyl/methylation_report.html" ]] && echo "<p>MicrobeMod report: <code>07_methyl/methylation_report.html</code></p>"
    else
        echo "<p class='meta'>No methylation data (MOD_DEMUX_DIR not set, or step 06 not run).</p>"
    fi
    echo "</div>"
    echo "</body></html>"
    } > "$RPT"
    log "  Wrote ${RPT}"
done

# =============================================================================
# B) Global report
# =============================================================================
GLOBAL="${OUT_DIR}/run_final_summary.html"
{
echo "<!DOCTYPE html><html lang='en'><head><meta charset='UTF-8'><title>Final summary — ${RUN_NAME}</title>${CSS}</head><body>"
echo "<h1>Phage pipeline — final report</h1><div class='meta'>Run ${RUN_NAME} · ${GEN}</div>"
echo "<p class='box'>Run-level metrics table: <a href='../09_summary/run_summary.html'>09_summary/run_summary.html</a> (machine-readable <code>run_summary.tsv</code>).</p>"
echo "<h2>Samples</h2><table><tr><th>Barcode</th><th>Sample</th><th>Mode</th><th>Contigs</th><th>Total bp</th><th>CheckV</th><th>Report</th></tr>"
for BARCODE in "${BARCODES[@]}"; do
    SAMPLE_NAME=$(get_sample_field "$BARCODE" "sample_name")
    MODE=$(get_sample_field "$BARCODE" "mode")
    LABEL="${BARCODE}_${SAMPLE_NAME}"
    POL="${RESULTS_DIR}/03_polished/${LABEL}/polished.fasta"
    AQC="${RESULTS_DIR}/04_assembly_qc/${LABEL}/asm_qc_summary.tsv"
    if [[ -s "$POL" ]]; then nctg=$(grep -c '^>' "$POL"); tot=$(grep -v '^>' "$POL" | tr -d '\n' | wc -c); else nctg="?"; tot="?"; fi
    ckq="?"; [[ -f "$AQC" ]] && ckq=$(tail -n +2 "$AQC" | head -1 | cut -f6)
    echo "<tr><td>${BARCODE}</td><td>${SAMPLE_NAME}</td><td>${MODE}</td><td>${nctg}</td><td>${tot}</td><td>${ckq}</td><td><a href='${LABEL}_report.html'>view</a></td></tr>"
done
echo "</table></body></html>"
} > "$GLOBAL"

log "  Global report: ${GLOBAL}"
log "  Step 10 done."
