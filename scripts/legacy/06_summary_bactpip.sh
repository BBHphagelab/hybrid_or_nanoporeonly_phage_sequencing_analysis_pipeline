#!/usr/bin/env bash
# =============================================================================
# Step 06 — Summary report across all samples
# =============================================================================
# Collecte les stats de tous les steps précédents et produit :
#   - run_summary.tsv   : tableau tabulé (toutes colonnes)
#   - run_summary.txt   : rapport texte lisible
#   - run_summary.html  : rapport HTML interactif
#
# Colonnes TSV (33) :
#   barcode | sample_name | mode
#   raw_reads | filtered_reads | coverage_est
#   ont_raw_mean_len | illumina_reads | illumina_mean_len
#   n_contigs | total_bp | N50 | circular
#   uni_long_bridges | uni_short_bridges
#   polypolish_changes | pypolca_sub_err | pypolca_indel_err | pypolca_qv_before | pypolca_qv_after
#   bakta_cds | bakta_hypothetical
#   compare_to | snp_count | snp_near_edge | snippy_mean_cov | snippy_pct_cov
#   breseq_int_mutations | breseq_ext_mutations
#   snippy_ext_snps | snippy_ext_near_edge | snippy_ext_mean_cov | snippy_ext_pct_cov
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/config.sh"

activate_env

step_banner "06 — Summary"

SAMPLE_SHEET="${RESULTS_DIR}/sample_sheet.tsv"
OUT="${RESULTS_DIR}/06_summary"
mkdir -p "$OUT"

SUMMARY_TSV="${OUT}/run_summary.tsv"
SUMMARY_TXT="${OUT}/run_summary.txt"

# ---------------------------------------------------------------------------
# Helper: compute N50 from a FASTA file
# ---------------------------------------------------------------------------
compute_n50() {
    local fasta="$1"
    python3 - "$fasta" <<'EOF'
import sys
lens = []
with open(sys.argv[1]) as f:
    cur = 0
    for line in f:
        if line.startswith('>'):
            if cur: lens.append(cur)
            cur = 0
        else:
            cur += len(line.strip())
    if cur: lens.append(cur)
lens.sort(reverse=True)
total = sum(lens)
running = 0
for l in lens:
    running += l
    if running >= total / 2:
        print(l)
        break
EOF
}

# ---------------------------------------------------------------------------
# Helper: parse pypolca.report
# ---------------------------------------------------------------------------
parse_pypolca_report() {
    local report="$1"
    local sub="?" indel="?" qv_before="?" qv_after="?"
    if [[ -f "$report" ]]; then
        sub=$(grep -iE 'Substitution.*(error|corrected)' "$report" \
              | grep -oP '\d+' | tail -1 || true)
        sub="${sub:-?}"
        indel=$(grep -iE 'Insertion.*(error|corrected|deletion)' "$report" \
                | grep -oP '\d+' | tail -1 || true)
        indel="${indel:-?}"
        qv_before=$(grep -iE 'QV.before' "$report" \
                    | grep -oP '\d+(\.\d+)?' | tail -1 || true)
        if [[ -z "$qv_before" ]]; then
            qv_before=$(grep -iE 'Consensus Quality Before' "$report" \
                        | grep -oP '\d+(\.\d+)?' | tail -1 || true)
        fi
        qv_before="${qv_before:-?}"
        qv_after=$(grep -iE 'QV.after' "$report" \
                   | grep -oP '\d+(\.\d+)?' | tail -1 || true)
        qv_after="${qv_after:-}"
        if [[ -z "$qv_after" ]]; then
            local asm_size
            asm_size=$(grep -iE 'Assembly Size' "$report" \
                       | grep -oP '\d+' | tail -1 || true)
            if [[ "$sub" =~ ^[0-9]+$ && "$indel" =~ ^[0-9]+$ && "$asm_size" =~ ^[0-9]+$ ]]; then
                local total_err=$(( sub + indel ))
                if [[ "$total_err" -eq 0 ]]; then
                    qv_after="${qv_before}"
                else
                    qv_after=$(python3 -c "
import math
s = int('${asm_size}')
print('>{:.1f}'.format(-10 * math.log10(1.0 / s)))
" 2>/dev/null || echo "?")
                fi
            else
                qv_after="?"
            fi
        fi
    fi
    sub="${sub//[$'\n\r\t']/}"
    indel="${indel//[$'\n\r\t']/}"
    qv_before="${qv_before//[$'\n\r\t']/}"
    qv_after="${qv_after//[$'\n\r\t']/}"
    echo "${sub}|${indel}|${qv_before}|${qv_after}"
}

# ---------------------------------------------------------------------------
# Helper: count bases changed by polypolish
# ---------------------------------------------------------------------------
count_fasta_changes() {
    local fasta_before="$1" fasta_after="$2" log="${3:-}"
    if [[ -n "$log" && -f "$log" ]]; then
        local n
        n=$(grep -oP '\d+(?=\s+changes?)' "$log" 2>/dev/null \
            | awk '{s+=$1} END{print (NR>0 ? s : "")}')
        if [[ -n "$n" ]]; then echo "$n"; return; fi
    fi
    if [[ -f "$fasta_before" && -f "$fasta_after" ]]; then
        python3 - "$fasta_before" "$fasta_after" <<'PYEOF'
import sys
def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name: seqs[name] = ''.join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
    if name: seqs[name] = ''.join(buf)
    return seqs
s1, s2 = read_fasta(sys.argv[1]), read_fasta(sys.argv[2])
total = sum(
    sum(a != b for a, b in zip(s1[n], s2[n])) + abs(len(s1[n]) - len(s2[n]))
    for n in s1 if n in s2
)
print(total)
PYEOF
        return
    fi
    echo "?"
}

# ---------------------------------------------------------------------------
# Helper: parse Unicycler log for bridge statistics
# ---------------------------------------------------------------------------
parse_unicycler_bridges() {
    local log="$1"
    local long_bridges="0" short_bridges="0"
    if [[ ! -f "$log" ]]; then echo "?|?"; return; fi
    local bridge_rows
    bridge_rows=$(awk '
        /Applying bridges \(/ { in_sec=1; next }
        in_sec && /Saving.*\.gfa/  { exit }
        in_sec && / -> /           { print }
    ' "$log" 2>/dev/null)
    if [[ -n "$bridge_rows" ]]; then
        long_bridges=$(echo "$bridge_rows" \
            | grep -ciE 'long_read|miniasm|racon' 2>/dev/null || echo "0")
        short_bridges=$(echo "$bridge_rows" \
            | grep -ciE 'spades|loop_unroll|short_read' 2>/dev/null || echo "0")
    fi
    long_bridges="${long_bridges//[$'\n\r\t ']/}"
    short_bridges="${short_bridges//[$'\n\r\t ']/}"
    [[ "$long_bridges"  =~ ^[0-9]+$ ]] || long_bridges="?"
    [[ "$short_bridges" =~ ^[0-9]+$ ]] || short_bridges="?"
    echo "${long_bridges}|${short_bridges}"
}

# ---------------------------------------------------------------------------
# TSV header (29 colonnes)
# ---------------------------------------------------------------------------
{
    printf "barcode\tsample_name\tmode\t"
    printf "raw_reads\tfiltered_reads\tcoverage_est\t"
    printf "ont_raw_mean_len\tillumina_reads\tillumina_mean_len\t"
    printf "n_contigs\ttotal_bp\tN50\tcircular\t"
    printf "uni_long_bridges\tuni_short_bridges\t"
    printf "polypolish_changes\tpypolca_sub_err\tpypolca_indel_err\tpypolca_qv_before\tpypolca_qv_after\t"
    printf "bakta_cds\tbakta_hypothetical\t"
    printf "compare_to\tsnp_count\tsnp_near_edge\tsnippy_mean_cov\tsnippy_pct_cov\t"
    printf "breseq_int_mutations\tbreseq_ext_mutations\t"
    printf "snippy_ext_snps\tsnippy_ext_near_edge\tsnippy_ext_mean_cov\tsnippy_ext_pct_cov\n"
} > "$SUMMARY_TSV"

mapfile -t BARCODES < <(list_barcodes)

for BARCODE in "${BARCODES[@]}"; do
    SAMPLE_NAME=$(get_sample_field "$BARCODE" "sample_name")
    MODE=$(get_sample_field         "$BARCODE" "mode")
    COMPARE_TO=$(get_sample_field     "$BARCODE" "compare_to")
    BRESEQ_REF=$(get_sample_field_opt "$BARCODE" "breseq_ref")
    BRESEQ_EXT=$(get_sample_field_opt "$BARCODE" "breseq_ext_ref")

    SAMPLE_DIR_QC="${RESULTS_DIR}/01_qc/${BARCODE}_${SAMPLE_NAME}"
    SAMPLE_DIR_ASM="${RESULTS_DIR}/02_assembly/${BARCODE}_${SAMPLE_NAME}"
    SAMPLE_DIR_POL="${RESULTS_DIR}/03_polished/${BARCODE}_${SAMPLE_NAME}"
    SAMPLE_DIR_ANN="${RESULTS_DIR}/04_annotated/${BARCODE}_${SAMPLE_NAME}"
    SAMPLE_DIR_VAR="${RESULTS_DIR}/05_variants/${BARCODE}_${SAMPLE_NAME}"

    # ── Stats reads ONT ───────────────────────────────────────────────────
    raw_reads="?"; filtered_reads="?"; coverage_est="?"; ont_raw_mean_len="?"
    cov_file="${SAMPLE_DIR_QC}/coverage_estimate.txt"
    if [[ -f "$cov_file" ]]; then
        raw_reads=$(grep       'raw_reads='             "$cov_file" | cut -d= -f2)
        filtered_reads=$(grep  'filtered_reads='        "$cov_file" | cut -d= -f2)
        coverage_est=$(grep    'filtered_coverage_est=' "$cov_file" | cut -d= -f2)
        ont_raw_mean_len=$(grep 'ont_raw_mean_len='     "$cov_file" | cut -d= -f2)
    fi

    # ── Stats reads Illumina ──────────────────────────────────────────────
    illumina_reads="N/A"; illumina_mean_len="N/A"
    illum_file="${SAMPLE_DIR_QC}/illumina_stats.txt"
    if [[ -f "$illum_file" ]]; then
        # Préférer le comptage réel si disponible, sinon estimation
        illumina_reads=$(grep 'illumina_total_real=' "$illum_file" | cut -d= -f2)
        [[ -z "$illumina_reads" ]] && \
            illumina_reads=$(grep 'illumina_total_est=' "$illum_file" | cut -d= -f2)
        illumina_mean_len=$(grep 'illumina_mean_len=' "$illum_file" | cut -d= -f2)
    fi

    # ── Assembly stats ────────────────────────────────────────────────────
    polished="${SAMPLE_DIR_POL}/polished.fasta"
    n_contigs="?"; total_bp="?"; n50="?"; circular="?"
    if [[ -f "$polished" ]]; then
        n_contigs=$(grep -c '^>' "$polished" || true); n_contigs="${n_contigs:-0}"
        total_bp=$(grep -v '^>' "$polished" | tr -d '\n' | wc -c)
        n50=$(compute_n50 "$polished" 2>/dev/null || echo "?")
        circ_count=$(grep '>' "$polished" | grep -c 'circular' 2>/dev/null || true)
        circ_count="${circ_count:-0}"
        [[ "$circ_count" -gt 0 ]] && circular="${circ_count}x" || circular="linear"
    fi

    # ── Unicycler bridge stats ─────────────────────────────────────────────
    uni_long_bridges="N/A"; uni_short_bridges="N/A"
    if [[ "$MODE" == "hybrid" ]]; then
        uni_log="${SAMPLE_DIR_ASM}/unicycler/unicycler.log"
        bridges=$(parse_unicycler_bridges "$uni_log")
        uni_long_bridges="${bridges%%|*}"
        uni_short_bridges="${bridges##*|}"
    fi

    # ── Polishing stats ───────────────────────────────────────────────────
    polypolish_changes="N/A"
    pypolca_sub="N/A"; pypolca_indel="N/A"; pypolca_qv_before="N/A"; pypolca_qv_after="N/A"
    if [[ "$MODE" == "hybrid" ]]; then
        polypolish_changes=$(count_fasta_changes \
            "${SAMPLE_DIR_ASM}/assembly.fasta" \
            "${SAMPLE_DIR_POL}/polypolish/polypolish_polished.fasta" \
            "${SAMPLE_DIR_POL}/polypolish/polish.log")
        pypolca_report="${SAMPLE_DIR_POL}/pypolca/pypolca.report"
        report_vals=$(parse_pypolca_report "$pypolca_report")
        pypolca_sub="${report_vals%%|*}";       rest="${report_vals#*|}"
        pypolca_indel="${rest%%|*}";            rest="${rest#*|}"
        pypolca_qv_before="${rest%%|*}"; pypolca_qv_after="${rest##*|}"
    elif [[ "$MODE" == "ont" ]]; then
        polypolish_changes="-"; pypolca_sub="-"; pypolca_indel="-"
        pypolca_qv_before="-"; pypolca_qv_after="-"
    fi

    # ── Annotation stats ──────────────────────────────────────────────────
    bakta_cds="?"; bakta_hypo="?"
    gff="${SAMPLE_DIR_ANN}/${SAMPLE_NAME}.gff3"
    bakta_tsv="${SAMPLE_DIR_ANN}/${SAMPLE_NAME}.tsv"
    if [[ -f "$gff" ]]; then
        bakta_cds=$(grep -c $'\tCDS\t' "$gff" 2>/dev/null || true); bakta_cds="${bakta_cds:-0}"
    fi
    if [[ -f "$bakta_tsv" ]]; then
        bakta_hypo=$(grep -ci 'hypothetical' "$bakta_tsv" 2>/dev/null || true); bakta_hypo="${bakta_hypo:-0}"
    fi

    # ── Snippy variant stats ──────────────────────────────────────────────
    snp_count="N/A"; snp_near_edge="N/A"
    snippy_mean_cov="N/A"; snippy_pct_cov="N/A"
    if [[ "${COMPARE_TO}" != "-" && -n "${COMPARE_TO}" ]]; then
        vcf="${SAMPLE_DIR_VAR}/snps.vcf"
        if [[ -f "$vcf" ]]; then
            snp_count=$(awk '!/^#/{n++} END{print n+0}' "$vcf" 2>/dev/null || echo 0)

            # Variants proches des bords de contigs
            ref_fasta="${RESULTS_DIR}/03_polished/${COMPARE_TO}_$(get_sample_field "$COMPARE_TO" "sample_name")/polished.fasta"
            if [[ -f "$ref_fasta" ]]; then
                edge_hits=$(python3 - "$vcf" "$ref_fasta" "${SNIPPY_CONTIG_EDGE_FILTER:-100}" <<'PYEOF'
import sys
vcf_path, fa_path, edge = sys.argv[1], sys.argv[2], int(sys.argv[3])
lengths = {}
cur, name = 0, None
with open(fa_path) as f:
    for line in f:
        if line.startswith('>'):
            if name: lengths[name] = cur
            name = line[1:].split()[0]; cur = 0
        else:
            cur += len(line.strip())
if name: lengths[name] = cur
n = 0
with open(vcf_path) as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.split('\t')
        chrom, pos = parts[0], int(parts[1])
        clen = lengths.get(chrom, 0)
        if clen > 0 and (pos <= edge or pos >= clen - edge):
            n += 1
print(n)
PYEOF
                )
                snp_near_edge="${edge_hits:-?}"
            fi

            # Couverture Snippy
            cov_file_snp="${SAMPLE_DIR_VAR}/coverage_stats.txt"
            if [[ -f "$cov_file_snp" ]]; then
                snippy_mean_cov=$(grep 'mean_depth='  "$cov_file_snp" | cut -d= -f2)
                snippy_pct_cov=$(grep  'pct_covered=' "$cov_file_snp" | cut -d= -f2)
            fi
        else
            snp_count="pending"
        fi
    fi

    # ── Breseq interne ────────────────────────────────────────────────────
    # Cherche d'abord breseq_internal/ (nouveau nom), puis breseq/ (ancien nom,
    # résultats générés par une version précédente du pipeline).
    breseq_int_mutations="N/A"
    if [[ "${BRESEQ_REF}" != "-" && -n "${BRESEQ_REF}" ]]; then
        vcf_int="${SAMPLE_DIR_VAR}/breseq_internal/output/output.vcf"
        vcf_int_old="${SAMPLE_DIR_VAR}/breseq/output/output.vcf"
        if [[ -f "$vcf_int" ]]; then
            breseq_int_mutations=$(grep -vc '^#' "$vcf_int" 2>/dev/null || echo 0)
        elif [[ -f "$vcf_int_old" ]]; then
            # Compatibilité ascendante : ancienne structure breseq/
            breseq_int_mutations=$(grep -vc '^#' "$vcf_int_old" 2>/dev/null || echo 0)
        elif [[ -d "${SAMPLE_DIR_VAR}/breseq_internal" || -d "${SAMPLE_DIR_VAR}/breseq" ]]; then
            breseq_int_mutations="error"
        fi
    # Fallback sans breseq_ref déclaré : si un dossier breseq/ existe (ancien pipeline),
    # on le remonte quand même dans le résumé pour ne pas perdre l'info.
    elif [[ -z "${BRESEQ_REF}" ]]; then
        vcf_int_old="${SAMPLE_DIR_VAR}/breseq/output/output.vcf"
        if [[ -f "$vcf_int_old" ]]; then
            breseq_int_mutations=$(grep -vc '^#' "$vcf_int_old" 2>/dev/null || echo 0)
        elif [[ -d "${SAMPLE_DIR_VAR}/breseq" ]]; then
            breseq_int_mutations="error"
        fi
    fi

    # ── Breseq externe ────────────────────────────────────────────────────
    breseq_ext_mutations="N/A"
    if [[ "${BRESEQ_EXT}" != "-" && -n "${BRESEQ_EXT}" ]]; then
        vcf_ext="${SAMPLE_DIR_VAR}/breseq_external/output/output.vcf"
        if [[ -f "$vcf_ext" ]]; then
            breseq_ext_mutations=$(grep -vc '^#' "$vcf_ext" 2>/dev/null || echo 0)
        elif [[ -d "${SAMPLE_DIR_VAR}/breseq_external" ]]; then
            breseq_ext_mutations="error"
        fi
    fi

    # ── Snippy externe ────────────────────────────────────────────────────
    SNIPPY_EXT=$(get_sample_field_opt "$BARCODE" "snippy_ext_ref")
    snippy_ext_snps="N/A"; snippy_ext_near_edge="N/A"
    snippy_ext_mean_cov="N/A"; snippy_ext_pct_cov="N/A"
    if [[ "${SNIPPY_EXT}" != "-" && -n "${SNIPPY_EXT}" ]]; then
        vcf_sext="${SAMPLE_DIR_VAR}/snippy_external/snps.vcf"
        if [[ -f "$vcf_sext" ]]; then
            snippy_ext_snps=$(awk '!/^#/{n++} END{print n+0}' "$vcf_sext" 2>/dev/null || echo 0)

            # Variants proches des bords — pas de FASTA de référence externe disponible
            # facilement, donc on utilise une estimation simple depuis le VCF (longueur max contig)
            sext_ref_fasta=""
            if [[ "${SNIPPY_EXT}" == "ext" && -n "${EXTERNAL_REF:-}" && -f "${EXTERNAL_REF}" ]]; then
                # Tenter d'extraire une FASTA de la référence externe si c'est un gbk
                sext_ref_fasta="${EXTERNAL_REF}"
            elif [[ "${SNIPPY_EXT}" != "ext" && -f "${SNIPPY_EXT}" ]]; then
                sext_ref_fasta="${SNIPPY_EXT}"
            fi
            if [[ -n "$sext_ref_fasta" && -f "$sext_ref_fasta" ]]; then
                sext_edge_hits=$(python3 - "$vcf_sext" "$sext_ref_fasta" "${SNIPPY_CONTIG_EDGE_FILTER:-100}" <<'PYEOF'
import sys, re
vcf_path, fa_path, edge = sys.argv[1], sys.argv[2], int(sys.argv[3])
lengths = {}
cur, name = 0, None
with open(fa_path) as f:
    for line in f:
        if line.startswith('>'):
            if name: lengths[name] = cur
            name = line[1:].split()[0]; cur = 0
        elif not line.startswith('#'):
            cur += len(line.strip())
if name: lengths[name] = cur
n = 0
with open(vcf_path) as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.split('\t')
        if len(parts) < 2: continue
        chrom, pos = parts[0], int(parts[1])
        clen = lengths.get(chrom, 0)
        if clen > 0 and (pos <= edge or pos >= clen - edge):
            n += 1
print(n)
PYEOF
                )
                snippy_ext_near_edge="${sext_edge_hits:-?}"
            else
                snippy_ext_near_edge="?"
            fi

            cov_file_sext="${SAMPLE_DIR_VAR}/snippy_external/coverage_stats.txt"
            if [[ -f "$cov_file_sext" ]]; then
                snippy_ext_mean_cov=$(grep 'mean_depth='  "$cov_file_sext" | cut -d= -f2)
                snippy_ext_pct_cov=$(grep  'pct_covered=' "$cov_file_sext" | cut -d= -f2)
            fi
        elif [[ -d "${SAMPLE_DIR_VAR}/snippy_external" ]]; then
            snippy_ext_snps="pending"
        fi
    fi

    # ── Sanitize + écriture ligne TSV ─────────────────────────────────────
    _s() { local v="${1//[$'\n\r\t']/}"; printf '%s' "${v:-?}"; }

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$(_s "$BARCODE")"              "$(_s "$SAMPLE_NAME")"       "$(_s "$MODE")" \
        "$(_s "$raw_reads")"            "$(_s "$filtered_reads")"    "$(_s "$coverage_est")" \
        "$(_s "$ont_raw_mean_len")"     "$(_s "$illumina_reads")"    "$(_s "$illumina_mean_len")" \
        "$(_s "$n_contigs")"            "$(_s "$total_bp")"          "$(_s "$n50")"          "$(_s "$circular")" \
        "$(_s "$uni_long_bridges")"     "$(_s "$uni_short_bridges")" \
        "$(_s "$polypolish_changes")" \
        "$(_s "$pypolca_sub")"          "$(_s "$pypolca_indel")"     "$(_s "$pypolca_qv_before")" "$(_s "$pypolca_qv_after")" \
        "$(_s "$bakta_cds")"            "$(_s "$bakta_hypo")" \
        "$(_s "$COMPARE_TO")"           "$(_s "$snp_count")"         "$(_s "$snp_near_edge")" \
        "$(_s "$snippy_mean_cov")"      "$(_s "$snippy_pct_cov")" \
        "$(_s "$breseq_int_mutations")" "$(_s "$breseq_ext_mutations")" \
        "$(_s "$snippy_ext_snps")"      "$(_s "$snippy_ext_near_edge")" \
        "$(_s "$snippy_ext_mean_cov")"  "$(_s "$snippy_ext_pct_cov")" \
        >> "$SUMMARY_TSV"
done

# ---------------------------------------------------------------------------
# Helper: annotate snps.tab against a GFF3 (inchangé)
# ---------------------------------------------------------------------------
_annotate_snps() {
    local snptab="$1" ref_gff="$2" ref_fa="$3" ref_fna="$4"
    local b2p=""
    if [[ -f "$ref_fa" && -f "$ref_fna" ]]; then
        b2p=$(paste \
            <(grep "^>" "$ref_fna" | sed 's/^>//;s/ .*//') \
            <(grep "^>" "$ref_fa"  | sed 's/^>//;s/ .*//'))
    fi
    awk -F'\t' -v b2p_str="$b2p" '
    BEGIN {
        n = split(b2p_str, lines, "\n")
        for (i=1; i<=n; i++) {
            split(lines[i], p, "\t")
            if (p[1]!="" && p[2]!="") b2p[p[1]] = p[2]
        }
    }
    FNR==NR && $3=="CDS" {
        ctg = ($1 in b2p) ? b2p[$1] : $1
        prod=$9; gsub(/.*product=/,"",prod); gsub(/;.*/,"",prod)
        gsub(/%2C/,",",prod); gsub(/%25/,"%",prod)
        if (prod=="") prod="hypothetical protein"
        nc++
        cctg[nc]=ctg; cs[nc]=$4+0; ce[nc]=$5+0; cprod[nc]=prod
        next
    }
    FNR==NR { next }
    NR>1 {
        chr=$1; pos=$2+0
        for (i=1; i<=nc; i++) {
            if (cctg[i]==chr && pos>=cs[i] && pos<=ce[i]) {
                print $0"\t"cprod[i]; break
            }
        }
    }' "$ref_gff" "$snptab"
}

# ---------------------------------------------------------------------------
# Rapport texte
# ---------------------------------------------------------------------------
{
    echo "================================================================"
    echo "  Nanopore Bacterial Pipeline — Run Summary"
    echo "  Generated: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "================================================================"
    echo ""
    echo "Results directory: ${RESULTS_DIR}"
    echo ""

    # ── Stats reads ───────────────────────────────────────────────────────
    echo "── Read statistics ──────────────────────────────────────────────"
    echo ""
    printf "%-12s %-12s %-7s %-10s %-10s %-10s %-14s %-14s\n" \
        "barcode" "sample" "mode" "ont_reads" "mean_len" "cov_est" "illumina_reads" "illum_len"
    printf "%-12s %-12s %-7s %-10s %-10s %-10s %-14s %-14s\n" \
        "-------" "------" "----" "---------" "--------" "-------" "--------------" "---------"
    while IFS=$'\t' read -r barcode sample mode raw filt cov ont_mean ill_reads ill_mean \
        nctg tbp n50 circ ulb usb ppc psub pind pqvb pqva cds hypo \
        cmpTo snps edge snp_cov snp_pct b_int b_ext \
        sext_snps sext_edge sext_cov sext_pct; do
        [[ "$barcode" == "barcode" ]] && continue
        printf "%-12s %-12s %-7s %-10s %-10s %-10s %-14s %-14s\n" \
            "$barcode" "$sample" "$mode" "$raw" "${ont_mean}bp" "${cov}x" "$ill_reads" "${ill_mean}bp"
    done < "$SUMMARY_TSV"
    echo ""

    # ── Assembly overview ─────────────────────────────────────────────────
    echo "── Per-sample assembly overview ─────────────────────────────────"
    echo ""
    printf "%-12s %-12s %-7s %-10s %-8s %-10s %-10s %-10s %-8s %-10s\n" \
        "barcode" "sample" "mode" "cov_est" "contigs" "total_bp" "N50" "circular" "SNPs" "edge_SNPs"
    printf "%-12s %-12s %-7s %-10s %-8s %-10s %-10s %-10s %-8s %-10s\n" \
        "-------" "------" "----" "-------" "-------" "--------" "---" "--------" "----" "---------"
    while IFS=$'\t' read -r barcode sample mode raw filt cov ont_mean ill_reads ill_mean \
        nctg tbp n50 circ ulb usb ppc psub pind pqvb pqva cds hypo \
        cmpTo snps edge snp_cov snp_pct b_int b_ext \
        sext_snps sext_edge sext_cov sext_pct; do
        [[ "$barcode" == "barcode" ]] && continue
        flag=""
        [[ "$nctg" =~ ^[0-9]+$ ]] && [[ "$nctg" -gt 50 ]]    && flag="${flag}[!contigs>50]"
        [[ "$n50"  =~ ^[0-9]+$ ]] && [[ "$n50" -lt 100000 ]]  && flag="${flag}[!N50<100k]"
        [[ "$pqva" =~ ^[0-9] ]]  && awk "BEGIN{exit !($pqva < 40)}" && flag="${flag}[!QV<40]"
        printf "%-12s %-12s %-7s %-10s %-8s %-10s %-10s %-10s %-8s %-10s %s\n" \
            "$barcode" "$sample" "$mode" "${cov}x" "$nctg" "$tbp" "$n50" "$circ" "$snps" "$edge" "$flag"
    done < "$SUMMARY_TSV"
    echo ""

    # ── Hybrid contribution ───────────────────────────────────────────────
    echo "── Hybrid assembly & Illumina polishing contribution ────────────"
    echo ""
    printf "%-12s %-12s %-10s %-10s %-12s %-12s %-10s %-10s\n" \
        "barcode" "sample" "ONT_bridg" "Ill_bridg" "polyp_chg" "pypolca_s+i" "qv_before" "qv_after"
    printf "%-12s %-12s %-10s %-10s %-12s %-12s %-10s %-10s\n" \
        "-------" "------" "---------" "---------" "---------" "-----------" "---------" "--------"
    while IFS=$'\t' read -r barcode sample mode raw filt cov ont_mean ill_reads ill_mean \
        nctg tbp n50 circ ulb usb ppc psub pind pqvb pqva cds hypo \
        cmpTo snps edge snp_cov snp_pct b_int b_ext \
        sext_snps sext_edge sext_cov sext_pct; do
        [[ "$barcode" == "barcode" ]] && continue
        [[ "$mode" != "hybrid" ]] && continue
        pypolca_display="?"
        if [[ "$psub" =~ ^[0-9]+$ && "$pind" =~ ^[0-9]+$ ]]; then
            pypolca_total=$(( psub + pind ))
            [[ "$pypolca_total" -gt 0 ]] && pypolca_display="${psub}+${pind}" || pypolca_display="0"
        fi
        printf "%-12s %-12s %-10s %-10s %-12s %-12s %-10s %-10s\n" \
            "$barcode" "$sample" "$ulb" "$usb" "$ppc" "$pypolca_display" "$pqvb" "$pqva"
    done < "$SUMMARY_TSV"
    echo ""

    # ── Couverture vs référence ───────────────────────────────────────────
    echo "── Couverture vs référence (Snippy) ─────────────────────────────"
    echo ""
    printf "%-12s %-12s %-12s %-14s %-12s\n" \
        "barcode" "sample" "compare_to" "mean_depth" "pct_covered"
    printf "%-12s %-12s %-12s %-14s %-12s\n" \
        "-------" "------" "----------" "----------" "-----------"
    while IFS=$'\t' read -r barcode sample mode raw filt cov ont_mean ill_reads ill_mean \
        nctg tbp n50 circ ulb usb ppc psub pind pqvb pqva cds hypo \
        cmpTo snps edge snp_cov snp_pct b_int b_ext \
        sext_snps sext_edge sext_cov sext_pct; do
        [[ "$barcode" == "barcode" ]] && continue
        [[ "$snp_cov" == "N/A" ]]     && continue
        printf "%-12s %-12s %-12s %-14s %-12s\n" \
            "$barcode" "$sample" "$cmpTo" "${snp_cov}x" "${snp_pct}%"
    done < "$SUMMARY_TSV"
    echo ""

    # ── Annotation ────────────────────────────────────────────────────────
    echo "── Annotation (Bakta) ───────────────────────────────────────────"
    echo ""
    printf "%-12s %-12s %-8s %-12s\n" "barcode" "sample" "CDS" "hypothetical"
    printf "%-12s %-12s %-8s %-12s\n" "-------" "------" "---" "------------"
    while IFS=$'\t' read -r barcode sample mode raw filt cov ont_mean ill_reads ill_mean \
        nctg tbp n50 circ ulb usb ppc psub pind pqvb pqva cds hypo \
        cmpTo snps edge snp_cov snp_pct b_int b_ext \
        sext_snps sext_edge sext_cov sext_pct; do
        [[ "$barcode" == "barcode" ]] && continue
        printf "%-12s %-12s %-8s %-12s\n" "$barcode" "$sample" "$cds" "$hypo"
    done < "$SUMMARY_TSV"
    echo ""

    # ── Variants dans des gènes annotés ──────────────────────────────────
    echo "── Variants dans des gènes annotés (Snippy) ────────────────────────"
    echo ""
    SNIPPY_ANY=false
    for BARCODE in "${BARCODES[@]}"; do
        sn=$(get_sample_field "$BARCODE" "sample_name" 2>/dev/null || echo "?")
        cmpTo=$(get_sample_field "$BARCODE" "compare_to" 2>/dev/null || echo "-")
        [[ -z "$cmpTo" || "$cmpTo" == "-" ]] && continue
        snptab="${RESULTS_DIR}/05_variants/${BARCODE}_${sn}/snps.tab"
        [[ ! -f "$snptab" ]] && continue
        ref_sn=$(get_sample_field "$cmpTo" "sample_name" 2>/dev/null || echo "")
        ref_gff="${RESULTS_DIR}/04_annotated/${cmpTo}_${ref_sn}/${ref_sn}.gff3"
        ref_fa="${RESULTS_DIR}/03_polished/${cmpTo}_${ref_sn}/polished.fasta"
        ref_fna="${RESULTS_DIR}/04_annotated/${cmpTo}_${ref_sn}/${ref_sn}.fna"
        [[ ! -f "$ref_gff" ]] && { echo "  ${BARCODE}/${sn}: [WARN] GFF3 introuvable"; continue; }
        hits=$(_annotate_snps "$snptab" "$ref_gff" "$ref_fa" "$ref_fna")
        [[ -z "$hits" ]] && continue
        SNIPPY_ANY=true
        echo "  ${BARCODE} / ${sn} (vs ${cmpTo}):"
        echo "$hits" | awk -F'\t' '{
            printf "    %-10s %-8s %-8s -> %-8s  %s\n", $2, $3, $4, $5, $NF
        }'
        echo ""
    done
    [[ "$SNIPPY_ANY" == "false" ]] && echo "  (aucun variant annoté, ou step 05 non terminé)"
    echo ""

    # ── Breseq comparaison ────────────────────────────────────────────────
    echo "── Variant Calls : Snippy / Breseq interne / Breseq externe ─────"
    echo ""
    printf "%-12s %-12s %-12s %-12s %-12s %-12s\n" \
        "barcode" "sample" "compare_to" "snippy_SNPs" "breseq_int" "breseq_ext"
    printf "%-12s %-12s %-12s %-12s %-12s %-12s\n" \
        "-------" "------" "----------" "-----------" "----------" "----------"
    while IFS=$'\t' read -r barcode sample mode raw filt cov ont_mean ill_reads ill_mean \
        nctg tbp n50 circ ulb usb ppc psub pind pqvb pqva cds hypo \
        cmpTo snps edge snp_cov snp_pct b_int b_ext \
        sext_snps sext_edge sext_cov sext_pct; do
        [[ "$barcode" == "barcode" ]] && continue
        [[ "$snps" == "N/A" && "$b_int" == "N/A" && "$b_ext" == "N/A" ]] && continue
        printf "%-12s %-12s %-12s %-12s %-12s %-12s\n" \
            "$barcode" "$sample" "$cmpTo" "$snps" "$b_int" "$b_ext"
    done < "$SUMMARY_TSV"
    echo ""
    echo "  N/A = non lancé  |  error = VCF absent (voir log)"
    echo "  Rapports breseq interne  : ${RESULTS_DIR}/05_variants/*/breseq_internal/output/index.html"
    echo "  Rapports breseq externe  : ${RESULTS_DIR}/05_variants/*/breseq_external/output/index.html"
    echo ""

    # ── Méthylation ───────────────────────────────────────────────────────
    echo "── Modified base methylation (step 07) ─────────────────────────"
    echo ""
    MODBASE_HEADER_PRINTED=false
    for BARCODE in "${BARCODES[@]}"; do
        sn=$(get_sample_field "$BARCODE" "sample_name" 2>/dev/null || echo "?")
        bed="${RESULTS_DIR}/07_modbase/${BARCODE}_${sn}/methylation.bed"
        if [[ -f "$bed" ]]; then
            if [[ "$MODBASE_HEADER_PRINTED" == "false" ]]; then
                printf "%-12s %-12s %-12s %-12s %-12s\n" \
                    "barcode" "sample" "6mA_hi-conf" "5mC_hi-conf" "mapped_reads"
                printf "%-12s %-12s %-12s %-12s %-12s\n" \
                    "-------" "------" "-----------" "-----------" "------------"
                MODBASE_HEADER_PRINTED=true
            fi
            n6mA=$(awk '$4=="a"&&$10>=8&&$11>=50{n++}END{print n+0}' "$bed")
            n5mC=$(awk '$4=="m"&&$10>=8&&$11>=50{n++}END{print n+0}' "$bed")
            bam="${RESULTS_DIR}/07_modbase/${BARCODE}_${sn}/aligned_mods.bam"
            mapped="?"
            [[ -f "$bam" ]] && mapped=$(samtools view -c -F 4 "$bam" 2>/dev/null || echo "?")
            printf "%-12s %-12s %-12s %-12s %-12s\n" \
                "$BARCODE" "$sn" "$n6mA" "$n5mC" "$mapped"
        fi
    done
    [[ "$MODBASE_HEADER_PRINTED" == "false" ]] && echo "  (step 07 non lancé)"
    echo ""

    echo "── Key output locations ─────────────────────────────────────────"
    echo "  Polished assemblies  : ${RESULTS_DIR}/03_polished/"
    echo "  Bakta annotations    : ${RESULTS_DIR}/04_annotated/"
    echo "  Snippy variants      : ${RESULTS_DIR}/05_variants/*/snps.vcf"
    echo "  Core SNP alignment   : ${RESULTS_DIR}/05_variants/core_vs_*/core.tab"
    echo "  Breseq interne       : ${RESULTS_DIR}/05_variants/*/breseq_internal/output/index.html"
    echo "  Breseq externe       : ${RESULTS_DIR}/05_variants/*/breseq_external/output/index.html"
    echo "  Modified bases       : ${RESULTS_DIR}/07_modbase/"
    echo ""
    echo "── Full TSV table ───────────────────────────────────────────────"
    column -t -s $'\t' "$SUMMARY_TSV"
} > "$SUMMARY_TXT"

# ---------------------------------------------------------------------------
# HTML summary report
# ---------------------------------------------------------------------------
HTML_REPORT="${OUT}/run_summary.html"

{
printf '<!DOCTYPE html>\n<html lang="en">\n<head>\n'
printf '<meta charset="UTF-8">\n'
printf '<title>Nanopore Pipeline Summary — %s</title>\n' "${RUN_NAME}"
cat << 'HTMLCSS'
<style>
  body  { font-family: Arial, sans-serif; margin: 2em; background: #f8f8f8; }
  h1    { color: #2c3e50; }
  h2    { color: #34495e; border-bottom: 2px solid #bdc3c7; padding-bottom: 4px; margin-top: 2em; }
  h3    { color: #555; font-size: 0.95em; margin: 1em 0 0.3em; }
  table { border-collapse: collapse; width: 100%; margin-bottom: 1.5em; background: white; }
  th    { background: #2c3e50; color: white; padding: 8px 12px; text-align: left;
          cursor: pointer; white-space: nowrap; user-select: none; }
  th:hover       { background: #3a5f8a; }
  th::after      { content: " ⇅"; font-size: 10px; opacity: .5; }
  td    { padding: 6px 12px; border-bottom: 1px solid #ecf0f1; white-space: nowrap; }
  tr:nth-child(even) { background: #f2f3f4; }
  tr:hover td   { background: #e8f4fd; }
  .flag { color: #e74c3c; font-weight: bold; }
  .ok   { color: #27ae60; }
  .na   { color: #aaa; }
  .warn { color: #e67e22; font-weight: bold; }
  a     { color: #2980b9; }
  .note { font-size: 0.85em; color: #666; margin-top: -0.8em; margin-bottom: 1em; }
  .filter { padding: 4px 8px; border: 1px solid #ccc; border-radius: 4px;
            font-size: 12px; width: 240px; margin-bottom: 8px; }
  .cards { display: flex; flex-wrap: wrap; gap: 12px; margin-bottom: 1.5em; }
  .card  { background: white; border: 1px solid #dde; border-radius: 6px;
           padding: 12px 20px; min-width: 130px; text-align: center; }
  .card .v { font-size: 24px; font-weight: 700; color: #2c7bb6; }
  .card .l { font-size: 11px; color: #666; margin-top: 2px; }
</style>
<script>
function sortTbl(th) {
  var tbl = th.closest('table'), col = th.cellIndex;
  var asc = (tbl.dataset.sc == col) ? (tbl.dataset.sa != '1') : true;
  var rows = Array.from(tbl.tBodies[0].rows);
  rows.sort(function(a,b){
    var va = a.cells[col] ? a.cells[col].innerText.trim() : '';
    var vb = b.cells[col] ? b.cells[col].innerText.trim() : '';
    var na = parseFloat(va), nb = parseFloat(vb);
    if (!isNaN(na) && !isNaN(nb)) return asc ? na-nb : nb-na;
    return asc ? va.localeCompare(vb) : vb.localeCompare(va);
  });
  rows.forEach(function(r){ tbl.tBodies[0].appendChild(r); });
  tbl.dataset.sc = col; tbl.dataset.sa = asc ? '1' : '0';
}
function filterTbl(inp, tblId) {
  var q = inp.value.toLowerCase();
  document.getElementById(tblId).tBodies[0].querySelectorAll('tr').forEach(function(r){
    r.style.display = r.innerText.toLowerCase().includes(q) ? '' : 'none';
  });
}
</script>
</head>
<body>
HTMLCSS

printf '<h1>Nanopore Pipeline Summary</h1>\n'
printf '<p><b>Run:</b> %s &nbsp;|&nbsp; <b>Generated:</b> %s &nbsp;|&nbsp; <b>Results dir:</b> %s</p>\n' \
    "${RUN_NAME}" "$(date '+%Y-%m-%d %H:%M:%S')" "${RESULTS_DIR}"

# ── Cartes résumé ────────────────────────────────────────────────────────────
N_SAMPLES=$(list_barcodes | wc -l | tr -d ' ')
N_HYBRID=$(grep -v '^#' "${RESULTS_DIR}/sample_sheet.tsv" | tail -n +2 | awk -F'\t' '$3=="hybrid"' | wc -l || echo 0)
N_ONT=$(grep -v '^#' "${RESULTS_DIR}/sample_sheet.tsv"    | tail -n +2 | awk -F'\t' '$3=="ont"'    | wc -l || echo 0)
N_COMPLETE=$(find "${RESULTS_DIR}/03_polished" -name "polished.fasta" 2>/dev/null | wc -l | tr -d ' ' || echo 0)

printf '<div class="cards">\n'
printf '  <div class="card"><div class="v">%s</div><div class="l">Samples</div></div>\n' "$N_SAMPLES"
printf '  <div class="card"><div class="v">%s</div><div class="l">Hybrid mode</div></div>\n'  "$N_HYBRID"
printf '  <div class="card"><div class="v">%s</div><div class="l">ONT mode</div></div>\n'     "$N_ONT"
printf '  <div class="card"><div class="v" style="color:#27ae60">%s</div><div class="l">Complete assemblies</div></div>\n' "$N_COMPLETE"
printf '</div>\n'

# ── Section 0 : Methods ───────────────────────────────────────────────────────
printf '<h2>0. Methods</h2>\n'
printf '<div style="background:white;border:1px solid #dde;border-radius:6px;padding:14px 18px;font-size:12px;line-height:1.7;margin-bottom:1.5em;">\n'

printf '<h3 style="color:#2c3e50;margin-bottom:6px">Read quality control</h3>\n'
printf '<p><strong>Tool:</strong> <code>chopper</code> &nbsp;|&nbsp; '
printf '<strong>Min quality:</strong> Q%s &nbsp;|&nbsp; ' "${MIN_READ_QUALITY}"
printf '<strong>Min length:</strong> %s bp</p>\n' "${MIN_READ_LENGTH}"

printf '<h3 style="color:#2c3e50;margin-top:10px;margin-bottom:6px">Assembly</h3>\n'
printf '<p>'
printf '<strong>ONT-only:</strong> <code>Flye --nano-hq</code>, genome size estimate: %s, followed by polishing with <code>Medaka</code> (model: <code>%s</code>).<br>' \
    "${GENOME_SIZE_ESTIMATE}" "${MEDAKA_MODEL}"
printf '<strong>Hybrid:</strong> <code>Unicycler</code> (short + long reads), followed by <code>Polypolish</code> (Illumina read polishing) and <code>Pypolca</code> (Illumina consensus correction).<br>'
printf 'All assemblies reoriented with <code>Dnaapler</code> (chromosome origin: <em>dnaA</em>; plasmid origin: <em>repA</em>).'
printf '</p>\n'
printf '<p><strong>Assembly quality flags in report:</strong> [!contigs&gt;50] = fragmented &nbsp;|&nbsp; [!N50&lt;100k] = low contiguity &nbsp;|&nbsp; [!QV&lt;40] = poor post-polishing quality.</p>\n'

printf '<h3 style="color:#2c3e50;margin-top:10px;margin-bottom:6px">Annotation</h3>\n'
printf '<p><strong>Tool:</strong> <code>Bakta</code> &nbsp;|&nbsp; <strong>Database:</strong> <code>%s</code></p>\n' "${BAKTA_DB}"

printf '<h3 style="color:#2c3e50;margin-top:10px;margin-bottom:6px">Variant calling — Snippy (internal reference)</h3>\n'
printf '<p><strong>Tool:</strong> <code>Snippy</code> &nbsp;|&nbsp; '
printf 'Reference: polished assembly of <code>compare_to</code> barcode.<br>'
printf '<strong>Mode hybrid:</strong> <code>snippy --R1/--R2</code> (Illumina reads) &nbsp;|&nbsp; '
printf '<strong>Mode ONT:</strong> <code>snippy --se --nanopore</code><br>'
printf '<strong>Parameters:</strong> <code>--minfrac %s</code> &nbsp;|&nbsp; <code>--mincov %s</code> &nbsp;|&nbsp; <code>--mapqual %s</code><br>' \
    "${SNIPPY_MIN_FRAC}" "${SNIPPY_MIN_COVERAGE}" "${SNIPPY_MAPQ}"
printf '<strong>Edge filter:</strong> variants within %s bp of a contig end are flagged as edge SNPs.</p>\n' "${SNIPPY_CONTIG_EDGE_FILTER:-100}"

printf '<h3 style="color:#2c3e50;margin-top:10px;margin-bottom:6px">Variant calling — Snippy (external reference)</h3>\n'
printf '<p><strong>Tool:</strong> <code>Snippy</code> &nbsp;|&nbsp; '
printf 'Reference: <code>EXTERNAL_REF</code> from config (<code>%s</code>), or absolute path in <code>snippy_ext_ref</code> column.<br>' \
    "$(basename "${EXTERNAL_REF:-not set}")"
printf 'Same parameters as internal Snippy. Output: <code>snippy_external/</code>.</p>\n'

printf '<h3 style="color:#2c3e50;margin-top:10px;margin-bottom:6px">Variant calling — Breseq</h3>\n'
printf '<p><strong>Tool:</strong> <code>Breseq</code> (hybrid mode only; requires Illumina reads from the sample).<br>'
printf '<strong>Internal reference:</strong> Bakta-annotated <code>.gbff</code> of the <code>breseq_ref</code> barcode (step 04 output).<br>'
printf '<strong>External reference:</strong> <code>EXTERNAL_REF</code> (<code>%s</code>), or absolute path in <code>breseq_ext_ref</code>.<br>' \
    "$(basename "${EXTERNAL_REF:-not set}")"
printf '<strong>ONT read handling:</strong> ONT filtered reads are fragmented into non-overlapping %s bp chunks (<code>seqkit sliding</code>, step %s bp) '  \
    "${ONT_BRESEQ_FRAGMENT_SIZE}" "${ONT_BRESEQ_FRAGMENT_STEP}"
printf 'and passed alongside Illumina R1+R2 reads. Fragments are deleted after use.</p>\n'

printf '<h3 style="color:#2c3e50;margin-top:10px;margin-bottom:6px">Coverage statistics</h3>\n'
printf '<p>Computed with <code>samtools coverage</code> on each alignment BAM, weighted across all reference contigs.<br>'
printf '<strong>Warning threshold:</strong> mean depth &lt; %s× &nbsp;|&nbsp; '  "${MIN_COVERAGE_WARN}"
printf '<strong>Failure threshold:</strong> mean depth &lt; %s×.</p>\n' "${MIN_COVERAGE_FAIL}"

printf '<h3 style="color:#2c3e50;margin-top:10px;margin-bottom:6px">Modified base methylation (step 07)</h3>\n'
printf '<p><strong>Basecalling:</strong> Dorado with <code>--modified-bases 6mA 5mCG_5hmCG</code> (mod-tagged BAMs).<br>'
printf '<strong>Alignment:</strong> <code>minimap2 -ax map-ont -y</code> vs polished assembly.<br>'
printf '<strong>Pileup:</strong> <code>modkit pileup</code> — reports per-position methylation frequency (0–100%%).<br>'
printf '<strong>Thresholds:</strong> high-confidence sites: coverage &ge;8× AND modification frequency &ge;50%%.<br>'
printf '<strong>Note:</strong> methylation analysis is based exclusively on Nanopore reads (no Illumina signal).</p>\n'

printf '</div>\n'

# ── Section 1 : Stats reads ──────────────────────────────────────────────────
printf '<h2>1. Read Statistics</h2>\n'

printf '<p class="note">Reads as they enter the pipeline (before ONT filtering). '
printf 'Illumina reads: exact count (R1+R2) — all reads are used in full by assembly (Unicycler) and polishing (Polypolish/Pypolca).</p>\n'
printf '<input class="filter" placeholder="Filter…" oninput="filterTbl(this,&apos;tbl-reads&apos;)">\n'
printf '<table id="tbl-reads">\n'
printf '<tr><th onclick="sortTbl(this)">Barcode</th>'
printf '<th onclick="sortTbl(this)">Sample</th>'
printf '<th onclick="sortTbl(this)">Mode</th>'
printf '<th onclick="sortTbl(this)">Raw ONT reads</th>'
printf '<th onclick="sortTbl(this)">Mean length ONT (bp)</th>'
printf '<th onclick="sortTbl(this)">Estimated coverage</th>'
printf '<th onclick="sortTbl(this)">Illumina reads (R1+R2)</th>'
printf '<th onclick="sortTbl(this)">Mean length Illumina (bp)</th></tr>\n'

while IFS=$'\t' read -r barcode sample mode raw filt cov ont_mean ill_reads ill_mean \
    nctg tbp n50 circ ulb usb ppc psub pind pqvb pqva cds hypo \
    cmpTo snps edge snp_cov snp_pct b_int b_ext; do
    [[ "$barcode" == "barcode" ]] && continue
    printf '<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%sx</td><td>%s</td><td>%s</td></tr>\n' \
        "$barcode" "$sample" "$mode" "$raw" "$ont_mean" "$cov" "$ill_reads" "$ill_mean"
done < "$SUMMARY_TSV"
printf '</table>\n'

# ── Section 2 : Assemblage & Polissage ───────────────────────────────────────
printf '<h2>2. Assembly &amp; Polishing</h2>\n'
printf '<p class="note">Flags: [!contigs&gt;50] fragmented assembly &nbsp;|&nbsp; [!N50&lt;100k] low contiguity &nbsp;|&nbsp; [!QV&lt;40] poor post-polishing quality</p>\n'
printf '<input class="filter" placeholder="Filter…" oninput="filterTbl(this,&apos;tbl-asm&apos;)">\n'
printf '<table id="tbl-asm">\n'
printf '<tr><th onclick="sortTbl(this)">Barcode</th>'
printf '<th onclick="sortTbl(this)">Sample</th>'
printf '<th onclick="sortTbl(this)">Mode</th>'
printf '<th onclick="sortTbl(this)">Contigs</th>'
printf '<th onclick="sortTbl(this)">Total bp</th>'
printf '<th onclick="sortTbl(this)">N50</th>'
printf '<th onclick="sortTbl(this)">Circular</th>'
printf '<th onclick="sortTbl(this)">Polypolish</th>'
printf '<th onclick="sortTbl(this)">Pypolca S+I</th>'
printf '<th onclick="sortTbl(this)">QV before</th>'
printf '<th onclick="sortTbl(this)">QV after</th>'
printf '<th>Flags</th></tr>\n'

while IFS=$'\t' read -r barcode sample mode raw filt cov ont_mean ill_reads ill_mean \
    nctg tbp n50 circ ulb usb ppc psub pind pqvb pqva cds hypo \
    cmpTo snps edge snp_cov snp_pct b_int b_ext; do
    [[ "$barcode" == "barcode" ]] && continue
    flags=""
    [[ "$nctg" =~ ^[0-9]+$ ]] && [[ "$nctg" -gt 50 ]] && flags="${flags}<span class='flag'>[!contigs&gt;50]</span> "
    [[ "$n50"  =~ ^[0-9]+$ ]] && [[ "$n50" -lt 100000 ]] && flags="${flags}<span class='flag'>[!N50&lt;100k]</span> "
    [[ "$pqva" =~ ^[0-9] ]] && awk "BEGIN{exit !($pqva < 40)}" && flags="${flags}<span class='flag'>[!QV&lt;40]</span> "
    [[ -z "$flags" ]] && flags="<span class='ok'>OK</span>"
    if [[ "${psub}" =~ ^[0-9]+$ && "${pind}" =~ ^[0-9]+$ ]]; then
        psum=$(( psub + pind ))
        [[ "$psum" -gt 0 ]] && py_cell="<span title='${psub} sub + ${pind} indels'>${psub}+${pind}</span>" || py_cell="0"
    else
        py_cell="${psub:-N/A}"
    fi
    printf '<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\n' \
        "$barcode" "$sample" "$mode" "$nctg" "$tbp" "$n50" "$circ" "$ppc" \
        "$py_cell" "$pqvb" "$pqva" "$flags"
done < "$SUMMARY_TSV"
printf '</table>\n'

# ── Section 3 : Couverture vs référence (Snippy & Breseq) ────────────────────
# Lue directement depuis les coverage_stats.txt de chaque outil (pas du TSV)
# pour pouvoir afficher Snippy + Breseq interne + Breseq externe sur le même tableau.
printf '<h2>3. Reference Coverage</h2>\n'
printf '<p class="note">'
printf 'Computed via <code>samtools coverage</code> on the alignment BAM, '
printf 'weighted across all contigs. '
printf '<span style="display:inline-block;background:#1565c0;color:#fff;padding:1px 7px;border-radius:10px;font-size:11px">Snippy</span> '
printf '= BAM <code>snps.bam</code> (internal reference, compare_to) &mdash; '
printf '<span style="display:inline-block;background:#6a1b9a;color:#fff;padding:1px 7px;border-radius:10px;font-size:11px">Snippy ext.</span> '
printf '= BAM <code>snippy_external/snps.bam</code> (external reference, snippy_ext_ref) &mdash; '
printf '<span style="display:inline-block;background:#2e7d32;color:#fff;padding:1px 7px;border-radius:10px;font-size:11px">Breseq int.</span> '
printf '= internal BAM <code>breseq_internal/data/reference.bam</code> &mdash; '
printf '<span style="display:inline-block;background:#e65100;color:#fff;padding:1px 7px;border-radius:10px;font-size:11px">Breseq ext.</span> '
printf '= external BAM <code>breseq_external/data/reference.bam</code>.</p>\n'
printf '<table id="tbl-cov">\n'
printf '<tr>'
printf '<th onclick="sortTbl(this)">Barcode</th>'
printf '<th onclick="sortTbl(this)">Sample</th>'
printf '<th onclick="sortTbl(this)">Tool</th>'
printf '<th onclick="sortTbl(this)">Reference used</th>'
printf '<th onclick="sortTbl(this)">Mean depth (&times;)</th>'
printf '<th onclick="sortTbl(this)">Genome covered (%%)</th>'
printf '</tr>\n'

COV_ANY=false
for _bc in "${BARCODES[@]}"; do
    _sn=$(get_sample_field "$_bc" "sample_name")
    _cmp=$(get_sample_field "$_bc" "compare_to")
    _bref=$(get_sample_field_opt "$_bc" "breseq_ref")
    _bext=$(get_sample_field_opt "$_bc" "breseq_ext_ref")
    _vdir="${RESULTS_DIR}/05_variants/${_bc}_${_sn}"

    # Fonction locale : émet une ligne si coverage_stats.txt existe et est complet
    # Args: badge_html ref_label cov_file
    _emit_cov_row() {
        local badge="$1" ref_lbl="$2" cov_f="$3"
        [[ ! -f "$cov_f" ]] && return 0
        local md pc
        md=$(grep 'mean_depth='  "$cov_f" | cut -d= -f2)
        pc=$(grep 'pct_covered=' "$cov_f" | cut -d= -f2)
        [[ -z "$md" || -z "$pc" ]] && return 0
        COV_ANY=true
        local cov_class="ok" pct_class="ok"
        if [[ "$md" =~ ^[0-9] ]]; then
            awk "BEGIN{exit !($md < ${MIN_COVERAGE_WARN})}" 2>/dev/null && cov_class="flag" || true
        fi
        if [[ "$pc" =~ ^[0-9] ]]; then
            if awk "BEGIN{exit !($pc < 50)}" 2>/dev/null; then
                pct_class="flag"
            elif awk "BEGIN{exit !($pc < 90)}" 2>/dev/null; then
                pct_class="warn"
            fi
        fi
        printf '<tr><td>%s</td><td>%s</td><td>%s</td><td style="font-size:11px">%s</td><td><span class="%s">%s&times;</span></td><td><span class="%s">%s%%</span></td></tr>\n' \
            "$_bc" "$_sn" "$badge" "$ref_lbl" "$cov_class" "$md" "$pct_class" "$pc"
    }

    # Snippy
    if [[ "${_cmp}" != "-" && -n "${_cmp}" ]]; then
        _cmp_sn=$(get_sample_field "$_cmp" "sample_name" 2>/dev/null || echo "$_cmp")
        _emit_cov_row \
            "<span style='background:#1565c0;color:#fff;padding:1px 7px;border-radius:10px;font-size:11px'>Snippy</span>" \
            "${_cmp_sn} (${_cmp})" \
            "${_vdir}/coverage_stats.txt"
    fi


    # Breseq interne
    if [[ "${_bref}" != "-" && -n "${_bref}" ]]; then
        _bref_sn=$(get_sample_field "$_bref" "sample_name" 2>/dev/null || echo "$_bref")
        _emit_cov_row \
            "<span style='background:#2e7d32;color:#fff;padding:1px 7px;border-radius:10px;font-size:11px'>Breseq int.</span>" \
            "${_bref_sn} (${_bref})" \
            "${_vdir}/breseq_internal/coverage_stats.txt"
    fi

    # Breseq externe
    if [[ "${_bext}" != "-" && -n "${_bext}" ]]; then
        if [[ "${_bext}" == "ext" ]]; then
            _ext_lbl="$(basename "${EXTERNAL_REF:-external reference}")"
        else
            _ext_lbl="$(basename "${_bext}")"
        fi
        _emit_cov_row \
            "<span style='background:#e65100;color:#fff;padding:1px 7px;border-radius:10px;font-size:11px'>Breseq ext.</span>" \
            "${_ext_lbl}" \
            "${_vdir}/breseq_external/coverage_stats.txt"
    fi

    # Snippy externe
    _sext=$(get_sample_field_opt "$_bc" "snippy_ext_ref" 2>/dev/null || true)
    if [[ "${_sext}" != "-" && -n "${_sext}" ]]; then
        if [[ "${_sext}" == "ext" ]]; then
            _sext_lbl="$(basename "${EXTERNAL_REF:-external reference}")"
        else
            _sext_lbl="$(basename "${_sext}")"
        fi
        _emit_cov_row \
            "<span style='background:#6a1b9a;color:#fff;padding:1px 7px;border-radius:10px;font-size:11px'>Snippy ext.</span>" \
            "${_sext_lbl}" \
            "${_vdir}/snippy_external/coverage_stats.txt"
    fi
done

if [[ "$COV_ANY" == "false" ]]; then
    printf '<tr><td colspan="6"><span class="na">(step 05 not complete or no reference defined)</span></td></tr>\n'
fi
printf '</table>\n'

# ── Section 4 : Variants Snippy & Breseq ─────────────────────────────────────
printf '<h2>4. Variant Calls — Snippy, Breseq internal &amp; Breseq external</h2>\n'
printf '<p class="note">'
printf 'edge_SNPs = number of mutations in the last %s bp of a contig end. If you find mutations here, check manually their distribution — this could be an assembly artefact. ' "${SNIPPY_CONTIG_EDGE_FILTER:-100}"
printf 'breseq_int = mutations vs Bakta gbff of the internal reference sample. '
printf 'breseq_ext = mutations vs external reference (.gbk/.gbff). '
printf 'snippy_ext = Snippy vs external reference (snippy_ext_ref in sample sheet; same reference file as breseq_ext by default). '
printf 'N/A = not run &nbsp;|&nbsp; error = VCF absent (check log).</p>\n'
printf '<input class="filter" placeholder="Filter…" oninput="filterTbl(this,&apos;tbl-var&apos;)">\n'
printf '<table id="tbl-var">\n'
printf '<tr><th onclick="sortTbl(this)">Barcode</th>'
printf '<th onclick="sortTbl(this)">Sample</th>'
printf '<th onclick="sortTbl(this)">compare_to</th>'
printf '<th onclick="sortTbl(this)">Snippy SNPs</th>'
printf '<th onclick="sortTbl(this)">Edge SNPs</th>'
printf '<th>Snippy report</th>'
printf '<th onclick="sortTbl(this)">Breseq int.</th>'
printf '<th>Breseq int. report</th>'
printf '<th onclick="sortTbl(this)">Breseq ext.</th>'
printf '<th>Breseq ext. report</th>'
printf '<th onclick="sortTbl(this)">Snippy ext. SNPs</th>'
printf '<th onclick="sortTbl(this)">Edge SNPs (ext.)</th>'
printf '<th>Snippy ext. report</th></tr>\n'

while IFS=$'\t' read -r barcode sample mode raw filt cov ont_mean ill_reads ill_mean \
    nctg tbp n50 circ ulb usb ppc psub pind pqvb pqva cds hypo \
    cmpTo snps edge snp_cov snp_pct b_int b_ext \
    sext_snps sext_edge sext_cov sext_pct; do
    [[ "$barcode" == "barcode" ]] && continue
    [[ "$snps" == "N/A" && "$b_int" == "N/A" && "$b_ext" == "N/A" && "$sext_snps" == "N/A" ]] && continue

    snippy_html_abs="${RESULTS_DIR}/05_variants/${barcode}_${sample}/snps.html"
    snippy_html_rel="../05_variants/${barcode}_${sample}/snps.html"
    snippy_link="<span class='na'>n/a</span>"
    [[ -f "$snippy_html_abs" ]] && snippy_link="<a href='${snippy_html_rel}'>snps.html</a>"

    # Lien breseq interne : cherche d'abord breseq_internal/ puis breseq/ (ancien nom)
    bint_html_abs="${RESULTS_DIR}/05_variants/${barcode}_${sample}/breseq_internal/output/index.html"
    bint_html_rel="../05_variants/${barcode}_${sample}/breseq_internal/output/index.html"
    if [[ ! -f "$bint_html_abs" ]]; then
        bint_html_abs="${RESULTS_DIR}/05_variants/${barcode}_${sample}/breseq/output/index.html"
        bint_html_rel="../05_variants/${barcode}_${sample}/breseq/output/index.html"
    fi
    bint_link="<span class='na'>n/a</span>"
    [[ -f "$bint_html_abs" ]] && bint_link="<a href='${bint_html_rel}'>index.html</a>"
    [[ "$b_int" == "error" ]] && bint_link="<span class='flag'>error</span>"

    bext_html_abs="${RESULTS_DIR}/05_variants/${barcode}_${sample}/breseq_external/output/index.html"
    bext_html_rel="../05_variants/${barcode}_${sample}/breseq_external/output/index.html"
    bext_link="<span class='na'>n/a</span>"
    [[ -f "$bext_html_abs" ]] && bext_link="<a href='${bext_html_rel}'>index.html</a>"
    [[ "$b_ext" == "error" ]] && bext_link="<span class='flag'>error</span>"

    sext_html_abs="${RESULTS_DIR}/05_variants/${barcode}_${sample}/snippy_external/snps.html"
    sext_html_rel="../05_variants/${barcode}_${sample}/snippy_external/snps.html"
    sext_link="<span class='na'>n/a</span>"
    [[ -f "$sext_html_abs" ]] && sext_link="<a href='${sext_html_rel}'>snps.html</a>"

    printf '<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\n' \
        "$barcode" "$sample" "$cmpTo" "$snps" "$edge" "$snippy_link" \
        "$b_int" "$bint_link" "$b_ext" "$bext_link" \
        "$sext_snps" "$sext_edge" "$sext_link"
done < "$SUMMARY_TSV"
printf '</table>\n'


# ── Section 5 : Annotation Bakta ─────────────────────────────────────────────
printf '<h2>5. Annotation (Bakta)</h2>\n'
printf '<table>\n'
printf '<tr><th onclick="sortTbl(this)">Barcode</th>'
printf '<th onclick="sortTbl(this)">Sample</th>'
printf '<th onclick="sortTbl(this)">CDS</th>'
printf '<th onclick="sortTbl(this)">Hypothetical</th>'
printf '<th>GenBank (.gbff)</th>'
printf '<th>Genome FASTA (.fna)</th>'
printf '<th>Proteins FASTA (.faa)</th></tr>\n'

while IFS=$'\t' read -r barcode sample mode raw filt cov ont_mean ill_reads ill_mean \
    nctg tbp n50 circ ulb usb ppc psub pind pqvb pqva cds hypo \
    cmpTo snps edge snp_cov snp_pct b_int b_ext; do
    [[ "$barcode" == "barcode" ]] && continue
    _ann_dir="${RESULTS_DIR}/04_annotated/${barcode}_${sample}"
    _ann_rel="../04_annotated/${barcode}_${sample}"
    gbk_link="<span class='na'>n/a</span>"
    fna_link="<span class='na'>n/a</span>"
    faa_link="<span class='na'>n/a</span>"
    [[ -f "${_ann_dir}/${sample}.gbff" ]] && gbk_link="<a href='${_ann_rel}/${sample}.gbff'>${sample}.gbff</a>"
    [[ -f "${_ann_dir}/${sample}.fna"  ]] && fna_link="<a href='${_ann_rel}/${sample}.fna'>${sample}.fna</a>"
    [[ -f "${_ann_dir}/${sample}.faa"  ]] && faa_link="<a href='${_ann_rel}/${sample}.faa'>${sample}.faa</a>"
    printf '<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\n' \
        "$barcode" "$sample" "$cds" "$hypo" "$gbk_link" "$fna_link" "$faa_link"
done < "$SUMMARY_TSV"
printf '</table>\n'

# ── Section 6 : Méthylation ───────────────────────────────────────────────────
printf '<h2>6. Modified Base Methylation (step 07)</h2>\n'
printf '<p class="note">Hi-conf = &ge;8&times; coverage AND &ge;50%% modification frequency.<br>'
printf '<strong>Note:</strong> This analysis is based exclusively on Nanopore reads (Dorado mod-tagged BAMs). No Illumina signal is used for methylation detection.</p>\n'
printf '<table>\n'
printf '<tr><th onclick="sortTbl(this)">Barcode</th>'
printf '<th onclick="sortTbl(this)">Sample</th>'
printf '<th onclick="sortTbl(this)">6mA hi-conf</th>'
printf '<th onclick="sortTbl(this)">5mC hi-conf</th>'
printf '<th onclick="sortTbl(this)">Mapped reads</th></tr>\n'

MODBASE_HTML_ANY=false
for BARCODE in "${BARCODES[@]}"; do
    sn=$(get_sample_field "$BARCODE" "sample_name" 2>/dev/null || echo "?")
    bed="${RESULTS_DIR}/07_modbase/${BARCODE}_${sn}/methylation.bed"
    if [[ -f "$bed" ]]; then
        MODBASE_HTML_ANY=true
        n6mA=$(awk '$4=="a"&&$10>=8&&$11>=50{n++}END{print n+0}' "$bed")
        n5mC=$(awk '$4=="m"&&$10>=8&&$11>=50{n++}END{print n+0}' "$bed")
        bam="${RESULTS_DIR}/07_modbase/${BARCODE}_${sn}/aligned_mods.bam"
        mapped="?"; mapped_int=0
        if [[ -f "$bam" ]]; then
            mapped=$(samtools view -c -F 4 "$bam" 2>/dev/null || echo "?")
            [[ "$mapped" =~ ^[0-9]+$ ]] && mapped_int=$mapped
        fi
        # Vérification cohérence Dorado vs reads FASTQ (step 01 QC)
        dorado_warn=""
        _qc_cov="${RESULTS_DIR}/01_qc/${BARCODE}_${sn}/coverage_estimate.txt"
        if [[ -f "$_qc_cov" ]]; then
            _fastq_reads=$(grep 'filtered_reads=' "$_qc_cov" | cut -d= -f2)
            if [[ "$_fastq_reads" =~ ^[0-9]+$ && "$_fastq_reads" -gt 0 && "$mapped_int" -gt 0 ]]; then
                _ratio=$(awk "BEGIN{printf \"%.3f\", ${mapped_int}/${_fastq_reads}}")
                if awk "BEGIN{exit !($_ratio < 0.10)}"; then
                    dorado_warn="<br><span style='color:#c0392b;font-weight:bold'>⚠ Seulement ${mapped_int}/${_fastq_reads} reads alignes ($(awk "BEGIN{printf \"%.1f\", $_ratio*100}")%%) — BAM Dorado probablement mal demultiplexe. Relancer <code>submit_modbase_prep.sh</code>.</span>"
                fi
            fi
        fi
        printf '<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s%s</td></tr>\n' \
            "$BARCODE" "$sn" "$n6mA" "$n5mC" "$mapped" "$dorado_warn"
    fi
done
[[ "$MODBASE_HTML_ANY" == "false" ]] && \
    printf '<tr><td colspan="5"><span class="na">(step 07 not run)</span></td></tr>\n'
printf '</table>\n'

printf '</body>\n</html>\n'
} > "$HTML_REPORT"

cat "$SUMMARY_TXT"
log "Summary written to  : ${SUMMARY_TSV}"
log "Text report         : ${SUMMARY_TXT}"
log "HTML report         : ${HTML_REPORT}"
log "Step 06 done."
