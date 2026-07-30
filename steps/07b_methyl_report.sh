#!/usr/bin/env bash
# =============================================================================
# Step 07 — Methylation analysis (MicrobeMod-centric)
# =============================================================================
# What this step does (each barcode analysed independently):
#   1. Global methylation statistics per barcode (modkit pileup BED)
#   2. MicrobeMod motif-level R-M system analysis per barcode
#   3. Functional annotation of motif sites via Pharokka GFF (step 05)
#   4. Self-contained HTML summary report
#
# Output: $RESULTS_DIR/07_methyl/
#   01_stats/         global modkit stats per barcode (TSV + per-sample TXT)
#   02_microbemod/    MicrobeMod call_methylation output per barcode
#   04_annotation/    motif sites annotated with nearest Pharokka gene
#   methylation_report.html
#
# Config parameters consumed (Section 4b of config.sh):
#   COV_MIN_METHYL       minimum coverage for modkit stats (default 8)
#
# Usage:
#   bash steps/07b_methyl_report.sh [--resume]
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/config.sh"

# --resume: skip blocs whose outputs already exist
RESUME=false
for _arg in "$@"; do [[ "$_arg" == "--resume" ]] && RESUME=true; done

step_banner "07 — Methylation analysis (MicrobeMod)${RESUME:+ (--resume)}"

# =============================================================================
# 0. Output paths and log
# =============================================================================
METHYL_OUT="${RESULTS_DIR}/07_methyl"
STATS_DIR="${METHYL_OUT}/01_stats"
MM_DIR="${METHYL_OUT}/02_microbemod"
ANN_DIR="${METHYL_OUT}/04_annotation"
HTML="${METHYL_OUT}/methylation_report.html"
LOGFILE="${METHYL_OUT}/07_methyl.log"

mkdir -p "${STATS_DIR}" "${MM_DIR}" "${ANN_DIR}"

if [[ "${RESUME}" == "true" && -f "${LOGFILE}" ]]; then
    echo "" >> "${LOGFILE}"
    echo "[$(date '+%H:%M:%S')] ── RESUME ─────────────────────────────────────" >> "${LOGFILE}"
else
    : > "${LOGFILE}"
fi

log_m()  { echo "[$(date '+%H:%M:%S')] $*" | tee -a "${LOGFILE}"; }
warn_m() { echo "[WARN]  $*" | tee -a "${LOGFILE}" >&2; }
die_m()  { echo "[ERROR] $*" | tee -a "${LOGFILE}" >&2; exit 1; }

# resume_skip FILE [FILE…] → 0 (skip) if --resume and all files non-empty
resume_skip() {
    [[ "${RESUME}" == "false" ]] && return 1
    for _rf in "$@"; do [[ ! -s "${_rf}" ]] && return 1; done
    return 0
}

log_m "Results dir : ${METHYL_OUT}"
log_m "Config      : COV_MIN_METHYL=${COV_MIN_METHYL:-8}"

# =============================================================================
# 1. Tool checks
# =============================================================================
activate_env
require_tool bedtools
require_tool python3

# MicrobeMod and dependencies: cluster modules
if command -v module &>/dev/null; then
    log_m "Loading MicrobeMod cluster modules..."

    # Load cluster modules including modkit.
    # MicrobeMod 1.1.0 requires modkit 0.2.x (uses --only-tabs flag removed in 0.3.0+).
    # The cluster modkit is typically 0.3.0+ which breaks MicrobeMod.
    # SOLUTION: load all cluster modules (including modkit) so that any module metadata /
    # environment variables are set up correctly, then PREPEND the conda bin to PATH below
    # so the conda modkit 0.2.6 takes priority over the cluster binary.
    # (Simply not loading modkit from modules can cause MicrobeMod to fail if it checks
    #  loaded module dependencies or expects modkit in a specific module-set environment.)
    module load prodigal blast+ hmmer cath-tools modkit meme MicrobeMod 2>>"${LOGFILE}" \
        && log_m "  Cluster modules loaded." \
        || warn_m "  Some cluster modules failed to load — check module availability."

    # PYTHONPATH patch: MicrobeMod 1.1.0 has min_coverage=10 hardcoded in read_modkit().
    # A patched copy (min_coverage=1) lives in ~/MicrobeMod_patched.
    if [[ -d "${HOME}/MicrobeMod_patched" ]]; then
        export PYTHONPATH="${HOME}/MicrobeMod_patched:${PYTHONPATH:-}"
        log_m "  PYTHONPATH patch: using ~/MicrobeMod_patched (min_coverage=1)"
    else
        log_m "  [NOTE] ~/MicrobeMod_patched not found — using default min_coverage=10"
    fi

    # Prepend conda env bin to PATH so the env's modkit 0.2.6 takes priority
    # over any cluster modkit that might have slipped through.
    export PATH="${MAMBA_ROOT_PREFIX}/envs/${MAMBA_ENV}/bin:${PATH}"
    log_m "  PATH: conda bin prepended (modkit=$(command -v modkit 2>/dev/null || echo 'not found'))"
    log_m "  MicrobeMod PATH: $(module show MicrobeMod 2>&1 | grep -i 'PATH\|bin' | head -2 || echo 'could not inspect')"
else
    warn_m "'module' command not found — are you on a Maestro node?"
fi

MICROBEMOD_BIN=""
for _c in MicrobeMod microbemod; do
    command -v "${_c}" &>/dev/null && { MICROBEMOD_BIN="${_c}"; break; }
done
[[ -z "${MICROBEMOD_BIN}" ]] && die_m "MicrobeMod not found in PATH after module load."
log_m "MicrobeMod : ${MICROBEMOD_BIN} — $("${MICROBEMOD_BIN}" --version 2>&1 | head -1)"

# STREME path (needed by MicrobeMod for motif finding)
STREME_PATH=""
if command -v streme &>/dev/null; then
    STREME_PATH="$(command -v streme)"
else
    STREME_PATH=$(find /opt/gensoft/exe/meme /opt/gensoft/exe/MEME 2>/dev/null \
                  -name "streme" -type f 2>/dev/null | head -1 || true)
    [[ -n "${STREME_PATH}" ]] \
        && log_m "  STREME: ${STREME_PATH}" \
        || warn_m "  STREME not found — MicrobeMod may produce partial results"
fi

COV_MIN="${COV_MIN_METHYL:-8}"

# =============================================================================
# 2. Parse sample sheet → list of barcodes (each analysed independently)
# =============================================================================
SAMPLE_SHEET="${RESULTS_DIR}/sample_sheet.tsv"
[[ -f "${SAMPLE_SHEET}" ]] || die_m "sample_sheet.tsv not found: ${SAMPLE_SHEET}"
log_m "Parsing sample sheet..."

declare -a ALL_BCS=()
declare -A BC_NAME=()

while IFS=$'\t' read -r bc name mode _rest; do
    [[ -z "$bc" || "$bc" == "barcode" ]] && continue
    BC_NAME["$bc"]="$name"
    METHYL_BED="${RESULTS_DIR}/06_modbase/${bc}_${name}/methylation.bed"
    [[ -f "$METHYL_BED" ]] || { warn_m "  No methylation.bed for ${bc} — skipping"; continue; }
    ALL_BCS+=("$bc")
done < <(grep -v '^#' "${SAMPLE_SHEET}" | tail -n +2)

if [[ ${#ALL_BCS[@]} -eq 0 ]]; then
    die_m "No barcodes with methylation.bed found under ${RESULTS_DIR}/06_modbase/<bc>_<name>/.
       → Did step 06 (modbase) run and write there? Note the pipeline was renumbered:
         older runs wrote to 07_modbase/ (modbase) and 08_methyl/ (methyl). If your data
         is under those old names, re-run 'submit_all.sh --step 06' so it lands in 06_modbase/."
fi

log_m "  Barcodes with data : ${#ALL_BCS[@]}"

# =============================================================================
# 3. Bloc 1 — Global methylation statistics
# =============================================================================
log_m ""
log_m "══ Bloc 1: Global methylation statistics ══════════════════"

STATS_SUMMARY_TSV="${STATS_DIR}/01_summary.tsv"
_b1_files=("${STATS_SUMMARY_TSV}")
for bc in "${ALL_BCS[@]}"; do _b1_files+=("${STATS_DIR}/${bc}_stats.txt"); done

if resume_skip "${_b1_files[@]}"; then
    log_m "  → already complete — skipped (--resume)"
else
    printf "barcode\tsample\tmod\tdescription\ttotal_sites\tcovered\t" \
        > "${STATS_SUMMARY_TSV}"
    printf "mean_pct_all\tmean_pct_methylated\tzero\thypo\tinter\thyper\tpct_hyper\n" \
        >> "${STATS_SUMMARY_TSV}"

    for bc in "${ALL_BCS[@]}"; do
        name="${BC_NAME[$bc]}"
        BED="${RESULTS_DIR}/06_modbase/${bc}_${name}/methylation.bed"
        OUT_PER="${STATS_DIR}/${bc}_stats.txt"
        log_m "  → ${bc} (${name})..."

        RAW=$(awk -v c="${COV_MIN}" '
        $10 >= c {
            mod=$4; frac=$11
            count_all[mod]++; sum_frac[mod]+=frac
            if (frac>0) { sum_frac_meth[mod]+=frac; count_meth[mod]++ }
            if      (frac==0)  zero[mod]++
            else if (frac<25)  hypo[mod]++
            else if (frac<=75) inter[mod]++
            else               hyper[mod]++
        }
        $10 < c { uncov[$4]++ }
        END {
            for (m in count_all) {
                cov=count_all[m]; tot=cov+uncov[m]+0
                m_all=(cov>0?sum_frac[m]/cov:0)
                m_met=(count_meth[m]>0?sum_frac_meth[m]/count_meth[m]:0)
                p_hyp=(cov>0?(hyper[m]+0)*100/cov:0)
                printf "%s\t%d\t%d\t%.2f\t%.2f\t%d\t%d\t%d\t%d\t%.2f\n",
                       m,tot,cov,m_all,m_met,zero[m]+0,hypo[m]+0,inter[m]+0,hyper[m]+0,p_hyp
            }
        }' "${BED}" | sort -k1,1)

        {
            printf "Barcode : %s\nSample  : %s\nCovMin  : %d\n\n" "${bc}" "${name}" "${COV_MIN}"
            printf "%-5s  %-14s  %-12s  %-10s  %-10s  %-6s  %-6s  %-6s  %-6s  %-8s\n" \
                   "Mod" "Description" "Covered" "MeanAll%" "MeanMeth%" "Zero" "Hypo" "Inter" "Hyper" "%Hyper"
            echo "----------------------------------------------------------------------"
            while IFS=$'\t' read -r mod tot cov m_all m_met zero hypo inter hyper p_hyp; do
                case "$mod" in a) desc="6mA";; m) desc="5mC";; 21839) desc="4mC";; h) desc="5hmC";; *) desc="$mod";; esac
                printf "%-5s  %-14s  %-12d  %-10.2f  %-10.2f  %-6d  %-6d  %-6d  %-6d  %-8.2f\n" \
                       "$mod" "$desc" "$cov" "$m_all" "$m_met" "$zero" "$hypo" "$inter" "$hyper" "$p_hyp"
            done <<< "$RAW"
        } > "${OUT_PER}"

        while IFS=$'\t' read -r mod tot cov m_all m_met zero hypo inter hyper p_hyp; do
            case "$mod" in a) desc="6mA";; m) desc="5mC";; 21839) desc="4mC";; h) desc="5hmC";; *) desc="$mod";; esac
            printf "%s\t%s\t%s\t%s\t%d\t%d\t%.2f\t%.2f\t%d\t%d\t%d\t%d\t%.2f\n" \
                   "$bc" "$name" "$mod" "$desc" "$tot" "$cov" \
                   "$m_all" "$m_met" "$zero" "$hypo" "$inter" "$hyper" "$p_hyp" \
                   >> "${STATS_SUMMARY_TSV}"
        done <<< "$RAW"
        log_m "    ✓ ${OUT_PER}"
    done
fi  # end Bloc 1

# =============================================================================
# 4. Bloc 2 — Reuse MicrobeMod outputs produced by 07a (per barcode)
# =============================================================================
# IMPORTANT: this gather step MUST NOT run `MicrobeMod call_methylation` itself.
# That call (modkit pileup + STREME + blast/hmmer on a very high-coverage BAM) is
# the slow part and is deliberately fanned out to the 07a job array — one barcode
# per task, in parallel, each with its OWN hard timeout and (07a) subsampling.
# Running it here, serially and on the FULL BAM with no timeout, is exactly what
# made 07b hit the 4 h SLURM wall (TIMEOUT). So here we ONLY reuse 07a's TSV:
# rebuild the BED if needed, and record barcodes whose 07a call is missing so the
# failure is surfaced loudly (never silently omitted) in the final summary.
log_m ""
log_m "══ Bloc 2: Reuse MicrobeMod outputs from 07a ══════════════"

# Two distinct no-TSV outcomes, kept separate so a legitimate empty result is NOT
# treated as a failure:
#   _MM_EMPTY   → 07a completed OK but MicrobeMod found no methylation above
#                 threshold (rc=0, no TSV — a valid biological result, e.g. a
#                 phage with no active methyltransferase). NOT a failure.
#   _MM_MISSING → 07a genuinely failed/timed out or was killed (STATUS FAILED/
#                 TIMEOUT, or no sentinel at all). A real problem to surface.
declare -a _MM_MISSING=()
declare -a _MM_EMPTY=()

for bc in "${ALL_BCS[@]}"; do
    name="${BC_NAME[$bc]}"
    MM_OUT="${MM_DIR}/${bc}_${name}"
    # BED: contig  start  end  motif  mod_type  strand
    MM_BED="${MM_OUT}/microbemod.bed"

    mkdir -p "${MM_OUT}"

    if resume_skip "${MM_BED}"; then
        log_m "  ${bc}: already done (--resume) — $(wc -l < "${MM_BED}") motif sites"
        continue
    fi

    # 07a is the ONLY producer of the MicrobeMod TSV. If it is absent, decide
    # whether that is a valid empty result or a real failure from the sentinel.
    MM_TSV=$(find "${MM_OUT}" -name "*.tsv" 2>/dev/null | head -1 || true)
    if [[ -z "${MM_TSV}" || ! -s "${MM_TSV}" ]]; then
        _st="${MM_OUT}/microbemod.STATUS"
        _st_status=""
        [[ -f "${_st}" ]] && _st_status=$(grep -m1 $'^status\t' "${_st}" | cut -f2)
        if [[ "${_st_status}" == "OK" ]]; then
            # MicrobeMod ran to completion but found no methylation above threshold.
            log_m "  ${bc}: MicrobeMod completed with NO methylation above threshold — 0 motifs (valid empty result)."
            : > "${MM_BED}"   # empty BED so downstream blocs treat it as 0 sites, not missing
            _MM_EMPTY+=("${bc}")
        else
            warn_m "  ${bc}: no MicrobeMod TSV and 07a did not complete OK — motif analysis MISSING."
            if [[ -n "${_st_status}" ]]; then
                warn_m "    → 07a sentinel: ${_st_status} (rc=$(grep -m1 $'^rc\t' "${_st}" | cut -f2))"
            else
                warn_m "    → no 07a sentinel — task probably hit its SLURM walltime (TIME_METHYL_CALL). Check logs/07a_methyl_call.*.out"
            fi
            warn_m "    → NOT re-running call_methylation here (that is 07a's job). Global modkit stats stay in the report."
            _MM_MISSING+=("${bc}")
        fi
        continue
    fi

    # TSV exists → (re)build the BED from it if 07a did not already, then continue.
    log_m "  ${bc}: reusing 07a TSV → $(basename "${MM_TSV}")"
    if [[ ! -s "${MM_BED}" ]]; then
        # Build BED from TSV (col3=Contig, col4=Pos, col5=Mod, col6=Strand, col19=Motif)
        awk 'NR>1 && $6!="." && $5!="" && $19!="" {
            printf "%s\t%d\t%d\t%s\t%s\t%s\n", $3, $4, $4+1, $19, $5, $6
        }' "${MM_TSV}" > "${MM_BED}"
    fi
    if [[ -s "${MM_BED}" ]]; then
        log_m "    → $(wc -l < "${MM_BED}") methylated motif sites  →  ${MM_BED}"
    else
        warn_m "  ${bc}: TSV present but no motif sites passed the BED filter."
    fi
done

# =============================================================================
# 5. Bloc 3 — Functional annotation of MicrobeMod motif sites (Pharokka GFF)
# =============================================================================
log_m ""
log_m "══ Bloc 3: Functional annotation ══════════════════════════"

# R-M / modification keywords — defined here (outside resume guard) because also used in HTML
# Phage-relevant terms: methyltransferases, restriction enzymes, anti-restriction,
# glucosyl modification (T4-like phages), and general defense/modification machinery.
RM_KEYWORDS="methyltransferase|restriction|modification|defense|dam|dcm|\
adenine.methylase|cytosine.methylase|type.I|type.II|type.III|type.IV|\
antirestriction|anti.restriction|glucosyl|glucosyltransferase|hydroxymethyl|\
dgt|alc|arn|ocr|kinase"

ANN_SUMMARY_TSV="${ANN_DIR}/annotated_motifs.tsv"

if resume_skip "${ANN_SUMMARY_TSV}"; then
    log_m "  → already complete — skipped (--resume)"
else
    printf "barcode\tsample\tcontig\tpos\tmotif\tmod_type\tstrand\t" > "${ANN_SUMMARY_TSV}"
    printf "g_contig\tg_start\tg_end\tg_strand\tg_locus\tg_product\tdist\n" >> "${ANN_SUMMARY_TSV}"

    # Cache: build gene BED for each barcode that has a Pharokka GFF
    declare -A _GENE_BED_DONE=()

    for bc in "${ALL_BCS[@]}"; do
        name="${BC_NAME[$bc]}"
        MM_BED="${MM_DIR}/${bc}_${name}/microbemod.bed"
        [[ -s "${MM_BED}" ]] || { warn_m "  No microbemod.bed for ${bc} — skip annotation"; continue; }

        # Pharokka produces <prefix>.gff (not .gff3)
        PHAROKKA_GFF="${RESULTS_DIR}/05_annotated/${bc}_${name}/${name}.gff"
        if [[ ! -f "${PHAROKKA_GFF}" ]]; then
            warn_m "  No Pharokka GFF for ${bc} — skipping annotation"
            warn_m "  Expected: ${PHAROKKA_GFF}"
            continue
        fi

        # Build gene BED once per barcode
        # Pharokka GFF attributes example:
        #   ID=CDS_0001;phrog=phrog_1234;product=DNA methyltransferase;function=DNA, RNA and nucleotide metabolism
        # We extract: ID (locus tag), product, and optionally function for R-M context.
        GENE_BED="${ANN_DIR}/${bc}_genes.bed"
        if [[ -z "${_GENE_BED_DONE[$bc]:-}" || ! -s "${GENE_BED}" ]]; then
            awk 'BEGIN{OFS="\t"} $3=="CDS" {
                locus=""; product=""; func=""
                n=split($9,a,";")
                for(i=1;i<=n;i++){
                    if(a[i]~/^ID=/)       { gsub(/^ID=/,"",a[i]);       locus=a[i] }
                    if(a[i]~/^product=/)  { gsub(/^product=/,"",a[i]);  product=a[i] }
                    if(a[i]~/^function=/) { gsub(/^function=/,"",a[i]); func=a[i] }
                }
                # Combine product + function so both are searchable in grep -iE step
                combined = (func != "" && func != "unknown function") ? product " [" func "]" : product
                gsub(/\t|\r/," ",locus);    if(locus=="")    locus="."
                gsub(/\t|\r/," ",combined); if(combined=="") combined="."
                print $1,$4-1,$5,locus,".",$7,combined
            }' "${PHAROKKA_GFF}" > "${GENE_BED}"
            _GENE_BED_DONE["$bc"]=1
        fi

        # bedtools closest — MM_BED (6 cols) vs GENE_BED (7 cols)
        # Output cols: 1-6 from MM_BED, 7-13 from GENE_BED, 14=distance
        bedtools sort -i "${MM_BED}" \
        | bedtools closest -a stdin -b "${GENE_BED}" -D ref -t first \
        | awk -v bc="${bc}" -v nm="${name}" 'BEGIN{OFS="\t"} {
            # $1=contig $2=start $3=end $4=motif $5=mod_type $6=strand
            # $7=g_contig $8=g_start $9=g_end $10=g_locus $11=. $12=g_strand $13=g_product $14=dist
            print bc, nm, $1, $2, $4, $5, $6, $7, $8, $9, $12, $10, $13, $14
        }' >> "${ANN_SUMMARY_TSV}"

        log_m "  ${bc}: annotation done"
    done

    # R-M candidate subset (rows whose gene product matches R-M keywords)
    RM_HITS="${ANN_DIR}/rm_candidates.tsv"
    head -1 "${ANN_SUMMARY_TSV}" > "${RM_HITS}"
    grep -iE "${RM_KEYWORDS}" "${ANN_SUMMARY_TSV}" >> "${RM_HITS}" || true
    N_RM=$(( $(wc -l < "${RM_HITS}") - 1 ))
    N_ANN=$(( $(wc -l < "${ANN_SUMMARY_TSV}") - 1 ))
    log_m "  Annotation: ${N_ANN} total sites → ${ANN_SUMMARY_TSV}"
    log_m "  R-M candidates: ${N_RM} sites → ${RM_HITS}"
fi  # end Bloc 3

# =============================================================================
# 6. Bloc 4 — HTML summary report
# =============================================================================
log_m ""
log_m "══ Bloc 4: Generating HTML report ═════════════════════════"

# Compute summary stats from output files (always, for final log + HTML)
TOTAL_BCS="${#ALL_BCS[@]}"

# Total MicrobeMod motif sites across all barcodes
N_MOTIF_TOTAL=0
declare -A BC_N_MOTIF=()
for bc in "${ALL_BCS[@]}"; do
    name="${BC_NAME[$bc]}"
    MM_BED="${MM_DIR}/${bc}_${name}/microbemod.bed"
    if [[ -s "${MM_BED}" ]]; then
        n=$(wc -l < "${MM_BED}")
        BC_N_MOTIF["$bc"]="$n"
        N_MOTIF_TOTAL=$(( N_MOTIF_TOTAL + n ))
    else
        BC_N_MOTIF["$bc"]=0
    fi
done

N_RM=0
RM_HITS="${ANN_DIR}/rm_candidates.tsv"
[[ -s "${RM_HITS}" ]] && N_RM=$(( $(wc -l < "${RM_HITS}") - 1 ))

RUN_DATE=$(date '+%Y-%m-%d %H:%M:%S')

# ── Verdict inputs: peak %Hyper across all barcode×mod, and mod types screened ──
# Drives the data-driven verdict box so the report never just shows empty cards.
MAX_HYPER=0; MAX_HYPER_BC="-"; MAX_HYPER_MOD="-"; MODS_PRESENT="none"
if [[ -s "${STATS_SUMMARY_TSV}" ]]; then
    # NB: awk must terminate the line with \n, else `read` hits EOF without a
    # delimiter and returns non-zero — under `set -e` that silently kills the
    # whole report step (regression 2026-07-09). `|| true` is a second guard so
    # a future missing-newline never blanks the report again.
    read -r MAX_HYPER MAX_HYPER_BC MAX_HYPER_MOD < <(awk -F'\t' '
        NR>1 && $13+0>m { m=$13+0; bc=$1; md=$3 }
        END {
            lbl=(md=="a"?"6mA":md=="m"?"5mC":md=="21839"?"4mC":md=="h"?"5hmC":md)
            printf "%.2f %s %s\n", m+0, (bc==""?"-":bc), (lbl==""?"-":lbl)
        }' "${STATS_SUMMARY_TSV}") || true
    MODS_PRESENT=$(awk -F'\t' '
        NR>1 { m=$3; lbl=(m=="a"?"6mA":m=="m"?"5mC":m=="21839"?"4mC":m=="h"?"5hmC":m)
               if(!(lbl in seen)){ seen[lbl]=1; order[++n]=lbl } }
        END { for(i=1;i<=n;i++) printf "%s%s",(i>1?", ":""),order[i] }' "${STATS_SUMMARY_TSV}")
    [[ -z "${MODS_PRESENT}" ]] && MODS_PRESENT="none"
fi
# %Hyper (% of covered sites hyper-methylated) above which a mod counts as real signal
VERDICT_MIN_HYPER="${VERDICT_MIN_HYPER:-5}"

if resume_skip "${HTML}"; then
    log_m "  → already complete — skipped (--resume)"
else

{
cat << 'HTML_STYLE'
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Methylation Report</title>
<style>
:root{--acc:#2c7bb6;--grn:#2e7d32;--red:#c62828;--org:#e65100;--bg:#f4f6f8;--card:#fff;}
*{box-sizing:border-box;margin:0;padding:0;}
body{font-family:Arial,Helvetica,sans-serif;font-size:13px;background:var(--bg);color:#222;}
header{background:#1a2638;color:#fff;padding:18px 28px;}
header h1{font-size:18px;font-weight:700;}
header .meta{font-size:11px;color:#9bb;margin-top:4px;}
main{max-width:1200px;margin:0 auto;padding:24px 20px;}
h2{font-size:14px;color:var(--acc);border-bottom:2px solid var(--acc);
   padding-bottom:4px;margin:28px 0 12px;}
h3{font-size:13px;color:#333;margin:16px 0 6px;}
.cards{display:flex;flex-wrap:wrap;gap:12px;margin-bottom:24px;}
.card{background:var(--card);border:1px solid #dde;border-radius:6px;
      padding:12px 18px;min-width:140px;text-align:center;}
.card .v{font-size:26px;font-weight:700;color:var(--acc);}
.card .l{font-size:11px;color:#666;margin-top:2px;}
.note{font-size:11px;color:#666;font-style:italic;margin-bottom:10px;}
.tbl-wrap{overflow-x:auto;margin-bottom:20px;}
table{border-collapse:collapse;width:100%;font-size:12px;}
th{background:#2c4a6e;color:#fff;padding:6px 10px;text-align:left;
   cursor:pointer;white-space:nowrap;user-select:none;}
th:hover{background:#3a5f8a;}
th::after{content:" ⇅";font-size:10px;opacity:.5;}
td{padding:5px 10px;border-bottom:1px solid #eee;white-space:nowrap;}
tr:hover td{background:#f0f6ff;}
.gained{color:var(--grn);font-weight:600;}
.lost  {color:var(--red);font-weight:600;}
.shared{color:#555;}
.hi {background:#e8f5e9;color:var(--grn);padding:1px 7px;border-radius:10px;font-size:11px;}
.me {background:#fff3e0;color:var(--org);padding:1px 7px;border-radius:10px;font-size:11px;}
.lo {background:#fce4ec;color:var(--red);padding:1px 7px;border-radius:10px;font-size:11px;}
.rm {background:#fce4ec;}
.filter{padding:4px 8px;border:1px solid #ccc;border-radius:4px;
        font-size:12px;width:240px;margin-bottom:8px;}
.pos-delta{color:#2e7d32;font-weight:700;}
.neg-delta{color:#c62828;font-weight:700;}
.intro-box{background:#fff;border-left:4px solid var(--acc);border-radius:4px;
           padding:14px 18px;margin-bottom:24px;font-size:12px;line-height:1.7;}
.intro-box h3{font-size:13px;color:var(--acc);margin-bottom:6px;}
.mod-tag{display:inline-block;padding:1px 8px;border-radius:10px;font-size:11px;
         font-weight:700;margin-right:4px;}
.tag-6mA{background:#e3f2fd;color:#1565c0;}
.tag-5mC{background:#e8f5e9;color:#2e7d32;}
.tag-4mC{background:#f3e5f5;color:#6a1b9a;}
.tag-5hmC{background:#fff3e0;color:#e65100;}
.verdict{border-radius:6px;padding:14px 18px;margin-bottom:24px;font-size:13px;line-height:1.6;}
.verdict.none{background:#e8f5e9;border-left:4px solid var(--grn);}
.verdict.signal{background:#fff3e0;border-left:4px solid var(--org);}
.verdict h3{font-size:14px;margin-bottom:6px;}
.verdict.none h3{color:var(--grn);}
.verdict.signal h3{color:var(--org);}
.distrib-bar{display:flex;height:14px;border-radius:3px;overflow:hidden;
             border:1px solid #ddd;width:140px;}
footer{margin-top:40px;padding:10px 0;border-top:1px solid #ddd;
       color:#aaa;font-size:11px;text-align:center;}
</style>
<script>
function sort_tbl(th){
  var tbl=th.closest('table'), col=th.cellIndex;
  var asc=(tbl.dataset.sc==col)?(tbl.dataset.sa!='1'):'1';
  var rows=Array.from(tbl.tBodies[0].rows);
  rows.sort(function(a,b){
    var va=a.cells[col]?a.cells[col].innerText.trim():'';
    var vb=b.cells[col]?b.cells[col].innerText.trim():'';
    var na=parseFloat(va),nb=parseFloat(vb);
    if(!isNaN(na)&&!isNaN(nb)) return asc=='1'?na-nb:nb-na;
    return asc=='1'?va.localeCompare(vb):vb.localeCompare(va);
  });
  rows.forEach(function(r){tbl.tBodies[0].appendChild(r);});
  tbl.dataset.sc=col; tbl.dataset.sa=asc;
}
function filter_tbl(inp,tblId){
  var q=inp.value.toLowerCase();
  document.getElementById(tblId).tBodies[0].querySelectorAll('tr').forEach(function(r){
    r.style.display=r.innerText.toLowerCase().includes(q)?'':'none';
  });
}
</script>
</head>
<body>
HTML_STYLE

cat << HTML_HEAD
<header>
  <h1>DNA Methylation Analysis — ${RUN_NAME}</h1>
  <div class="meta">Generated: ${RUN_DATE} &nbsp;|&nbsp; Coverage threshold: ≥${COV_MIN} reads &nbsp;|&nbsp; Tool: MicrobeMod ${MICROBEMOD_BIN:+v}$("${MICROBEMOD_BIN}" --version 2>&1 | head -1)</div>
</header>
<main>
HTML_HEAD

# ── Biological intro box ──────────────────────────────────────────────────────
cat << 'HTML_INTRO'
<div class="intro-box">
<h3>Biological context — methylations detected by this pipeline</h3>
<p>This pipeline detects three main types of bacterial DNA methylation by direct Nanopore sequencing (Dorado), called in <strong>all sequence contexts</strong> (not restricted to CpG):</p>
<ul style="margin:8px 0 8px 18px;line-height:1.8">
  <li><span class="mod-tag tag-6mA">6mA</span> <strong>N6-methyladenine</strong> (modkit code <code>a</code>) — the most frequent modification in bacteria. Deposited mainly at <strong>GATC</strong> motifs by the Dam methylase, and at other motifs by type I, II and III R-M systems. Key role in mismatch repair (MutHLS system), replication control and gene expression regulation. Loss of 6mA at a given motif in an evolved clone may indicate an inactivating mutation in the corresponding methylase.</li>
  <li><span class="mod-tag tag-5mC">5mC</span> <strong>5-methylcytosine</strong> (modkit code <code>m</code>) — less frequent than 6mA in bacteria. Associated with <strong>CCWGG</strong> motifs (Dcm methylase) and some type II R-M systems. Called here in all contexts (not only CpG), so non-CpG 5mC is visible. Inter-clone 5mC variation can indicate methylase mutations or R-M system rearrangements.</li>
  <li><span class="mod-tag tag-4mC">4mC</span> <strong>N4-methylcytosine</strong> (modkit ChEBI code <code>21839</code>, shown as <code>4mC</code>) — a cytosine modification typical of some restriction-modification systems. Scattered low-frequency 4mC calls are usually background rather than a genuine motif.</li>
  <li><span class="mod-tag tag-5hmC">5hmC</span> <strong>5-hydroxymethylcytosine</strong> (modkit code <code>h</code>) — detected sporadically, often associated with phages or atypical modification pathways. Not called by the current model complex unless enabled; interpret in the relevant biological context.</li>
</ul>
<p>MicrobeMod identifies recognition motifs of active methylases from detected methylated positions. <strong>Section 3 (R-M gene overlaps)</strong> cross-references these motifs with the Pharokka annotation to locate methylation, restriction and defence genes — the key players of R-M systems.</p>
</div>
HTML_INTRO

# ── Summary cards ─────────────────────────────────────────────────────────────
cat << HTML_CARDS
<h2>Summary</h2>
<div class="cards">
  <div class="card"><div class="v">${TOTAL_BCS}</div><div class="l">Barcodes analysed</div></div>
  <div class="card"><div class="v">${N_MOTIF_TOTAL}</div><div class="l">Total motif sites</div></div>
  <div class="card"><div class="v" style="color:var(--red)">${N_RM}</div><div class="l">R-M gene overlaps</div></div>
</div>
HTML_CARDS

# ── Verdict box (data-driven: no significant methylation vs signal present) ─────
if awk "BEGIN{exit !(${MAX_HYPER}+0 < ${VERDICT_MIN_HYPER}+0 && ${N_MOTIF_TOTAL}+0==0)}"; then
cat << HTML_VERDICT
<div class="verdict none">
  <h3>Verdict — no significant methylation detected</h3>
  <p>Across all ${TOTAL_BCS} barcode(s), no methylation motif passed the MicrobeMod threshold
     (<strong>${N_MOTIF_TOTAL}</strong> motif site(s)), and global modkit statistics peak at only
     <strong>${MAX_HYPER}%</strong> hyper-methylated sites (${MAX_HYPER_BC}, ${MAX_HYPER_MOD}) —
     far below the ${VERDICT_MIN_HYPER}% signal threshold.
     Modifications screened in all sequence contexts: <strong>${MODS_PRESENT}</strong>.</p>
  <p style="margin-top:6px">The near-zero <em>% Hyper</em> values in Section 2 are background noise,
     not constitutive methylation. See Section 2 for the per-barcode breakdown.</p>
</div>
HTML_VERDICT
else
cat << HTML_VERDICT
<div class="verdict signal">
  <h3>Verdict — methylation signal present</h3>
  <p>At least one barcode shows methylation worth review: <strong>${N_MOTIF_TOTAL}</strong> MicrobeMod
     motif site(s); peak <em>% Hyper</em> <strong>${MAX_HYPER}%</strong> (${MAX_HYPER_BC}, ${MAX_HYPER_MOD}).
     Modifications screened in all sequence contexts: <strong>${MODS_PRESENT}</strong>.
     See Sections 1 and 3 for details.</p>
</div>
HTML_VERDICT
fi

# ── Section 0: Methods ────────────────────────────────────────────────────────
cat << HTML_METHODS
<h2>0. Methods</h2>
<div style="background:var(--card);border:1px solid #dde;border-radius:6px;padding:14px 18px;font-size:12px;line-height:1.7;margin-bottom:24px;">

<p><strong>Data source:</strong> All methylation analyses are based <em>exclusively on Nanopore reads</em>
(Dorado mod-tagged BAMs). Illumina reads are not used for methylation detection.</p>

<p><strong>Basecalling prerequisite (step 06):</strong>
  Dorado super-accuracy (SUP) re-basecalling, model complex
  <code>${DORADO_MODEL_COMPLEX:-sup@v5.0.0,6mA,4mC_5mC}</code> — calls
  <code>${DORADO_MOD_BASES:-6mA 4mC_5mC}</code> in <strong>all sequence contexts</strong> (not CpG-restricted):
  N6-methyladenine (6mA), N4-methylcytosine (4mC) and 5-methylcytosine (5mC).
  Calling 5mC/4mC in every context (rather than the older CpG-only <code>5mCG_5hmCG</code> model)
  removes the blind spot on non-CpG cytosine methylation (e.g. Dcm CCWGG, type II R-M motifs).
  Demultiplexed per-barcode BAMs placed in <code>MOD_DEMUX_DIR</code>.</p>

<p><strong>Read alignment (step 06):</strong>
  <code>minimap2 -ax map-ont -y</code> — aligns mod-tagged reads to the per-sample polished assembly (step 03).
  The <code>-y</code> flag propagates modification tags into the output BAM (<code>aligned_mods.bam</code>).</p>

<p><strong>Methylation pileup (step 06):</strong>
  <code>modkit pileup</code> with the polished assembly as reference.
  Outputs a BED file (<code>methylation.bed</code>) with per-position modification frequency (0–100% of reads).
  Note: modkit reports 4mC under its ChEBI code <code>21839</code>; this report relabels it <code>4mC</code>.
  Coverage threshold applied in this report: ≥${COV_MIN}× (positions below this are excluded from statistics).</p>

<p><strong>Global methylation statistics (Section 2):</strong>
  Per barcode and modification type, computed from the pileup BED.
  Each covered site is classified by modification frequency:
  absent (0%), low (1–24%), moderate (25–75%), hyper-methylated (&gt;75%).
  <em>% Hyper</em>: proportion of covered sites with &gt;75% modification frequency —
  main indicator of constitutive methylase activity.</p>

<p><strong>MicrobeMod motif identification (Section 1):</strong>
  <code>MicrobeMod call_methylation</code> — identifies recognition motifs of active methylases
  directly from the mod-tagged BAM (using STREME for motif finding). Threads: ${THREADS:-8}.
  Note: MicrobeMod internal <code>min_coverage</code> is patched to 1 (copy in <code>~/MicrobeMod_patched</code>)
  to avoid discarding low-coverage motifs; the ≥${COV_MIN}× threshold is enforced upstream in the pileup BED.</p>

<p><strong>Functional annotation (Section 3):</strong>
  Nearest CDS from per-sample Pharokka GFF (step 05) via <code>bedtools closest</code>.
  R-M candidates flagged by product name matching keywords:
  methyltransferase | restriction | modification | defense | dam | dcm |
  adenine/cytosine methylase | type I/II/III/IV | antirestriction | glucosyl | dgt | alc | arn | ocr | kinase.
  Distance 0 = methylated site is inside the CDS.</p>

</div>
HTML_METHODS

# ── Section 1: MicrobeMod per barcode ─────────────────────────────────────────
cat << HTML_S1
<h2>1. Methylated motifs per sample (MicrobeMod)</h2>
<p class="note">
  One row per barcode &times; motif.
  <strong>Modification:</strong>
  <span class="mod-tag tag-6mA">a = 6mA</span>
  <span class="mod-tag tag-5mC">m = 5mC</span>
  <span class="mod-tag tag-4mC">21839 = 4mC</span>
  <span class="mod-tag tag-5hmC">h = 5hmC</span>.<br>
  <strong>% of total methylated</strong>: proportion of this motif among <em>all</em> methylated sites in the sample.
  A motif at 60% means 60 out of every 100 methylated positions carry that recognition sequence.<br>
  <strong>Methylase activity at motif</strong>: % of recognition sites fully methylated (&gt;75% of reads).
  High = methylase constitutively active; low = partial activity or mutated/absent methylase.
  <span class="hi">High</span> ≥50% · <span class="me">Moderate</span> ≥10% · <span class="lo">Low</span> &lt;10%.
</p>
<input class="filter" placeholder="Filter…" oninput="filter_tbl(this,'tbl-mm')">
<div class="tbl-wrap"><table id="tbl-mm">
<thead><tr>
  <th onclick="sort_tbl(this)">Barcode</th>
  <th onclick="sort_tbl(this)">Sample</th>
  <th onclick="sort_tbl(this)">Motif</th>
  <th onclick="sort_tbl(this)">Mod type</th>
  <th onclick="sort_tbl(this)">N sites</th>
  <th onclick="sort_tbl(this)">% sites / total</th>
  <th onclick="sort_tbl(this)">Confidence (modkit %Hyper)</th>
</tr></thead><tbody>
HTML_S1

for bc in "${ALL_BCS[@]}"; do
    name="${BC_NAME[$bc]}"
    MM_BED="${MM_DIR}/${bc}_${name}/microbemod.bed"

    if [[ ! -s "${MM_BED}" ]]; then
        echo "<tr><td>${bc}</td><td>${name}</td><td colspan='5'><em>No motif above MicrobeMod threshold (0 sites &gt;66%) — see Section 2 for global stats</em></td></tr>"
        continue
    fi

    # Count sites per motif×mod from BED (col4=motif, col5=mod_type)
    TOTAL_SITES=$(wc -l < "${MM_BED}")

    awk -v bc="${bc}" -v name="${name}" \
        -v tot="${TOTAL_SITES}" '{cnt[$4"|"$5]++}
    END {
        for(k in cnt){
            n=split(k,a,"|"); motif=a[1]; mod=a[2]
            pct=sprintf("%.1f%%",cnt[k]*100/tot)
            printf "%s\t%s\t%s\t%s\t%d\t%s\n",
                   bc,name,motif,mod,cnt[k],pct
        }
    }' "${MM_BED}" | sort -t$'\t' -k5,5rn \
    | while IFS=$'\t' read -r bc2 nm2 motif mod nsites pct_sites; do

        # Look up pct_hyper from modkit stats for this barcode × mod type
        pct_hyp=$(awk -v bc="${bc}" -v mod="${mod}" \
                  '$1==bc && $3==mod {print $13; exit}' \
                  "${STATS_SUMMARY_TSV}" 2>/dev/null || echo ".")
        # Assign confidence
        if [[ "${pct_hyp}" != "." ]]; then
            h=$(awk "BEGIN{print (${pct_hyp}+0 >= 50) ? 1 : 0}" 2>/dev/null || echo 0)
            m=$(awk "BEGIN{print (${pct_hyp}+0 >= 10) ? 1 : 0}" 2>/dev/null || echo 0)
            if   [[ "$h" == "1" ]]; then conf_cls="hi"; conf_lbl="High"
            elif [[ "$m" == "1" ]]; then conf_cls="me"; conf_lbl="Medium"
            else                         conf_cls="lo"; conf_lbl="Low"
            fi
            conf_txt="${pct_hyp}%"
        else
            conf_cls="lo"; conf_lbl="N/A"; conf_txt="."
        fi
        echo "<tr><td>${bc}</td><td>${nm2}</td><td><code>${motif}</code></td><td>${mod}</td><td>${nsites}</td><td>${pct_sites}</td><td><span class=\"${conf_cls}\">${conf_lbl}</span> (${conf_txt})</td></tr>"
    done
done

echo "</tbody></table></div>"

# ── Section 2: Global methylation stats (modkit) ──────────────────────────────
cat << HTML_S3
<h2>2. Global methylation statistics (modkit)</h2>
<p class="note">
  Per barcode and modification type. Coverage threshold: ≥${COV_MIN} reads.<br>
  <strong>% Hyper (&gt;75%)</strong>: proportion of covered sites methylated at &gt;75% — indicates constitutive methylation (methylase active at saturation).<br>
  <strong>Profile</strong>: visual distribution — hover over segments for percentages.
  <span style="display:inline-flex;gap:6px;vertical-align:middle;margin-left:4px">
    <span style="display:flex;align-items:center;gap:3px"><span style="display:inline-block;width:12px;height:10px;background:#9e9e9e;border-radius:2px"></span>0% (absent)</span>
    <span style="display:flex;align-items:center;gap:3px"><span style="display:inline-block;width:12px;height:10px;background:#ffee58;border-radius:2px"></span>1-24% (low)</span>
    <span style="display:flex;align-items:center;gap:3px"><span style="display:inline-block;width:12px;height:10px;background:#ff9800;border-radius:2px"></span>25-75% (moderate)</span>
    <span style="display:flex;align-items:center;gap:3px"><span style="display:inline-block;width:12px;height:10px;background:#c62828;border-radius:2px"></span>&gt;75% (high)</span>
  </span>
</p>
<input class="filter" placeholder="Filter…" oninput="filter_tbl(this,'tbl-stats')">
<div class="tbl-wrap"><table id="tbl-stats">
<thead><tr>
  <th onclick="sort_tbl(this)">Barcode</th>
  <th onclick="sort_tbl(this)">Sample</th>
  <th onclick="sort_tbl(this)" title="a=6mA  m=5mC  21839=4mC  h=5hmC">Mod</th>
  <th onclick="sort_tbl(this)">Description</th>
  <th onclick="sort_tbl(this)" title="Positions with ≥${COV_MIN} aligned reads">Covered sites</th>
  <th onclick="sort_tbl(this)" title="Mean frequency over all covered sites">Mean % (all)</th>
  <th onclick="sort_tbl(this)" title="Mean frequency among non-zero sites only">Mean % (methylated only)</th>
  <th onclick="sort_tbl(this)" title="% of covered sites with >75% modified reads">% Hyper (&gt;75%)</th>
  <th title="Distribution: Absent | Low | Moderate | High — hover for percentages">Profile</th>
</tr></thead><tbody>
HTML_S3

while IFS=$'\t' read -r bc nm mod desc tot cov m_all m_met zero hypo inter hyper p_hyp; do
    [[ "$bc" == "barcode" ]] && continue
    distrib_bar=$(awk -v z="$zero" -v o="$hypo" -v i="$inter" -v h="$hyper" 'BEGIN{
        tot=z+o+i+h
        if (tot==0) { print "<span style=\"color:#aaa;font-size:11px\">n/a</span>"; exit }
        W=140
        wz=int(z/tot*W); wo=int(o/tot*W); wi=int(i/tot*W); wh=W-wz-wo-wi
        pz=z/tot*100; po=o/tot*100; pi=i/tot*100; ph=h/tot*100
        printf "<div class=\"distrib-bar\">"
        if (wz>0) printf "<div style=\"width:%dpx;background:#9e9e9e\" title=\"Absent (0%%): %.1f%%\"></div>",wz,pz
        if (wo>0) printf "<div style=\"width:%dpx;background:#ffee58\" title=\"Low (1-24%%): %.1f%%\"></div>",wo,po
        if (wi>0) printf "<div style=\"width:%dpx;background:#ff9800\" title=\"Moderate (25-75%%): %.1f%%\"></div>",wi,pi
        if (wh>0) printf "<div style=\"width:%dpx;background:#c62828\" title=\"High (>75%%): %.1f%%\"></div>",wh,ph
        printf "</div>"
    }')
    echo "<tr><td>${bc}</td><td>${nm}</td><td><code>${mod}</code></td><td>${desc}</td>"
    echo "<td>${cov}</td><td>${m_all}</td><td>${m_met}</td><td>${p_hyp}</td>"
    echo "<td>${distrib_bar}</td></tr>"
done < "${STATS_SUMMARY_TSV}"

echo "</tbody></table></div>"

# ── Section 3: R-M gene annotation ────────────────────────────────────────────
cat << HTML_S4
<h2>3. Overlaps with R-M genes (restriction-modification)</h2>
<p class="note">
  Methylated sites (MicrobeMod) whose nearest gene carries a Pharokka annotation associated with R-M systems
  (keywords: methyltransferase, restriction, modification, defense, dam, dcm, antirestriction, glucosyl, dgt, alc, arn, ocr, kinase).<br>
  <strong>Distance</strong>: bp between the methylated site and the nearest CDS
  (0 = site is <em>inside</em> the gene; positive = downstream; negative = upstream).<br>
  Full table: <code>04_annotation/annotated_motifs.tsv</code>.
</p>
HTML_S4

if [[ ! -s "${RM_HITS}" ]] || [[ $(wc -l < "${RM_HITS}") -le 1 ]]; then
    echo "<p class=\"note\">No R-M system gene overlaps found.</p>"
else
    echo "<input class=\"filter\" placeholder=\"Filter…\" oninput=\"filter_tbl(this,'tbl-rm')\">"
    echo "<div class=\"tbl-wrap\"><table id=\"tbl-rm\"><thead><tr>"
    echo "<th onclick=\"sort_tbl(this)\">Barcode</th>"
    echo "<th onclick=\"sort_tbl(this)\">Sample</th>"
    echo "<th onclick=\"sort_tbl(this)\">Motif</th>"
    echo "<th onclick=\"sort_tbl(this)\">Mod</th>"
    echo "<th onclick=\"sort_tbl(this)\">Gene locus</th>"
    echo "<th onclick=\"sort_tbl(this)\">Product</th>"
    echo "<th onclick=\"sort_tbl(this)\">Dist (bp)</th>"
    echo "</tr></thead><tbody>"

    tail -n +2 "${RM_HITS}" \
    | while IFS=$'\t' read -r bc nm contig pos motif mod strand \
          g_contig g_start g_end g_strand g_locus g_product dist; do
        echo "<tr class=\"rm\"><td>${bc}</td><td>${nm}</td>"
        echo "<td><code>${motif}</code></td><td>${mod}</td>"
        echo "<td><code>${g_locus}</code></td><td>${g_product}</td>"
        echo "<td>${dist}</td></tr>"
    done
    echo "</tbody></table></div>"
fi

# ── Footer ────────────────────────────────────────────────────────────────────
cat << HTML_FOOT
<footer>
  Step 07 — Institut Pasteur / Maestro — Run: ${RUN_NAME} | ${RUN_DATE}
</footer>
</main>
</body>
</html>
HTML_FOOT

} > "${HTML}"

log_m "  HTML report → ${HTML}"

fi  # end Bloc 4 resume guard

# =============================================================================
# 8. Final summary
# =============================================================================
# ── Surface MicrobeMod outcomes (non-silent), distinguishing failure from empty ─
# _MM_MISSING (Bloc 2) = REAL failures only: 07a FAILED/TIMEOUT or was killed with
#   no sentinel. Decorate each with its sentinel detail (status/rc) if present.
# _MM_EMPTY (Bloc 2) = 07a completed OK but MicrobeMod found no methylation above
#   threshold — a valid biological result, NOT a failure and NEVER fatal.
declare -a _MM_FAILED=()
for _mbc in "${_MM_MISSING[@]:-}"; do
    [[ -z "${_mbc}" ]] && continue
    _st=$(find "${MM_DIR}" -maxdepth 2 -name microbemod.STATUS -path "*/${_mbc}_*/*" 2>/dev/null | head -1)
    if [[ -n "${_st}" && -f "${_st}" ]]; then
        _MM_FAILED+=("${_mbc} ($(grep -m1 $'^status\t' "${_st}" | cut -f2), rc=$(grep -m1 $'^rc\t' "${_st}" | cut -f2))")
    else
        _MM_FAILED+=("${_mbc} (NO OUTPUT — 07a killed/absent, no sentinel)")
    fi
done

log_m ""
log_m "════════════════════════════════════════════════════════"
log_m "  Step 07 complete"
log_m "  Barcodes analysed    : ${TOTAL_BCS}"
log_m "  Motif sites total    : ${N_MOTIF_TOTAL}"
log_m "  R-M gene overlaps    : ${N_RM}"
log_m "  Report               : ${HTML}"
if [[ ${#_MM_EMPTY[@]} -gt 0 ]]; then
    log_m  "  No methylation found : ${#_MM_EMPTY[@]}/${TOTAL_BCS} barcode(s) — ${_MM_EMPTY[*]}"
    log_m  "  → MicrobeMod completed but found no sites above threshold (valid empty result, 0 motifs)."
fi
if [[ ${#_MM_FAILED[@]} -gt 0 ]]; then
    warn_m "  MicrobeMod FAILED    : ${#_MM_FAILED[@]}/${TOTAL_BCS} barcode(s) — ${_MM_FAILED[*]}"
    warn_m "  → global modkit stats are still in the report; motif/R-M analysis is missing for those."
fi
log_m "  Copy to desktop   : scp ${EMAIL%%@*}@maestro-submit:${HTML} ~/Desktop/"
log_m "════════════════════════════════════════════════════════"

# Non-silent exit: fail ONLY if every analysed barcode had a REAL MicrobeMod
# failure. A run where MicrobeMod completed everywhere but found no methylation
# (all _MM_EMPTY) is a SUCCESS — the report legitimately shows 0 motifs.
if [[ ${TOTAL_BCS} -gt 0 && ${#_MM_FAILED[@]} -ge ${TOTAL_BCS} ]]; then
    die_m "MicrobeMod FAILED (error/timeout/killed) for ALL ${TOTAL_BCS} barcode(s) — see sentinels in ${MM_DIR}/*/microbemod.STATUS and the modkit self-check in logs/07a_methyl_call.*.out."
fi
