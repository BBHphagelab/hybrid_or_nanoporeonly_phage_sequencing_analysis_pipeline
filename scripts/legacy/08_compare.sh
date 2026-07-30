#!/usr/bin/env bash
# =============================================================================
# Step 08 — Assembly-vs-assembly comparison (nucdiff / MUMmer)  [SLURM array]
# =============================================================================
# Replaces the former read-based variant calling (breseq / snippy on ONT-
# fragmented reads), which produced methylation- and homopolymer-driven false
# positives and could not represent the large structural events expected
# between engineered recombinants. Here we compare FINISHED ASSEMBLIES directly:
# nucdiff (built on MUMmer/nucmer) reports SNPs, small indels AND structural
# rearrangements (insertions/deletions of blocks, relocations, translocations,
# inversions, duplications), with an optional coordinate dotplot.
#
# Per sample, up to two comparisons (driven by the sample sheet):
#   • EXTERNAL  (cmp_ext_ref column): polished assembly vs a published reference
#     (absolute .gbk/.gbff/.fasta path, or shortcut ext/ext1…ext9 → config.sh).
#   • INTERNAL  (cmp_int_bc  column): polished assembly vs a PARENT sample's own
#     polished assembly from this same run (e.g. recombinant vs its backbone).
#
# QUERY  = this sample's 03_polished/<label>/polished.fasta
# Output : 08_compare/<label>/
#            vs_external_<refname>/   nucdiff results + dotplot + summary.tsv
#            vs_internal_<parent>/    nucdiff results + dotplot + summary.tsv
#            compare_index.tsv        one row per comparison (for steps 09/10)
#
# Self-skips per comparison when the relevant column is '-'.
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/config.sh"

activate_env
require_tool nucdiff
require_tool nucmer

SAMPLE_SHEET="${RESULTS_DIR}/sample_sheet.tsv"
[[ -f "${SAMPLE_SHEET}" ]] || die "Sample sheet not found: ${SAMPLE_SHEET}"

IDX="${ARRAY_IDX:-${SLURM_ARRAY_TASK_ID:-0}}"
mapfile -t BARCODES < <(list_barcodes)
if [[ "${IDX}" -ge "${#BARCODES[@]}" ]]; then
    echo "[SKIP] ARRAY_IDX=${IDX} >= barcodes (${#BARCODES[@]}) — nothing to do."
    exit 0
fi
BARCODE="${BARCODES[${IDX}]}"
SAMPLE_NAME=$(get_sample_field "$BARCODE" "sample_name")
INT_REF_BC=$(get_sample_field "$BARCODE" "cmp_int_bc")
EXT_REF_RAW=$(get_sample_field "$BARCODE" "cmp_ext_ref")
LABEL="${BARCODE}_${SAMPLE_NAME}"

step_banner "08 — Compare (nucdiff)  [${BARCODE} / ${SAMPLE_NAME}]"

QUERY="${RESULTS_DIR}/03_polished/${LABEL}/polished.fasta"
OUT="${RESULTS_DIR}/08_compare/${LABEL}"
mkdir -p "$OUT"
INDEX="${OUT}/compare_index.tsv"

# ── Input guard ──────────────────────────────────────────────────────────────
if [[ ! -s "$QUERY" ]]; then
    echo "[SKIP] ${LABEL}: polished.fasta missing — step 03 did not complete. Nothing to do."
    exit 0
fi

if [[ "${INT_REF_BC}" == "-" || -z "${INT_REF_BC}" ]] \
   && [[ "${EXT_REF_RAW}" == "-" || -z "${EXT_REF_RAW}" ]]; then
    echo "[SKIP] ${LABEL}: cmp_int_bc and cmp_ext_ref both '-' — no comparison requested."
    exit 0
fi

printf "comparison\treference\tref_label\tsubstitutions\tinsertions\tdeletions\tstructural\tquery_unaligned_bp\n" > "$INDEX"

# ---------------------------------------------------------------------------
# Resolve an external-ref shortcut (ext / ext1 … ext9) to a path
# ---------------------------------------------------------------------------
resolve_ext_ref() {
    case "$1" in
        ext)  echo "${EXTERNAL_REF:-}"   ;;
        ext1) echo "${EXTERNAL_REF_1:-}" ;;
        ext2) echo "${EXTERNAL_REF_2:-}" ;;
        ext3) echo "${EXTERNAL_REF_3:-}" ;;
        ext4) echo "${EXTERNAL_REF_4:-}" ;;
        ext5) echo "${EXTERNAL_REF_5:-}" ;;
        ext6) echo "${EXTERNAL_REF_6:-}" ;;
        ext7) echo "${EXTERNAL_REF_7:-}" ;;
        ext8) echo "${EXTERNAL_REF_8:-}" ;;
        ext9) echo "${EXTERNAL_REF_9:-}" ;;
        *)    echo "$1" ;;   # treat as a literal path
    esac
}

# ---------------------------------------------------------------------------
# Ensure a FASTA reference: GenBank (.gbk/.gbff) → FASTA via any2fasta.
# Echoes the path to a usable FASTA (in $OUT), or empty on failure.
# ---------------------------------------------------------------------------
ensure_fasta() {
    local src="$1" tag="$2" dst
    case "$src" in
        *.gbk|*.gbff|*.gb)
            dst="${OUT}/ref_${tag}.fasta"
            if command -v any2fasta &>/dev/null && any2fasta "$src" > "$dst" 2>/dev/null && [[ -s "$dst" ]]; then
                echo "$dst"
            else
                echo ""
            fi ;;
        *.fasta|*.fa|*.fna)
            echo "$src" ;;
        *)
            # Unknown extension — try any2fasta, else assume FASTA
            dst="${OUT}/ref_${tag}.fasta"
            if command -v any2fasta &>/dev/null && any2fasta "$src" > "$dst" 2>/dev/null && [[ -s "$dst" ]]; then
                echo "$dst"
            else
                echo "$src"
            fi ;;
    esac
}

# ---------------------------------------------------------------------------
# run_compare <ref_fasta> <ref_label> <comparison_kind> <outdir> [ref_gff]
#   Runs nucdiff + a mummerplot dotplot, parses counts, appends a row to INDEX.
# ---------------------------------------------------------------------------
run_compare() {
    local ref_fasta="$1" ref_label="$2" kind="$3" cdir="$4" ref_gff="${5:-}"
    mkdir -p "$cdir"
    local prefix="cmp"

    log "  [${kind}] ${LABEL}  vs  ${ref_label}"

    local nd_ann=()
    [[ -n "$ref_gff" && -f "$ref_gff" ]] && nd_ann=(--ref_annotation "$ref_gff")

    set +e
    nucdiff --vcf yes "${nd_ann[@]}" "$ref_fasta" "$QUERY" "$cdir" "$prefix" \
        2>&1 | tee "${cdir}/nucdiff.log"
    local nd_exit=${PIPESTATUS[0]}
    set -e
    if [[ ${nd_exit} -ne 0 ]]; then
        warn "    nucdiff failed (exit=${nd_exit}) — see ${cdir}/nucdiff.log"
        printf "%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\n" "$kind" "$ref_fasta" "$ref_label" >> "$INDEX"
        return
    fi

    local stat="${cdir}/results/${prefix}_stat.out"
    local subs ins del struct unaln
    subs="NA"; ins="NA"; del="NA"; struct="NA"; unaln="NA"
    if [[ -f "$stat" ]]; then
        subs=$(grep -iE 'Substitutions'          "$stat" | grep -oE '[0-9]+' | head -1)
        ins=$(grep -iE  'Insertions'             "$stat" | grep -oE '[0-9]+' | head -1)
        del=$(grep -iE  'Deletions'              "$stat" | grep -oE '[0-9]+' | head -1)
        # Structural events = relocations + translocations + inversions + duplications
        struct=$(grep -iE 'Relocations|Translocations|Inversions|duplications' "$stat" \
                  | grep -oE '[0-9]+' | awk '{s+=$1} END{print s+0}')
        unaln=$(grep -iE 'Unaligned' "$stat" | grep -oE '[0-9]+' | head -1)
        : "${subs:=0}" "${ins:=0}" "${del:=0}" "${struct:=0}" "${unaln:=NA}"
    fi
    log "    → SNP=${subs}  ins=${ins}  del=${del}  structural=${struct}"

    # Per-comparison summary
    {
        printf "metric\tvalue\n"
        printf "comparison\t%s\n" "$kind"
        printf "query\t%s\n" "$LABEL"
        printf "reference\t%s\n" "$ref_label"
        printf "substitutions\t%s\n" "$subs"
        printf "insertions\t%s\n" "$ins"
        printf "deletions\t%s\n" "$del"
        printf "structural_events\t%s\n" "$struct"
        printf "query_unaligned_bp\t%s\n" "$unaln"
    } > "${cdir}/summary.tsv"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$kind" "$ref_fasta" "$ref_label" "$subs" "$ins" "$del" "$struct" "$unaln" >> "$INDEX"

    # ── Dotplot (coordinate alignment) — non-fatal ────────────────────────────
    if [[ "${COMPARE_DOTPLOT:-true}" == "true" ]] && command -v mummerplot &>/dev/null; then
        set +e
        ( cd "$cdir" \
          && nucmer --maxmatch -p dot "$ref_fasta" "$QUERY" 2>dot_nucmer.log \
          && mummerplot --png --large -p dot dot.delta 2>dot_mummerplot.log )
        set -e
        [[ -f "${cdir}/dot.png" ]] && log "    dotplot: ${cdir}/dot.png" \
            || warn "    dotplot not produced (gnuplot missing?) — see ${cdir}/dot_*.log"
    fi
}

# =============================================================================
# 1. External comparison
# =============================================================================
if [[ "${EXT_REF_RAW}" != "-" && -n "${EXT_REF_RAW}" ]]; then
    EXT_PATH=$(resolve_ext_ref "${EXT_REF_RAW}")
    if [[ -z "${EXT_PATH}" || ! -f "${EXT_PATH}" ]]; then
        warn "  External reference '${EXT_REF_RAW}' → '${EXT_PATH:-<empty>}' not found — skipping external comparison."
    else
        EXT_NAME=$(basename "${EXT_PATH}"); EXT_NAME="${EXT_NAME%.*}"
        EXT_FASTA=$(ensure_fasta "${EXT_PATH}" "ext")
        # Gene-level annotation of variants if the reference is GenBank-derived
        EXT_GFF=""
        case "${EXT_PATH}" in *.gbk|*.gbff|*.gb)
            if command -v any2fasta &>/dev/null; then : ; fi ;; esac
        if [[ -n "${EXT_FASTA}" && -s "${EXT_FASTA}" ]]; then
            run_compare "${EXT_FASTA}" "${EXT_NAME}" "external" "${OUT}/vs_external_${EXT_NAME}" "${EXT_GFF}"
        else
            warn "  Could not obtain a FASTA for external reference ${EXT_PATH} — skipping."
        fi
    fi
fi

# =============================================================================
# 2. Internal comparison (vs a parent barcode assembled in this run)
# =============================================================================
if [[ "${INT_REF_BC}" != "-" && -n "${INT_REF_BC}" ]]; then
    REF_NAME=$(get_sample_field "${INT_REF_BC}" "sample_name" 2>/dev/null || echo "")
    if [[ -z "${REF_NAME}" ]]; then
        warn "  cmp_int_bc='${INT_REF_BC}' not found in sample sheet — skipping internal comparison."
    else
        REF_LABEL="${INT_REF_BC}_${REF_NAME}"
        REF_FASTA="${RESULTS_DIR}/03_polished/${REF_LABEL}/polished.fasta"
        # Prefer the parent's Pharokka GFF for gene-level variant annotation
        REF_GFF="${RESULTS_DIR}/05_annotated/${REF_LABEL}/${REF_NAME}.gff"
        [[ -f "${REF_GFF}" ]] || REF_GFF=""
        if [[ ! -s "${REF_FASTA}" ]]; then
            warn "  Internal reference assembly missing (${REF_FASTA}) — skipping internal comparison."
        else
            run_compare "${REF_FASTA}" "${REF_LABEL}" "internal" "${OUT}/vs_internal_${REF_LABEL}" "${REF_GFF}"
        fi
    fi
fi

log "  Comparison index: ${INDEX}"
log "  Step 08 done."
