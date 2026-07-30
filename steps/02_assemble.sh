#!/usr/bin/env bash
# =============================================================================
# Step 02 — Assembly: Autocycler (default) or Flye  [SLURM array]
# =============================================================================
# ASSEMBLER=autocycler (config.sh, recommended): multi-assembler consensus
#   (Ryan Wick's successor to Trycycler) — the best assemblies for accurate
#   terminus detection (PhageTerm, step 03). Autocycler does its OWN read
#   subsampling, so very high ONT coverage is handled internally.
#
# ASSEMBLER=flye  (or automatic fallback if Autocycler / its sub-assemblers are
#   unavailable): a single Flye --nano-hq assembly. The input is FIRST capped to
#   FILTLONG_TARGET_COV with Filtlong, because very high coverage (5 000–13 000×
#   on these runs) degrades Flye (heavier graph, over-fragmentation). ~100× is
#   the recommended cap.
#
# Illumina reads are NEVER used here (structure is ONT-defined); they only refine
# bases during short-read polishing in step 03. A single contig is a REPORTED QC
# metric, not a target — if a genome stays fragmented, the assembly graph is
# exposed as for_bandage.gfa (see MANUAL_FINISHING.md).
#
# Output per sample: $RESULTS_DIR/02_assembly/<barcode>_<sample>/
#   assembly.fasta | assembly_prefilter.fasta | assembler_used.txt
#   finishing_report.tsv | for_bandage.gfa (if fragmented) | assembly.log
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/config.sh"

activate_env
load_samtools_cluster
require_tool minimap2

SAMPLE_SHEET="${RESULTS_DIR}/sample_sheet.tsv"
mapfile -t BARCODES < <(list_barcodes)
if [[ "${ARRAY_IDX}" -ge "${#BARCODES[@]}" ]]; then
    echo "[SKIP] ARRAY_IDX=${ARRAY_IDX} >= barcodes in sample_sheet (${#BARCODES[@]}) — rien a faire."
    exit 0
fi
BARCODE="${BARCODES[${ARRAY_IDX}]}"
SAMPLE_NAME=$(get_sample_field "$BARCODE" "sample_name")
LABEL="${BARCODE}_${SAMPLE_NAME}"

step_banner "02 — Assembly  [${BARCODE} / ${SAMPLE_NAME}]"

ONT_READS="${RESULTS_DIR}/01_qc/${LABEL}/filtered.fastq.gz"
OUT="${RESULTS_DIR}/02_assembly/${LABEL}"
ASSEMBLY_FASTA="${OUT}/assembly.fasta"
mkdir -p "$OUT"

# ── Input guard ──────────────────────────────────────────────────────────────
if [[ ! -s "$ONT_READS" ]]; then
    echo "[SKIP] ${LABEL}: filtered reads missing (${ONT_READS}) — step 01 did not complete. Nothing to do."
    exit 0
fi

# Prefer host-depleted reads (step 01b) when HOST_DEPLETION is enabled and present
if [[ "${HOST_DEPLETION:-false}" == "true" ]]; then
    DEPLETED="${RESULTS_DIR}/01b_hostdepleted/${LABEL}/filtered.host_depleted.fastq.gz"
    if [[ -s "$DEPLETED" ]]; then
        ONT_READS="$DEPLETED"
        log "  Using host-depleted reads (step 01b)."
    else
        warn "  HOST_DEPLETION=true but no depleted reads for ${LABEL} — using full QC reads."
    fi
fi
if [[ -f "${ASSEMBLY_FASTA}" && "${FORCE:-false}" != "true" ]]; then
    log "  Assembly exists — skipping (set FORCE=true to rerun)"
    exit 0
fi

# Genome size string ("100k"/"1.5m") → integer bp
to_bp() {
    local v="${1,,}"
    case "$v" in
        *k) awk -v x="${v%k}" 'BEGIN{printf "%d", x*1000}' ;;
        *m) awk -v x="${v%m}" 'BEGIN{printf "%d", x*1000000}' ;;
        *)  awk -v x="$v"     'BEGIN{printf "%d", x}' ;;
    esac
}
GS_NUM=$(to_bp "${GENOME_SIZE_ESTIMATE}")
GFA_SRC=""
ASM_EFF="${ASSEMBLER:-autocycler}"

# ---------------------------------------------------------------------------
# Load Autocycler + sub-assembler tools (env preferred, else Lmod). Non-fatal.
# ---------------------------------------------------------------------------
if [[ "$ASM_EFF" == "autocycler" ]]; then
    if [[ -n "${AUTOCYCLER_ENV:-}" ]]; then
        log "  [env] activating micromamba env '${AUTOCYCLER_ENV}' for assembly tools"
        micromamba activate "${AUTOCYCLER_ENV}" 2>/dev/null \
            || warn "  [env] could not activate '${AUTOCYCLER_ENV}'."
    elif command -v module &>/dev/null; then
        log "  [modules] module load ${AUTOCYCLER_MODULES}"
        # ORDER MATTERS: assemblers first, Autocycler last (modulefile prerequisite).
        module load ${AUTOCYCLER_MODULES} 2>/dev/null \
            || warn "  [modules] 'module load' returned an error — check module names/order."
    else
        warn "  [autocycler] no AUTOCYCLER_ENV and no Lmod — autocycler must already be in PATH."
    fi
    if ! command -v autocycler &>/dev/null; then
        warn "  [assembler] 'autocycler' not found — falling back to Flye + Filtlong."
        ASM_EFF="flye"
    fi
fi

# Resolve which Autocycler sub-assemblers are actually installed
autocycler_binary_for() {
    case "$1" in
        flye) echo flye;; raven) echo raven;; miniasm) echo miniasm;;
        plassembler) echo plassembler;; canu) echo canu;; necat) echo necat;;
        nextdenovo) echo nextDenovo;; metamdbg) echo metaMDBG;; *) echo "$1";;
    esac
}
ACTIVE_ASSEMBLERS=""
if [[ "$ASM_EFF" == "autocycler" ]]; then
    for a in ${AUTOCYCLER_ASSEMBLERS:-flye}; do
        command -v "$(autocycler_binary_for "$a")" &>/dev/null && ACTIVE_ASSEMBLERS+="${a} "
    done
    ACTIVE_ASSEMBLERS="$(echo "$ACTIVE_ASSEMBLERS" | xargs || true)"
    if [[ -z "$ACTIVE_ASSEMBLERS" ]]; then
        warn "  [autocycler] none of (${AUTOCYCLER_ASSEMBLERS}) installed — falling back to Flye."
        ASM_EFF="flye"
    else
        log "  [autocycler] sub-assemblers available: ${ACTIVE_ASSEMBLERS}"
    fi
fi

# ---------------------------------------------------------------------------
# run_autocycler <reads> <outdir> <gs> -> echoes consensus fasta path; rc!=0 on fail
# ---------------------------------------------------------------------------
run_autocycler() {
    local reads="$1" outdir="$2" gs="$3"
    rm -rf "$outdir"; mkdir -p "$outdir/assemblies"
    autocycler subsample --reads "$reads" --out_dir "$outdir/subsampled" \
        --genome_size "$gs" >>"$outdir/autocycler.log" 2>&1 || return 1
    local got=0 f i a
    for f in "$outdir"/subsampled/sample_*.fastq; do
        [[ -e "$f" ]] || continue
        i="$(basename "$f" .fastq)"
        for a in $ACTIVE_ASSEMBLERS; do
            _hrc=0
            timeout "${AUTOCYCLER_HELPER_TIMEOUT:-2700}" autocycler helper "$a" --reads "$f" \
                   --out_prefix "$outdir/assemblies/${a}_${i}" \
                   --threads "${THREADS}" --genome_size "$gs" \
                   >>"$outdir/autocycler.log" 2>&1 || _hrc=$?
            if [[ ${_hrc} -eq 0 ]]; then
                got=$((got+1))
            elif [[ ${_hrc} -eq 124 ]]; then
                warn "    autocycler helper ${a} on ${i} TIMED OUT (${AUTOCYCLER_HELPER_TIMEOUT:-2700}s) — skipped."
            else
                warn "    autocycler helper ${a} failed on ${i} (rc=${_hrc}, continuing)."
            fi
        done
    done
    [[ "$got" -gt 0 ]] || return 1
    autocycler compress -i "$outdir/assemblies" -a "$outdir/autocycler_out" >>"$outdir/autocycler.log" 2>&1 || return 1
    autocycler cluster  -a "$outdir/autocycler_out" >>"$outdir/autocycler.log" 2>&1 || return 1
    local c had=0
    for c in "$outdir"/autocycler_out/clustering/qc_pass/cluster_*; do
        [[ -d "$c" ]] || continue
        had=1
        autocycler trim    -c "$c" >>"$outdir/autocycler.log" 2>&1 || true
        autocycler resolve -c "$c" >>"$outdir/autocycler.log" 2>&1 || true
    done
    [[ "$had" -eq 1 ]] || return 1
    autocycler combine -a "$outdir/autocycler_out" \
        -i "$outdir"/autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa \
        >>"$outdir/autocycler.log" 2>&1 || return 1
    [[ -s "$outdir/autocycler_out/consensus_assembly.fasta" ]] || return 1
    echo "$outdir/autocycler_out/consensus_assembly.fasta"
}

# ---------------------------------------------------------------------------
# Assemble
# ---------------------------------------------------------------------------
if [[ "$ASM_EFF" == "autocycler" ]]; then
    log "  Running Autocycler consensus (${ACTIVE_ASSEMBLERS})..."
    if CONS=$(run_autocycler "$ONT_READS" "${OUT}/autocycler" "$GS_NUM"); then
        cp "$CONS" "$ASSEMBLY_FASTA"
        GFA_SRC="${OUT}/autocycler/autocycler_out/consensus_assembly.gfa"
        echo "autocycler (${ACTIVE_ASSEMBLERS})" > "${OUT}/assembler_used.txt"
        log "  Autocycler consensus: $(basename "$CONS")"
    else
        warn "  Autocycler failed (see ${OUT}/autocycler/autocycler.log) — falling back to Flye."
        ASM_EFF="flye"
    fi
fi

if [[ "$ASM_EFF" == "flye" ]]; then
    require_tool flye

    # Cap input coverage with Filtlong (recommended for Flye on very high coverage)
    FLYE_INPUT="$ONT_READS"
    if [[ "${FILTLONG_TARGET_COV:-0}" -gt 0 ]] && command -v filtlong &>/dev/null; then
        _tgt=$(( FILTLONG_TARGET_COV * GS_NUM ))
        log "  Filtlong: capping input to ~${FILTLONG_TARGET_COV}× (target_bases=${_tgt})..."
        if filtlong --target_bases "$_tgt" "$ONT_READS" 2>"${OUT}/filtlong.log" | gzip > "${OUT}/reads_capped.fastq.gz" \
             && [[ -s "${OUT}/reads_capped.fastq.gz" ]]; then
            FLYE_INPUT="${OUT}/reads_capped.fastq.gz"
        else
            warn "  Filtlong failed — using full read set for Flye."
        fi
    elif [[ "${FILTLONG_TARGET_COV:-0}" -gt 0 ]]; then
        warn "  filtlong not found — using full read set (high coverage may degrade Flye)."
    fi

    FLYE_COV_ARG=()
    if [[ "${FLYE_ASM_COVERAGE:-0}" -gt 0 ]]; then
        FLYE_COV_ARG=(--asm-coverage "${FLYE_ASM_COVERAGE}")
    fi
    FLYE_META_ARG=()
    for _meta_bc in ${FLYE_META_BARCODES:-}; do
        if [[ "$_meta_bc" == "$BARCODE" ]]; then
            FLYE_META_ARG=(--meta)
            log "  Flye --meta enabled for ${BARCODE} (DTR / repeat-graph workaround)"
            [[ "${#FLYE_COV_ARG[@]}" -gt 0 ]] && { FLYE_COV_ARG=(); log "  --asm-coverage disabled (incompatible with --meta)"; }
            break
        fi
    done

    log "  Running Flye (--nano-hq, genome ~${GENOME_SIZE_ESTIMATE})..."
    flye --nano-hq "$FLYE_INPUT" --genome-size "${GENOME_SIZE_ESTIMATE}" \
        "${FLYE_COV_ARG[@]}" "${FLYE_META_ARG[@]}" \
        --out-dir "${OUT}/flye" --threads "${THREADS}" \
        2>&1 | tee "${OUT}/assembly.log"
    [[ -f "${OUT}/flye/assembly.fasta" ]] \
        || die "Flye failed — assembly.fasta not found. See ${OUT}/assembly.log"
    cp "${OUT}/flye/assembly.fasta" "$ASSEMBLY_FASTA"
    cp "${OUT}/flye/assembly_info.txt" "${OUT}/assembly_info.txt" 2>/dev/null || true
    GFA_SRC="${OUT}/flye/assembly_graph.gfa"
    echo "flye${FILTLONG_TARGET_COV:+ (filtlong ~${FILTLONG_TARGET_COV}x)}" > "${OUT}/assembler_used.txt"
fi

n_contigs=$(grep -c '^>' "$ASSEMBLY_FASTA" || echo 0)
total_len=$(grep -v '^>' "$ASSEMBLY_FASTA" | tr -d '\n' | wc -c)
log "  Assembly done: ${n_contigs} contig(s), total ${total_len} bp (assembler: $(cat "${OUT}/assembler_used.txt"))"
[[ "$n_contigs" -gt 5 ]] && warn "  ${n_contigs} contigs — check for contamination / mixed infection."

# Preserve raw assembly before parasite filter / finishing
cp "$ASSEMBLY_FASTA" "${OUT}/assembly_prefilter.fasta"

# ---------------------------------------------------------------------------
# Finishing: DTR / circular-overhang detection + parasite filter (never merges)
# ---------------------------------------------------------------------------
if command -v python3 &>/dev/null; then
    FIN_REPORT="${OUT}/finishing_report.tsv"
    log "  Finishing check (DTR / circular-overhang detection)..."
    FIN_TRIM_ARG=()
    [[ "${FINISH_TRIM_OVERHANG:-false}" == "true" ]] && FIN_TRIM_ARG=(--trim)
    if python3 "${SCRIPT_DIR}/finish_assembly.py" \
        --in "$ASSEMBLY_FASTA" --out "${ASSEMBLY_FASTA}.fin" --report "$FIN_REPORT" \
        --min-contig-len "${MIN_CONTIG_LEN:-0}" --circ-max "${FINISH_CIRC_MAX:-1000}" \
        "${FIN_TRIM_ARG[@]}" 2>&1 | tee -a "${OUT}/assembly.log"; then
        [[ -s "${ASSEMBLY_FASTA}.fin" ]] && { mv "${ASSEMBLY_FASTA}.fin" "$ASSEMBLY_FASTA"; n_contigs=$(grep -c '^>' "$ASSEMBLY_FASTA" || echo 0); }
    else
        warn "  Finishing helper failed — keeping assembly unchanged."
        rm -f "${ASSEMBLY_FASTA}.fin"
    fi
fi

# Expose the assembly graph for manual Bandage finishing if still fragmented
if [[ "${n_contigs}" -gt 1 ]]; then
    if [[ -n "$GFA_SRC" && -f "$GFA_SRC" ]]; then
        cp "$GFA_SRC" "${OUT}/for_bandage.gfa"
        warn "  ${n_contigs} contigs remain — open ${OUT}/for_bandage.gfa in Bandage (see MANUAL_FINISHING.md);"
        warn "  after manual circularisation, replace assembly.fasta and rerun from step 03 with FORCE=true."
    else
        warn "  ${n_contigs} contigs remain and no assembly graph (GFA) found for Bandage."
    fi
fi

# Reorientation happens in step 03 (PhageTerm/dnaapler on polished.fasta).
log "  Step 02 done."
