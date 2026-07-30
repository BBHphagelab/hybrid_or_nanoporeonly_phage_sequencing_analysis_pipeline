#!/usr/bin/env bash
# =============================================================================
# rebasecall_dorado.sh — Re-basecall POD5 with BROADER modified-base models
# =============================================================================
# Why: MinKNOW called 6mA + 5mCG_5hmCG for this run — 5mC only in CpG context and
# NO 4mC. To VERIFY the three canonical prokaryotic methylations in ALL contexts
# (6mA, 4mC, 5mC), re-basecall the raw POD5 with the SUP model + all-context mods.
#
# Pipeline (2 SLURM jobs):
#   Job A (GPU) : dorado basecaller  <model> --modified-bases ${DORADO_MOD_BASES}
#                 --kit-name ${BARCODE_KIT}   → one unaligned mod-tagged BAM
#   Job B (CPU) : dorado demux --no-classify → split on the BC tags Job A wrote.
#                 (--no-classify and --kit-name are mutually exclusive; passing
#                 the kit again re-classifies from scratch on trimmed reads.
#                 Omitting BOTH → dorado exits 1 immediately.)
#                 → per-barcode BAMs in ${REBC_DEMUX_DIR}
#
# It does NOT touch the MinKNOW bam_pass: outputs go to REBC_DEMUX_DIR. After it
# finishes you point MOD_DEMUX_DIR at REBC_DEMUX_DIR and rerun steps 06→07.
#
# Usage:
#   bash rebasecall_dorado.sh --check       # validate paths/tools/GPU, submit NOTHING
#   bash rebasecall_dorado.sh --dry-run     # print the sbatch commands, submit nothing
#   bash rebasecall_dorado.sh --demux-only  # skip Job A, resubmit demux on existing BAM
#   bash rebasecall_dorado.sh               # submit Job A + Job B
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

MODE="submit"
DEMUX_ONLY=false
for a in "$@"; do
    case "$a" in
        --check)      MODE="check" ;;
        --dry-run)    MODE="dryrun" ;;
        --demux-only) MODE="submit"; DEMUX_ONLY=true ;;  # skip Job A, resubmit demux only
        -h|--help) sed -n '3,30p' "$0"; exit 0 ;;
        *) die "Unknown option: $a (use --check, --dry-run, --demux-only, or no flag)" ;;
    esac
done

MAESTRO_TMPDIR="/local/scratch/tmp"
DORADO_BIN="${DORADO_BIN_DIR}/bin/dorado"
MOD_CALLS_BAM="${REBC_DEMUX_DIR}/../rebasecall_mod_calls.bam"
mkdir -p "${RESULTS_DIR}/logs"

echo "════════════════════════════════════════════════════════"
echo "  Re-basecalling (Dorado) — methylation verification"
echo "════════════════════════════════════════════════════════"

# ---------------------------------------------------------------------------
# Resolve POD5 directory (config override, else auto-detect under RUN_DIR)
# ---------------------------------------------------------------------------
_pod5_dir="${POD5_DIR:-}"
if [[ -z "${_pod5_dir}" ]]; then
    for cand in "${RUN_DIR}/pod5_pass" "${RUN_DIR}/pod5" "${RUN_DIR}/fast5_pass"; do
        [[ -d "${cand}" ]] && { _pod5_dir="${cand}"; break; }
    done
fi

PREFLIGHT_OK=true
_ok()   { echo "  [OK]   $*"; }
_warn() { echo "  [WARN] $*"; }
_fail() { echo "  [FAIL] $*"; PREFLIGHT_OK=false; }

echo ""
echo "  Configuration:"
echo "    RUN_DIR          : ${RUN_DIR}"
echo "    POD5_DIR         : ${_pod5_dir:-'(not found)'}"
echo "    Simplex model    : ${DORADO_SIMPLEX_MODEL}"
echo "    Modified bases   : ${DORADO_MOD_BASES}"
echo "    Barcode kit      : ${BARCODE_KIT}"
echo "    Output (demux)   : ${REBC_DEMUX_DIR}"
echo "    GPU partition    : ${GPU_PARTITION}"
echo ""

# ── Checks ──────────────────────────────────────────────────────────────────
# Job A (GPU basecalling) needs POD5. In --demux-only we skip Job A, so POD5 is
# irrelevant; instead the existing mod-calls BAM from a prior Job A must be present.
if [[ "${DEMUX_ONLY}" == true ]]; then
    if [[ -s "${MOD_CALLS_BAM}" ]]; then
        _ok "mod-calls BAM present ($(ls -lh "${MOD_CALLS_BAM}" | awk '{print $5}')): ${MOD_CALLS_BAM}"
    else
        _fail "mod-calls BAM missing/empty: ${MOD_CALLS_BAM} — run Job A first (no --demux-only)."
    fi
elif [[ -n "${_pod5_dir}" && -d "${_pod5_dir}" ]]; then
    n_pod5=$(find "${_pod5_dir}" -name '*.pod5' 2>/dev/null | wc -l | tr -d ' ')
    n_fast5=$(find "${_pod5_dir}" -name '*.fast5' 2>/dev/null | wc -l | tr -d ' ')
    if [[ "${n_pod5}" -gt 0 ]]; then
        _ok "POD5 dir has ${n_pod5} .pod5 file(s)"
    elif [[ "${n_fast5}" -gt 0 ]]; then
        _warn "No .pod5 but ${n_fast5} .fast5 found — dorado reads .pod5; convert with 'pod5 convert fast5' first."
    else
        _fail "No .pod5/.fast5 under ${_pod5_dir}"
    fi
else
    _fail "POD5 directory not found. Set POD5_DIR in config.sh (looked for pod5_pass/pod5/fast5_pass under RUN_DIR)."
fi

if [[ -x "${DORADO_BIN}" ]]; then
    _ok "dorado binary: ${DORADO_BIN}  ($("${DORADO_BIN}" --version 2>&1 | head -1))"
else
    _warn "dorado not found at ${DORADO_BIN} — Job A will attempt to download it (needs internet on the GPU node)."
fi

# BARCODE_KIT sanity: RBK (rapid) vs NBD (native) — a wrong kit → all-unclassified demux.
# This run = native barcoding (confirmed): SQK-NBD114-24.
case "${BARCODE_KIT}" in
    SQK-NBD114-24) _ok "Barcode kit: ${BARCODE_KIT} (native barcoding)" ;;
    SQK-RBK114-24) _warn "Barcode kit is ${BARCODE_KIT} (rapid) — but this run is NATIVE (SQK-NBD114-24). Fix BARCODE_KIT in config.sh." ;;
    *) _warn "Barcode kit '${BARCODE_KIT}' unexpected — this run should be SQK-NBD114-24 (native)." ;;
esac

if [[ "${DEMUX_ONLY}" != true ]] && command -v sinfo &>/dev/null; then
    if sinfo -p "${GPU_PARTITION}" -h &>/dev/null && [[ -n "$(sinfo -p "${GPU_PARTITION}" -h 2>/dev/null)" ]]; then
        _ok "GPU partition '${GPU_PARTITION}' exists"
    else
        _warn "GPU partition '${GPU_PARTITION}' not visible via sinfo — set GPU_PARTITION in config.sh (try: sinfo -o '%P')."
    fi
fi

echo ""
if [[ "${PREFLIGHT_OK}" == "true" ]]; then
    echo "  ✅  Preflight OK."
else
    echo "  ❌  Preflight FAILED — fix the [FAIL] items above."
    [[ "${MODE}" == "submit" ]] && exit 1
fi

if [[ "${MODE}" == "check" ]]; then
    echo ""
    echo "  (--check) Nothing submitted. Commands that WOULD run:"
    echo "    pre-dl: ${DORADO_BIN} download --model ${DORADO_MODEL_COMPLEX} --data ${_pod5_dir} --models-directory ${DORADO_MODELS_DIR}   (login node)"
    echo "    Job A: ${DORADO_BIN} basecaller ${DORADO_MODEL_COMPLEX} ${_pod5_dir} \\"
    echo "             --kit-name ${BARCODE_KIT} --recursive --device cuda:all \\"
    echo "             --models-directory ${DORADO_MODELS_DIR} > ${MOD_CALLS_BAM}"
    echo "    Job B: ${DORADO_BIN} demux --output-dir ${REBC_DEMUX_DIR} --no-classify --emit-summary ${MOD_CALLS_BAM}"
    exit 0
fi

# ---------------------------------------------------------------------------
# Pre-download the model on the LOGIN node (has internet). Maestro GPU nodes may
# be offline; caching the model complex here means Job A finds it locally in
# --models-directory and needs no network. Best-effort — Job A also auto-resolves.
# ---------------------------------------------------------------------------
if [[ "${MODE}" == "submit" && "${DEMUX_ONLY}" != true ]]; then
    if [[ -x "${DORADO_BIN}" ]]; then
        echo ""
        echo "  Pre-downloading model complex on login node → ${DORADO_MODELS_DIR}"
        mkdir -p "${DORADO_MODELS_DIR}"
        # A model COMPLEX triggers dorado's automatic model detection, which needs
        # --data <pod5> to inspect the chemistry (else: "Must set --data when using
        # automatic model detection"). Pointing at the POD5 dir resolves it.
        "${DORADO_BIN}" download --model "${DORADO_MODEL_COMPLEX}" --data "${_pod5_dir}" --models-directory "${DORADO_MODELS_DIR}" \
            && _ok "model(s) cached locally" \
            || _warn "pre-download failed — Job A will try to fetch on the GPU node (needs internet there)."
    else
        _warn "dorado not on login node — skipping pre-download; Job A must fetch on the GPU node."
    fi
fi

# ---------------------------------------------------------------------------
# sbatch (with dry-run override)
# ---------------------------------------------------------------------------
if [[ "${MODE}" == "dryrun" ]]; then
    sbatch() { echo "  [DRY-RUN] sbatch $*" | fold -s -w 120 >&2; echo "DRY"; }
fi

gpu_flags() {
    echo -n "--parsable --partition=${GPU_PARTITION}"
    [[ -n "${ACCOUNT}" ]] && echo -n " --account=${ACCOUNT}"
    [[ -n "${EMAIL}"   ]] && echo -n " --mail-user=${EMAIL} --mail-type=begin,end,fail"
}
cpu_flags() {
    echo -n "--parsable"
    [[ -n "${PARTITION}" ]] && echo -n " --partition=${PARTITION}"
    [[ -n "${ACCOUNT}"   ]] && echo -n " --account=${ACCOUNT}"
    [[ -n "${EMAIL}"     ]] && echo -n " --mail-user=${EMAIL} --mail-type=begin,end,fail"
}

mkdir -p "${REBC_DEMUX_DIR}"

# ── Job A — GPU basecaller with all-context modified bases ───────────────────
# Skipped in --demux-only (Job A already produced MOD_CALLS_BAM); Job B then runs
# standalone with no afterok dependency.
DEP=""
if [[ "${DEMUX_ONLY}" != true ]]; then
JOB_BC=$(sbatch \
    $(gpu_flags) \
    --job-name="phage_rebc_dorado" \
    --time="16:00:00" \
    --ntasks=1 --cpus-per-task=16 --mem=64G \
    --gres=gpu:1 \
    --output="${RESULTS_DIR}/logs/rebc_dorado.%j.out" \
    --wrap="set -euo pipefail
DORADO='${DORADO_BIN}'
if [[ ! -x \"\${DORADO}\" ]]; then
    echo '[INFO] dorado missing — downloading to ${DORADO_BIN_DIR}'
    mkdir -p '${DORADO_BIN_DIR}'
    curl -fsSL -L 'https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.1-linux-x64.tar.gz' \
        | tar -xz -C '${DORADO_BIN_DIR}' --strip-components=1
fi
echo '[INFO] dorado: '\$(\"\${DORADO}\" --version 2>&1 | head -1)
mkdir -p '${DORADO_MODELS_DIR}'
n=\$(find '${_pod5_dir}' -name '*.pod5' | wc -l); echo \"[INFO] \${n} POD5 files\"
# Model COMPLEX (speed@version,mods) → dorado auto-resolves/downloads. Passing the
# full model NAME instead makes dorado treat it as a path and fail.
\"\${DORADO}\" basecaller '${DORADO_MODEL_COMPLEX}' '${_pod5_dir}' \
    --kit-name '${BARCODE_KIT}' \
    --recursive \
    --device cuda:all \
    --models-directory '${DORADO_MODELS_DIR}' \
    > '${MOD_CALLS_BAM}'
echo '[INFO] basecalling done: '\$(ls -lh '${MOD_CALLS_BAM}')")
[[ "${MODE}" != "dryrun" ]] && echo "  Job A (basecaller): ${JOB_BC}"
DEP="--dependency=afterok:${JOB_BC}"
else
    echo "  (--demux-only) Job A skipped — reusing ${MOD_CALLS_BAM}"
fi

# ── Job B — demux per barcode ────────────────────────────────────────────────
# --no-classify (NOT --kit-name) splits on the BC tags Job A already wrote. dorado
# demux requires exactly one of --no-classify / --kit-name; omitting both → it
# searches for barcodes with no kit and exits 1 immediately (the 59031899 failure).
JOB_DX=$(sbatch \
    $(cpu_flags) \
    --job-name="phage_rebc_demux" \
    --time="03:00:00" \
    ${DEP} \
    --ntasks=1 --cpus-per-task=8 --mem=16G \
    --output="${RESULTS_DIR}/logs/rebc_demux.%j.out" \
    --wrap="set -euo pipefail
DORADO='${DORADO_BIN}'
mkdir -p '${REBC_DEMUX_DIR}'
\"\${DORADO}\" demux --output-dir '${REBC_DEMUX_DIR}' --no-classify --emit-summary '${MOD_CALLS_BAM}'
echo '[INFO] demux done. Per-barcode BAMs:'
ls -lh '${REBC_DEMUX_DIR}'/*.bam 2>/dev/null | head
echo ''
echo '================ NEXT STEPS ================'
echo '1. Confirm one BAM per barcode above (names contain barcodeNN).'
echo '2. In config.sh set:  MOD_DEMUX_DIR=\"${REBC_DEMUX_DIR}\"'
echo '3. Rerun modbase + methylation (FORCE re-does the cached step 06):'
echo '     FORCE=true bash submit_all.sh --step 06'
echo '     bash submit_all.sh --step 07a   # then, once done:'
echo '     bash submit_all.sh --step 07b'
echo '==========================================='")
[[ "${MODE}" != "dryrun" ]] && echo "  Job B (demux):      ${JOB_DX}"

echo ""
echo "════════════════════════════════════════════════════════"
if [[ "${MODE}" == "dryrun" ]]; then
    echo "  DRY-RUN — nothing submitted."
elif [[ "${DEMUX_ONLY}" == true ]]; then
    echo "  Submitted:  B(demux only)=${JOB_DX}  (Job A skipped, reused existing BAM)"
else
    echo "  Submitted:  A(GPU basecall)=${JOB_BC} → B(demux)=${JOB_DX}"
fi
echo "  Monitor : squeue -u \$(whoami)"
echo "  Logs    : ${RESULTS_DIR}/logs/rebc_*.out"
echo "════════════════════════════════════════════════════════"
