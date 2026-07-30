#!/usr/bin/env bash
# =============================================================================
# submit_all.sh — Submit the full nanopore PHAGE pipeline via SLURM
# =============================================================================
# Run from the pipeline directory:
#   bash submit_all.sh                    # full run
#   bash submit_all.sh --dry-run          # print sbatch commands, submit nothing
#   bash submit_all.sh --check            # pre-flight checks only
#   bash submit_all.sh --step 04          # submit only step 04 (no dependency)
#   bash submit_all.sh --step 04 --resume # submit step 04 only for unfinished samples
#   bash submit_all.sh --resume           # full run, skip already-completed samples
#
# Job dependency chain:
#
#   [00_setup]         ← global (builds sample_sheet.tsv, merges FASTQs)
#        ↓  (array per sample)
#   [01_qc]            ← chopper + NanoStat per sample
#        ↓
#   [02_assemble]      ← Autocycler (default) or Flye per sample (ONT defines structure)
#        ↓
#   [03_polish]        ← Medaka [+ Polypolish/Pypolca if Illumina]
#       ↓               ↓               ↓
#  [04_assembly_qc]  [05_annotate]   [06_modbase*]   ← parallel after polishing
#                                         ↓
#                                    [07a_methyl*]
#                                         ↓
#                                    [07b_methyl_report*]   ← * only if MOD_DEMUX_DIR set
#
#  [09_summary]  ← 04+05+06           [10_final_report]  ← 09+07b
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# =============================================================================
# Parse arguments
# =============================================================================
DRY_RUN=false
CHECK_ONLY=false
RESUME=false
STEP_ONLY=""   # empty = run all steps; "04" = run only that step

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)  DRY_RUN=true;   shift ;;
        --check)    CHECK_ONLY=true; shift ;;
        --resume)   RESUME=true;     shift ;;
        --step)
            [[ -z "${2:-}" ]] && { echo "[ERROR] --step requires a step number (e.g. --step 04)" >&2; exit 1; }
            STEP_ONLY="$2"; shift 2 ;;
        -h|--help)
            sed -n '3,19p' "$0"
            exit 0 ;;
        *) echo "[ERROR] Unknown option: $1  (use --dry-run, --check, --resume, --step N, or --help)" >&2; exit 1 ;;
    esac
done

# =============================================================================
# Dry-run override
# =============================================================================
_DRY_COUNTER_FILE=""
if [[ "${DRY_RUN}" == "true" ]]; then
    _DRY_COUNTER_FILE=$(mktemp)
    echo 0 > "${_DRY_COUNTER_FILE}"
    trap 'rm -f "${_DRY_COUNTER_FILE}"' EXIT

    sbatch() {
        local n=$(( $(cat "${_DRY_COUNTER_FILE}") + 1 ))
        echo "${n}" > "${_DRY_COUNTER_FILE}"
        echo "  [DRY-RUN] sbatch $*" | fold -s -w 120 >&2
        echo "DRY${n}"
    }
    echo ""
    echo "════════════════════════════════════════════════════════"
    echo "  DRY-RUN MODE — no jobs will be submitted"
    echo "════════════════════════════════════════════════════════"
fi

# =============================================================================
# Pre-flight checks
# =============================================================================
echo ""
echo "════════════════════════════════════════════════════════"
echo "  Pre-flight checks"
echo "════════════════════════════════════════════════════════"

PREFLIGHT_OK=true
_ok()   { echo "  [OK]   $*"; }
_warn() { echo "  [WARN] $*"; }
_fail() { echo "  [FAIL] $*"; PREFLIGHT_OK=false; }

echo ""
echo "  Configuration:"
echo "    NANOPORE_BASE  : ${NANOPORE_BASE}"
echo "    RUN_NAME       : ${RUN_NAME}"
echo "    RUN_DIR        : ${RUN_DIR}"
echo "    RESULTS_DIR    : ${RESULTS_DIR}"
echo "    EMAIL          : ${EMAIL}"
echo "    PARTITION      : ${PARTITION:-'(default)'}"
echo "    THREADS        : ${THREADS}"
echo "    MAMBA_ENV      : ${MAMBA_ENV}"
echo "    GENOME_SIZE    : ${GENOME_SIZE_ESTIMATE}"
echo "    MOD_DEMUX_DIR  : ${MOD_DEMUX_DIR:-'(not set — step 07 skipped)'}"
echo ""

echo "  Input directories:"
[[ -d "${NANOPORE_BASE}" ]] \
    && _ok  "NANOPORE_BASE exists" \
    || _fail "NANOPORE_BASE not found: ${NANOPORE_BASE}"

[[ -d "${RUN_DIR}" ]] \
    && _ok  "RUN_DIR exists" \
    || _fail "RUN_DIR not found: ${RUN_DIR}  ← is RUN_NAME correct in config.sh?"

[[ -d "${FASTQ_PASS_DIR}" ]] \
    && _ok  "FASTQ_PASS_DIR exists  ($(find "${FASTQ_PASS_DIR}" -maxdepth 1 -type d -name 'barcode*' 2>/dev/null | wc -l) barcode folders)" \
    || _fail "FASTQ_PASS_DIR not found: ${FASTQ_PASS_DIR}"

[[ -d "${ILLUMINA_DIR:-}" ]] \
    && _ok  "ILLUMINA_DIR exists  ($(find "${ILLUMINA_DIR}" -maxdepth 1 -type d -name 'barcode*' 2>/dev/null | wc -l) barcode folders)" \
    || _warn "ILLUMINA_DIR not found: ${ILLUMINA_DIR:-unset}  (OK — ONT-only mode)"

echo ""
echo "  Micromamba environment:"
MAMBA_BIN="${MAMBA_ROOT_PREFIX}/bin/micromamba"
[[ -x "${MAMBA_BIN}" ]] \
    && _ok  "micromamba binary found: ${MAMBA_BIN}" \
    || _fail "micromamba not found — run: bash ${SCRIPT_DIR}/00_install_env.sh"

[[ -d "${MAMBA_ROOT_PREFIX}/envs/${MAMBA_ENV}" ]] \
    && _ok  "env '${MAMBA_ENV}' exists" \
    || _fail "env '${MAMBA_ENV}' not found — run: bash ${SCRIPT_DIR}/00_install_env.sh"

echo ""
echo "  Databases:"
[[ -d "${PHAROKKA_DB}" ]] \
    && _ok  "Pharokka DB found: ${PHAROKKA_DB}" \
    || _fail "Pharokka DB not found at ${PHAROKKA_DB}  ← run: pharokka_db_install.py --db_dir <path>"

echo ""
echo "  Step scripts:"
for s in 00_setup 01_qc 01b_hostdeplete 02_assemble 03_polish 04_assembly_qc 05_annotate 06_modbase 07a_methyl_call 07b_methyl_report \
         09_summary 10_final_report; do
    f="${STEPS_DIR}/${s}.sh"
    [[ -f "$f" ]] && _ok "$s.sh" || _fail "$s.sh MISSING from ${STEPS_DIR}/"
done

echo ""
echo "  QC tools (step 04):"
command -v checkv  &>/dev/null && _ok "checkv found" || _warn "checkv not in PATH — step 04 CheckV will be skipped"

if [[ "${HOST_DEPLETION:-false}" == "true" ]]; then
    echo ""
    echo "  Host depletion (step 01b — ENABLED):"
    if [[ -n "${HOST_REF:-}" && -f "${HOST_REF}" ]]; then
        _ok "HOST_REF found: ${HOST_REF}"
    else
        _warn "HOST_DEPLETION=true but HOST_REF missing ('${HOST_REF:-}') — 01b will pass reads through unchanged"
    fi
fi

if [[ -n "${MOD_DEMUX_DIR:-}" ]]; then
    echo ""
    echo "  Methylation (steps 06+07):"
    if [[ -d "${MOD_DEMUX_DIR}" ]]; then
        n_bams=$(find "${MOD_DEMUX_DIR}" -name "*.bam" -maxdepth 2 2>/dev/null | wc -l)
        _ok "MOD_DEMUX_DIR exists  (${n_bams} BAM file(s) found)"
    else
        _warn "MOD_DEMUX_DIR set but not found: ${MOD_DEMUX_DIR}  (will be needed before step 07 runs)"
    fi

    # modkit 0.2.x is REQUIRED by MicrobeMod 1.1.0 (step 07a): it uses --only-tabs,
    # removed in modkit >= 0.3. Without a pinned 0.2.x, every barcode fails rc=1.
    if [[ -n "${MODKIT_02X_BIN:-}" && -x "${MODKIT_02X_BIN}" ]]; then
        _mk_ver=$("${MODKIT_02X_BIN}" --version 2>&1 | head -1)
        case "${_mk_ver}" in
            *0.2.*) _ok  "modkit 0.2.x for MicrobeMod: ${_mk_ver}  (${MODKIT_02X_BIN})" ;;
            *)      _warn "MODKIT_02X_BIN is ${_mk_ver} (expected 0.2.x) — MicrobeMod (07a) will fail on --only-tabs" ;;
        esac
    else
        _warn "modkit 0.2.x not found (MODKIT_02X_BIN='${MODKIT_02X_BIN:-unset}') — step 07a MicrobeMod will FAIL (rc=1)."
        _warn "  Install once: micromamba create -n modkit_0.2 -c bioconda -c conda-forge ont-modkit=0.2.6"
    fi
fi

echo ""
if [[ "${PREFLIGHT_OK}" == "true" ]]; then
    echo "  ✅  All checks passed."
else
    echo "  ❌  One or more checks FAILED — fix the issues above before submitting."
    [[ "${DRY_RUN}" == "false" && "${CHECK_ONLY}" == "false" ]] && exit 1
fi
echo ""

[[ "${CHECK_ONLY}" == "true" ]] && { echo "  (--check mode: stopping here)"; exit 0; }

# =============================================================================
# Helpers
# =============================================================================
sbatch_flags() {
    local time="$1"
    echo -n "--parsable"
    echo -n " --time=${time}"
    [[ -n "${PARTITION}" ]] && echo -n " --partition=${PARTITION}"
    [[ -n "${ACCOUNT}"   ]] && echo -n " --account=${ACCOUNT}"
    [[ -n "${EMAIL}"     ]] && echo -n " --mail-user=${EMAIL} --mail-type=begin,end,fail"
}

MAESTRO_TMPDIR="/local/scratch/tmp"

# =============================================================================
# Array sizing
# =============================================================================
mkdir -p "${RESULTS_DIR}/logs"

SAMPLE_SHEET="${RESULTS_DIR}/sample_sheet.tsv"
if [[ -f "$SAMPLE_SHEET" ]]; then
    N_SAMPLES=$(count_barcodes)
    echo "  Array sizing: ${N_SAMPLES} barcodes from existing sample_sheet.tsv"
else
    N_SAMPLES=$(find "${FASTQ_PASS_DIR}" -maxdepth 1 -mindepth 1 -type d -name 'barcode*' 2>/dev/null | wc -l | tr -d ' ')
    [[ "${N_SAMPLES}" -eq 0 ]] && { echo "  WARNING: cannot detect sample count — defaulting to 1."; N_SAMPLES=1; }
    echo "  Array sizing: ${N_SAMPLES} barcodes detected from FASTQ_PASS_DIR"
fi
ARRAY_RANGE="0-$((N_SAMPLES - 1))"

# =============================================================================
# --resume helper: compute which array indices still need to run for a step.
#
# Sentinel template uses %BARCODE% and %SAMPLE% as placeholders.
#   e.g. "${RESULTS_DIR}/05_annotated/%BARCODE%_%SAMPLE%/%SAMPLE%.gff"
# Returns (stdout): comma-separated indices, "none" if all done, or full ARRAY_RANGE.
#
# BUG NOTES:
#   - Use tail -n +2 + process substitution to skip the header WITHOUT a subshell
#     (a pipeline `cmd | while` runs the while in a subshell; array assignments
#     inside would be lost). Process substitution keeps the loop in the parent shell.
#   - Placeholder syntax %BARCODE%/%SAMPLE% avoids ambiguity with __BARCODE____SAMPLE__
#     where the double-underscore boundary is unclear when the two tokens are adjacent.
# =============================================================================
resume_array() {
    local step="$1"
    local sentinel_tpl="$2"

    # Without --resume, always return the full range
    if [[ "${RESUME}" != "true" ]]; then
        echo "${ARRAY_RANGE}"
        return
    fi

    # Need sample_sheet to resolve barcode→sample_name
    if [[ ! -f "${SAMPLE_SHEET}" ]]; then
        echo "  [RESUME] sample_sheet.tsv not found — cannot check completion, submitting full array" >&2
        echo "${ARRAY_RANGE}"
        return
    fi

    local pending=()
    local idx=0
    # Use process substitution (< <(...)) so the while loop runs in the current
    # shell — array variable 'pending' is visible after the loop.
    # tail -n +2 skips the header row without an extra [[ ]] test inside the loop.
    while IFS=$'\t' read -r barcode sample_name _rest; do
        local sentinel="${sentinel_tpl//%BARCODE%/${barcode}}"
        sentinel="${sentinel//%SAMPLE%/${sample_name}}"
        if [[ -f "${sentinel}" ]]; then
            echo "  [RESUME] step ${step} idx ${idx} (${barcode}/${sample_name}) — already done, skipping" >&2
        else
            echo "  [RESUME] step ${step} idx ${idx} (${barcode}/${sample_name}) — pending" >&2
            pending+=("${idx}")
        fi
        idx=$(( idx + 1 ))
    done < <(tail -n +2 "${SAMPLE_SHEET}")

    if [[ ${#pending[@]} -eq 0 ]]; then
        echo "none"
    else
        local IFS=','
        echo "${pending[*]}"
    fi
}

# =============================================================================
# Wrapper: submit one array step, honouring --step and --resume.
#
# Stdout returns ONLY the SLURM job ID (or "skipped_step"/"skipped_resume").
# All human-readable messages go to stderr so they don't pollute the job-ID
# capture done by callers (JOB_XX=$(submit_array_step ...)).
# =============================================================================
submit_array_step() {
    local label="$1";    shift
    local step_num="$1"; shift   # e.g. "04"
    local sentinel="$1"; shift   # sentinel template for --resume (%BARCODE%, %SAMPLE%)
    # Remaining args are passed verbatim to sbatch

    # If --step is set and this isn't the requested step, skip silently
    if [[ -n "${STEP_ONLY}" && "${STEP_ONLY}" != "${step_num}" ]]; then
        echo "skipped_step"
        return
    fi

    local arr
    arr=$(resume_array "${step_num}" "${sentinel}")

    if [[ "${arr}" == "none" ]]; then
        echo "  [RESUME] step ${step_num} — all samples already done, skipping submission" >&2
        echo "skipped_resume"
        return
    fi

    # Progress message → stderr only (so stdout stays clean for job-ID capture)
    echo "Step ${step_num} — ${label}  (array ${arr})" >&2
    sbatch --array="${arr}" "$@"
    # sbatch --parsable prints only the job ID to stdout → caller captures it
}

# =============================================================================
# Job submission
# =============================================================================
echo ""
echo "════════════════════════════════════════════════════════"
[[ "${DRY_RUN}" == "true" ]] \
    && echo "  Jobs that WOULD be submitted:" \
    || echo "  Submitting jobs..."
if [[ -n "${STEP_ONLY}" ]]; then
    echo "  --step ${STEP_ONLY}: submitting only step ${STEP_ONLY} (no upstream dependency)"
fi
[[ "${RESUME}" == "true" ]] && echo "  --resume: completed samples will be skipped"
echo "════════════════════════════════════════════════════════"
echo ""

# When --step is used, we don't build a dependency chain — each step submitted
# standalone so it can run immediately.
_dep() {
    local job_id="$1"
    # Default afterany: a downstream step STARTS once upstream finishes in ANY
    # state. Combined with the per-barcode input guards at the top of each step
    # (clean [SKIP] when an input is missing), one failed/empty barcode never
    # blocks the others, and a non-zero-but-output-producing task still lets the
    # pipeline advance. Pass "afterok" explicitly where strict success is wanted.
    local dep_type="${2:-afterany}"
    # If --step is set, suppress dependencies entirely
    if [[ -n "${STEP_ONLY}" ]]; then
        return
    fi
    # If the upstream job was skipped (resume or step filter), no dependency needed
    if [[ "${job_id}" == "skipped_"* ]]; then
        return
    fi
    echo -n " --dependency=${dep_type}:${job_id}"
}

# Build a combined dependency from several upstream job IDs (all must reach a
# terminal state — afterany). Skipped/standalone jobs are ignored. Prints e.g.
#   --dependency=afterany:123:456
# or nothing if --step is set or no real upstream jobs remain.
_dep_multi() {
    [[ -n "${STEP_ONLY}" ]] && return
    local ids=()
    local j
    for j in "$@"; do
        [[ "${j}" == "skipped_"* || -z "${j}" ]] || ids+=("${j}")
    done
    [[ ${#ids[@]} -eq 0 ]] && return
    local IFS=':'
    echo -n " --dependency=afterany:${ids[*]}"
}

# ── Step 00 ───────────────────────────────────────────────────────────────────
JOB_SETUP="skipped_step"
if [[ -z "${STEP_ONLY}" || "${STEP_ONLY}" == "00" ]]; then
    # --resume: skip if sample_sheet already exists
    if [[ "${RESUME}" == "true" && -f "${SAMPLE_SHEET}" ]]; then
        echo "  [RESUME] step 00 — sample_sheet.tsv already exists, skipping setup"
        JOB_SETUP="skipped_resume"
    else
        echo "Step 00 — setup (global)"
        JOB_SETUP=$(sbatch \
            $(sbatch_flags "${TIME_SETUP}") \
            --job-name="phage_00_setup" \
            --ntasks=1 --cpus-per-task=2 --mem=4G \
            --output="${RESULTS_DIR}/logs/00_setup.%j.out" \
            --wrap="export TMPDIR=${MAESTRO_TMPDIR} && source ${SCRIPT_DIR}/config.sh && activate_env && bash ${STEPS_DIR}/00_setup.sh")
        echo "  Job ID: ${JOB_SETUP}"
    fi
fi

# ── Step 01 ───────────────────────────────────────────────────────────────────
SENTINEL_01="${RESULTS_DIR}/01_qc/%BARCODE%_%SAMPLE%/filtered.fastq.gz"
JOB_QC=$(submit_array_step "QC" "01" "${SENTINEL_01}" \
    $(sbatch_flags "${TIME_QC}") \
    $(_dep "${JOB_SETUP}" "afterok") \
    --job-name="phage_01_qc" \
    --ntasks=1 --cpus-per-task="${THREADS}" --mem=8G \
    --output="${RESULTS_DIR}/logs/01_qc.%a.%j.out" \
    --wrap="export TMPDIR=${MAESTRO_TMPDIR} && source ${SCRIPT_DIR}/config.sh && activate_env && \
            ARRAY_IDX=\${SLURM_ARRAY_TASK_ID} bash ${STEPS_DIR}/01_qc.sh")
[[ "${JOB_QC}"       != "skipped_"* ]] && echo "  Job ID: ${JOB_QC}"

# ── Step 01b — OPTIONAL host depletion (only if HOST_DEPLETION=true) ──────────
JOB_HOSTDEP="skipped_step"
if [[ "${HOST_DEPLETION:-false}" == "true" ]]; then
    SENTINEL_01B="${RESULTS_DIR}/01b_hostdepleted/%BARCODE%_%SAMPLE%/filtered.host_depleted.fastq.gz"
    JOB_HOSTDEP=$(submit_array_step "host depletion" "01b" "${SENTINEL_01B}" \
        $(sbatch_flags "${TIME_HOSTDEPLETE}") \
        $(_dep "${JOB_QC}") \
        --job-name="phage_01b_hostdeplete" \
        --ntasks=1 --cpus-per-task="${THREADS}" --mem=24G \
        --output="${RESULTS_DIR}/logs/01b_hostdeplete.%a.%j.out" \
        --wrap="export TMPDIR=${MAESTRO_TMPDIR} && source ${SCRIPT_DIR}/config.sh && activate_env && \
                ARRAY_IDX=\${SLURM_ARRAY_TASK_ID} bash ${STEPS_DIR}/01b_hostdeplete.sh")
    [[ "${JOB_HOSTDEP}" != "skipped_"* ]] && echo "  Job ID: ${JOB_HOSTDEP}"
fi

# Step 02 depends on 01b when host depletion ran, else directly on 01
DEP_FOR_02="${JOB_HOSTDEP}"
[[ "${JOB_HOSTDEP}" == "skipped_step" ]] && DEP_FOR_02="${JOB_QC}"

# ── Step 02 ───────────────────────────────────────────────────────────────────
SENTINEL_02="${RESULTS_DIR}/02_assembly/%BARCODE%_%SAMPLE%/assembly.fasta"
JOB_ASSEMBLE=$(submit_array_step "assembly" "02" "${SENTINEL_02}" \
    $(sbatch_flags "${TIME_ASSEMBLE}") \
    $(_dep "${DEP_FOR_02}") \
    --job-name="phage_02_assemble" \
    --ntasks=1 --cpus-per-task="${THREADS}" --mem="${MAX_MEMORY_GB}G" \
    --gres=disk:20000 \
    --output="${RESULTS_DIR}/logs/02_assemble.%a.%j.out" \
    --wrap="export TMPDIR=${MAESTRO_TMPDIR} && source ${SCRIPT_DIR}/config.sh && activate_env && \
            ARRAY_IDX=\${SLURM_ARRAY_TASK_ID} bash ${STEPS_DIR}/02_assemble.sh")
[[ "${JOB_ASSEMBLE}"  != "skipped_"* ]] && echo "  Job ID: ${JOB_ASSEMBLE}"

# ── Step 03 ───────────────────────────────────────────────────────────────────
SENTINEL_03="${RESULTS_DIR}/03_polished/%BARCODE%_%SAMPLE%/polished.fasta"
JOB_POLISH=$(submit_array_step "polishing" "03" "${SENTINEL_03}" \
    $(sbatch_flags "${TIME_POLISH}") \
    $(_dep "${JOB_ASSEMBLE}") \
    --job-name="phage_03_polish" \
    --ntasks=1 --cpus-per-task="${THREADS}" --mem=24G \
    --gres=disk:10000 \
    --output="${RESULTS_DIR}/logs/03_polish.%a.%j.out" \
    --wrap="export TMPDIR=${MAESTRO_TMPDIR} && source ${SCRIPT_DIR}/config.sh && activate_env && \
            ARRAY_IDX=\${SLURM_ARRAY_TASK_ID} bash ${STEPS_DIR}/03_polish.sh")
[[ "${JOB_POLISH}"    != "skipped_"* ]] && echo "  Job ID: ${JOB_POLISH}"

# ── Step 04 ───────────────────────────────────────────────────────────────────
# Sentinel: asm_qc_summary.tsv (produced by assembly QC)
SENTINEL_04="${RESULTS_DIR}/04_assembly_qc/%BARCODE%_%SAMPLE%/asm_qc_summary.tsv"
JOB_ASMQC=$(submit_array_step "assembly QC (CheckV + coverage)" "04" "${SENTINEL_04}" \
    $(sbatch_flags "${TIME_ASMQC}") \
    $(_dep "${JOB_POLISH}") \
    --job-name="phage_04_assembly_qc" \
    --ntasks=1 --cpus-per-task="${THREADS}" --mem=24G \
    --output="${RESULTS_DIR}/logs/04_assembly_qc.%a.%j.out" \
    --wrap="export TMPDIR=${MAESTRO_TMPDIR} && source ${SCRIPT_DIR}/config.sh && activate_env && \
            ARRAY_IDX=\${SLURM_ARRAY_TASK_ID} bash ${STEPS_DIR}/04_assembly_qc.sh")
[[ "${JOB_ASMQC}"     != "skipped_"* ]] && echo "  Job ID: ${JOB_ASMQC}"

# ── Step 05 — Annotation (Pharokka) ──────────────────────────────────────────
SENTINEL_05="${RESULTS_DIR}/05_annotated/%BARCODE%_%SAMPLE%/%SAMPLE%.gff"
JOB_ANNOTATE=$(submit_array_step "annotation" "05" "${SENTINEL_05}" \
    $(sbatch_flags "${TIME_ANNOTATE}") \
    $(_dep "${JOB_POLISH}") \
    --job-name="phage_05_annotate" \
    --ntasks=1 --cpus-per-task="${THREADS}" --mem=32G \
    --output="${RESULTS_DIR}/logs/05_annotate.%a.%j.out" \
    --wrap="export TMPDIR=${MAESTRO_TMPDIR} && source ${SCRIPT_DIR}/config.sh && activate_env && \
            ARRAY_IDX=\${SLURM_ARRAY_TASK_ID} bash ${STEPS_DIR}/05_annotate.sh")
[[ "${JOB_ANNOTATE}"  != "skipped_"* ]] && echo "  Job ID: ${JOB_ANNOTATE}"

# ── Steps 06 + 07 (methylation — only if MOD_DEMUX_DIR is set) ───────────────
JOB_MODBASE="skipped_step"
JOB_METHYL_MM="skipped_step"
JOB_METHYL="skipped_step"
if [[ -n "${MOD_DEMUX_DIR:-}" ]]; then
    SENTINEL_06="${RESULTS_DIR}/06_modbase/%BARCODE%_%SAMPLE%/methylation.bed"
    JOB_MODBASE=$(submit_array_step "modified base detection" "06" "${SENTINEL_06}" \
        $(sbatch_flags "${TIME_MODBASE}") \
        $(_dep "${JOB_POLISH}") \
        --job-name="phage_06_modbase" \
        --ntasks=1 --cpus-per-task="${THREADS}" --mem=16G \
        --output="${RESULTS_DIR}/logs/06_modbase.%a.%j.out" \
        --wrap="export TMPDIR=${MAESTRO_TMPDIR} && source ${SCRIPT_DIR}/config.sh && activate_env && \
                ARRAY_IDX=\${SLURM_ARRAY_TASK_ID} bash ${STEPS_DIR}/06_modbase.sh")
    [[ "${JOB_MODBASE}"   != "skipped_"* ]] && echo "  Job ID: ${JOB_MODBASE}"

    # Step 07a — MicrobeMod per barcode (ARRAY, parallel) — the slow part.
    # The gather (07b) then reuses these per-barcode TSVs and builds the report.
    SENTINEL_07A="${RESULTS_DIR}/07_methyl/02_microbemod/%BARCODE%_%SAMPLE%/microbemod.bed"
    JOB_METHYL_MM=$(submit_array_step "MicrobeMod (per barcode)" "07a" "${SENTINEL_07A}" \
        $(sbatch_flags "${TIME_METHYL_CALL}") \
        $(_dep "${JOB_MODBASE}") \
        --job-name="phage_07a_methyl_call" \
        --ntasks=1 --cpus-per-task="${THREADS}" --mem=24G \
        --output="${RESULTS_DIR}/logs/07a_methyl_call.%a.%j.out" \
        --wrap="export TMPDIR=${MAESTRO_TMPDIR} && source ${SCRIPT_DIR}/config.sh && activate_env && \
                ARRAY_IDX=\${SLURM_ARRAY_TASK_ID} bash ${STEPS_DIR}/07a_methyl_call.sh")
    [[ "${JOB_METHYL_MM}" != "skipped_"* ]] && echo "  Job ID: ${JOB_METHYL_MM}"

    if [[ -z "${STEP_ONLY}" || "${STEP_ONLY}" == "07b" ]]; then
        echo "Step 07b — methylation gather (stats + annotation + HTML)"
        DEP_07B=$(_dep "${JOB_METHYL_MM}")
        JOB_METHYL=$(sbatch \
            $(sbatch_flags "${TIME_METHYL}") \
            --job-name="phage_07b_methyl_report" \
            ${DEP_07B} \
            --ntasks=1 --cpus-per-task=8 --mem=24G \
            --output="${RESULTS_DIR}/logs/07b_methyl_report.%j.out" \
            --wrap="export TMPDIR=${MAESTRO_TMPDIR} && source ${SCRIPT_DIR}/config.sh && activate_env && \
                    bash ${STEPS_DIR}/07b_methyl_report.sh")
        echo "  Job ID: ${JOB_METHYL}"
    fi
else
    echo "Steps 06+07 — SKIPPED  (MOD_DEMUX_DIR not set in config.sh)"
fi

# ── Step 09 — Run summary (single global job) ────────────────────────────────
JOB_SUMMARY="skipped_step"
if [[ -z "${STEP_ONLY}" || "${STEP_ONLY}" == "09" ]]; then
    echo "Step 09 — summary  (depends on 04, 05, 06)"
    DEP_09=$(_dep_multi "${JOB_ASMQC}" "${JOB_ANNOTATE}" "${JOB_MODBASE}")
    JOB_SUMMARY=$(sbatch \
        $(sbatch_flags "${TIME_SUMMARY}") \
        --job-name="phage_09_summary" \
        ${DEP_09} \
        --ntasks=1 --cpus-per-task=2 --mem=4G \
        --output="${RESULTS_DIR}/logs/09_summary.%j.out" \
        --wrap="export TMPDIR=${MAESTRO_TMPDIR} && source ${SCRIPT_DIR}/config.sh && activate_env && bash ${STEPS_DIR}/09_summary.sh")
    echo "  Job ID: ${JOB_SUMMARY}"
fi

# ── Step 10 — Final HTML report (single global job, after all analysis steps) ─
JOB_FINAL="skipped_step"
if [[ -z "${STEP_ONLY}" || "${STEP_ONLY}" == "10" ]]; then
    echo "Step 10 — final HTML report  (depends on 09, 07b)"
    DEP_10=$(_dep_multi "${JOB_SUMMARY}" "${JOB_METHYL}")
    JOB_FINAL=$(sbatch \
        $(sbatch_flags "${TIME_FINAL_REPORT}") \
        --job-name="phage_10_report" \
        ${DEP_10} \
        --ntasks=1 --cpus-per-task=2 --mem=4G \
        --output="${RESULTS_DIR}/logs/10_final_report.%j.out" \
        --wrap="export TMPDIR=${MAESTRO_TMPDIR} && source ${SCRIPT_DIR}/config.sh && activate_env && \
                bash ${STEPS_DIR}/10_final_report.sh")
    echo "  Job ID: ${JOB_FINAL}"
fi

# =============================================================================
# Final summary
# =============================================================================
_show_jid() {
    local id="$1"
    case "$id" in
        skipped_step)   echo "—" ;;
        skipped_resume) echo "done" ;;
        *)              echo "$id" ;;
    esac
}

echo ""
echo "════════════════════════════════════════════════════════"
if [[ "${DRY_RUN}" == "true" ]]; then
    echo "  DRY-RUN complete — no jobs were submitted."
    echo "  If the output above looks correct, run without --dry-run."
else
    echo "  Submission complete."
fi

echo ""
echo "  Job chain:"
echo "    00=$(_show_jid ${JOB_SETUP}) → 01=$(_show_jid ${JOB_QC}) → 01b=$(_show_jid ${JOB_HOSTDEP}) → 02=$(_show_jid ${JOB_ASSEMBLE}) → 03=$(_show_jid ${JOB_POLISH})"
echo "    03 → 04 asm_qc=$(_show_jid ${JOB_ASMQC})   05 annotate=$(_show_jid ${JOB_ANNOTATE})   06 modbase=$(_show_jid ${JOB_MODBASE})"
echo "    06 → 07a methyl=$(_show_jid ${JOB_METHYL_MM}) → 07b report=$(_show_jid ${JOB_METHYL})"
echo "    04+05+06 → 09 summary=$(_show_jid ${JOB_SUMMARY})"
echo "    09+07b → 10 report=$(_show_jid ${JOB_FINAL})"
echo ""
echo "  Results : ${RESULTS_DIR}"
echo "  Logs    : ${RESULTS_DIR}/logs/"
echo ""
if [[ "${DRY_RUN}" == "false" ]]; then
    echo "  Monitor : squeue -u \$(whoami)"
fi
echo "════════════════════════════════════════════════════════"
