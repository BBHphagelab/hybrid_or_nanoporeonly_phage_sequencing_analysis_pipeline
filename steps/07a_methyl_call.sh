#!/usr/bin/env bash
# =============================================================================
# Step 07a — MicrobeMod call_methylation, ONE barcode per SLURM array task
# =============================================================================
# The slow part of methylation analysis (MicrobeMod = prodigal + blast + hmmer +
# STREME) is split here into a job array so all barcodes run in PARALLEL.
# The gather step (07b_methyl_report.sh) then reuses these per-barcode outputs:
# it finds the TSV already present and only rebuilds the BED, then produces the
# stats, comparison, annotation and HTML report quickly.
#
# Output (per barcode): $RESULTS_DIR/07_methyl/02_microbemod/<bc>_<name>/
#   <bc>_<name>*.tsv     MicrobeMod methylation table
#   microbemod.bed       contig,start,end,motif,mod_type,strand (built here)
#
# Usage (SLURM array):  managed by submit_all.sh
# Manual (single):      ARRAY_IDX=0 bash steps/07a_methyl_call.sh
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/config.sh"

MM_DIR="${RESULTS_DIR}/07_methyl/02_microbemod"
mkdir -p "${MM_DIR}"

# ---------------------------------------------------------------------------
# Tools / MicrobeMod cluster modules (same setup as the gather step)
# MicrobeMod 1.1.0 needs modkit 0.2.x: load cluster modules, then PREPEND the
# conda bin so the env's modkit 0.2.x wins over any cluster 0.3.x.
# ---------------------------------------------------------------------------
activate_env
require_tool samtools
if command -v module &>/dev/null; then
    module load prodigal blast+ hmmer cath-tools modkit meme MicrobeMod 2>/dev/null \
        || echo "[WARN] some cluster modules failed to load" >&2
    [[ -d "${HOME}/MicrobeMod_patched" ]] && export PYTHONPATH="${HOME}/MicrobeMod_patched:${PYTHONPATH:-}"
    export PATH="${MAMBA_ROOT_PREFIX}/envs/${MAMBA_ENV}/bin:${PATH}"
fi

MICROBEMOD_BIN=""
for _c in MicrobeMod microbemod; do
    command -v "${_c}" &>/dev/null && { MICROBEMOD_BIN="${_c}"; break; }
done
[[ -z "${MICROBEMOD_BIN}" ]] && die "MicrobeMod not found in PATH after module load."

STREME_PATH=""
if command -v streme &>/dev/null; then
    STREME_PATH="$(command -v streme)"
else
    STREME_PATH=$(find /opt/gensoft/exe/meme /opt/gensoft/exe/MEME 2>/dev/null \
                  -name streme -type f 2>/dev/null | head -1 || true)
fi

# ---------------------------------------------------------------------------
# Resolve barcode from array index
# ---------------------------------------------------------------------------
mapfile -t BARCODES < <(list_barcodes)
IDX="${ARRAY_IDX:-${SLURM_ARRAY_TASK_ID:-0}}"
if [[ "${IDX}" -ge "${#BARCODES[@]}" ]]; then
    echo "[SKIP] ARRAY_IDX=${IDX} >= barcodes (${#BARCODES[@]}) — nothing to do."
    exit 0
fi
bc="${BARCODES[${IDX}]}"
name=$(get_sample_field "${bc}" "sample_name")

step_banner "07a — MicrobeMod  [${bc} / ${name}]"

ALIGNED_BAM="${RESULTS_DIR}/06_modbase/${bc}_${name}/aligned_mods.bam"
POLISHED="${RESULTS_DIR}/03_polished/${bc}_${name}/polished.fasta"
MM_OUT="${MM_DIR}/${bc}_${name}"
MM_BED="${MM_OUT}/microbemod.bed"

[[ -f "${ALIGNED_BAM}" ]] || { echo "[SKIP] No aligned_mods.bam for ${bc} — run step 06 first."; exit 0; }
[[ -f "${POLISHED}"   ]] || { echo "[SKIP] No polished.fasta for ${bc}."; exit 0; }

if [[ -s "${MM_BED}" && "${FORCE:-false}" != "true" ]]; then
    echo "[SKIP] microbemod.bed already exists for ${bc} (set FORCE=true to rerun)."
    exit 0
fi

# Load the modern cluster samtools FIRST. The micromamba env samtools is
# prehistoric (< 0.1.19, see load_samtools_cluster) and reading a modern
# high-depth mod-tagged BAM with it is the prime suspect for the multi-hour hang
# observed before MicrobeMod ever launched (no microbemod_run.log was created).
load_samtools_cluster

# Force modkit 0.2.x to the front of PATH. MicrobeMod 1.1.0 calls
# `modkit pileup … --only-tabs`; modkit ≥ 0.3 removed that flag, so the env's
# modkit 0.6.2 makes MicrobeMod fail with exit 2 (→ rc=1). Non-fatal if missing:
# MicrobeMod will fail and the sentinel + self-check will record why.
use_modkit_02x || true

# Pre-check: BAM must carry MM/Mm modification tags. Bounded by a timeout so a
# pathological read (wrong samtools, unindexed/huge BAM) can never eat the wall.
# `grep -c ... || true` inside makes the inner shell always exit 0 (grep exits 1
# on zero matches), so the ONLY way the outer command fails is the timeout firing.
# A timeout is treated as "assume tagged" — MicrobeMod itself validates and the
# sentinel logic surfaces a real absence of tags.
_mm=$(timeout 600 bash -c 'samtools view "$1" 2>/dev/null | head -2000 | grep -c "MM:Z\|Mm:Z" || true' _ "${ALIGNED_BAM}") || _mm="TIMEOUT"
if [[ "${_mm}" == "TIMEOUT" ]]; then
    echo "[WARN] ${bc}: MM/Mm tag pre-check timed out (600s) — proceeding anyway (MicrobeMod will validate)." >&2
elif [[ "${_mm}" -eq 0 ]]; then
    echo "[WARN] ${bc}: aligned_mods.bam has no MM/Mm tags in first 2000 reads — skipping MicrobeMod." >&2
    exit 0
else
    echo "[INFO] ${bc}: MM tag check OK (${_mm} tagged reads in first 2000)."
fi

mkdir -p "${MM_OUT}"

# ── Subsample very-high-coverage mod-BAM for MicrobeMod speed ─────────────────
# Motif discovery (STREME) and the gene searches scale with data volume; at
# 5 000–13 000x MicrobeMod can run for many hours. Subsample to ~MICROBEMOD_MAX_COV
# (>> COV_MIN_METHYL, so calls stay confident). Full BAM is left untouched.
MM_BAM="${ALIGNED_BAM}"
# Subsampling runs BEFORE the MicrobeMod timeout, so bound it separately: on a
# multi-GB high-depth BAM `samtools coverage`/`view -s` could otherwise eat the
# whole walltime and starve the MicrobeMod timeout (the failure mode that made
# the task get SLURM-SIGKILLed at the wall with no sentinel).
SUBSAMPLE_TO="${MICROBEMOD_SUBSAMPLE_TIMEOUT:-3600}"
if [[ "${MICROBEMOD_MAX_COV:-0}" -gt 0 ]]; then
    mean_dp=$(timeout "${SUBSAMPLE_TO}" samtools coverage "${ALIGNED_BAM}" 2>/dev/null \
              | awk 'NR>1{s+=$7*$3; w+=$3} END{if(w>0) printf "%.0f", s/w; else print 0}')
    if [[ "${mean_dp:-0}" -gt "${MICROBEMOD_MAX_COV}" ]]; then
        frac=$(awk -v t="${MICROBEMOD_MAX_COV}" -v m="${mean_dp}" 'BEGIN{printf "%.4f", t/m}')
        echo "[INFO] ${bc}: mod-BAM mean depth ~${mean_dp}x -> subsampling to ~${MICROBEMOD_MAX_COV}x (frac ${frac})."
        SUB="${MM_OUT}/aligned_mods.subsampled.bam"
        if timeout "${SUBSAMPLE_TO}" samtools view -@ "${THREADS:-8}" -b -s "${frac}" "${ALIGNED_BAM}" > "${SUB}" 2>/dev/null \
             && timeout "${SUBSAMPLE_TO}" samtools index "${SUB}" 2>/dev/null; then
            MM_BAM="${SUB}"
        else
            echo "[WARN] ${bc}: subsampling failed or timed out (${SUBSAMPLE_TO}s) — using full BAM." >&2
            rm -f "${SUB}" "${SUB}.bai" 2>/dev/null || true
        fi
    else
        echo "[INFO] ${bc}: mod-BAM mean depth ~${mean_dp:-?}x ≤ ${MICROBEMOD_MAX_COV}x — no subsampling."
    fi
fi
MM_STREME_ARG=()
[[ -n "${STREME_PATH:-}" ]] && ! command -v streme &>/dev/null && MM_STREME_ARG=(--streme_path "${STREME_PATH}")
MM_RUN_LOG="${MM_OUT}/microbemod_run.log"

echo "[INFO] ${bc}: running ${MICROBEMOD_BIN} call_methylation (threads ${THREADS:-8}, timeout ${MICROBEMOD_TIMEOUT:-28800}s)..."
# Non-fatal on failure: exit 0 so the gather's afterok dependency still runs.
# `timeout -k 120`: if MicrobeMod ignores the SIGTERM at the deadline, escalate to
# SIGKILL 120 s later so the task cannot hang past its budget. rc=124 on timeout.
_mmrc=0
timeout -k 120 "${MICROBEMOD_TIMEOUT:-28800}" "${MICROBEMOD_BIN}" call_methylation \
    --bam              "${MM_BAM}" \
    --reference        "${POLISHED}" \
    --output_directory "${MM_OUT}" \
    --output_prefix    "${bc}_${name}" \
    --threads          "${THREADS:-8}" \
    "${MM_STREME_ARG[@]}" \
    > "${MM_RUN_LOG}" 2>&1 || _mmrc=$?
# --------------------------------------------------------------------------
# Record a machine-readable failure sentinel so the gather step (07b) and the
# run summary (09) can SURFACE the failure instead of silently omitting the
# sample. We still exit 0 so the array task's afterok dependency keeps running,
# but the failure is no longer invisible: STATUS file + loud banner + a modkit
# self-diagnostic (MicrobeMod hides modkit's real stderr behind exit 2).
# --------------------------------------------------------------------------
STATUS_FILE="${MM_OUT}/microbemod.STATUS"
write_status() {  # write_status <OK|FAILED|TIMEOUT> <rc>
    { printf "barcode\t%s\n" "${bc}"
      printf "sample\t%s\n"  "${name}"
      printf "status\t%s\n"  "$1"
      printf "rc\t%s\n"      "$2"
      printf "modkit\t%s (%s)\n" "$(command -v modkit 2>/dev/null || echo 'not found')" \
                                 "$(modkit --version 2>&1 | head -1 || echo '?')"
    } > "${STATUS_FILE}"
}

if [[ ${_mmrc} -eq 0 ]]; then
    echo "[INFO] ${bc}: MicrobeMod done."
    write_status OK 0
elif [[ ${_mmrc} -eq 124 ]]; then
    write_status TIMEOUT 124
    echo "[WARN] ${bc}: MicrobeMod TIMED OUT (${MICROBEMOD_TIMEOUT:-28800}s) — sentinel written, non-fatal." >&2
    tail -25 "${MM_RUN_LOG}" >&2 || true
    exit 0
else
    write_status FAILED "${_mmrc}"
    echo "════════════════════════════════════════════════════════" >&2
    echo "[FAIL] ${bc} (${name}): MicrobeMod call_methylation FAILED (rc=${_mmrc})" >&2
    echo "       Sentinel: ${STATUS_FILE}" >&2
    echo "       Log     : ${MM_RUN_LOG} (last lines):" >&2
    echo "════════════════════════════════════════════════════════" >&2
    tail -25 "${MM_RUN_LOG}" >&2 || true
    # MicrobeMod swallows modkit's stderr; probe modkit directly for the real cause.
    echo "[DIAG] modkit self-check ($(modkit --version 2>&1 | head -1)):" >&2
    modkit pileup -t 2 "${MM_BAM}" "${MM_OUT}/selfcheck_low.bed" \
        -r "${POLISHED}" --only-tabs --filter-threshold 0.66 >&2 2>&1 || \
        echo "[DIAG] ↑ this is the real modkit error MicrobeMod hid (rc=$?)." >&2
    echo "[WARN] ${bc}: continuing (non-fatal) so the gather step can run." >&2
    exit 0
fi

# Build BED from the MicrobeMod TSV (col3=Contig col4=Pos col5=Mod col6=Strand col19=Motif)
MM_TSV=$(find "${MM_OUT}" -name "*.tsv" 2>/dev/null | head -1 || true)
if [[ -n "${MM_TSV}" && -s "${MM_TSV}" ]]; then
    awk 'NR>1 && $6!="." && $5!="" && $19!="" {
        printf "%s\t%d\t%d\t%s\t%s\t%s\n", $3, $4, $4+1, $19, $5, $6
    }' "${MM_TSV}" > "${MM_BED}"
    echo "[INFO] ${bc}: $(wc -l < "${MM_BED}") motif sites → ${MM_BED}"
else
    echo "[WARN] ${bc}: MicrobeMod produced no TSV." >&2
fi
echo "[INFO] Step 07a done for ${bc}."
