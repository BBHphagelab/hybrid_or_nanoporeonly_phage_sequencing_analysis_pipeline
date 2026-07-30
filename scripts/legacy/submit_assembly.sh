#!/usr/bin/env bash
# =============================================================================
# submit_assembly.sh
# -----------------------------------------------------------------------------
# Submits assemble_phages_hostdepleted.sh to SLURM as a SINGLE compute-node job
# (so it doesn't hog the submit node). All requested barcodes run sequentially
# inside one allocation — no concurrency, no race conditions on the shared phage
# reference / summary file. Resources mirror the pipeline's step 02
# (cpus=THREADS, mem=MAX_MEMORY_GB, --gres=disk, TMPDIR=/local/scratch/tmp).
#
# USAGE
#   bash submit_assembly.sh --host-ref /path/host_combined.fasta [options]
#
# OPTIONS (passed through to assemble_phages_hostdepleted.sh)
#   --host-ref PATH    REQUIRED. Combined host genome FASTA.
#   --assembler NAME   autocycler (default) or flye.
#   --phage-ref PATH   Optional. Phage reference FASTA (else built from config refs).
#   --barcodes LIST    Optional. Comma list; default = all barcodes.
#   --force            Re-run even if assembly.fasta exists.
# SUBMISSION-ONLY OPTIONS
#   --time HH:MM:SS    SLURM walltime (default 12:00:00).
#   --dry-run          Print the sbatch command, submit nothing.
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

MAESTRO_TMPDIR="/local/scratch/tmp"
WALLTIME="12:00:00"   # Autocycler runs several assemblers x subsamples x barcodes
DRY_RUN=false
HOST_REF=""
PASS_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --host-ref)  HOST_REF="$2"; PASS_ARGS+=( --host-ref "$2" ); shift 2 ;;
        --phage-ref) PASS_ARGS+=( --phage-ref "$2" ); shift 2 ;;
        --assembler) PASS_ARGS+=( --assembler "$2" ); shift 2 ;;
        --barcodes)  PASS_ARGS+=( --barcodes "$2" ); shift 2 ;;
        --force)     PASS_ARGS+=( --force ); shift ;;
        --time)      WALLTIME="$2"; shift 2 ;;
        --dry-run)   DRY_RUN=true; shift ;;
        -h|--help)   sed -n '2,30p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *)           echo "[ERROR] Unknown argument: $1 (use --help)" >&2; exit 1 ;;
    esac
done

[[ -n "$HOST_REF" ]] || { echo "[ERROR] --host-ref is required." >&2; exit 1; }
[[ -f "$HOST_REF" ]] || { echo "[ERROR] host-ref not found: $HOST_REF" >&2; exit 1; }

mkdir -p "${RESULTS_DIR}/logs"

# Quote the pass-through args safely for embedding in --wrap
PASS_STR=""
for a in "${PASS_ARGS[@]}"; do PASS_STR+=" $(printf '%q' "$a")"; done

# Forward selected env vars into the job (so e.g. AUTOCYCLER_ENV=... reaches it
# regardless of SLURM --export settings).
ENVPREFIX=""
for _v in AUTOCYCLER_ENV AUTOCYCLER_ASSEMBLERS AUTOCYCLER_MODULES HOST_COV_FRAC FILTLONG_TARGET_COV; do
    [[ -n "${!_v:-}" ]] && ENVPREFIX+="${_v}=$(printf '%q' "${!_v}") "
done

INNER="export TMPDIR=${MAESTRO_TMPDIR} && source ${SCRIPT_DIR}/config.sh && activate_env && ${ENVPREFIX}bash ${SCRIPT_DIR}/assemble_phages_hostdepleted.sh${PASS_STR}"

# Build sbatch flags (mirror step 02: cpus=THREADS, mem=MAX_MEMORY_GB, disk scratch)
SBATCH_ARGS=(
    --job-name="phage_02b_hostdepleted"
    --ntasks=1
    --cpus-per-task="${THREADS}"
    --mem="${MAX_MEMORY_GB}G"
    --gres=disk:20000
    --time="${WALLTIME}"
    --output="${RESULTS_DIR}/logs/02b_hostdepleted.%j.out"
)
[[ -n "${PARTITION:-}" ]] && SBATCH_ARGS+=( --partition="${PARTITION}" )
[[ -n "${ACCOUNT:-}"   ]] && SBATCH_ARGS+=( --account="${ACCOUNT}" )
[[ -n "${EMAIL:-}"     ]] && SBATCH_ARGS+=( --mail-user="${EMAIL}" --mail-type=begin,end,fail )

echo "Submitting host-depleted assembly as a single SLURM job:"
echo "  cpus=${THREADS}  mem=${MAX_MEMORY_GB}G  disk=20000  time=${WALLTIME}"
echo "  host-ref: ${HOST_REF}"
echo "  log: ${RESULTS_DIR}/logs/02b_hostdepleted.%j.out"
echo ""

if [[ "$DRY_RUN" == true ]]; then
    echo "[DRY-RUN] sbatch ${SBATCH_ARGS[*]} --wrap=\"${INNER}\""
    exit 0
fi

JOBID=$(sbatch --parsable "${SBATCH_ARGS[@]}" --wrap="${INNER}")
echo "Submitted batch job ${JOBID}"
echo ""
echo "Monitor:   squeue -u \$(whoami)"
echo "Live log:  tail -f ${RESULTS_DIR}/logs/02b_hostdepleted.${JOBID}.out"
echo ""
echo "When it finishes, check ${RESULTS_DIR}/02b_hostdepleted/decontamination_summary.tsv,"
echo "then:  bash clean_for_reassembly.sh --yes  &&  bash submit_all.sh"
