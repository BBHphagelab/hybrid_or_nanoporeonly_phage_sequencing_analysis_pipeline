#!/usr/bin/env bash
# =============================================================================
# submit_modbase_prep.sh — Extract POD5 archive + run Dorado mod-base calling
# =============================================================================
# Run ONCE before step 07 (modified base detection).
#
# What this does:
#   Job 1 (CPU): Extract the TAR archive containing POD5 files to NFS
#   Job 2 (GPU): Run Dorado basecaller with 6mA + 5mCG_5hmCG models on the PODs
#   Job 3 (CPU): Run Dorado demux to split the BAM per barcode
#
# After this completes:
#   1. Set MOD_DEMUX_DIR in config.sh to the demux output directory
#   2. Re-submit step 07: bash submit_all.sh   (it chains correctly from step 03)
#      Or submit just step 07:
#        source config.sh
#        sbatch ... --dependency=afterok:<polish_job> ... bash steps/07_modbase.sh
#
# Prerequisites:
#   - The TAR file must already be on NFS (accessible from compute nodes)
#   - Dorado binary: downloaded automatically to DORADO_BIN_DIR if missing
#   - Dorado model: downloaded automatically by dorado basecaller if missing
#   - GPU partition available on Maestro (Maestro has a 'gpu' partition)
#     Check with: sinfo -p gpu
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# ---------------------------------------------------------------------------
# ★ CONFIGURE THESE ★
# ---------------------------------------------------------------------------

# Full path to the TAR file on NFS (Z:\ maps to /pasteur/helix/projects/...)
POD5_TAR="/pasteur/helix/projects/Phage3C/_projects/Nanopore/LS_16042026_1stnanoporebact_AT_csivelle__Autonomous__Sat_Apr_18_06h21m10_2026.tar"

# Where to extract the POD5 files (needs ~50-200 GB free on NFS)
POD5_EXTRACT_DIR="/pasteur/helix/projects/Phage3C/_projects/Nanopore/pod5_pass_extracted"

# Where Dorado will write its output
MOD_CALLS_BAM="/pasteur/helix/projects/Phage3C/_projects/Nanopore/modbase/mod_calls.bam"
MOD_DEMUX_OUTPUT="/pasteur/helix/projects/Phage3C/_projects/Nanopore/modbase/demuxed"

# Dorado binary directory — downloaded here if not already present
DORADO_BIN_DIR="/pasteur/helix/projects/Phage3C/_projects/tools/dorado"

# Dorado model name — must match the chemistry used during sequencing
# This run: R10.4.1 e8.2 400bps HAC v5
# NOTE: Dorado 0.9.x does NOT auto-download models from a name string.
# The model is downloaded via 'dorado download' into DORADO_MODELS_DIR,
# then the downloaded directory path is passed to basecaller.
DORADO_MODEL="dna_r10.4.1_e8.2_400bps_hac@v5.0.0"
DORADO_MODELS_DIR="/pasteur/helix/projects/Phage3C/_projects/tools/dorado_models"

# Kit used for barcode demultiplexing (check your MinKNOW run config)
BARCODE_KIT="SQK-RBK114-24"

# GPU partition name on Maestro (check with: sinfo -p gpu)
GPU_PARTITION="gpu"
GPU_ACCOUNT="${ACCOUNT:-}"  # leave empty if not needed

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
gpu_flags() {
    echo -n "--parsable"
    echo -n " --partition=${GPU_PARTITION}"
    [[ -n "${GPU_ACCOUNT}" ]] && echo -n " --account=${GPU_ACCOUNT}"
    [[ -n "${EMAIL}" ]] && echo -n " --mail-user=${EMAIL} --mail-type=begin,end,fail"
}
cpu_flags() {
    echo -n "--parsable"
    [[ -n "${PARTITION}" ]] && echo -n " --partition=${PARTITION}"
    [[ -n "${ACCOUNT}" ]] && echo -n " --account=${ACCOUNT}"
    [[ -n "${EMAIL}" ]] && echo -n " --mail-user=${EMAIL} --mail-type=begin,end,fail"
}

mkdir -p "${RESULTS_DIR}/logs" "$(dirname "$MOD_CALLS_BAM")"

# ---------------------------------------------------------------------------
# Job 1 — Extract POD5 TAR archive
# ---------------------------------------------------------------------------
# TAR extraction is CPU + I/O bound. Single node, no GPU needed.
# The archive likely contains pod5_pass/barcode01/ ... structure.
# We extract into POD5_EXTRACT_DIR; the actual pod5_pass dir will be inside.
echo "Submitting Job 1 — Extract POD5 archive..."
echo "  TAR: ${POD5_TAR}"
echo "  Destination: ${POD5_EXTRACT_DIR}"

JOB_EXTRACT=$(sbatch \
    $(cpu_flags) \
    --job-name="modbase_01_extract" \
    --time="04:00:00" \
    --ntasks=1 --cpus-per-task=4 --mem=8G \
    --output="${RESULTS_DIR}/logs/modbase_01_extract.%j.out" \
    --wrap="
set -euo pipefail
echo '['\$(date '+%H:%M:%S')'] Extracting ${POD5_TAR} ...'
mkdir -p '${POD5_EXTRACT_DIR}'
# Use pv if available for progress, otherwise plain tar
if command -v pv &>/dev/null; then
    pv '${POD5_TAR}' | tar -x -C '${POD5_EXTRACT_DIR}'
else
    tar -xvf '${POD5_TAR}' -C '${POD5_EXTRACT_DIR}' 2>&1 | tail -1
    echo '['\$(date '+%H:%M:%S')'] Extraction complete.'
fi
echo 'Contents of ${POD5_EXTRACT_DIR}:'
ls -lh '${POD5_EXTRACT_DIR}' | head -20
echo '['\$(date '+%H:%M:%S')'] Job 1 done.'
")
echo "  Job 1 submitted: ${JOB_EXTRACT}"

# ---------------------------------------------------------------------------
# Job 2 — Dorado basecaller with modified bases
# ---------------------------------------------------------------------------
# Dorado needs a GPU. On Maestro, request 1 GPU on the gpu partition.
# The input is the extracted pod5_pass directory.
# --modified-bases 6mA 5mCG_5hmCG tells Dorado to output MM/ML BAM tags.
#
# NOTE on POD5 path: the TAR likely extracts to a structure like:
#   ${POD5_EXTRACT_DIR}/pod5_pass/   (containing per-barcode subdirs or flat .pod5)
# Dorado basecaller can handle nested directories recursively.
# Adjust POD5_INPUT if needed after extraction.

POD5_INPUT="${POD5_EXTRACT_DIR}"   # dorado will recurse to find all .pod5 files

echo "Submitting Job 2 — Dorado mod-base basecalling (GPU)..."
JOB_DORADO=$(sbatch \
    $(gpu_flags) \
    --job-name="modbase_02_dorado" \
    --time="12:00:00" \
    --dependency=afterok:${JOB_EXTRACT} \
    --ntasks=1 --cpus-per-task=16 --mem=64G \
    --gres=gpu:1 \
    --output="${RESULTS_DIR}/logs/modbase_02_dorado.%j.out" \
    --wrap="
set -euo pipefail
DORADO_BIN='${DORADO_BIN_DIR}/bin/dorado'

# Download Dorado if not installed
if [[ ! -x \"\${DORADO_BIN}\" ]]; then
    echo '['\$(date '+%H:%M:%S')'] Dorado not found — downloading...'
    mkdir -p '${DORADO_BIN_DIR}'
    # Latest Dorado for Linux x86_64 — check https://github.com/nanoporetech/dorado/releases
    DORADO_URL='https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.1-linux-x64.tar.gz'
    curl -fsSL -L \"\${DORADO_URL}\" | tar -xz -C '${DORADO_BIN_DIR}' --strip-components=1
    echo '['\$(date '+%H:%M:%S')'] Dorado installed: '\$(\${DORADO_BIN} --version)
fi

echo '['\$(date '+%H:%M:%S')'] Dorado version: '\$(\${DORADO_BIN} --version)
echo '['\$(date '+%H:%M:%S')'] Input POD5 dir: ${POD5_INPUT}'
echo '['\$(date '+%H:%M:%S')'] Model: ${DORADO_MODEL}'
echo '['\$(date '+%H:%M:%S')'] Modifications: 6mA 5mCG_5hmCG'

# List POD5 files found
n_pod5=\$(find '${POD5_INPUT}' -name '*.pod5' | wc -l)
echo '['\$(date '+%H:%M:%S')'] Found '\${n_pod5}' POD5 files'

# Download the model if not already present.
# Dorado 0.9.x requires a local model directory — it does NOT auto-download
# from a model name string. 'dorado download' fetches the model weights once.
DORADO_MODELS_DIR='${DORADO_MODELS_DIR}'
mkdir -p \"\${DORADO_MODELS_DIR}\"
MODEL_PATH=\"\${DORADO_MODELS_DIR}/${DORADO_MODEL}\"
if [[ ! -d \"\${MODEL_PATH}\" ]]; then
    echo '['\$(date '+%H:%M:%S')'] Downloading Dorado model: ${DORADO_MODEL}'
    \${DORADO_BIN} download \
        --model '${DORADO_MODEL}' \
        --directory \"\${DORADO_MODELS_DIR}\"
    echo '['\$(date '+%H:%M:%S')'] Model downloaded to '\"\${MODEL_PATH}\"
else
    echo '['\$(date '+%H:%M:%S')'] Model already present: '\"\${MODEL_PATH}\"
fi

\${DORADO_BIN} basecaller \
    \"\${MODEL_PATH}\" \
    '${POD5_INPUT}' \
    --modified-bases 6mA 5mCG_5hmCG \
    --kit-name '${BARCODE_KIT}' \
    --recursive \
    --device cuda:all \
    > '${MOD_CALLS_BAM}'

echo '['\$(date '+%H:%M:%S')'] Dorado basecalling done.'
echo 'Output BAM: '\$(ls -lh '${MOD_CALLS_BAM}')
echo '['\$(date '+%H:%M:%S')'] Job 2 done.'
")
echo "  Job 2 submitted: ${JOB_DORADO}"

# ---------------------------------------------------------------------------
# Job 3 — Dorado demux (split BAM per barcode)
# ---------------------------------------------------------------------------
echo "Submitting Job 3 — Dorado demux (per-barcode BAM split)..."
JOB_DEMUX=$(sbatch \
    $(cpu_flags) \
    --job-name="modbase_03_demux" \
    --time="02:00:00" \
    --dependency=afterok:${JOB_DORADO} \
    --ntasks=1 --cpus-per-task=8 --mem=16G \
    --output="${RESULTS_DIR}/logs/modbase_03_demux.%j.out" \
    --wrap="
set -euo pipefail
DORADO_BIN='${DORADO_BIN_DIR}/bin/dorado'
mkdir -p '${MOD_DEMUX_OUTPUT}'

echo '['\$(date '+%H:%M:%S')'] Running Dorado demux...'
# IMPORTANT: do NOT pass --kit-name here.
# The basecaller (Job 2) already wrote BC:Z: barcode tags into mod_calls.bam.
# Passing --kit-name again causes demux to re-classify from scratch (ignoring
# those tags), which produces ~96% unclassified. Without --kit-name, demux
# simply reads the existing BC:Z: tags and splits by barcode — recovering
# the full classification rate from basecalling.
\${DORADO_BIN} demux \
    --output-dir '${MOD_DEMUX_OUTPUT}' \
    --emit-summary \
    '${MOD_CALLS_BAM}'

echo '['\$(date '+%H:%M:%S')'] Demux done. Per-barcode BAM files:'
ls -lh '${MOD_DEMUX_OUTPUT}'/*.bam 2>/dev/null | head -20

echo ''
echo '================================================================'
echo 'NEXT STEPS:'
echo '1. Check that barcode01.bam ... barcode07.bam exist in:'
echo '   ${MOD_DEMUX_OUTPUT}'
echo ''
echo '2. Edit config.sh and set:'
echo '   MOD_DEMUX_DIR=\"${MOD_DEMUX_OUTPUT}\"'
echo ''
echo '3. Submit step 07:'
echo '   cd /pasteur/helix/projects/Phage3C/_projects/Nanopore/cluster'
echo '   source config.sh'
echo '   sbatch --time=03:00:00 --array=0-6 --ntasks=1 --cpus-per-task=16 --mem=16G \\'
echo '       --job-name=nano_07_modbase \\'
echo '       --output=\${RESULTS_DIR}/logs/07_modbase.%a.%j.out \\'
echo '       --wrap=\"export TMPDIR=/local/scratch/tmp && source \$(pwd)/config.sh && activate_env && ARRAY_IDX=\\\${SLURM_ARRAY_TASK_ID} bash \$(pwd)/steps/07_modbase.sh\"'
echo '================================================================'
echo '['\$(date '+%H:%M:%S')'] Job 3 done.'
")
echo "  Job 3 submitted: ${JOB_DEMUX}"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "All modbase prep jobs submitted."
echo ""
echo "Dependency chain:"
echo "  [Extract TAR: ${JOB_EXTRACT}]"
echo "     → [Dorado basecall: ${JOB_DORADO}]"
echo "         → [Dorado demux: ${JOB_DEMUX}]"
echo "             → then manually submit step 07"
echo ""
echo "Monitor: squeue -u $(whoami)"
echo "Logs:    ${RESULTS_DIR}/logs/modbase_*.out"
echo ""
echo "NOTE: Job 2 (Dorado) requests a GPU node (partition: ${GPU_PARTITION})."
echo "  Check GPU queue: squeue -p ${GPU_PARTITION}"
echo "  If Maestro uses a different GPU partition name, edit GPU_PARTITION in this script."
echo ""
echo "After Job 3 completes, set in config.sh:"
echo "  MOD_DEMUX_DIR=\"${MOD_DEMUX_OUTPUT}\""
