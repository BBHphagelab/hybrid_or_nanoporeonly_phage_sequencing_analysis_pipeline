#!/usr/bin/env bash
# =============================================================================
# submit.sh — Phage nanopore pipeline — CLI entry point
# =============================================================================
# Soumet tous les jobs SLURM du pipeline et rend immédiatement la main.
# Remplace l'édition manuelle de config.sh pour les paramètres de run.
#
# Usage:
#   bash submit.sh --run-dir /path/to/run [OPTIONS]
#
# ── Obligatoire ───────────────────────────────────────────────────────────────
#   --run-dir DIR        Répertoire de run MinKNOW (contient MinKnow_output/)
#
# ── Chemins (optionnels, dérivés de --run-dir par défaut) ────────────────────
#   --out-dir DIR        Répertoire de résultats
#                        (défaut: <run-dir>/results_nanopore)
#   --mod-bam-dir DIR    Dossier bam_pass de MinKNOW (BAMs mod-tagged par barcode)
#                        (défaut: <run-dir>/MinKnow_output/bam_pass)
#
# ── Paramètres d'analyse (valeurs par défaut dans config.sh) ─────────────────
#   --genome-size SIZE   Taille estimée du génome phagique, ex: 50k, 100k, 200k
#   --medaka-model STR   Modèle Medaka (ex: r1041_e82_400bps_hac_v5.0.0)
#
# ── Paramètres SLURM ─────────────────────────────────────────────────────────
#   --email ADDR         Email pour notifications SLURM
#   --partition NAME     Partition SLURM (ex: fast, common, long)
#   --account NAME       Compte SLURM
#   --threads N          CPUs par job (défaut: config.sh)
#   --memory-gb N        RAM en GB (défaut: config.sh)
#
# ── Contrôle d'exécution ─────────────────────────────────────────────────────
#   --from STEP          Reprendre depuis l'étape N (0-8).
#                        Les étapes < N sont supposées déjà terminées.
#                        Exemple: --from 3  →  soumet les étapes 3,4,6,7,8
#   --dry-run            Affiche les commandes sbatch sans soumettre
#   --check              Vérifications pre-flight uniquement, pas de soumission
#
# ── Exemples ─────────────────────────────────────────────────────────────────
#   # Run complet
#   bash submit.sh --run-dir /pasteur/helix/.../20260427_ONT_run1
#
#   # Run avec paramètres personnalisés
#   bash submit.sh \
#       --run-dir    /pasteur/helix/.../20260427_ONT_run1 \
#       --genome-size 80k \
#       --partition  fast \
#       --account    phage3c \
#       --threads    12
#
#   # Reprendre depuis l'étape de méthylation
#   bash submit.sh --run-dir /pasteur/.../run1 --from 7
#
#   # Vérifier sans soumettre
#   bash submit.sh --run-dir /pasteur/.../run1 --dry-run
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# =============================================================================
# Valeurs par défaut des arguments CLI
# =============================================================================
CLI_RUN_DIR=""
CLI_OUT_DIR=""
CLI_MOD_BAM_DIR=""
CLI_GENOME_SIZE=""
CLI_MEDAKA_MODEL=""
CLI_EMAIL=""
CLI_PARTITION=""
CLI_ACCOUNT=""
CLI_THREADS=""
CLI_MEMORY_GB=""
FROM_STEP=0
DRY_RUN=false
CHECK_ONLY=false

# =============================================================================
# Parsing des arguments
# =============================================================================
usage() { grep '^#' "$0" | grep -v '^#!/' | sed 's/^# \{0,3\}//' | head -60; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-dir)      CLI_RUN_DIR="$2";      shift 2 ;;
        --out-dir)      CLI_OUT_DIR="$2";      shift 2 ;;
        --mod-bam-dir)  CLI_MOD_BAM_DIR="$2";  shift 2 ;;
        --genome-size)  CLI_GENOME_SIZE="$2";  shift 2 ;;
        --medaka-model) CLI_MEDAKA_MODEL="$2"; shift 2 ;;
        --email)        CLI_EMAIL="$2";        shift 2 ;;
        --partition)    CLI_PARTITION="$2";    shift 2 ;;
        --account)      CLI_ACCOUNT="$2";      shift 2 ;;
        --threads)      CLI_THREADS="$2";      shift 2 ;;
        --memory-gb)    CLI_MEMORY_GB="$2";    shift 2 ;;
        --from)         FROM_STEP="$2";        shift 2 ;;
        --dry-run)      DRY_RUN=true;          shift ;;
        --check)        CHECK_ONLY=true;       shift ;;
        -h|--help)      usage; exit 0 ;;
        *) echo "[ERROR] Option inconnue : $1" >&2; usage >&2; exit 1 ;;
    esac
done

[[ -n "$CLI_RUN_DIR" ]] || {
    echo "[ERROR] --run-dir est obligatoire." >&2
    echo "        Exemple : bash submit.sh --run-dir /pasteur/.../20260427_ONT_run1" >&2
    exit 1
}

# =============================================================================
# Chargement de config.sh (paramètres statiques : bases de données, env, outils)
# =============================================================================
source "${SCRIPT_DIR}/config.sh"

# =============================================================================
# Application des overrides CLI
# =============================================================================
RUN_DIR="$CLI_RUN_DIR"
RESULTS_DIR="${CLI_OUT_DIR:-${RUN_DIR}/results_nanopore}"
FASTQ_PASS_DIR="${RUN_DIR}/MinKnow_output/fastq_pass"
SEQ_SUMMARY_DIR="${RUN_DIR}/MinKnow_output"
MOD_DEMUX_DIR="${CLI_MOD_BAM_DIR:-${RUN_DIR}/MinKnow_output/bam_pass}"

[[ -n "$CLI_GENOME_SIZE"  ]] && GENOME_SIZE_ESTIMATE="$CLI_GENOME_SIZE"
[[ -n "$CLI_MEDAKA_MODEL" ]] && MEDAKA_MODEL="$CLI_MEDAKA_MODEL"
[[ -n "$CLI_EMAIL"        ]] && EMAIL="$CLI_EMAIL"
[[ -n "$CLI_PARTITION"    ]] && PARTITION="$CLI_PARTITION"
[[ -n "$CLI_ACCOUNT"      ]] && ACCOUNT="$CLI_ACCOUNT"
[[ -n "$CLI_THREADS"      ]] && THREADS="$CLI_THREADS"
[[ -n "$CLI_MEMORY_GB"    ]] && MAX_MEMORY_GB="$CLI_MEMORY_GB"

# =============================================================================
# Mode dry-run : sbatch devient une fonction d'affichage
# =============================================================================
_DRY_COUNTER=0
if [[ "$DRY_RUN" == "true" ]]; then
    sbatch() {
        _DRY_COUNTER=$(( _DRY_COUNTER + 1 ))
        echo "  [DRY-RUN] sbatch $*" | fold -s -w 110 >&2
        echo "DRY${_DRY_COUNTER}"
    }
    echo ""
    echo "════════════════════════════════════════════════════════"
    echo "  DRY-RUN — aucun job ne sera soumis"
    echo "════════════════════════════════════════════════════════"
fi

# =============================================================================
# Pre-flight checks
# =============================================================================
echo ""
echo "════════════════════════════════════════════════════════"
echo "  Vérifications pre-flight"
echo "════════════════════════════════════════════════════════"

PREFLIGHT_OK=true
_ok()   { echo "  [OK]   $*"; }
_warn() { echo "  [WARN] $*"; }
_fail() { echo "  [FAIL] $*"; PREFLIGHT_OK=false; }

echo ""
echo "  Paramètres du run :"
echo "    RUN_DIR        : ${RUN_DIR}"
echo "    RESULTS_DIR    : ${RESULTS_DIR}"
echo "    GENOME_SIZE    : ${GENOME_SIZE_ESTIMATE}"
echo "    MOD_DEMUX_DIR  : ${MOD_DEMUX_DIR}"
echo "    THREADS        : ${THREADS}"
echo "    MEMORY_GB      : ${MAX_MEMORY_GB}"
echo "    PARTITION      : ${PARTITION:-'(défaut cluster)'}"
echo "    FROM_STEP      : ${FROM_STEP}"
echo ""

[[ -d "${RUN_DIR}" ]] \
    && _ok  "RUN_DIR existe" \
    || _fail "RUN_DIR introuvable : ${RUN_DIR}"

[[ -d "${FASTQ_PASS_DIR}" ]] \
    && _ok  "FASTQ_PASS_DIR existe  ($(find "${FASTQ_PASS_DIR}" -maxdepth 1 -type d -name 'barcode*' 2>/dev/null | wc -l) dossiers barcode)" \
    || _fail "FASTQ_PASS_DIR introuvable : ${FASTQ_PASS_DIR}"

MAMBA_BIN="${MAMBA_ROOT_PREFIX}/bin/micromamba"
[[ -x "${MAMBA_BIN}" ]] \
    && _ok  "micromamba : ${MAMBA_BIN}" \
    || _fail "micromamba introuvable — exécuter : bash ${SCRIPT_DIR}/00_install_env.sh"

[[ -d "${MAMBA_ROOT_PREFIX}/envs/${MAMBA_ENV}" ]] \
    && _ok  "Environnement '${MAMBA_ENV}' présent" \
    || _fail "Environnement '${MAMBA_ENV}' absent — exécuter : bash ${SCRIPT_DIR}/00_install_env.sh"

[[ -d "${PHAROKKA_DB}" ]] \
    && _ok  "Pharokka DB : ${PHAROKKA_DB}" \
    || _fail "Pharokka DB introuvable : ${PHAROKKA_DB}"

if [[ -d "${MOD_DEMUX_DIR}" ]]; then
    n_bams=$(find "${MOD_DEMUX_DIR}" -name "*.bam" -maxdepth 2 2>/dev/null | wc -l)
    _ok "MOD_DEMUX_DIR existe  (${n_bams} fichier(s) BAM)"
else
    _warn "MOD_DEMUX_DIR absent : ${MOD_DEMUX_DIR}  (steps 07+08 soumis mais échoueront si le dossier n'est pas créé avant leur exécution)"
fi

echo ""
if [[ "${PREFLIGHT_OK}" == "true" ]]; then
    echo "  ✅  Toutes les vérifications sont OK."
else
    echo "  ❌  Des vérifications ont échoué — corriger avant de soumettre."
    [[ "${DRY_RUN}" == "false" && "${CHECK_ONLY}" == "false" ]] && exit 1
fi
echo ""

[[ "${CHECK_ONLY}" == "true" ]] && { echo "  (mode --check : arrêt ici)"; exit 0; }

# =============================================================================
# Préparation
# =============================================================================
mkdir -p "${RESULTS_DIR}/logs"

# ── Dimensionnement du tableau SLURM ─────────────────────────────────────────
SAMPLE_SHEET="${RESULTS_DIR}/sample_sheet.tsv"
if [[ -f "$SAMPLE_SHEET" ]]; then
    N_SAMPLES=$(count_barcodes)
    log "Array sizing : ${N_SAMPLES} barcodes depuis sample_sheet.tsv existant"
else
    N_SAMPLES=$(find "${FASTQ_PASS_DIR}" -maxdepth 1 -mindepth 1 -type d -name 'barcode*' 2>/dev/null | wc -l | tr -d ' ')
    [[ "${N_SAMPLES}" -eq 0 ]] && { _warn "Impossible de détecter le nombre de samples — fallback à 1"; N_SAMPLES=1; }
    log "Array sizing : ${N_SAMPLES} barcodes détectés dans FASTQ_PASS_DIR"
fi
ARRAY_RANGE="0-$(( N_SAMPLES - 1 ))"

# ── Fichier d'overrides pour les jobs SLURM ──────────────────────────────────
# Les jobs SLURM sourcent config.sh (paramètres statiques) puis ce fichier
# (overrides CLI), de sorte que les valeurs CLI prennent toujours la priorité.
OVERRIDES_FILE="${RESULTS_DIR}/.run_overrides.sh"
cat > "${OVERRIDES_FILE}" << OVERRIDES
# Auto-généré par submit.sh le $(date '+%Y-%m-%d %H:%M:%S') — ne pas éditer
RUN_DIR='${RUN_DIR}'
RESULTS_DIR='${RESULTS_DIR}'
FASTQ_PASS_DIR='${FASTQ_PASS_DIR}'
SEQ_SUMMARY_DIR='${SEQ_SUMMARY_DIR}'
MOD_DEMUX_DIR='${MOD_DEMUX_DIR}'
GENOME_SIZE_ESTIMATE='${GENOME_SIZE_ESTIMATE}'
MEDAKA_MODEL='${MEDAKA_MODEL}'
THREADS='${THREADS}'
MAX_MEMORY_GB='${MAX_MEMORY_GB}'
OVERRIDES

# =============================================================================
# Helpers de soumission
# =============================================================================
MAESTRO_TMPDIR="/local/scratch/tmp"
declare -A JOB_IDS=()

# dep_of STEP → job ID si soumis, "" sinon
dep_of() {
    local id="${JOB_IDS[$1]:-}"
    [[ "$id" == "skipped" || -z "$id" ]] && echo "" || echo "$id"
}

# submit STEP LABEL DEP IS_ARRAY TIME CPUS MEM EXTRA_FLAGS STEP_SCRIPT
submit() {
    local step="$1" label="$2" dep="${3:-}" is_array="${4:-false}"
    local time="$5" cpus="$6" mem="${7}G" extra="${8:-}" step_script="$9"

    # Étape ignorée si --from demande de démarrer plus loin
    if [[ "$step" -lt "$FROM_STEP" ]]; then
        JOB_IDS["$step"]="skipped"
        printf '  Étape %02d (%s) — ignorée (--from %d)\n' "$step" "$label" "$FROM_STEP"
        return
    fi

    # Construction de la commande wrap pour le job SLURM
    local wrap="export TMPDIR=${MAESTRO_TMPDIR}"
    wrap+=" && source ${SCRIPT_DIR}/config.sh"
    wrap+=" && source ${OVERRIDES_FILE}"
    wrap+=" && activate_env"
    if [[ "$is_array" == "true" ]]; then
        wrap+=" && ARRAY_IDX=\${SLURM_ARRAY_TASK_ID} bash ${STEPS_DIR}/${step_script}"
    else
        wrap+=" && bash ${STEPS_DIR}/${step_script}"
    fi

    # Construction des flags sbatch
    local log_tag
    log_tag="${RESULTS_DIR}/logs/$(printf '%02d' $step)_${label}"
    [[ "$is_array" == "true" ]] && log_tag+=".%a.%j.out" || log_tag+=".%j.out"

    local flags="--parsable"
    flags+=" --job-name=phage_$(printf '%02d' $step)_${label}"
    flags+=" --time=${time}"
    flags+=" --ntasks=1 --cpus-per-task=${cpus} --mem=${mem}"
    flags+=" --output=${log_tag}"
    [[ "$is_array" == "true" ]] && flags+=" --array=${ARRAY_RANGE}"
    [[ -n "${dep}" ]]            && flags+=" --dependency=afterok:${dep}"
    [[ -n "${PARTITION:-}" ]]    && flags+=" --partition=${PARTITION}"
    [[ -n "${ACCOUNT:-}" ]]      && flags+=" --account=${ACCOUNT}"
    [[ -n "${EMAIL:-}" ]]        && flags+=" --mail-user=${EMAIL} --mail-type=begin,end,fail"
    [[ -n "$extra" ]]            && flags+=" ${extra}"

    local job_id
    job_id=$(sbatch $flags --wrap="$wrap")
    JOB_IDS["$step"]="$job_id"

    printf '  Étape %02d (%s) → Job %-12s%s\n' \
        "$step" "$label" "$job_id" \
        "${dep:+  ← attend ${dep}}"
}

# =============================================================================
# Soumission des jobs
# =============================================================================
echo "════════════════════════════════════════════════════════"
[[ "$DRY_RUN" == "true" ]] \
    && echo "  Jobs qui SERAIENT soumis (dry-run) :" \
    || echo "  Soumission des jobs..."
echo "════════════════════════════════════════════════════════"
echo ""

#        STEP  LABEL        DEP              ARR    TIME                  CPUS      MEM              EXTRA                SCRIPT
submit   0  "setup"     ""               false  "${TIME_SETUP}"       2         4                ""                   "00_setup.sh"
submit   1  "qc"        "$(dep_of 0)"   true   "${TIME_QC}"          "${THREADS}" 8             ""                   "01_qc.sh"
submit   2  "assemble"  "$(dep_of 1)"   true   "${TIME_ASSEMBLE}"    "${THREADS}" "${MAX_MEMORY_GB}" "--gres=disk:20000" "02_assemble.sh"
submit   3  "polish"    "$(dep_of 2)"   true   "${TIME_POLISH}"      "${THREADS}" 16             "--gres=disk:10000"  "03_polish.sh"
submit   4  "annotate"  "$(dep_of 3)"   true   "${TIME_ANNOTATE}"    "${THREADS}" 8              ""                   "04_annotate.sh"
submit   6  "summary"   "$(dep_of 4)"   false  "${TIME_SUMMARY}"     2         4                ""                   "06_summary.sh"
submit   7  "modbase"   "$(dep_of 3)"   true   "${TIME_MODBASE}"     "${THREADS}" 16             ""                   "07_modbase.sh"
submit   8  "methyl"    "$(dep_of 7)"   false  "${TIME_METHYL}"      8         24               ""                   "08_methyl_analysis.sh"

# =============================================================================
# Résumé final
# =============================================================================
echo ""
echo "════════════════════════════════════════════════════════"
if [[ "$DRY_RUN" == "true" ]]; then
    echo "  DRY-RUN terminé — aucun job soumis."
    echo "  Si la sortie ci-dessus est correcte : relancer sans --dry-run."
else
    echo "  Jobs soumis. Le terminal est libre."
fi
echo ""
echo "  Chaîne :"
echo "    00=${JOB_IDS[0]:-—} → 01=${JOB_IDS[1]:-—} → 02=${JOB_IDS[2]:-—} → 03=${JOB_IDS[3]:-—}"
echo "    03 → 04=${JOB_IDS[4]:-—} → 06=${JOB_IDS[6]:-—}"
echo "    03 → 07=${JOB_IDS[7]:-—} → 08=${JOB_IDS[8]:-—}"
echo ""
echo "  Résultats : ${RESULTS_DIR}"
echo "  Logs      : ${RESULTS_DIR}/logs/"
echo ""
if [[ "$DRY_RUN" == "false" ]]; then
    echo "  Suivi   : squeue -u $(whoami)"
    # Collect non-skipped job IDs for cancel command
    ALL_JOBS=()
    for step in 0 1 2 3 4 6 7 8; do
        id="${JOB_IDS[$step]:-}"
        [[ -n "$id" && "$id" != "skipped" ]] && ALL_JOBS+=("$id")
    done
    [[ ${#ALL_JOBS[@]} -gt 0 ]] && \
        echo "  Annuler : scancel ${ALL_JOBS[*]}"
fi
echo "════════════════════════════════════════════════════════"
