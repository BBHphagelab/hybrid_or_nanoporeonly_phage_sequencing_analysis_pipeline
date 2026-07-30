#!/usr/bin/env bash
# =============================================================================
# 00_install_env.sh — One-time cluster environment setup (Maestro / Institut Pasteur)
# =============================================================================
# !! RUN THIS INTERACTIVELY ON A LOGIN NODE — NOT via sbatch !!
#
#   ssh csivelle@maestro.pasteur.fr
#   cd /pasteur/helix/projects/Phage3C/_projects/Nanopore/CS_phage_nanopore/
#   bash 00_install_env.sh
#
# Why micromamba and NOT conda:
#   Anaconda/conda is license-prohibited on Institut Pasteur's Maestro cluster.
#   Pasteur policy: use virtualenv+pip for pure Python, or micromamba for Bioconda packages.
#   micromamba is a standalone C++ reimplementation of conda — fully compatible with
#   conda-forge and bioconda channels, but requires no Anaconda license.
#
# What this does:
#   1. Installs micromamba binary to ~/.mamba/bin/ (no sudo needed)
#   2. Creates the 'nanopore_phage' environment with all pipeline tools
#   3. Installs Medaka via pip inside the env (avoids dependency conflicts)
#   4. Downloads modkit binary from GitHub (nanoporetech/modkit)
#   5. Downloads the Pharokka database (~1.5 GB)
#   6. Verifies all tools and prints a post-install checklist
#
# Takes ~20–40 minutes depending on cluster internet speed.
# Rerunning is safe — each step checks if already done.
# =============================================================================
set -euo pipefail

MAMBA_ENV="nanopore_phage"
MAMBA_ROOT_PREFIX="${HOME}/.mamba"           # micromamba root — must match config.sh

# modkit is NOT in bioconda — released by ONT as a pre-compiled binary only.
# Check https://github.com/nanoporetech/modkit/releases for newer versions.
MODKIT_VERSION="0.6.2"

# Pharokka database download destination.
# Must match PHAROKKA_DB in config.sh.
PHAROKKA_DB_DIR="/pasteur/helix/projects/Phage3C/_projects/databases/pharokka_db"  # ★

MAMBA_BIN="${MAMBA_ROOT_PREFIX}/bin/micromamba"

# ---------------------------------------------------------------------------
# Colours
# ---------------------------------------------------------------------------
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
ok()   { echo -e "${GREEN}[OK]${NC}  $*"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $*"; }
info() { echo -e "       $*"; }
die()  { echo -e "${RED}[ERROR]${NC} $*" >&2; exit 1; }

echo ""
echo "════════════════════════════════════════════════════════"
echo "  Nanopore PHAGE pipeline — Maestro environment setup"
echo "  $(date '+%Y-%m-%d %H:%M:%S')"
echo "════════════════════════════════════════════════════════"
echo ""

# ---------------------------------------------------------------------------
# Step 1 — Install micromamba
# ---------------------------------------------------------------------------
echo "── Step 1: micromamba ──────────────────────────────────"

if [[ -x "${MAMBA_BIN}" ]]; then
    ok "micromamba already installed: $(${MAMBA_BIN} --version)"
else
    info "Downloading micromamba binary to ${MAMBA_ROOT_PREFIX}/bin/ ..."
    mkdir -p "${MAMBA_ROOT_PREFIX}/bin"
    curl -fsSL "https://micro.mamba.pm/api/micromamba/linux-64/latest" \
        | tar -xvj -C "${MAMBA_ROOT_PREFIX}/bin" --strip-components=1 bin/micromamba
    chmod +x "${MAMBA_BIN}"
    ok "micromamba installed: $(${MAMBA_BIN} --version)"
    "${MAMBA_BIN}" shell init --shell bash --root-prefix "${MAMBA_ROOT_PREFIX}"
    warn "Shell hook added to ~/.bashrc (effective for new interactive shells)."
    warn "For this current shell, run: source ~/.bashrc"
fi

export MAMBA_ROOT_PREFIX
eval "$("${MAMBA_BIN}" shell hook --shell bash)"
ok "micromamba shell hook loaded"

# ---------------------------------------------------------------------------
# Step 2 — Create the micromamba environment
# ---------------------------------------------------------------------------
echo ""
echo "── Step 2: nanopore_phage environment ──────────────────"

if [[ -d "${MAMBA_ROOT_PREFIX}/envs/${MAMBA_ENV}" ]]; then
    ok "Environment '${MAMBA_ENV}' already exists"
    warn "To reinstall from scratch:"
    warn "  micromamba env remove -n ${MAMBA_ENV} && bash 00_install_env.sh"
else
    info "Creating '${MAMBA_ENV}' — this will take 10–20 minutes..."

    info "  Pass 1/2 — core utilities + alignment tools..."
    micromamba create -y -n "${MAMBA_ENV}" \
        --channel-priority strict \
        -c conda-forge -c bioconda \
        python=3.10 \
        \
        `# ── Read QC ──────────────────────────────────────────` \
        chopper \
        nanostat \
        filtlong \
        \
        `# ── Alignment / mapping ──────────────────────────────` \
        minimap2 \
        samtools \
        \
        `# ── General utilities ────────────────────────────────` \
        seqkit \
        pigz \
        any2fasta \
        bedtools

    ok "Pass 1 done."

    info "  Pass 2/2 — assemblers, polishers, annotators (longer)..."
    micromamba install -y -n "${MAMBA_ENV}" \
        --channel-priority strict \
        -c conda-forge -c bioconda \
        \
        `# ── Assembly ─────────────────────────────────────────` \
        flye \
        \
        `# ── Assembly post-processing ─────────────────────────` \
        dnaapler \
        \
        `# ── Annotation ───────────────────────────────────────` \
        pharokka \
        \
        `# ── Assembly QC (step 04) ────────────────────────────` \
        checkv \
        \
        `# ── Illumina short-read polishing (step 03, optional) ` \
        bwa \
        polypolish \
        pypolca

    ok "Environment '${MAMBA_ENV}' created (both passes complete)."
fi

# ---------------------------------------------------------------------------
# Step 3 — Install Medaka via pip (inside the micromamba env)
# ---------------------------------------------------------------------------
echo ""
echo "── Step 3: Medaka (via pip) ─────────────────────────────"

micromamba activate "${MAMBA_ENV}"

if python -c "import medaka" 2>/dev/null; then
    ok "Medaka already installed: $(medaka --version 2>&1 | head -1)"
else
    info "Installing medaka via pip..."
    pip install medaka
    ok "Medaka installed: $(medaka --version 2>&1 | head -1)"
fi

# ---------------------------------------------------------------------------
# Optional tools (NOT installed here — Autocycler & sub-assemblers come from
# Maestro Lmod modules or a dedicated env; PhageTerm is optional). Informational.
# ---------------------------------------------------------------------------
echo ""
echo "── Optional: Autocycler / PhageTerm ─────────────────────"
if command -v autocycler &>/dev/null; then
    ok "autocycler: $(autocycler --version 2>&1 | head -1)"
else
    warn "autocycler not in PATH — step 02 prefers it if available (Maestro: 'module load Autocycler'"
    warn "  with assembler modules loaded first, or set AUTOCYCLER_ENV). Falls back to Flye+Filtlong."
fi
for _a in raven miniasm plassembler; do
    if command -v "$_a" &>/dev/null; then ok "sub-assembler ${_a} present"; else warn "sub-assembler ${_a} not found (optional)"; fi
done
if command -v PhageTerm.py &>/dev/null || command -v phageterm &>/dev/null; then
    ok "PhageTerm present (terminus detection, step 03)"
else
    warn "PhageTerm not in PATH — step 03 falls back to dnaapler. Set PHAGETERM_BIN in config.sh to enable."
fi

micromamba deactivate

# ---------------------------------------------------------------------------
# Step 4 — Modkit binary installation
# ---------------------------------------------------------------------------
# modkit is released by ONT as a pre-compiled binary only — NOT in bioconda.
# We place it inside the micromamba env's bin/ so it behaves like any env tool.
# ---------------------------------------------------------------------------
echo ""
echo "── Step 4: Modkit binary ────────────────────────────────"

MODKIT_BIN="${MAMBA_ROOT_PREFIX}/envs/${MAMBA_ENV}/bin/modkit"

if [[ -x "${MODKIT_BIN}" ]]; then
    ok "modkit already installed: $(${MODKIT_BIN} --version 2>&1 | head -1)"
else
    info "Downloading modkit v${MODKIT_VERSION} binary..."
    # u16 = Ubuntu 16 build — most compatible across cluster node OS versions
    MODKIT_URL="https://github.com/nanoporetech/modkit/releases/download/v${MODKIT_VERSION}/modkit_v${MODKIT_VERSION}_u16_x86_64.tar.gz"
    TMP_MODKIT=$(mktemp -d)
    curl -fsSL -L "${MODKIT_URL}" | tar -xz -C "${TMP_MODKIT}"
    MODKIT_EXTRACTED=$(find "${TMP_MODKIT}" -name "modkit" -type f | head -1)
    [[ -z "${MODKIT_EXTRACTED}" ]] && die "modkit binary not found in tarball — check MODKIT_VERSION and URL"
    mv "${MODKIT_EXTRACTED}" "${MODKIT_BIN}"
    chmod +x "${MODKIT_BIN}"
    rm -rf "${TMP_MODKIT}"
    ok "modkit installed: $(${MODKIT_BIN} --version 2>&1 | head -1)"
fi

# ---------------------------------------------------------------------------
# Step 5 — Pharokka database download
# ---------------------------------------------------------------------------
echo ""
echo "── Step 5: Pharokka database ────────────────────────────"

if [[ -d "${PHAROKKA_DB_DIR}" ]] && [[ -n "$(ls -A "${PHAROKKA_DB_DIR}" 2>/dev/null)" ]]; then
    ok "Pharokka database already exists at ${PHAROKKA_DB_DIR}"
else
    info "Downloading Pharokka database → ${PHAROKKA_DB_DIR}"
    info "(~1.5 GB, uses PHROGs HMM profiles)"
    mkdir -p "${PHAROKKA_DB_DIR}"
    "${MAMBA_ROOT_PREFIX}/envs/${MAMBA_ENV}/bin/pharokka_db_install.py" \
        --db_dir "${PHAROKKA_DB_DIR}"
    ok "Pharokka database downloaded."
fi

# ---------------------------------------------------------------------------
# Step 5b — CheckV database download (step 04 assembly QC)
# ---------------------------------------------------------------------------
echo ""
echo "── Step 5b: CheckV database ─────────────────────────────"
CHECKV_DB_DIR="/pasteur/helix/projects/Phage3C/_projects/databases/checkv-db"  # ★ match CHECKV_DB in config.sh
if compgen -G "${CHECKV_DB_DIR}/checkv-db-*" >/dev/null 2>&1; then
    ok "CheckV database already present under ${CHECKV_DB_DIR}"
else
    info "Downloading CheckV database → ${CHECKV_DB_DIR} (~1.5 GB)"
    mkdir -p "${CHECKV_DB_DIR}"
    if checkv download_database "${CHECKV_DB_DIR}" 2>/dev/null; then
        ok "CheckV database downloaded."
    else
        warn "checkv download_database failed (network?). Step 04 will skip CheckV until CHECKV_DB is populated."
    fi
fi

# ---------------------------------------------------------------------------
# Step 6 — MicrobeMod + methylation analysis dependencies
# ---------------------------------------------------------------------------
# MicrobeMod/1.1.0 is a CLUSTER MODULE on Maestro — do NOT install via pip or conda.
# Its dependencies (prodigal, blast+, hmmer, cath-tools, meme) are also cluster modules,
# loaded automatically at runtime in step 07a.
#
# bedtools is already installed in pass 1 above.
# ---------------------------------------------------------------------------
echo ""
echo "── Step 6: MicrobeMod (Maestro module — no install needed) ─"
warn "MicrobeMod and its deps (prodigal, blast+, hmmer, cath-tools, meme) are"
warn "Maestro cluster modules — loaded automatically at runtime by step 08."
warn "No conda installation needed."

# ---------------------------------------------------------------------------
# Step 7 — Verify all tools
# ---------------------------------------------------------------------------
echo ""
echo "── Step 7: Verification ─────────────────────────────────"
micromamba activate "${MAMBA_ENV}"

TOOLS_OK=true
check_tool() {
    local tool="$1" cmd="${2:-$1}"
    if command -v "$cmd" &>/dev/null; then
        local ver; ver=$($cmd --version 2>&1 | head -1 | tr -d '\n')
        ok "${tool}: ${ver}"
    else
        echo -e "${RED}[MISSING]${NC} ${tool} — check installation"
        TOOLS_OK=false
    fi
}

check_tool "chopper"
check_tool "NanoStat"
check_tool "filtlong"
check_tool "Flye"          "flye"
check_tool "Dnaapler"      "dnaapler"
check_tool "Pharokka"      "pharokka.py"
check_tool "Medaka"        "medaka"
check_tool "CheckV"        "checkv"
check_tool "Polypolish"    "polypolish"
check_tool "Pypolca"       "pypolca"
check_tool "bwa"
check_tool "Samtools"      "samtools"
check_tool "minimap2"
check_tool "Modkit"        "modkit"
check_tool "bedtools"
check_tool "seqkit"
check_tool "any2fasta"
check_tool "MicrobeMod"    "MicrobeMod"

micromamba deactivate

# ---------------------------------------------------------------------------
# Post-install summary
# ---------------------------------------------------------------------------
echo ""
if [[ "${TOOLS_OK}" == "true" ]]; then
    echo -e "${GREEN}All tools verified OK.${NC}"
else
    echo -e "${RED}Some tools are MISSING — see [MISSING] lines above.${NC}"
    echo -e "  call_methylation missing? → load module: module load MicrobeMod"
fi
echo ""
echo "════════════════════════════════════════════════════════"
echo "  Post-install checklist"
echo "════════════════════════════════════════════════════════"
echo "  1. Edit config.sh — Section 1 (RUN_NAME, EMAIL)"
echo "  2. Edit config.sh — Section 3 (CHECKV_DB, PHAROKKA_DB)"
echo "  3. Vérifier MOD_DEMUX_DIR dans config.sh (activé par défaut → bam_pass/)"
echo "  4. (option) Autocycler: module load Autocycler (Maestro) OU AUTOCYCLER_ENV ; PhageTerm: PHAGETERM_BIN"
echo "  5. Lancer : bash submit_all.sh --check  puis  bash submit_all.sh"
echo "════════════════════════════════════════════════════════"
