#!/usr/bin/env bash
# =============================================================================
# Step 03 — Polishing: Medaka (ONT-only) or Polypolish+Pypolca (hybrid)
# =============================================================================
# ONT-only:  Medaka consensus using the HAC v5 model.
#            Corrects systematic ONT errors (especially homopolymers).
#            Model must match the basecalling model set in config.sh.
#
# Hybrid:    Two-stage Illumina polishing, following Ryan Wick's Perfect Bacterial
#            Genome Tutorial (2024):
#              1. Polypolish --careful  (read-pair based, great for indels/SNPs)
#              2. Pypolca --careful     (pileup-based, catches errors Polypolish misses,
#                                       especially in homopolymer runs)
#            No Medaka needed after Unicycler hybrid assembly.
#
# Output: $RESULTS_DIR/03_polished/<barcode>_<sample>/polished.fasta
#         This is the final assembly used by annotation and methylation steps.
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/config.sh"

activate_env
load_samtools_cluster   # samtools/1.21 — conda env's samtools is too old for sort/index

SAMPLE_SHEET="${RESULTS_DIR}/sample_sheet.tsv"
mapfile -t BARCODES < <(list_barcodes)
if [[ "${ARRAY_IDX}" -ge "${#BARCODES[@]}" ]]; then
    echo "[SKIP] ARRAY_IDX=${ARRAY_IDX} >= barcodes in sample_sheet (${#BARCODES[@]}) — rien a faire."
    exit 0
fi
BARCODE="${BARCODES[${ARRAY_IDX}]}"
SAMPLE_NAME=$(get_sample_field "$BARCODE" "sample_name")
MODE=$(get_sample_field "$BARCODE" "mode")

step_banner "03 — Polish  [${BARCODE} / ${SAMPLE_NAME} / mode: ${MODE}]"

ASSEMBLY="${RESULTS_DIR}/02_assembly/${BARCODE}_${SAMPLE_NAME}/assembly.fasta"
ONT_READS="${RESULTS_DIR}/01_qc/${BARCODE}_${SAMPLE_NAME}/filtered.fastq.gz"
OUT="${RESULTS_DIR}/03_polished/${BARCODE}_${SAMPLE_NAME}"
POLISHED="${OUT}/polished.fasta"
mkdir -p "$OUT"

# ── Input guard: skip this barcode cleanly if assembly is missing ────────────
# (allows afterany dependency chaining — one failed barcode never blocks others)
if [[ ! -s "$ASSEMBLY" ]]; then
    echo "[SKIP] ${BARCODE}_${SAMPLE_NAME}: assembly missing (${ASSEMBLY}) — step 02 did not complete. Nothing to do."
    exit 0
fi

# Purge any stale FASTA index left by a previous run. If assembly.fasta was
# re-generated (e.g. re-assembly with Autocycler) while an old .fai describing
# DIFFERENT contigs remained next to it, medaka/pysam reads the stale index and
# dies with:
#   [E::fai_retrieve] Failed to retrieve block: error reading file
#   ValueError: failure when retrieving sequence on 'contig_NN'
# (root cause of the barcode02 medaka failure). Removing it forces a correct,
# fresh index when medaka/minimap2/samtools reopen the draft.
rm -f "${ASSEMBLY}.fai" "${ASSEMBLY}.mmi"

if [[ -f "$POLISHED" && "${FORCE:-false}" != "true" ]]; then
    log "  Polished assembly exists — skipping (set FORCE=true to rerun)"
    exit 0
fi

# ---------------------------------------------------------------------------
# Medaka read downsampling (speed). Medaka gains nothing above ~100-200x and is
# very slow at the 1000-13000x coverage of these phage runs (this caused the
# step-03 time-limit failures). Downsample the ONT reads for Medaka only
# (filtlong keeps the longest/highest-quality reads). The hybrid Illumina polish
# still uses the full short-read set. Falls back to full reads if filtlong is absent.
# ---------------------------------------------------------------------------
MEDAKA_READS="$ONT_READS"
if command -v filtlong &>/dev/null; then
    _mtgt=$(( ${MEDAKA_MAX_COV:-200} * 100000 ))   # ~MEDAKA_MAX_COV x of a ~100 kb phage
    if filtlong --target_bases "$_mtgt" "$ONT_READS" 2>/dev/null | gzip > "${OUT}/medaka_reads.fastq.gz" \
         && [[ -s "${OUT}/medaka_reads.fastq.gz" ]]; then
        MEDAKA_READS="${OUT}/medaka_reads.fastq.gz"
        log "  Medaka: downsampled ONT reads to ~${MEDAKA_MAX_COV:-200}x for speed."
    fi
fi

# ---------------------------------------------------------------------------
# ONT-only: Medaka
# ---------------------------------------------------------------------------
if [[ "$MODE" == "ont" ]]; then
    require_tool medaka
    require_tool minimap2
    require_tool samtools

    log "  Running Medaka pipeline (inference API: minimap2 → inference → sequence)..."
    log "    Model  : ${MEDAKA_MODEL}"
    log "    Assembly: $(basename $ASSEMBLY)"
    log "    Reads   : $(basename $ONT_READS)"

    # 3-step medaka pipeline (this installation exposes inference/sequence subcommands,
    # not the unified 'medaka consensus' command):
    #   1. minimap2 map-ont → sorted BAM
    #   2. medaka inference BAM → predictions HDF5
    #   3. medaka sequence HDF5 + draft → consensus FASTA
    #
    # Run on LOCAL SCRATCH (TMPDIR) to avoid NFS I/O issues with large BAM files.
    MEDAKA_LOCAL="${TMPDIR:-/tmp}/medaka_${BARCODE}_${SAMPLE_NAME}_$$"
    rm -rf "$MEDAKA_LOCAL"
    mkdir -p "$MEDAKA_LOCAL"

    log "  Step 1 — Aligning reads to assembly (minimap2)..."
    minimap2 -ax map-ont \
        -t "${THREADS}" \
        "$ASSEMBLY" \
        "$MEDAKA_READS" \
    | samtools sort \
        -@ "${THREADS}" \
        -o "${MEDAKA_LOCAL}/aligned.bam" \
        2>"${MEDAKA_LOCAL}/minimap2.log"
    samtools index "${MEDAKA_LOCAL}/aligned.bam"

    log "  Step 2 — medaka inference..."
    medaka inference \
        "${MEDAKA_LOCAL}/aligned.bam" \
        "${MEDAKA_LOCAL}/predictions.hdf" \
        --model   "${MEDAKA_MODEL}" \
        --threads "${THREADS}" \
        2>&1 | tee "${OUT}/medaka_inference.log"

    [[ -f "${MEDAKA_LOCAL}/predictions.hdf" ]] \
        || die "medaka inference failed — predictions.hdf not found. See ${OUT}/medaka_inference.log"

    log "  Step 3 — medaka sequence (consensus FASTA)..."
    medaka sequence \
        "${MEDAKA_LOCAL}/predictions.hdf" \
        "$ASSEMBLY" \
        "${MEDAKA_LOCAL}/consensus.fasta" \
        --threads "${THREADS}" \
        2>&1 | tee "${OUT}/medaka_sequence.log"

    [[ -s "${MEDAKA_LOCAL}/consensus.fasta" ]] \
        || die "medaka sequence failed — consensus.fasta not found. See ${OUT}/medaka_sequence.log"

    cp "${MEDAKA_LOCAL}/consensus.fasta" "$POLISHED"
    # Copy BAM + HDF5 to NFS for record-keeping (non-fatal)
    mkdir -p "${OUT}/medaka"
    cp -r "$MEDAKA_LOCAL"/. "${OUT}/medaka"/ \
        || warn "  Could not copy medaka intermediates to NFS — polished.fasta already saved."
    rm -rf "$MEDAKA_LOCAL"
    log "  Medaka done."

# ---------------------------------------------------------------------------
# Hybrid: Polypolish --careful  +  Pypolca --careful
# (Ryan Wick Perfect Bacterial Genome Tutorial, 2024 recommendations)
# ---------------------------------------------------------------------------
elif [[ "$MODE" == "hybrid" ]]; then
    require_tool bwa
    require_tool polypolish
    require_tool pypolca

    ILLUMINA_R1=$(get_sample_field "$BARCODE" "illumina_r1")
    ILLUMINA_R2=$(get_sample_field "$BARCODE" "illumina_r2")
    check_file "$ILLUMINA_R1"
    check_file "$ILLUMINA_R2"

    # ── Stage 0: Medaka (ONT consensus) BEFORE Illumina polishing ────────────
    # The hybrid assemblies are now Flye (ONT) drafts (host-depleted reads), not
    # Unicycler. Flye drafts carry systematic ONT errors (homopolymers, methylated
    # motifs), so we Medaka-polish with the ONT reads FIRST, then let Polypolish +
    # Pypolca correct residual errors with the more accurate Illumina bases LAST
    # (canonical "long-read then short-read" order). Controlled by HYBRID_MEDAKA
    # (config.sh, default true). Non-fatal: on failure we keep the raw assembly.
    HYBRID_ASM="$ASSEMBLY"
    if [[ "${HYBRID_MEDAKA:-true}" == "true" ]] && command -v medaka &>/dev/null; then
        log "  ── Stage 0: Medaka (ONT pre-polish of the Flye draft) ──"
        check_file "$ONT_READS"
        MEDAKA_LOCAL="${TMPDIR:-/tmp}/medaka_${BARCODE}_${SAMPLE_NAME}_$$"
        rm -rf "$MEDAKA_LOCAL"; mkdir -p "$MEDAKA_LOCAL"
        set +e
        minimap2 -ax map-ont -t "${THREADS}" "$ASSEMBLY" "$MEDAKA_READS" \
            2>"${MEDAKA_LOCAL}/minimap2.log" \
          | samtools sort -@ "${THREADS}" -o "${MEDAKA_LOCAL}/aligned.bam" 2>>"${MEDAKA_LOCAL}/minimap2.log"
        _M1=$?
        samtools index "${MEDAKA_LOCAL}/aligned.bam" 2>>"${MEDAKA_LOCAL}/minimap2.log"
        medaka inference "${MEDAKA_LOCAL}/aligned.bam" "${MEDAKA_LOCAL}/predictions.hdf" \
            --model "${MEDAKA_MODEL}" --threads "${THREADS}" \
            2>&1 | tee "${OUT}/medaka_inference.log"; _M2=${PIPESTATUS[0]}
        medaka sequence "${MEDAKA_LOCAL}/predictions.hdf" "$ASSEMBLY" \
            "${MEDAKA_LOCAL}/consensus.fasta" --threads "${THREADS}" \
            2>&1 | tee "${OUT}/medaka_sequence.log"; _M3=${PIPESTATUS[0]}
        set -e
        if [[ ${_M1} -eq 0 && ${_M2} -eq 0 && ${_M3} -eq 0 && -s "${MEDAKA_LOCAL}/consensus.fasta" ]]; then
            mkdir -p "${OUT}/medaka"
            cp "${MEDAKA_LOCAL}/consensus.fasta" "${OUT}/medaka/medaka_consensus.fasta"
            HYBRID_ASM="${OUT}/medaka/medaka_consensus.fasta"
            log "    Medaka pre-polish done → $(basename "$HYBRID_ASM")"
        else
            warn "    Medaka pre-polish failed (m1=${_M1} m2=${_M2} m3=${_M3}) — using raw Flye assembly for Illumina polishing."
        fi
        rm -rf "$MEDAKA_LOCAL"
    else
        log "  Medaka pre-polish skipped (HYBRID_MEDAKA=${HYBRID_MEDAKA:-true}, medaka present: $(command -v medaka &>/dev/null && echo yes || echo no))."
    fi

    POLY_DIR="${OUT}/polypolish"
    PYPOLCA_DIR="${OUT}/pypolca"
    mkdir -p "$POLY_DIR"
    # Note: do NOT pre-create PYPOLCA_DIR — pypolca creates it itself and refuses
    # to run if the directory already exists (even when empty).

    log "  ── Stage 1: Polypolish (Illumina read-pair based) ──"
    log "    Assembly: $(basename $HYBRID_ASM)"
    log "    R1: $(basename $ILLUMINA_R1)"
    log "    R2: $(basename $ILLUMINA_R2)"

    log "  Indexing assembly with BWA..."
    bwa index "$HYBRID_ASM" 2>"${POLY_DIR}/bwa_index.log"

    # Map R1 and R2 separately with -a flag (retain all alignments — Polypolish needs this)
    log "  Mapping R1 (all alignments)..."
    bwa mem \
        -t "${THREADS}" \
        -a "$HYBRID_ASM" \
        "$ILLUMINA_R1" \
        2>"${POLY_DIR}/bwa_R1.log" \
        > "${POLY_DIR}/R1.sam"

    log "  Mapping R2 (all alignments)..."
    bwa mem \
        -t "${THREADS}" \
        -a "$HYBRID_ASM" \
        "$ILLUMINA_R2" \
        2>"${POLY_DIR}/bwa_R2.log" \
        > "${POLY_DIR}/R2.sam"

    # Filter unpaired alignments (Polypolish insert filter step)
    # Some Illumina datasets cause polypolish filter to fail; fall back to unfiltered SAMs.
    log "  Filtering insert sizes..."
    R1_FOR_POLISH="${POLY_DIR}/R1_filtered.sam"
    R2_FOR_POLISH="${POLY_DIR}/R2_filtered.sam"
    if polypolish filter \
        --in1 "${POLY_DIR}/R1.sam" \
        --in2 "${POLY_DIR}/R2.sam" \
        --out1 "$R1_FOR_POLISH" \
        --out2 "$R2_FOR_POLISH" \
        2>"${POLY_DIR}/filter.log"; then
        log "  Insert-size filter succeeded."
    else
        warn "  polypolish filter failed (see ${POLY_DIR}/filter.log)."
        warn "  Falling back to unfiltered SAMs — this is safe but skips outlier removal."
        warn "  Likely cause: non-standard Illumina read-pair naming in this dataset."
        R1_FOR_POLISH="${POLY_DIR}/R1.sam"
        R2_FOR_POLISH="${POLY_DIR}/R2.sam"
    fi

    log "  Polishing with Polypolish --careful..."
    POLYPOLISH_OUT="${POLY_DIR}/polypolish_polished.fasta"
    polypolish polish \
        --careful \
        "$HYBRID_ASM" \
        "$R1_FOR_POLISH" \
        "$R2_FOR_POLISH" \
        > "$POLYPOLISH_OUT" \
        2>"${POLY_DIR}/polish.log"

    [[ -s "$POLYPOLISH_OUT" ]] \
        || die "Polypolish produced empty output — check ${POLY_DIR}/polish.log"

    # Clean up large SAM files
    rm -f "${POLY_DIR}/R1.sam" "${POLY_DIR}/R2.sam" \
          "${POLY_DIR}/R1_filtered.sam" "${POLY_DIR}/R2_filtered.sam"
    log "  Polypolish done."

    # ── Stage 2: Pypolca --careful (pileup-based, catches homopolymer errors) ──
    # Pypolca complements Polypolish: different algorithm (POLCA-style pileup)
    # catches errors that Polypolish's read-pair approach misses, especially in
    # homopolymer runs which are the primary error mode of R10.4.1 HAC.
    log "  ── Stage 2: Pypolca --careful (pileup-based correction) ──"
    # Remove any previous pypolca dir — pypolca's --force uses shutil.rmtree which
    # fails on NFS filesystems; shell rm -rf is reliable.
    rm -rf "$PYPOLCA_DIR"
    # Run pypolca with output on LOCAL SCRATCH (TMPDIR), not NFS.
    # Rationale: pypolca writes a large temp_bwa.sam inside its output directory.
    # Writing several GB of SAM to NFS under concurrent jobs causes bwa to return
    # exit code 1 (failed write). TMPDIR is /local/scratch/tmp on Maestro — fast,
    # node-local, and guaranteed to have space (--gres=disk:10000 was reserved).
    PYPOLCA_LOCAL="${TMPDIR:-/tmp}/pypolca_${BARCODE}_$$"
    rm -rf "$PYPOLCA_LOCAL"
    log "  Running pypolca (local scratch: ${PYPOLCA_LOCAL})..."
    if pypolca run \
        --careful \
        --assembly  "$POLYPOLISH_OUT" \
        --reads1    "$ILLUMINA_R1" \
        --reads2    "$ILLUMINA_R2" \
        --output    "$PYPOLCA_LOCAL" \
        --threads   "${THREADS}" \
        2>&1 | tee "${OUT}/pypolca_run.log"; then

        PYPOLCA_OUT="${PYPOLCA_LOCAL}/pypolca_corrected.fasta"
        if [[ -s "$PYPOLCA_OUT" ]]; then
            cp "$PYPOLCA_OUT" "$POLISHED"
            # Copy full results to NFS for record-keeping (report, VCF, logs).
            cp -r "$PYPOLCA_LOCAL" "$PYPOLCA_DIR" \
                || warn "Could not copy pypolca results to ${PYPOLCA_DIR} — polished.fasta already saved."
            rm -rf "$PYPOLCA_LOCAL"
            log "  Pypolca done."
        else
            warn "  Pypolca ran without error but produced no corrected FASTA."
            warn "  Falling back to Polypolish output as final assembly."
            cp "$POLYPOLISH_OUT" "$POLISHED"
            rm -rf "$PYPOLCA_LOCAL"
        fi

    else
        # Pypolca failed — dump bwa error files to SLURM log for diagnosis.
        warn "  Pypolca FAILED for ${BARCODE} (see ${OUT}/pypolca_run.log)."
        warn "  ── bwa error files from local scratch ──"
        for f in "${PYPOLCA_LOCAL}"/logs/bwa_*.err; do
            [[ -f "$f" ]] && { warn "  === $(basename "$f") ==="; cat "$f" >&2; }
        done
        warn "  ── end bwa errors ──"
        warn "  Falling back to Polypolish output as polished assembly."
        warn "  QV may be slightly lower (~QV50-55 vs ~QV60 with pypolca)."
        cp "$POLYPOLISH_OUT" "$POLISHED"
        rm -rf "$PYPOLCA_LOCAL"
    fi

else
    die "Unknown mode '${MODE}'"
fi

# ---------------------------------------------------------------------------
# Genome reorientation — PhageTerm FIRST, dnaapler as FALLBACK. Done HERE so
# polished.fasta is the SINGLE reoriented reference shared by steps 04/05/06/08
# (and avoids a frame mismatch vs the methylation BED in step 07).
#
#   PhageTerm  : characterises physical termini (DTR, cos, pac, headful, mu-like)
#                from the ONT read-mapping pattern and, on a clear call,
#                repositions the genome to its biological start. Optional
#                (PHAGETERM_BIN or PhageTerm in PATH); non-fatal.
#   dnaapler   : fallback — reorients to the terL gene when PhageTerm is absent
#                or inconclusive. 1 contig → dnaapler phage; N → dnaapler all.
# A termini_report.txt records which method was used and the terminus class.
# ---------------------------------------------------------------------------
REORIENTED=false
TERMINI_REPORT="${OUT}/termini_report.txt"
printf "method\tdetail\n" > "${TERMINI_REPORT}"

# Resolve PhageTerm executable (PHAGETERM_BIN may include args after the binary)
PT_BIN="${PHAGETERM_BIN:-}"
# If not preset and not already in PATH, try the Lmod module (Maestro).
if [[ -z "${PT_BIN}" ]] && ! command -v PhageTerm.py &>/dev/null && ! command -v phageterm &>/dev/null \
   && command -v module &>/dev/null && [[ -n "${PHAGETERM_MODULE:-}" ]]; then
    log "  [modules] module load ${PHAGETERM_MODULE}"
    module load ${PHAGETERM_MODULE} 2>/dev/null || warn "  [modules] could not load ${PHAGETERM_MODULE}"
fi
if [[ -z "${PT_BIN}" ]]; then
    for _c in PhageTerm.py phageterm PhageTerm; do
        command -v "$_c" &>/dev/null && { PT_BIN="$_c"; break; }
    done
fi

if [[ -n "${PT_BIN}" ]] && command -v "${PT_BIN%% *}" &>/dev/null; then
    log "  ── PhageTerm terminus detection / reorientation (${PT_BIN}) ──"
    PT_DIR="${OUT}/phageterm"; rm -rf "$PT_DIR"; mkdir -p "$PT_DIR"
    PT_READS="$ONT_READS"
    if command -v filtlong &>/dev/null; then
        if filtlong --target_bases 10000000 "$ONT_READS" 2>/dev/null | gzip > "${PT_DIR}/pt_reads.fastq.gz" \
             && [[ -s "${PT_DIR}/pt_reads.fastq.gz" ]]; then
            PT_READS="${PT_DIR}/pt_reads.fastq.gz"
        fi
    fi
    set +e
    ( cd "$PT_DIR" && "${PT_BIN%% *}" ${PT_BIN#* } \
        -f "$PT_READS" -r "$POLISHED" -p "${SAMPLE_NAME}" -c "${THREADS}" \
        > phageterm.log 2>&1 )
    _PTEXIT=$?
    set -e
    PT_SEQ=$(ls "${PT_DIR}/${SAMPLE_NAME}"*sequence.fasta 2>/dev/null | head -1 || true)
    [[ -z "$PT_SEQ" ]] && PT_SEQ=$(ls "${PT_DIR}"/*sequence.fasta 2>/dev/null | head -1 || true)
    if [[ ${_PTEXIT} -eq 0 && -s "${PT_SEQ}" \
          && "$(grep -c '^>' "${PT_SEQ}" || true)" -eq "$(grep -c '^>' "$POLISHED" || true)" ]]; then
        PT_STAT=$(ls "${PT_DIR}"/*statistics.csv "${PT_DIR}"/*_report.csv 2>/dev/null | head -1 || true)
        PT_CLASS="detected"
        [[ -n "$PT_STAT" ]] && PT_CLASS=$(grep -iE 'DTR|COS|PAC|Headful|Mu-like|terminal|class' "$PT_STAT" 2>/dev/null | head -1 | tr ',\t' '  ' | tr -s ' ' | cut -c1-80 || echo "detected")
        cp "${PT_SEQ}" "$POLISHED"
        REORIENTED=true
        printf "PhageTerm\t%s\n" "${PT_CLASS:-detected}" > "${TERMINI_REPORT}"
        log "    PhageTerm reoriented genome (termini: ${PT_CLASS:-detected})."
    else
        warn "    PhageTerm gave no usable reorientation (exit=${_PTEXIT}) — trying dnaapler."
        printf "PhageTerm\tinconclusive\n" > "${TERMINI_REPORT}"
    fi
else
    log "  PhageTerm not available (set PHAGETERM_BIN) — using dnaapler."
fi

if [[ "${REORIENTED}" != "true" ]] && command -v dnaapler &>/dev/null; then
    N_CTG=$(grep -c '^>' "$POLISHED" || true)
    REORIENT_DIR="${TMPDIR:-/tmp}/dnaapler_${BARCODE}_${SAMPLE_NAME}_$$"
    rm -rf "$REORIENT_DIR"
    if [[ "${N_CTG}" -eq 1 ]]; then _DMODE="phage"; else _DMODE="all"; fi
    log "  Reorienting with dnaapler (${_DMODE}, ${N_CTG} contig(s))..."
    set +e
    dnaapler "${_DMODE}" -i "$POLISHED" -o "$REORIENT_DIR" -p reoriented \
        -t "${THREADS}" --force 2>&1 | tee "${OUT}/dnaapler.log"
    _DEXIT=${PIPESTATUS[0]}
    set -e
    _DOUT="${REORIENT_DIR}/reoriented_reoriented.fasta"
    if [[ ${_DEXIT} -eq 0 && -s "${_DOUT}" \
          && "$(grep -c '^>' "${_DOUT}" || true)" -eq "${N_CTG}" ]]; then
        cp "${_DOUT}" "$POLISHED"
        REORIENTED=true
        printf "dnaapler\tterL start\n" > "${TERMINI_REPORT}"
        log "  Reorientation applied (dnaapler, terL start)."
    else
        warn "  dnaapler did not reorient (exit=${_DEXIT}) — keeping current orientation."
    fi
    rm -rf "$REORIENT_DIR"
fi

if [[ "${REORIENTED}" != "true" ]]; then
    grep -q . "${TERMINI_REPORT}" 2>/dev/null &&         printf "none\tno reorientation (PhageTerm/dnaapler unavailable or inconclusive)\n" >> "${TERMINI_REPORT}"
    warn "  No reorientation applied — keeping assembly orientation as-is."
fi

# ---------------------------------------------------------------------------
# Normalise contig names → <SAMPLE>_ctgNNN  (collision-proof for Pharokka)
# ---------------------------------------------------------------------------
# Pharokka's pyrodigal-gv names proteins "contig_N_geneI". If the assembler
# left contigs called "contig_1" / "1" / etc., they collide with these internal
# protein IDs and Pharokka crashes (ValueError: Duplicate key 'contig_N').
# Renaming here in step 03 means polished.fasta — and therefore EVERY downstream
# step (05 annotate, 06 modbase, 07 methyl) — shares identical, unique contig
# names. Required for the methylation x annotation bedtools overlap in step 07
# (methylation BED contig names must match the Pharokka GFF names).
log "  Normalising contig names to ${SAMPLE_NAME}_ctgNNN..."
awk -v s="${SAMPLE_NAME}" 'BEGIN{n=0}
    /^>/ { n++; printf(">%s_ctg%03d\n", s, n); next }
    { print }' "$POLISHED" > "${POLISHED}.tmp" \
    && mv "${POLISHED}.tmp" "$POLISHED"
log "  Contig names normalised ($(grep -c '^>' "$POLISHED") contig(s))."

# Report polished assembly stats
n_contigs=$(grep -c '^>' "$POLISHED" || echo 0)
total_len=$(grep -v '^>' "$POLISHED" | tr -d '\n' | wc -c)
log "  Polished: ${n_contigs} contig(s), ${total_len} bp total"
log "  Output: ${POLISHED}"
log "  Step 03 done."
