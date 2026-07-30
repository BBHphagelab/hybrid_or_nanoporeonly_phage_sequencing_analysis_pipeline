#!/usr/bin/env bash
# =============================================================================
# assemble_phages_hostdepleted.sh
# -----------------------------------------------------------------------------
# Standalone replacement for step 02 that fixes the contaminated-assembly
# problem (ONT reads carry host *Pseudomonas* gDNA -> de novo assembly produces
# hundreds of contigs / multi-Mb totals).
#
# Strategy (per barcode):
#   1. HOST DEPLETION  — competitive minimap2 (host vs phage refs): a read is
#      removed only if it aligns better to the host than to any phage reference
#      (and covers >= HOST_COV_FRAC of its length on host). Keeps recombinant
#      junctions + novel/engineered inserts; protects against host prophage.
#   2. ASSEMBLY        — Autocycler (default): multi-assembler consensus of
#      several independent read subsamples -> the most accurate single-contig
#      circular genome with reliable ends (best for downstream PhageTerm).
#      Autocycler does its own subsampling, so very high coverage is handled
#      internally. Fallback: Flye --nano-hq on reads downsampled to ~100x with
#      Filtlong (used if Autocycler / its sub-assemblers are unavailable or fail).
#   3. CONTIG SELECTION— drop residual host debris (coverage/length/circularity
#      for Flye; length window for Autocycler). Never empties the assembly.
#   4. FINISHING       — DTR / circular-overhang detection (finish_assembly.py).
#      Reorientation + contig renaming are LEFT to step 03 (unchanged), so
#      PhageTerm can later run on a non-reoriented assembly if desired.
#
# Output is written where the pipeline expects it:
#     ${RESULTS_DIR}/02_assembly/<barcode>_<sample>/assembly.fasta
# Then resume the normal pipeline:  bash submit_all.sh   (step 02 self-skips).
#
# USAGE
#   bash assemble_phages_hostdepleted.sh --host-ref host_combined.fasta \
#        [--assembler autocycler|flye] [--phage-ref phage_refs.fasta] \
#        [--barcodes barcode03,barcode04,barcode05] [--force]
#
#   --host-ref   REQUIRED. Combined host genome FASTA (PAK + CHA).
#   --assembler  autocycler (default) or flye.
#   --phage-ref  Optional phage reference FASTA; else built from EXTERNAL_REF_1/2/3.
#   --barcodes   Comma list; default = all barcodes in the sample sheet.
#   --force      Re-run even if assembly.fasta exists.
#
# Tunables (env):
#   HOST_COV_FRAC=0.5     COV_FRAC=0.10     MIN_PHAGE_LEN=20000  MAX_PHAGE_LEN=200000
#   AUTOCYCLER_ASSEMBLERS="flye raven miniasm plassembler"  (auto-filtered to installed)
#   FILTLONG_TARGET_COV=100   (x coverage Filtlong keeps for the Flye fallback)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# ---------------------------------------------------------------------------
# Defaults / argument parsing
# ---------------------------------------------------------------------------
HOST_REF=""
PHAGE_REF=""
BARCODES_ARG=""
ASSEMBLER="autocycler"
HOST_COV_FRAC="${HOST_COV_FRAC:-0.5}"
COV_FRAC="${COV_FRAC:-0.10}"
MIN_PHAGE_LEN="${MIN_PHAGE_LEN:-20000}"
MAX_PHAGE_LEN="${MAX_PHAGE_LEN:-200000}"
FILTLONG_TARGET_COV="${FILTLONG_TARGET_COV:-100}"
AUTOCYCLER_ASSEMBLERS="${AUTOCYCLER_ASSEMBLERS:-flye raven miniasm plassembler}"
# PREFERRED: a dedicated micromamba env with autocycler + the assemblers (fast,
# reproducible, no slow Lmod conda-modules). If AUTOCYCLER_ENV is set, it is
# activated for the assembly tools and the Lmod modules below are NOT used.
# Create it with (see README):
#   micromamba create -y -n autocycler -c conda-forge -c bioconda \
#       autocycler flye raven-assembler miniasm minipolish plassembler \
#       minimap2 samtools seqkit filtlong any2fasta
AUTOCYCLER_ENV="${AUTOCYCLER_ENV:-}"
# FALLBACK (no env): Lmod modules. MINIMAL set, and ORDER MATTERS — the Autocycler
# modulefile requires an assembler module to be loaded BEFORE it, so Autocycler is
# LAST. We avoid the heavy conda-based modules (canu, metaMDBG, Minipolish,
# Unicycler, dnaapler, plassembler, wtdbg2) that hung the job. The assembler list
# above is auto-filtered to whatever ends up in PATH (here: flye + raven).
AUTOCYCLER_MODULES="${AUTOCYCLER_MODULES:-Flye raven minimap2 samtools/1.21 Autocycler/0.5.2}"
FORCE="${FORCE:-false}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --host-ref)   HOST_REF="$2"; shift 2 ;;
        --phage-ref)  PHAGE_REF="$2"; shift 2 ;;
        --assembler)  ASSEMBLER="$2"; shift 2 ;;
        --barcodes)   BARCODES_ARG="$2"; shift 2 ;;
        --force)      FORCE="true"; shift ;;
        -h|--help)    sed -n '2,55p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *)            die "Unknown argument: $1 (use --help)" ;;
    esac
done

[[ "$ASSEMBLER" == "autocycler" || "$ASSEMBLER" == "flye" ]] \
    || die "--assembler must be 'autocycler' or 'flye' (got '${ASSEMBLER}')"
[[ -n "$HOST_REF" ]] || die "--host-ref is required (combined host genome FASTA)."
check_file "$HOST_REF"

activate_env
load_samtools_cluster
require_tool minimap2
require_tool samtools

# Load Autocycler + sub-assembler Lmod modules (Maestro). flye stays from the env.
# Module binaries are prepended to PATH so `autocycler`, raven, miniasm, plassembler
# become available inside the job. Non-fatal: missing modules just warn.
if [[ "$ASSEMBLER" == "autocycler" ]]; then
    if [[ -n "$AUTOCYCLER_ENV" ]]; then
        log "[env] activating micromamba env '${AUTOCYCLER_ENV}' for assembly tools"
        micromamba activate "${AUTOCYCLER_ENV}" \
            || warn "[env] could not activate '${AUTOCYCLER_ENV}' — check the env name."
    elif command -v module &>/dev/null; then
        log "[modules] module load ${AUTOCYCLER_MODULES}"
        # ORDER MATTERS: assemblers first, Autocycler last (modulefile prerequisite).
        module load ${AUTOCYCLER_MODULES} \
            || warn "[modules] 'module load' returned an error — check module names/order."
    else
        warn "[autocycler] no AUTOCYCLER_ENV and no Lmod — autocycler must already be in PATH."
    fi
fi

# Convert a genome-size string like "100k" / "1.5m" to integer bp
to_bp() {
    local v="${1,,}"
    case "$v" in
        *k) awk -v x="${v%k}" 'BEGIN{printf "%d", x*1000}' ;;
        *m) awk -v x="${v%m}" 'BEGIN{printf "%d", x*1000000}' ;;
        *)  awk -v x="$v"     'BEGIN{printf "%d", x}' ;;
    esac
}
GS_NUM=$(to_bp "${GENOME_SIZE_ESTIMATE}")

# Decide effective assembler now (fall back to flye if autocycler not installed)
if [[ "$ASSEMBLER" == "autocycler" ]] && ! command -v autocycler &>/dev/null; then
    warn "[assembler] 'autocycler' not found in PATH — falling back to Flye + Filtlong."
    ASSEMBLER="flye"
fi
if [[ "$ASSEMBLER" == "flye" ]]; then
    require_tool flye
fi
log "[assembler] Using: ${ASSEMBLER}  (genome size ~${GS_NUM} bp)"

SAMPLE_SHEET="${RESULTS_DIR}/sample_sheet.tsv"
check_file "$SAMPLE_SHEET"

# ---------------------------------------------------------------------------
# Build the phage reference FASTA (union of external refs) for competitive
# host-vs-phage read assignment, if not supplied.
# ---------------------------------------------------------------------------
WORK_REF_DIR="${RESULTS_DIR}/02b_hostdepleted/_refs"
mkdir -p "$WORK_REF_DIR"

gbk_to_fasta() {  # $1 in.gbk  $2 out.fasta (append)  -> 0 on success
    local in="$1" out="$2"
    [[ -f "$in" ]] || return 1
    if command -v any2fasta &>/dev/null; then any2fasta "$in" >> "$out" 2>/dev/null && return 0; fi
    if command -v seqkit   &>/dev/null; then seqkit seq "$in" >> "$out" 2>/dev/null && return 0; fi
    python3 - "$in" "$out" <<'PY' 2>/dev/null && return 0
import sys, re
inp, outp = sys.argv[1], sys.argv[2]
recs=[]; name=None; seq=[]; ino=False
for line in open(inp):
    if line.startswith("LOCUS"): name=line.split()[1]; seq=[]; ino=False
    elif line.startswith("ORIGIN"): ino=True
    elif line.startswith("//"):
        if name and seq: recs.append((name,"".join(seq)))
        name=None; ino=False
    elif ino: seq.append(re.sub(r"[^A-Za-z]","",line))
if not recs: sys.exit(1)
with open(outp,"a") as fh:
    for n,s in recs: fh.write(f">{n}\n{s.upper()}\n")
PY
    return 1
}

if [[ -z "$PHAGE_REF" ]]; then
    CANDIDATE="${WORK_REF_DIR}/phage_refs_union.fasta"
    : > "$CANDIDATE"; built=0
    for r in "${EXTERNAL_REF_1:-}" "${EXTERNAL_REF_2:-}" "${EXTERNAL_REF_3:-}"; do
        [[ -n "$r" && -f "$r" ]] || continue
        gbk_to_fasta "$r" "$CANDIDATE" && built=$((built+1))
    done
    if [[ "$built" -gt 0 && -s "$CANDIDATE" ]]; then
        PHAGE_REF="$CANDIDATE"
        log "[refs] Built phage reference union from ${built} external ref(s)."
    else
        warn "[refs] Could not build phage reference FASTA — using HOST-ONLY depletion."
        PHAGE_REF=""
    fi
else
    check_file "$PHAGE_REF"
fi

# ---------------------------------------------------------------------------
# Resolve which Autocycler sub-assemblers are actually installed
# ---------------------------------------------------------------------------
autocycler_binary_for() {
    case "$1" in
        flye) echo flye;; raven) echo raven;; miniasm) echo miniasm;;
        plassembler) echo plassembler;; canu) echo canu;; necat) echo necat;;
        nextdenovo) echo nextDenovo;; metamdbg) echo metaMDBG;; *) echo "$1";;
    esac
}
ACTIVE_ASSEMBLERS=""
if [[ "$ASSEMBLER" == "autocycler" ]]; then
    for a in ${AUTOCYCLER_ASSEMBLERS}; do
        command -v "$(autocycler_binary_for "$a")" &>/dev/null && ACTIVE_ASSEMBLERS+="${a} "
    done
    ACTIVE_ASSEMBLERS="$(echo "$ACTIVE_ASSEMBLERS" | xargs || true)"
    if [[ -z "$ACTIVE_ASSEMBLERS" ]]; then
        warn "[autocycler] none of (${AUTOCYCLER_ASSEMBLERS}) are installed — falling back to Flye."
        ASSEMBLER="flye"; require_tool flye
    else
        log "[autocycler] sub-assemblers available: ${ACTIVE_ASSEMBLERS}"
    fi
fi

# ---------------------------------------------------------------------------
# run_autocycler <reads> <outdir> <gs>  -> echoes consensus fasta path; rc!=0 on fail
# ---------------------------------------------------------------------------
run_autocycler() {
    local reads="$1" outdir="$2" gs="$3"
    rm -rf "$outdir"; mkdir -p "$outdir/assemblies"
    autocycler subsample --reads "$reads" --out_dir "$outdir/subsampled" \
        --genome_size "$gs" >>"$outdir/autocycler.log" 2>&1 || return 1
    local got=0 f i
    for f in "$outdir"/subsampled/sample_*.fastq; do
        [[ -e "$f" ]] || continue
        i="$(basename "$f" .fastq)"   # sample_01
        for a in $ACTIVE_ASSEMBLERS; do
            if autocycler helper "$a" --reads "$f" \
                   --out_prefix "$outdir/assemblies/${a}_${i}" \
                   --threads "${THREADS}" --genome_size "$gs" \
                   >>"$outdir/autocycler.log" 2>&1; then
                got=$((got+1))
            else
                warn "    autocycler helper ${a} failed on ${i} (continuing)."
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
# run_flye <reads> <outdir> <barcode>  -> sets globals FLYE_DRAFT, FLYE_INFO
# Downsamples to ~FILTLONG_TARGET_COV x with Filtlong first (if available).
# ---------------------------------------------------------------------------
run_flye() {
    local reads="$1" outdir="$2" barcode="$3"
    rm -rf "$outdir"; mkdir -p "$outdir"
    local asm_reads="$reads"
    if command -v filtlong &>/dev/null; then
        local target=$(( FILTLONG_TARGET_COV * GS_NUM ))
        log "    Filtlong downsample to ~${FILTLONG_TARGET_COV}x (target ${target} bp)..."
        filtlong --target_bases "$target" "$reads" 2>"$outdir/filtlong.log" \
            | gzip > "$outdir/downsampled.fastq.gz"
        [[ -s "$outdir/downsampled.fastq.gz" ]] && asm_reads="$outdir/downsampled.fastq.gz" \
            || warn "    Filtlong produced no output — using full read set."
    else
        warn "    filtlong not found — assembling without downsampling (high coverage may hurt Flye)."
    fi
    local meta_arg=()
    for m in ${FLYE_META_BARCODES:-}; do [[ "$m" == "$barcode" ]] && meta_arg=(--meta); done
    flye --nano-hq "$asm_reads" --genome-size "${GENOME_SIZE_ESTIMATE}" \
        "${meta_arg[@]}" --out-dir "$outdir/flye" --threads "${THREADS}" \
        2>&1 | tee "$outdir/flye.log"
    FLYE_DRAFT="$outdir/flye/assembly.fasta"
    FLYE_INFO="$outdir/flye/assembly_info.txt"
    [[ -s "$FLYE_DRAFT" ]]
}

# ---------------------------------------------------------------------------
# Resolve barcode list + summary header
# ---------------------------------------------------------------------------
if [[ -n "$BARCODES_ARG" ]]; then
    IFS=',' read -r -a BARCODES <<< "$BARCODES_ARG"
else
    mapfile -t BARCODES < <(list_barcodes)
fi
[[ "${#BARCODES[@]}" -gt 0 ]] || die "No barcodes to process."

SUMMARY="${RESULTS_DIR}/02b_hostdepleted/decontamination_summary.tsv"
mkdir -p "$(dirname "$SUMMARY")"
[[ -f "$SUMMARY" ]] || printf "barcode\tsample\tmode\tassembler\treads_in\treads_kept\tpct_kept\tcontigs_raw\tcontigs_kept\ttotal_bp_kept\tcircular\n" > "$SUMMARY"

# ===========================================================================
# Per-barcode processing
# ===========================================================================
for BARCODE in "${BARCODES[@]}"; do
    SAMPLE_NAME=$(get_sample_field "$BARCODE" "sample_name")
    MODE=$(get_sample_field "$BARCODE" "mode")
    step_banner "02b host-depleted assembly  [${BARCODE} / ${SAMPLE_NAME} / ${ASSEMBLER}]"

    ONT_READS="${RESULTS_DIR}/01_qc/${BARCODE}_${SAMPLE_NAME}/filtered.fastq.gz"
    OUT="${RESULTS_DIR}/02_assembly/${BARCODE}_${SAMPLE_NAME}"
    WORK="${RESULTS_DIR}/02b_hostdepleted/${BARCODE}_${SAMPLE_NAME}"
    ASSEMBLY_FASTA="${OUT}/assembly.fasta"
    mkdir -p "$OUT" "$WORK"

    if [[ ! -f "$ONT_READS" ]]; then
        warn "  ${BARCODE}: filtered ONT reads not found (${ONT_READS}) — run step 01 first. Skipping."
        continue
    fi
    if [[ -f "$ASSEMBLY_FASTA" && "$FORCE" != "true" ]]; then
        log "  assembly.fasta exists — skipping (use --force to rerun)."
        continue
    fi
    [[ -f "$ASSEMBLY_FASTA" ]] && cp -f "$ASSEMBLY_FASTA" "${OUT}/assembly.contaminated.bak.fasta" 2>/dev/null || true

    # -----------------------------------------------------------------------
    # 1. HOST DEPLETION (competitive host-vs-phage)
    # -----------------------------------------------------------------------
    log "  ── Host depletion ──"
    HOST_PAF="${WORK}/reads_vs_host.paf"
    PHAGE_PAF="${WORK}/reads_vs_phage.paf"
    HOST_IDS="${WORK}/host_read_ids.txt"
    CLEAN_READS="${WORK}/clean.fastq.gz"

    minimap2 -x map-ont -t "${THREADS}" --secondary=no "$HOST_REF" "$ONT_READS" \
        > "$HOST_PAF" 2>"${WORK}/minimap2_host.log"
    if [[ -n "$PHAGE_REF" ]]; then
        minimap2 -x map-ont -t "${THREADS}" --secondary=no "$PHAGE_REF" "$ONT_READS" \
            > "$PHAGE_PAF" 2>"${WORK}/minimap2_phage.log"
    else
        : > "$PHAGE_PAF"
    fi

    python3 - "$HOST_PAF" "$PHAGE_PAF" "$HOST_IDS" "$HOST_COV_FRAC" <<'PY'
import sys
host_paf, phage_paf, out_ids, frac = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
def parse(path):
    bs={}; bf={}
    try: f=open(path)
    except FileNotFoundError: return bs,bf
    for line in f:
        c=line.rstrip("\n").split("\t")
        if len(c)<11: continue
        q=c[0]
        try: qlen=int(c[1]); qs=int(c[2]); qe=int(c[3]); m=int(c[9])
        except ValueError: continue
        fr=(qe-qs)/qlen if qlen else 0.0
        if m>bs.get(q,-1): bs[q]=m
        if fr>bf.get(q,-1.0): bf[q]=fr
    f.close(); return bs,bf
hs,hf=parse(host_paf); ps,_=parse(phage_paf)
host=[q for q,sc in hs.items() if hf.get(q,0.0)>=frac and sc>ps.get(q,0)]
open(out_ids,"w").write("\n".join(host)+("\n" if host else ""))
print(f"[host-depletion] host-mapped={len(hs)} phage-mapped={len(ps)} flagged-host={len(host)}")
PY

    n_host=$(wc -l < "$HOST_IDS" | tr -d ' ')
    if command -v seqkit &>/dev/null; then
        if [[ "$n_host" -gt 0 ]]; then
            seqkit grep -v -f "$HOST_IDS" "$ONT_READS" -o "$CLEAN_READS" 2>"${WORK}/seqkit.log"
        else cp -f "$ONT_READS" "$CLEAN_READS"; fi
    else
        python3 - "$ONT_READS" "$HOST_IDS" "$CLEAN_READS" <<'PY'
import sys, gzip
reads, ids_file, out = sys.argv[1], sys.argv[2], sys.argv[3]
host=set(l.strip() for l in open(ids_file) if l.strip())
op=gzip.open if reads.endswith(".gz") else open
ow=gzip.open if out.endswith(".gz") else open
with op(reads,"rt") as f, ow(out,"wt") as g:
    while True:
        h=f.readline()
        if not h: break
        s=f.readline(); p=f.readline(); q=f.readline()
        if h[1:].split()[0] not in host: g.write(h); g.write(s); g.write(p); g.write(q)
PY
    fi

    reads_in=$(( $(zcat -f "$ONT_READS" | wc -l) / 4 ))
    reads_kept=$(( $(zcat -f "$CLEAN_READS" | wc -l) / 4 ))
    pct_kept=$(awk -v a="$reads_kept" -v b="$reads_in" 'BEGIN{printf (b>0)?"%.1f":"0", 100*a/b}')
    log "    Reads: ${reads_in} -> ${reads_kept} kept (${pct_kept}%); ${n_host} host removed."
    [[ "$reads_kept" -lt 100 ]] && warn "    Very few reads survived (${reads_kept}) — check HOST_COV_FRAC / host ref."

    # -----------------------------------------------------------------------
    # 2. ASSEMBLY (Autocycler, fallback Flye+Filtlong)
    # -----------------------------------------------------------------------
    DRAFT=""; DRAFT_INFO=""
    if [[ "$ASSEMBLER" == "autocycler" ]]; then
        log "  ── Autocycler consensus (${ACTIVE_ASSEMBLERS}) ──"
        if CONS=$(run_autocycler "$CLEAN_READS" "${WORK}/autocycler" "$GS_NUM"); then
            DRAFT="$CONS"
            log "    Autocycler consensus: $(basename "$DRAFT")"
        else
            warn "    Autocycler failed (see ${WORK}/autocycler/autocycler.log) — falling back to Flye for ${BARCODE}."
        fi
    fi
    if [[ -z "$DRAFT" ]]; then
        log "  ── Flye --nano-hq (Filtlong ${FILTLONG_TARGET_COV}x) ──"
        require_tool flye
        if run_flye "$CLEAN_READS" "${WORK}/flye_run" "$BARCODE"; then
            DRAFT="$FLYE_DRAFT"; DRAFT_INFO="$FLYE_INFO"
        else
            warn "  Flye produced no assembly for ${BARCODE} — see ${WORK}/flye_run/flye.log. Skipping."
            continue
        fi
    fi
    n_raw=$(grep -c '^>' "$DRAFT" || echo 0)
    log "    Draft: ${n_raw} contig(s)."

    # -----------------------------------------------------------------------
    # 3. CONTIG SELECTION (drop residual host debris; never empty)
    # -----------------------------------------------------------------------
    SELECTED="${WORK}/selected.fasta"
    if [[ -n "$DRAFT_INFO" && -f "$DRAFT_INFO" ]]; then
        # Flye: coverage + length + circular
        python3 - "$DRAFT_INFO" "$DRAFT" "$SELECTED" "$COV_FRAC" "$MIN_PHAGE_LEN" "$MAX_PHAGE_LEN" <<'PY'
import sys
info,fasta,out,cov_frac,min_len,max_len=sys.argv[1:7]
cov_frac=float(cov_frac); min_len=int(min_len); max_len=int(max_len)
rows=[]
with open(info) as f:
    f.readline()
    for line in f:
        p=line.split("\t")
        if len(p)<4: continue
        try: length=int(p[1]); cov=float(p[2])
        except ValueError: continue
        circ=p[3].strip().lower() in ("y","yes","+","true","1")
        rows.append((p[0],length,cov,circ))
keep=set()
if rows:
    mc=max(c for _,_,c,_ in rows); thr=cov_frac*mc
    for n,l,c,ci in rows:
        if ci and l>=min_len: keep.add(n)
        elif c>=thr and min_len<=l<=max_len: keep.add(n)
    if not keep: keep.add(max(rows,key=lambda r:(r[2],r[1]))[0])
out_recs=[]; name=None; seq=[]
def flush():
    if name is not None and name.split()[0] in keep: out_recs.append((name,"".join(seq)))
with open(fasta) as f:
    for line in f:
        if line.startswith(">"): flush(); name=line[1:].strip(); seq=[]
        else: seq.append(line.strip())
    flush()
with open(out,"w") as g:
    for h,s in out_recs: g.write(f">{h}\n{s}\n")
print(f"[select-flye] contigs={len(rows)} kept={len(out_recs)}")
PY
    else
        # Autocycler: length window only (host already depleted, consensus is clean)
        python3 - "$DRAFT" "$SELECTED" "$MIN_PHAGE_LEN" "$MAX_PHAGE_LEN" <<'PY'
import sys
fasta,out,min_len,max_len=sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4])
recs=[]; name=None; seq=[]
def flush():
    if name is not None: recs.append((name,"".join(seq)))
with open(fasta) as f:
    for line in f:
        if line.startswith(">"): flush(); name=line[1:].strip(); seq=[]
        else: seq.append(line.strip())
    flush()
keep=[(n,s) for n,s in recs if min_len<=len(s)<=max_len]
if not keep and recs: keep=[max(recs,key=lambda r:len(r[1]))]
with open(out,"w") as g:
    for n,s in keep: g.write(f">{n}\n{s}\n")
print(f"[select-autocycler] contigs={len(recs)} kept={len(keep)}")
PY
    fi

    cp "$SELECTED" "$ASSEMBLY_FASTA"
    # Drop any stale indexes from a previous assembly so downstream tools
    # (medaka/pysam, bwa, minimap2) never read an index that no longer matches.
    rm -f "${ASSEMBLY_FASTA}.fai" "${ASSEMBLY_FASTA}.mmi" \
          "${ASSEMBLY_FASTA}.amb" "${ASSEMBLY_FASTA}.ann" "${ASSEMBLY_FASTA}.bwt" \
          "${ASSEMBLY_FASTA}.pac" "${ASSEMBLY_FASTA}.sa"
    [[ -n "$DRAFT_INFO" && -f "$DRAFT_INFO" ]] && cp "$DRAFT_INFO" "${OUT}/assembly_info.txt" || true
    cp "$ASSEMBLY_FASTA" "${OUT}/assembly_prefilter.fasta"
    # expose an assembly graph for Bandage if present
    for g in "${WORK}/autocycler/autocycler_out/consensus_assembly.gfa" \
             "${WORK}/flye_run/flye/assembly_graph.gfa"; do
        [[ -f "$g" ]] && { cp "$g" "${OUT}/for_bandage.gfa"; break; }
    done

    n_kept=$(grep -c '^>' "$ASSEMBLY_FASTA" || echo 0)
    bp_kept=$(grep -v '^>' "$ASSEMBLY_FASTA" | tr -d '\n' | wc -c)
    log "    Selected: ${n_kept} contig(s), ${bp_kept} bp."

    # -----------------------------------------------------------------------
    # 4. FINISHING (DTR / circular) — reorientation left to step 03
    # -----------------------------------------------------------------------
    if [[ -f "${SCRIPT_DIR}/finish_assembly.py" ]]; then
        FIN_TRIM_ARG=(); [[ "${FINISH_TRIM_OVERHANG:-false}" == "true" ]] && FIN_TRIM_ARG=(--trim)
        if python3 "${SCRIPT_DIR}/finish_assembly.py" \
            --in "$ASSEMBLY_FASTA" --out "${ASSEMBLY_FASTA}.fin" \
            --report "${OUT}/finishing_report.tsv" \
            --min-contig-len "${MIN_CONTIG_LEN:-0}" \
            --circ-max "${FINISH_CIRC_MAX:-1000}" "${FIN_TRIM_ARG[@]}" \
            2>&1 | tee "${OUT}/finishing.log"; then
            [[ -s "${ASSEMBLY_FASTA}.fin" ]] && mv "${ASSEMBLY_FASTA}.fin" "$ASSEMBLY_FASTA"
        else
            warn "    Finishing helper failed — keeping assembly unchanged."; rm -f "${ASSEMBLY_FASTA}.fin"
        fi
        n_kept=$(grep -c '^>' "$ASSEMBLY_FASTA" || echo 0)
    fi

    n_circ=$(grep '^>' "$ASSEMBLY_FASTA" | grep -ic 'circular\|dtr' || true)
    n_circ="${n_circ:-0}"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$BARCODE" "$SAMPLE_NAME" "$MODE" "$ASSEMBLER" "$reads_in" "$reads_kept" "$pct_kept" \
        "$n_raw" "$n_kept" "$bp_kept" "$n_circ" >> "$SUMMARY"
    log "  Done: ${BARCODE} ${SAMPLE_NAME} -> ${ASSEMBLY_FASTA} (${n_kept} contig(s), ${bp_kept} bp)"
done

echo ""
log "All requested barcodes processed."
log "Decontamination/assembly summary: ${SUMMARY}"
echo ""
echo "NEXT STEPS:"
echo "  1. Remove stale downstream results (every step self-skips if output exists):"
echo "        bash clean_for_reassembly.sh            # dry-run"
echo "        bash clean_for_reassembly.sh --yes      # delete"
echo "  2. Relaunch:  bash submit_all.sh              # 02 self-skips, keeps new asm"
echo ""
echo "Step 03 then: Medaka (ONT) + Polypolish + Pypolca (Illumina) for hybrids,"
echo "Medaka for ONT-only, reorients (dnaapler) + renames -> polished.fasta."
echo "# end of script"
