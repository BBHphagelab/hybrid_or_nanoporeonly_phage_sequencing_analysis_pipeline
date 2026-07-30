#!/usr/bin/env bash
# =============================================================================
# fetch_host_refs.sh
# -----------------------------------------------------------------------------
# Downloads the host *Pseudomonas aeruginosa* genome(s) used to decontaminate
# the ONT reads, and concatenates them into a single host_combined.fasta.
#
# RUN THIS ON A MACHINE WITH INTERNET (your laptop or a login node) — Maestro
# compute nodes have no outbound network. Then transfer host_combined.fasta to
# the cluster and pass it to assemble_phages_hostdepleted.sh via --host-ref.
#
# Default references:
#   PAK   CP020659.1   P. aeruginosa PAK, complete chromosome (the host of
#                      PAK_P1 / PAK_P4 and their recombinants — bc01/02/04/05).
#   PAO1  NC_002516.2  P. aeruginosa PAO1 reference — universal backstop. The
#                      P. aeruginosa core genome is >99% identical across
#                      strains, so PAO1 mops up host reads from the CHA lineage
#                      even without the exact CHA assembly.
#
# BETTER, if you have it: add your in-house *P. aeruginosa CHA* genome. Drop its
# FASTA path in EXTRA_HOST_FASTA below (or pass it as $1) — it will be appended.
# An exact public "P. aeruginosa CHA" assembly is not unambiguous on NCBI, and
# the actual host you propagated the phages on is the most accurate decontaminant.
#
# USAGE
#   bash fetch_host_refs.sh [extra_host.fasta ...]
# Output: ./host_refs/host_combined.fasta
# =============================================================================
set -euo pipefail

OUTDIR="./host_refs"
mkdir -p "$OUTDIR"
COMBINED="${OUTDIR}/host_combined.fasta"
: > "$COMBINED"

# accession -> label
ACCESSIONS=( "CP020659.1:PAK" "NC_002516.2:PAO1" )

fetch_efetch() {  # $1 accession  $2 outfile
    local acc="$1" out="$2"
    local url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text"
    if command -v curl &>/dev/null; then
        curl -fsSL "$url" -o "$out"
    elif command -v wget &>/dev/null; then
        wget -q -O "$out" "$url"
    else
        echo "[ERROR] need curl or wget (or the NCBI 'datasets' CLI)." >&2
        return 1
    fi
    [[ -s "$out" ]] && head -1 "$out" | grep -q '^>'
}

echo "[fetch] Downloading host references from NCBI..."
for entry in "${ACCESSIONS[@]}"; do
    acc="${entry%%:*}"; label="${entry##*:}"
    out="${OUTDIR}/${label}_${acc}.fasta"
    echo "  - ${label}  (${acc})"
    if [[ -s "$out" ]]; then
        echo "    already present, reusing ${out}"
    else
        fetch_efetch "$acc" "$out" || { echo "[ERROR] failed to fetch ${acc}" >&2; exit 1; }
    fi
    # prefix headers so host contigs are unambiguous in the combined file
    awk -v p="host_${label}" '/^>/{sub(/^>/,">"p"|")}1' "$out" >> "$COMBINED"
done

# Optional in-house CHA (or any extra) host genome(s)
for extra in "$@"; do
    [[ -f "$extra" ]] || { echo "[WARN] extra host FASTA not found: $extra" >&2; continue; }
    echo "  - extra host genome: ${extra}"
    awk -v p="host_extra" '/^>/{sub(/^>/,">"p"|")}1' "$extra" >> "$COMBINED"
done

n=$(grep -c '^>' "$COMBINED" || echo 0)
bp=$(grep -v '^>' "$COMBINED" | tr -d '\n' | wc -c)
echo ""
echo "[done] ${COMBINED}"
echo "       ${n} sequence(s), ${bp} bp total"
echo ""
echo "Transfer it to the cluster, e.g.:"
echo "   scp ${COMBINED} <user>@maestro:/pasteur/.../20260520_ONT_phages/refs/"
echo "Then run on the cluster:"
echo "   bash assemble_phages_hostdepleted.sh --host-ref /pasteur/.../refs/host_combined.fasta \\"
echo "        --barcodes barcode03,barcode04,barcode05 --force"
