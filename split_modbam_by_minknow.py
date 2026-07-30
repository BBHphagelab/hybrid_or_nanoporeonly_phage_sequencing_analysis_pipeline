#!/usr/bin/env python3
"""
split_modbam_by_minknow.py
--------------------------
Split mod_calls.bam (Dorado-basecalled, unaligned) into per-barcode BAMs using
the barcode assignments from MinKNOW's sequencing_summary.txt.

This bypasses Dorado demux (which only classifies ~3.6% of reads on this run)
and instead uses MinKNOW's barcode classification (~94.1%) to assign reads.

The output BAMs are named to match the pattern expected by step 07:
    *_barcode01.bam, *_barcode02.bam, ..., *_barcode07.bam

Usage (on cluster, with nanopore_bact env active):
    python3 split_modbam_by_minknow.py \\
        --mod-bam   /path/to/mod_calls.bam \\
        --summary   /path/to/sequencing_summary_FBF13730_11d92855_f48f5fca.txt \\
        --output-dir /path/to/modbase/demuxed_minknow/ \\
        [--prefix minknow_SQK-RBK114-24] \\
        [--pass-only]   # default: include all reads (pass+fail)

Requirements:
    pip install pysam  (or: micromamba install -c bioconda pysam)
"""

import os
import sys
import argparse
from collections import defaultdict

try:
    import pysam
except ImportError:
    sys.exit("ERROR: pysam not found. Install with: micromamba install -c bioconda pysam")


# ---------------------------------------------------------------------------
# Load read_id -> barcode map from sequencing_summary.txt
# ---------------------------------------------------------------------------
def load_barcode_map(summary_path: str, pass_only: bool) -> dict:
    """
    Returns {read_id: barcode_arrangement} for all reads in the summary.
    barcode_arrangement is e.g. "barcode01", "barcode02", ..., "unclassified".
    """
    barcode_map = {}
    skipped_fail = 0

    with open(summary_path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        try:
            ri = header.index('read_id')
            bi = header.index('barcode_arrangement')
            pi = header.index('passes_filtering')
        except ValueError as e:
            sys.exit(f"ERROR: sequencing_summary.txt missing expected column: {e}")

        for lineno, line in enumerate(fh, start=2):
            parts = line.rstrip('\n').split('\t')
            if len(parts) <= max(ri, bi, pi):
                continue
            if pass_only and parts[pi].upper() != 'TRUE':
                skipped_fail += 1
                continue
            barcode_map[parts[ri]] = parts[bi]

    classified = sum(1 for v in barcode_map.values() if v != 'unclassified')
    print(f"  Loaded {len(barcode_map):,} reads from summary"
          f" ({classified:,} with barcode, "
          f"{len(barcode_map) - classified:,} unclassified"
          + (f", {skipped_fail:,} fail reads skipped" if pass_only else "") + ")")
    return barcode_map


# ---------------------------------------------------------------------------
# Split mod_calls.bam by barcode
# ---------------------------------------------------------------------------
def split_bam(mod_bam_path: str, barcode_map: dict,
              output_dir: str, prefix: str,
              target_barcodes: list) -> None:
    """
    Streams through mod_calls.bam once.
    For each read whose query_name is in barcode_map and maps to a target
    barcode, writes it to the appropriate per-barcode output BAM.
    """
    os.makedirs(output_dir, exist_ok=True)

    out_files: dict[str, pysam.AlignmentFile] = {}
    counts: dict[str, int] = defaultdict(int)
    total = 0
    classified = 0
    no_summary = 0

    with pysam.AlignmentFile(mod_bam_path, "rb", check_sq=False) as bam_in:
        header = bam_in.header.to_dict()

        for read in bam_in:
            total += 1
            read_id = read.query_name
            barcode = barcode_map.get(read_id)

            if barcode is None:
                no_summary += 1
                continue
            if barcode == 'unclassified' or barcode not in target_barcodes:
                continue

            classified += 1
            if barcode not in out_files:
                out_path = os.path.join(output_dir, f"{prefix}_{barcode}.bam")
                out_files[barcode] = pysam.AlignmentFile(out_path, "wb",
                                                          header=header)
                print(f"  Opened output: {os.path.basename(out_path)}", flush=True)

            out_files[barcode].write(read)
            counts[barcode] += 1

            if total % 200_000 == 0:
                pct = 100.0 * classified / total
                print(f"  ... {total:,} reads processed, "
                      f"{classified:,} classified ({pct:.1f}%)", flush=True)

    for f in out_files.values():
        f.close()

    # Summary
    print()
    if no_summary > 0:
        pct = 100.0 * no_summary / total
        print(f"  Note: {no_summary:,} reads ({pct:.1f}%) not found in "
              f"sequencing_summary (dorado re-called reads may have new IDs).")
    pct_class = 100.0 * classified / total if total > 0 else 0
    print(f"  Total: {total:,} reads processed, {classified:,} written ({pct_class:.1f}%)")
    print()
    for bc in sorted(counts):
        print(f"    {bc}: {counts[bc]:,} reads")

    # Index all output BAMs (required for step 07 if needed)
    print()
    for barcode in sorted(out_files):
        out_path = os.path.join(output_dir, f"{prefix}_{barcode}.bam")
        print(f"  Indexing {os.path.basename(out_path)} ...", flush=True)
        pysam.index(out_path)

    print("\nDone.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Split mod_calls.bam using MinKNOW sequencing_summary barcode assignments",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Usage")[1] if "Usage" in __doc__ else ""
    )
    parser.add_argument("--mod-bam",    required=True,
                        help="Dorado mod_calls.bam (unaligned, with MM/ML tags)")
    parser.add_argument("--summary",    required=True,
                        help="MinKNOW sequencing_summary.txt")
    parser.add_argument("--output-dir", required=True,
                        help="Directory for per-barcode output BAMs")
    parser.add_argument("--prefix",     default="minknow_SQK-RBK114-24",
                        help="Output filename prefix (default: minknow_SQK-RBK114-24)")
    parser.add_argument("--barcodes",   nargs="+",
                        default=[f"barcode{i:02d}" for i in range(1, 8)],
                        help="Barcodes to extract (default: barcode01..barcode07)")
    parser.add_argument("--pass-only",  action="store_true",
                        help="Only include reads where passes_filtering=TRUE in summary")
    args = parser.parse_args()

    print(f"=== MinKNOW barcode split ===")
    print(f"  Input BAM:   {args.mod_bam}")
    print(f"  Summary:     {args.summary}")
    print(f"  Output dir:  {args.output_dir}")
    print(f"  Prefix:      {args.prefix}")
    print(f"  Barcodes:    {', '.join(args.barcodes)}")
    print(f"  Pass-only:   {args.pass_only}")
    print()

    print("Loading barcode map ...", flush=True)
    barcode_map = load_barcode_map(args.summary, pass_only=args.pass_only)

    print(f"Splitting {os.path.basename(args.mod_bam)} ...", flush=True)
    split_bam(args.mod_bam, barcode_map, args.output_dir,
              args.prefix, args.barcodes)


if __name__ == "__main__":
    main()
