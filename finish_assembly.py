#!/usr/bin/env python3
"""
finish_assembly.py — semi-automatic single-contig finishing for phage assemblies.

Two conservative operations:

1. Parasite/contaminant length filter (--min-contig-len):
   Remove contigs shorter than the threshold (likely host fragments, small
   plasmids, or assembly debris). SAFETY: if NO contig reaches the threshold
   (e.g. a genuinely fragmented small genome), the filter is skipped so the whole
   assembly is never wiped.

2. Terminal-repeat (DTR / circular-overhang) detection:
   For each kept contig, compute the longest exact terminal repeat r (the longest
   proper prefix that is also a suffix, via the KMP border / prefix-function,
   O(n)). This is exactly the largest r with seq[:r] == seq[-r:].
   Classification:
     - r == 0                 → 'linear'
     - 0 < r <= --circ-max    → 'overhang?'        (circular overhang or short DTR)
     - r  > --circ-max        → 'complete(DTR)'    (Direct Terminal Repeat; kept)
     - r  > 0.45*L            → 'ambiguous(repeat)' (very repetitive; inspect)
   Trimming is OFF by default and only applied with --trim (removes the redundant
   suffix for terminal repeats <= --circ-max). Conservative because a real phage
   DTR (e.g. T7 ~160 bp) must NOT be trimmed.

Never merges separate contigs — use Bandage on the assembly GFA for that.
Writes the finished FASTA and a per-contig TSV report (kept + dropped).
"""
import argparse, sys


def read_fasta(path):
    name, seq, out = None, [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(seq)))
                name = line[1:]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        out.append((name, "".join(seq)))
    return out


def longest_border(s):
    """Longest proper prefix of s that is also a suffix (KMP prefix function). O(n)."""
    n = len(s)
    if n == 0:
        return 0
    pi = [0] * n
    k = 0
    for i in range(1, n):
        while k > 0 and s[i] != s[k]:
            k = pi[k - 1]
        if s[i] == s[k]:
            k += 1
        pi[i] = k
    return pi[-1]


def classify(r, L, circ_max):
    if r == 0:
        return "linear", ""
    if r > 0.45 * L:
        return "ambiguous(repeat)", f"terminal_repeat={r}bp (>45% of contig — inspect manually)"
    if r > circ_max:
        return "complete(DTR)", f"DTR={r}bp (genome complete; DTR kept)"
    return "overhang?", f"terminal_repeat={r}bp (circular overhang or short DTR)"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", required=True, help="input assembly FASTA")
    ap.add_argument("--out", required=True, help="output (finished) FASTA")
    ap.add_argument("--report", required=True, help="per-contig TSV report")
    ap.add_argument("--min-contig-len", type=int, default=0,
                    help="drop contigs shorter than this many bp as likely parasites/"
                         "contaminants; 0 disables. SAFETY: skipped if no contig reaches "
                         "the threshold, so the assembly is never fully removed (default 0)")
    ap.add_argument("--circ-max", type=int, default=1000,
                    help="terminal repeat <= this is treated as a circular overhang "
                         "(trimmed only with --trim); larger is reported as a DTR "
                         "(default 1000)")
    ap.add_argument("--trim", action="store_true",
                    help="trim circular overhangs (terminal repeat <= --circ-max). "
                         "OFF by default to protect real short DTRs.")
    a = ap.parse_args()

    contigs = read_fasta(a.inp)
    if not contigs:
        sys.stderr.write("[finish] ERROR: no sequences in %s\n" % a.inp)
        sys.exit(1)

    # ── 1. Parasite/contaminant length filter (with never-drop-all safety) ──
    kept, dropped = list(contigs), []
    if a.min_contig_len > 0:
        has_big = any(len(s) >= a.min_contig_len for _, s in contigs)
        if has_big:
            kept = [(n, s) for n, s in contigs if len(s) >= a.min_contig_len]
            dropped = [(n, s) for n, s in contigs if len(s) < a.min_contig_len]
            for n, s in dropped:
                sys.stderr.write(f"[finish] dropped {n} ({len(s)}bp < {a.min_contig_len}bp) "
                                 f"— likely parasite/contaminant\n")
        else:
            sys.stderr.write(f"[finish] WARNING: all {len(contigs)} contig(s) are shorter "
                             f"than {a.min_contig_len}bp — length filter SKIPPED so the "
                             f"assembly is not wiped (genome may be fragmented)\n")

    # ── 2. Terminal-repeat analysis on kept contigs ──
    rows, out_recs = [], []
    for name, seq in kept:
        L = len(seq)
        r = longest_border(seq)
        cls, note = classify(r, L, circ_max=a.circ_max)
        new = seq
        if a.trim and cls == "overhang?":
            new = seq[:-r]
            cls = "circular(trimmed)"
            note = f"overhang={r}bp trimmed"
        header = name
        if (cls.startswith("circular") or cls.startswith("complete")) and "circular" not in header.lower():
            header = name + " circular=true"
        out_recs.append((header, new))
        rows.append((name, str(L), str(len(new)), str(r), cls, note))

    for name, seq in dropped:
        rows.append((name, str(len(seq)), "0", "-", "DROPPED(short)",
                     f"< {a.min_contig_len}bp — removed as likely parasite/contaminant"))

    with open(a.out, "w") as fo:
        for header, seq in out_recs:
            fo.write(">" + header + "\n")
            for i in range(0, len(seq), 80):
                fo.write(seq[i:i + 80] + "\n")

    with open(a.report, "w") as fr:
        fr.write("contig\tlen_in\tlen_out\tterminal_repeat_bp\tclassification\tnote\n")
        for row in rows:
            fr.write("\t".join(row) + "\n")

    ncomplete = sum(1 for r in rows if r[4].startswith(("circular", "complete")))
    print(f"[finish] {len(contigs)} contig(s) in; {len(dropped)} dropped (<{a.min_contig_len}bp); "
          f"{len(kept)} kept; {ncomplete} circular/complete; report: {a.report}")
    for row in rows:
        print(f"[finish]   {row[0]}: len {row[1]}->{row[2]}, "
              f"terminal_repeat {row[3]}bp → {row[4]} {row[5]}")


if __name__ == "__main__":
    main()
