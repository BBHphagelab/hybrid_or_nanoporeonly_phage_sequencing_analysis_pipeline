# Manual assembly finishing (Bandage) — out-of-pipeline procedure

The pipeline **never forces a single contig**. The assembler output stands as-is,
and the number of contigs is reported as a QC metric (step 04 / step 09 summary).
When a genome stays fragmented, step 02 writes the assembly graph so you can
finish it by hand. This file is the exact procedure.

## When is this needed?

Only when step 02 leaves **more than one contig** for a sample. In that case it
writes:

```
<RESULTS_DIR>/02_assembly/<barcode>_<sample>/for_bandage.gfa
```

If that file is absent, the assembly is already single-contig — nothing to do.

## Tool

[Bandage](https://rrwick.github.io/Bandage/) — a free graphical assembly-graph
viewer. Install it **on your own computer** (not the cluster). Copy the
`for_bandage.gfa` file locally (e.g. `scp`), then:

## Steps

1. **Open the graph** — Bandage → *File ▸ Load graph* → select `for_bandage.gfa`.
2. **Draw it** — click *Draw graph*. Nodes are contigs, edges are candidate
   junctions. A phage genome that is circular or carries Direct Terminal Repeats
   (DTR) typically shows up as a single loop or a node repeated at both ends.
3. **Identify the real path** — find the unambiguous walk through the graph that
   reconstructs the genome (usually one circle for a phage). Use coverage and
   BLAST (*Output ▸ Web BLAST selected nodes*) to disambiguate if needed.
4. **Export the finished sequence** — select the path / the node(s) you want,
   then *Output ▸ Save selected sequences to FASTA* (or the whole graph if a
   single component). Save it as your finished genome FASTA.

## Re-inject into the pipeline

Replace the raw assembly with your finished FASTA and rerun from polishing:

```bash
# 1. drop your manually finished genome in place
cp my_finished_genome.fasta \
   "<RESULTS_DIR>/02_assembly/<barcode>_<sample>/assembly.fasta"

# 2. re-run polish for THIS sample only (ARRAY_IDX = its row index, 0-based,
#    in sample_sheet.tsv, header excluded)
FORCE=true ARRAY_IDX=<index> bash steps/03_polish.sh

# 3. let the rest of the chain catch up (skips already-done samples)
bash submit_all.sh --resume
```

Notes:
- `<RESULTS_DIR>` is the value set in `config.sh`
  (currently `…/<run>/result_nanopore_optimisation`).
- Step 03 re-polishes, reorients (dnaapler, terL) and renames contigs, so the
  finished genome flows correctly into QC (04), annotation (05), modbase (06),
  methylation (07), run summary (09) and final report (10).
- Trimming of circular overhangs is **off by default** (`FINISH_TRIM_OVERHANG`
  in `config.sh`) to protect short phage DTRs — keep it off unless you are sure.
