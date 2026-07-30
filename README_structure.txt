Nanopore phage pipeline — repository structure
==============================================

The pipeline is a single, flat set of scripts (no separate cluster/ and local/
copies). The same scripts run on Maestro (SLURM + Lmod modules) and locally;
helpers such as load_samtools_cluster() fall back gracefully when the cluster
module system is absent.

phage_nanopore/
├── config.sh                 ← ALL settings. Edit Section 1 per run.
│                               (paths, DBs, references, thresholds, walltimes,
│                                + internal helper functions)
├── submit_all.sh             ← THE entry point. SLURM job graph with per-sample
│                               input guards and afterany dependencies.
│                               Usage: bash submit_all.sh [--check|--dry-run|
│                                      --resume|--step NN]
├── 00_install_env.sh         ← one-time micromamba env + databases setup
├── finish_assembly.py        ← DTR/circular detection + GFA export (used by 02)
├── split_modbam_by_minknow.py← helper (modbase BAM splitting, if needed)
│
├── steps/                    ← pipeline steps, in execution order
│   ├── 00_setup.sh           ← detect samples, merge ONT FASTQs, sample_sheet.tsv
│   ├── 01_qc.sh              ← chopper (Q≥12, len≥200) + NanoStat + coverage
│   ├── 01b_hostdeplete.sh    ← OPTIONAL host-read depletion (HOST_DEPLETION=true)
│   ├── 02_assemble.sh        ← Autocycler consensus (default) | Flye+Filtlong~100x; ONT defines structure
│   ├── 03_polish.sh          ← Medaka [+ Polypolish + Pypolca if Illumina present]
│   │                            then PhageTerm termini + reorient (dnaapler fallback) + contig rename
│   ├── 04_assembly_qc.sh     ← CheckV (completeness/contamination/host) +
│   │                            read-remap coverage uniformity (depth CV, %low)
│   ├── 05_annotate.sh        ← Pharokka annotation
│   ├── 06_modbase.sh         ← minimap2 (mod BAM) + modkit pileup  [if MOD_DEMUX_DIR]
│   ├── 07a_methyl_call.sh    ← MicrobeMod per barcode (array)      [if MOD_DEMUX_DIR]
│   ├── 07b_methyl_report.sh  ← methylation stats + HTML report     [if MOD_DEMUX_DIR]
│   ├── 09_summary.sh         ← run-wide TSV + TXT + HTML summary (QC, CheckV,
│   │                            coverage, annotation, methylation)
│   └── 10_final_report.sh    ← per-sample + global HTML reports
│
├── scripts/legacy/           ← archived, non-canonical scripts (do not use)
├── MANUAL_FINISHING.md       ← Bandage manual-finishing procedure (out of pipeline)
├── PROJECT_STATUS.md         ← run status, biological context, design notes
├── TROUBLESHOOTING.md        ← known issues and fixes
└── README_hostdepleted_assembly.md  ← host-depletion side-workflow notes

Dependency graph (submit_all.sh)
--------------------------------
  00 setup
     → 01 qc → 02 assemble → 03 polish
                               ├─ 04 assembly_qc
                               ├─ 05 annotate
                               └─ 06 modbase* → 07a methyl_call* → 07b methyl_report*
        04 + 05 + 06 ───────► 09 summary
        09 + 07b ───────────► 10 final_report
  (* steps 06/07 run only if MOD_DEMUX_DIR is set in config.sh)

Note: assembly-vs-assembly comparison and read-based variant calling
(former step 08 / breseq / snippy / nucdiff) have been removed from the
pipeline. Each barcode is assembled and analysed independently. The archived
comparison script lives under scripts/legacy/.
