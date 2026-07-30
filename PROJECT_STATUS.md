# Nanopore Phage Pipeline — Project Status
*Run: `20260520_ONT_phages` — last updated: 2026-07-17*
*Pipeline: Autocycler assembly (Flye+Filtlong fallback), Medaka [+Illumina] polishing, PhageTerm/dnaapler reorientation, CheckV + coverage QC, Pharokka annotation, MicrobeMod methylation. Each barcode analysed independently. Results in `result_nanopore_optimisation/`.*

> **Comparaison retirée le 2026-07-17** — tout le volet comparaison (variant calling
> read-based breseq/snippy/gdtools, comparaison d'assemblage nucdiff, comparaison de
> méthylation clone-vs-référence) a été supprimé du pipeline. Le script archivé se
> trouve sous `scripts/legacy/` (voir aussi `08_compare.sh`).

---

## Quick orientation

| What | Where |
|---|---|
| Pipeline root | `phage_nanopore/` (synced to Maestro) |
| Config (edit Section 1 per run) | `config.sh` |
| Main orchestrator | `submit_all.sh` (steps 00→10) |
| Pre-flight only | `bash submit_all.sh --check` |
| Results on cluster | `…/Nanopore/20260520_ONT_phages/result_nanopore_optimisation/` |
| Logs | `${RESULTS_DIR}/logs/` |
| Micromamba env | `nanopore_phage` |

**Key files to read first:** `config.sh` (Section 1), `submit_all.sh`, `TROUBLESHOOTING.md`,
`result_nanopore_optimisation/sample_sheet.tsv`.

---

## Biological context

**Goal:** de novo assembly, QC, annotation and DNA-methylation profiling of
*Pseudomonas* phages — two parental phages (PAK_P1, PAK_P4) and three engineered /
recombinant derivatives. Each barcode is assembled and analysed independently
(no assembly-vs-assembly or read-based variant comparison).

**Samples (run 20260520, SQK-NBD114-24 native barcoding, R10.4.1):**

| Barcode | Sample | Role | Mode | Illumina |
|---|---|---|---|---|
| barcode01 | PAK_P1 | Parent phage | ont | — |
| barcode02 | PAK_P4 | Parent phage | ont | — |
| barcode03 | CHA_P1_rec | Recombinant (CHA_P1 lineage) | hybrid | yes (~50 bp reads) |
| barcode04 | PAK_P4_rec | Recombinant of PAK_P4 | hybrid | yes (~50 bp reads) |
| barcode05 | PAK_P4_CHA | PAK_P4 × CHA recombinant | ont (+Illumina) | yes (~50 bp reads) |

>  Illumina reads are ~50 bp mean length (native). They are used ONLY for short-read
> polishing in step 03 (Polypolish --careful + Pypolca --careful), never for assembly:
> high depth refines bases in unique regions but cannot resolve structure across
> repeats/DTRs, so genome structure is left to the ONT data.
>
> **Assembly:** Autocycler consensus by default (multi-assembler; best for terminus
> detection), with automatic fallback to Flye + Filtlong (~100x input cap). The former
> Unicycler hybrid and the
> forced single-contig SPAdes fallback have been removed). A single contig is a REPORTED
> QC metric, not a target. `finish_assembly.py` flags terminal repeats (DTR = complete
> genome) and circular overhangs and, if the assembly stays fragmented, copies the
> assembly graph to `for_bandage.gfa` for manual finishing (see MANUAL_FINISHING.md;
> overhang trimming off by default to protect short phage DTRs). It also drops contigs
> shorter than `MIN_CONTIG_LEN` (35 kb, config) as likely parasites/contaminants — with a
> safety that skips the filter if no contig reaches the threshold, so a fragmented genome
> is never wiped.

**Chemistry / basecalling:** R10.4.1 (FLO-MIN114), modified-base calling enabled in MinKNOW
(6mA + 5mCG_5hmCG). Per-barcode mod-tagged BAMs are in `bam_pass` (`MOD_DEMUX_DIR`).

**Observed raw ONT coverage (step 01):** bc01 ~5 250×, bc02 ~12 980×, bc03 ~9 270×,
bc04 ~8 660×, bc05 ~6 490× — all far above the fail threshold (20×).

---

## Pipeline overview

```
[00 setup]    -> builds sample_sheet.tsv, merges per-barcode FASTQ chunks
      |
[01 qc]       -> chopper (Q>=12, len>=200) + NanoStat + coverage estimate
[01b host*]   -> (optional) competitive host-read depletion before assembly  (* if HOST_DEPLETION=true)
      |
[02 assemble] -> Autocycler consensus (default) OR Flye+Filtlong(~100x); ONT defines structure
                 + parasite filter (<35 kb) + DTR/circular detection
                 + GFA for Bandage if still fragmented (see MANUAL_FINISHING.md)
      |
[03 polish]   -> Medaka; if Illumina present: + Polypolish --careful + Pypolca --careful
                 + PhageTerm termini detection -> reorient (dnaapler fallback) + contig rename
                 -> polished.fasta is the single reference for ALL downstream steps
      |
   +----------------+----------------+----------------+
[04 assembly_qc]  [05 annotate]   [06 modbase*]        (parallel after polishing)
 CheckV +          Pharokka        minimap2 -y +
 coverage CV                       modkit pileup
                                        |
                                  [07a methyl_call*] -> [07b methyl_report*]
                                   MicrobeMod            stats + HTML

[04]+[05]+[06]          -> [09 summary]  TSV / TXT / HTML
[09]+[07b]              -> [10 final_report]  per-sample + global HTML

  (* steps 06/07 run only if MOD_DEMUX_DIR is set)
```

All array->array dependencies use `afterany` plus per-barcode input guards, so a
failed or empty barcode is skipped cleanly and never blocks the others.

---

## Environment & cluster reference

| Item | Value |
|---|---|
| Conda **prohibited** on Maestro | use **micromamba** only (Institut Pasteur policy) |
| Main env | `nanopore_phage` — assembly, polishing, annotation, modbase |
| samtools | env copy is **broken** (stub: `Missing samtools executable`) → ALWAYS via `load_samtools_cluster()` (`samtools/1.21` Lmod). Reinstall pending: `micromamba install -n nanopore_phage samtools=1.21` |
| Pharokka DB | `…/databases/pharokka_db` (`PHAROKKA_DB`) |
| MicrobeMod (step 07a) | cluster modules `prodigal blast+ hmmer cath-tools modkit meme MicrobeMod`; needs modkit 0.2.x |
| Scratch | always `TMPDIR=/local/scratch/tmp` — never `/tmp` |
| Genome size estimate | `GENOME_SIZE_ESTIMATE=100k` (Flye `--genome-size`, coverage calc) |
| Flye --meta barcodes | `FLYE_META_BARCODES="barcode02"` (DTR over-fragmentation workaround) |

---

## Recent fixes (2026-07-07 session) — methylation TIMEOUT

Full re-run of the renumbered pipeline: 07a array ended `Mixed / MaxSignal[9]`, 07b `TIMEOUT` at
4 h. Root cause was NOT MicrobeMod — it never launched (`microbemod_run.log` absent everywhere). The
mod-tag pre-check `samtools view … | head | grep -c MM:Z` ran with the micromamba env samtools, which
is **corrupted** (stub `Missing samtools executable`); the cluster `samtools/1.21` was only loaded
*after* that pre-check. Separately, 07b re-ran `call_methylation` (full BAM, no timeout) as a fallback.

| Script | Issue | Fix |
|---|---|---|
| `steps/07a_methyl_call.sh` | pre-check used broken env samtools (before `load_samtools_cluster`); subsampling + MicrobeMod not effectively time-bounded | `load_samtools_cluster` moved BEFORE pre-check; pre-check `timeout 600`; subsampling `timeout ${MICROBEMOD_SUBSAMPLE_TIMEOUT}`; MicrobeMod `timeout -k 120 ${MICROBEMOD_TIMEOUT}` |
| `steps/07b_methyl_report.sh` | Bloc 2 re-ran MicrobeMod (no timeout, full BAM) → 4 h TIMEOUT | never runs MicrobeMod now — reuses 07a TSV only; missing barcodes surfaced loudly (incl. SIGKILL-without-sentinel case) |
| `config.sh` | single 4 h wall for both 07a & 07b; MicrobeMod timeout too close to wall | split `TIME_METHYL_CALL=10:00:00` (07a) / `TIME_METHYL=02:00:00` (07b); `MICROBEMOD_TIMEOUT` 3 h→8 h; new `MICROBEMOD_SUBSAMPLE_TIMEOUT=3600` |
| `submit_all.sh` | 07a used `TIME_METHYL` | 07a wired to `TIME_METHYL_CALL` |

All four scripts pass `bash -n`. Every 07 step now terminates with an honest sentinel — no silent SLURM
TIMEOUT. Relaunch: resync → `bash submit_all.sh --step 07a`, wait, then `--step 07b`.
See `TROUBLESHOOTING.md` → "Step 07a/07b — TIMEOUT SLURM de 4 h".

**Confirmed root cause after relaunch (2026-07-07):** with the hang gone, MicrobeMod finally launched and
the 07a self-check exposed the real error — `modkit 0.6.2: error: unexpected argument '--only-tabs'`.
MicrobeMod 1.1.0 needs **modkit 0.2.x** (`--only-tabs` was removed in modkit ≥ 0.3); the env modkit is
0.6.2. Fix: install pinned `ont-modkit=0.2.6` in its own env and force it to the front of PATH.

| Script | Change |
|---|---|
| `config.sh` | new `MODKIT_02X_BIN` (→ `envs/modkit_0.2/bin/modkit`) + helper `use_modkit_02x()` (per-process PATH wrapper) |
| `steps/07a_methyl_call.sh` | calls `use_modkit_02x` after `load_samtools_cluster` → MicrobeMod's `modkit` = 0.2.x |
| `submit_all.sh` | `--check` warns if modkit 0.2.x missing/wrong version |

**Prerequisite before relaunch:** `micromamba create -n modkit_0.2 -c bioconda -c conda-forge ont-modkit=0.2.6`.

## Recent fixes (2026-06-12 session)

These resolved a pipeline that was blocked at step 01 (each QC job finished its work then
exited on a trailing syntax error → SLURM marked it FAILED → the `afterok` chain never
advanced). Root cause: partial-write corruption that left stray text + trailing NUL bytes
at the end of several scripts.

| Script | Issue | Fix |
|---|---|---|
| `01_qc.sh` | stray `}"` + 3 NUL bytes at EOF | removed |
| `02_assemble.sh` | dnaapler if/else inverted; could kill step 02 if no terminase | if/else fixed + call made non-fatal |
| `03_polish.sh` | (new) contig names left as `contig_N` collide with pyrodigal-gv inside Pharokka | contig names normalised to `<SAMPLE>_ctgNNN` |
| `03_polish.sh` / `04` / `02` | reorientation duplicated and in the wrong place (frame mismatch vs methylation) | dnaapler reorientation consolidated into step 03; removed from 02/04 |

(Historical note: the 2026-06-12 session also touched the now-removed comparison
scripts — see `scripts/legacy/`.)

---

## Known issues / TODO

- **Genome reorientation:** consolidated into step 03 — PhageTerm/dnaapler now reorients
  `polished.fasta` once, so annotation (05), modbase (06) and methylation (07) all share a
  single reoriented, contig-renamed reference. **Fixed 2026-06-12.**
- **Sample sheet on the cluster:** step 00 preserves an existing `sample_sheet.tsv`
  (columns: `barcode sample_name mode illumina_r1 illumina_r2`). To regenerate, rerun
  step 00 with `FORCE=true`.
- See `TROUBLESHOOTING.md` for the samtools/Lmod, Pharokka duplicate-contig, MicrobeMod
  modkit-version, Flye DTR `--meta`, and step 07a/07b TIMEOUT issues.
- **Methylation (step 07) — 2026-07-07:** THREE issues fixed in sequence. (1) TIMEOUT: env samtools
  broken + missing time-bounds. (2) MicrobeMod rc=1: modkit 0.6.2 rejects `--only-tabs` (needs 0.2.x) →
  install `modkit_0.2` env. (3) After (1)+(2), MicrobeMod runs cleanly (rc=0) but for barcode01 found
  **0 methylated sites >66%** (only 36 6mA sites >33% out of 212 637; median 117× stranded) → no TSV.
  07b used to treat "no TSV" as failure and `die`; now it distinguishes **07a-OK-but-empty** (valid, 0
  motifs, report still produced, exit 0) from **07a-FAILED/killed** (real failure). **Resolved:** both
  independent methods agree across all 5 barcodes — MicrobeMod 0 sites >66%, Bloc 1 `pct_hyper` ≈ 0
  (0.03–0.27%) → **no significant methylation, in particular no 6mA** (all-context call → solid). Pipeline
  is trustworthy and works end-to-end. **Caveat:** basecalling used `5mCG_5hmCG`, which only detects 5mC
  in **CpG context** → non-CpG 5mC (Dcm CCWGG, R-M motifs) is invisible; the "0 5mC" is not conclusive.
  For non-CpG 5mC, re-basecall with an all-context 5mC model (step-06/basecalling change).
- **Broken env samtools (dette):** `nanopore_phage` env samtools is a non-functional stub. Not
  blocking (all steps use `load_samtools_cluster`), but reinstall when convenient:
  `micromamba install -n nanopore_phage samtools=1.21`.
- **Re-basecalling for broader methylation (2026-07-07):** `rebasecall_dorado.sh` (new, standalone)
  re-calls POD5 with SUP + all-context `6mA 4mC_5mC` (vs MinKNOW's `6mA 5mCG_5hmCG`), demuxes per
  barcode into `REBC_DEMUX_DIR`, leaving MinKNOW `bam_pass` untouched. Run `--check` first, then submit;
  then set `MOD_DEMUX_DIR=REBC_DEMUX_DIR` and `FORCE=true bash submit_all.sh --step 06` → 07a → 07b.
  Config: `DORADO_SIMPLEX_MODEL`, `DORADO_MOD_BASES`, `POD5_DIR`, `REBC_DEMUX_DIR`, `GPU_PARTITION`.
- **Barcode kit — RESOLVED (2026-07-07):** this run is **native barcoding → `SQK-NBD114-24`**
  (confirmed by user). `config.sh` set accordingly. The kit drives dorado demux; the earlier
  `SQK-RBK114-24` (rapid) value was wrong and has been corrected.

---

## Run the pipeline

```bash
# from the pipeline directory on Maestro:
bash submit_all.sh --check        # pre-flight (env, DB, refs, sample sheet)
bash submit_all.sh                # full chain 00→10
bash submit_all.sh --step 04      # single step, no dependencies
bash submit_all.sh --resume       # skip already-completed samples
squeue -u $(whoami)               # monitor
```
