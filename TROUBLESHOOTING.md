# Troubleshooting — phage_nanopore pipeline

Entrées taguées pour `grep` rapide :
- `[ENV]` cluster / SLURM / environnement
- `[INPUT]` fichiers d'entrée / sample sheet
- `[BIO]` problème biologique ou d'interprétation
- `[CODE]` bug ou limitation du pipeline

---

## [ENV] samtools trop ancien dans l'environnement micromamba

### Symptômes

```
[main] unrecognized command 'coverage'
[main] unrecognized command '--version'
[main] unrecognized command 'version'
sort: invalid option -- 'O'
```

Dans les logs des étapes qui appellent samtools (04 assembly_qc, 06 modbase).

### Cause

L'environnement micromamba `nanopore_phage` contient une copie de samtools
antérieure à la version 0.1.19 — si ancienne qu'elle ne reconnaît même pas
`samtools --version`. Les tentatives de mise à jour avec
`micromamba install "samtools>=1.17"` n'ont pas résolu le problème
(conflit de dépendances ou mauvais channel).

Le cluster Maestro dispose de `samtools/1.21` via son système de modules Lmod.

### Fix — charger samtools/1.21 depuis les modules Maestro

`load_samtools_cluster()` (dans `config.sh`) est appelée après `activate_env`
dans chaque étape qui utilise samtools (04, 06). Elle crée un wrapper
per-process (chemin absolu vers le samtools du module) et le préfixe à PATH —
concurrency-safe pour les tâches array. Le log confirme :
`[INFO] samtools active: samtools 1.21`.

---

## [INPUT] Doublons de contigs dans polished.fasta (barcode02)

### Symptôme

Step 04 (pharokka) plante avec :

```
ValueError: Duplicate key 'contig_270'
```

dans `post_processing.py` — pharokka vérifie les doublons dès le début
(`validate_fasta`).

### Cause — couche 1 : doublons dans polished.fasta

Flye a produit une assembly fragmentée pour barcode02 (14 contigs, couverture
12 817×). La fragmentation est due à la présence de DTR (Direct Terminal
Repeats) dans le génome phagique. Dans ce cas, Flye peut générer des IDs
non-uniques (ex. deux contigs nommés `contig_270`).

Cette première couche a été corrigée en renommant les doublons (ajout `_c2`...).
Pharokka confirme ensuite `All headers are unique`. **Mais le crash persiste.**

### Cause — couche 2 : collision de noms entre contig et protéines pyrodigal-gv

Même après correction des doublons, pharokka v1.9.1 plante dans
`post_processing.py` en lisant `prodigal-gv_aas_tmp.fasta`. Ce fichier est
généré par pyrodigal-gv qui nomme ses protéines `contig_N_geneI`. Quand les
contigs de polished.fasta s'appellent `contig_1`, `contig_2`, ..., `contig_270`,
ces noms entrent en collision avec les IDs de protéines produits en interne :
un ID de protéine peut valoir exactement `contig_270`, ce qui crée un doublon
dans le dict Python.

**Fix définitif** : renommer TOUS les contigs avec un préfixe spécifique à
l'échantillon (ex. `PAK_P4_ctg001`) qui ne peut pas entrer en collision avec
la nomenclature interne de pyrodigal-gv (`contig_N`).

### Commande cluster (fix définitif)

```bash
POLISHED="/pasteur/helix/projects/Phage3C/_projects/Nanopore/20260520_ONT_phages/results_nanopore/03_polished/barcode02_PAK_P4/polished.fasta"
SAMPLE="PAK_P4"

cp "${POLISHED}" "${POLISHED}.bak2"
python3 -c "
sample = '${SAMPLE}'
with open('${POLISHED}') as fi, open('${POLISHED}.fixed', 'w') as fo:
    n = 0
    for line in fi:
        if line.startswith('>'):
            n += 1
            fo.write(f'>{sample}_ctg{n:03d}\n')
        else:
            fo.write(line)
"
mv "${POLISHED}.fixed" "${POLISHED}"
grep '^>' "${POLISHED}"   # doit afficher PAK_P4_ctg001 ... PAK_P4_ctg014

RESULTS="/pasteur/helix/projects/Phage3C/_projects/Nanopore/20260520_ONT_phages/results_nanopore"
rm -rf "${RESULTS}/04_annotated/barcode02_PAK_P4"
cd /pasteur/helix/projects/Phage3C/_projects/Nanopore/phage_nanopore
bash submit_all.sh --step 04
```

# Supprimer l'ancien dossier d'annotation et relancer
RESULTS="/pasteur/helix/projects/Phage3C/_projects/Nanopore/20260520_ONT_phages/results_nanopore"
rm -rf "${RESULTS}/04_annotated/barcode02_PAK_P4"
cd /pasteur/helix/projects/Phage3C/_projects/Nanopore/phage_nanopore
bash submit_all.sh --step 04
```

### Note long terme

La fragmentation en 14 contigs est un problème d'assembly (DTR non résolus),
pas de couverture. Pour l'éliminer proprement, utiliser Bandage pour identifier
les jonctions terminales, puis `dnaapler` pour recirculariser. Sans reads
Illumina pour barcode02, Unicycler hybrid n'est pas possible.

---

## [INPUT] Désalignement des colonnes dans sample_sheet.tsv

### Symptôme

Une étape saute des barcodes, ou lit une mauvaise valeur (mode/illumina),
parce qu'une ligne n'a pas le bon nombre de colonnes.

### Cause

Ligne mal formée : tabulation manquante ou surnuméraire, colonne décalée.
Le sample sheet attend exactement 5 colonnes :
`barcode  sample_name  mode  illumina_r1  illumina_r2`.

### Fix

Régénérer avec `bash steps/00_setup.sh` (FORCE=true) ou vérifier avec :

```bash
awk -F'\t' '!/^#/{print NF, $1}' <RESULTS_DIR>/sample_sheet.tsv
```

Toutes les lignes de données doivent afficher `5 barcodeXX`.

---

## [ENV] Jobs SLURM array avec indices hors limite

### Symptôme

```
[SKIP] ARRAY_IDX=5 >= barcodes in sample_sheet (5) — rien a faire.
```

Pour les indices 5, 6, 7, 8.

### Cause

`submit_all.sh` calcule le nombre de barcodes en comptant les sous-dossiers
dans `fastq_pass/` (qui peut en avoir plus que le sample sheet). Les jobs
excédentaires sortent proprement avec exit 0.

### Fix

Pas un bug — comportement attendu. Ignorer ces messages dans les logs.

Amé
---

## [ASSEMBLY] Fragmentation Flye pour phages avec DTR — option --meta

### Symptôme

Flye produit 5–20 petits contigs au lieu d'un génome circulaire unique pour
un phage connu pour avoir des DTR (Direct Terminal Repeats, ex : barcode02
PAK_P4).  La couverture estimée est très élevée (~10 000×) mais l'assemblage
est fragmenté.

### Cause

Le graphe de répétitions de Flye confond les DTR avec des régions de répétition
internes et fragmente l'assemblage lors de la simplification du graphe.

### Fix

Utiliser le flag `--meta` de Flye, qui désactive la simplification agressive du
graphe de répétitions (conçue à l'origine pour les métagénomes, mais efficace
pour les génomes phagiques avec DTR).

Activer dans `config.sh` (Section 4) pour les barcodes concernés :

```bash
FLYE_META_BARCODES="barcode02"   # espace-séparé si plusieurs barcodes
```

`02_assemble.sh` applique automatiquement `--meta` pour chaque barcode listé.

### Notes

- `--meta` est **incompatible avec `--asm-coverage`** dans Flye (erreur fatale).
  `02_assemble.sh` désactive automatiquement `--asm-coverage` quand `--meta` est actif.
- La fragmentation se voit dans `06_summary.html` (colonne **Contigs** > 5,
  flag `[!contigs>5]`).
- Si `--meta` ne suffit pas, considérer : longueur minimale de reads plus
  élevée (`MIN_READ_LENGTH` dans `config.sh`), ou assembly hybride Unicycler
  avec des reads Illumina.

---

## [BUG] MicrobeMod (07a) — call_methylation échoue silencieusement pour tous les barcodes

### Symptômes

```
[WARN]    MicrobeMod failed for barcode01 (check log)
[WARN]    MicrobeMod failed for barcode02 (check log)
...
```

Chaque barcode échoue en ~1 seconde. `0 motif sites` dans le rapport HTML.

### Causes possibles et diagnostic

**1. BAM sans tags de modification (cause la plus fréquente)**

Si le BAM produit par l'étape 07 ne contient pas de tags `MM`/`Ml` (modification
base), MicrobeMod échoue immédiatement car aucune information de méthylation n'est
disponible.

Diagnostic :
```bash
samtools view /path/to/07_modbase/barcode01_PAK_P1/aligned_mods.bam \
  | head -2000 | grep -c "MM:Z\|Mm:Z"
# → 0 = pas de tags → problème de basecalling ou de l'étape 07
```

Solution : vérifier que MinKNOW a été configuré en **mod-base aware basecalling**
(HAC ou SUP avec le modèle incluant 6mA/5mC), et relancer l'étape 07.

**2. Version modkit incompatible**

MicrobeMod 1.1.0 requiert modkit **0.2.x**. Le module cluster est généralement
0.3.x qui a supprimé le flag `--only-tabs` utilisé par MicrobeMod.

Fix appliqué dans `08_methyl_analysis.sh` :
- Charger tous les modules cluster (`prodigal blast+ hmmer cath-tools modkit meme MicrobeMod`)
- Puis **préfixer** `PATH` avec le bin conda pour que le modkit 0.2.x du conda ait
  priorité sur celui du module cluster

Diagnostic version :
```bash
/pasteur/helix/users/csivelle/.mamba/envs/nanopore_phage/bin/modkit --version
# → doit afficher 0.2.x pour MicrobeMod 1.1.0
```

**3. Erreurs MicrobeMod — où les lire**

Depuis la correction de `08_methyl_analysis.sh`, chaque barcode a un log dédié :
```bash
cat results_nanopore/08_methyl/02_microbemod/barcode01_PAK_P1/microbemod_run.log
```

Le log global :
```bash
cat results_nanopore/08_methyl/08_methyl.log
```

**4. Streme / MEME non disponible**

MicrobeMod utilise STREME (MEME suite) pour la découverte de motifs. Si `streme`
n'est pas dans PATH, MicrobeMod plante.
```bash
module load meme
which streme
```

### Modules cluster requis (tous obligatoires)

```bash
module load prodigal blast+ hmmer cath-tools modkit meme MicrobeMod
```

Tous chargés automatiquement dans `08_methyl_analysis.sh`.

> **Note (2026-07-07) :** ce blocage `exit 2` immédiat est RÉSOLU. Le nouveau symptôme
> observé (TIMEOUT 4 h de 07a/07b) a une cause DIFFÉRENTE — voir la section suivante.

---

## [ENV][CODE] Step 07a/07b — TIMEOUT SLURM de 4 h (le pré-check samtools bloque avant MicrobeMod)

### Symptômes

```
Slurm Array  phage_07a_methyl_call  Ended, Mixed, MaxSignal [9]
Slurm Job    phage_07b_methyl_report  Failed, TIMEOUT, ExitCode 0  (run time 04:00:17)
```

- `07a_methyl_call.*.out` : n'affiche que la bannière, puis `CANCELLED ... DUE TO TIME LIMIT` à 4 h pile.
- `07_methyl.log` (07b) : s'arrête sur l'en-tête `══ Bloc 2 ══`, **sans** la ligne `MM tag check OK`.
- `02_microbemod/<bc>_<name>/microbemod_run.log` : **absent** → `MicrobeMod call_methylation`
  n'a JAMAIS démarré (ce fichier est créé par la redirection au lancement de MicrobeMod).

### Cause racine — le samtools de l'env micromamba est CASSÉ

Le hang de 4 h n'est PAS MicrobeMod : c'est le pré-check des tags de modification exécuté juste avant,
`samtools view aligned_mods.bam | head -2000 | grep -c "MM:Z\|Mm:Z"`, avec le **samtools de l'env
micromamba** (07a ne chargeait le samtools cluster qu'APRÈS ce pré-check ; l'ancien 07b ne le chargeait
jamais). Test cluster 2026-07-07 :

```
$ /pasteur/helix/users/csivelle/.mamba/envs/nanopore_phage/bin/samtools --version
Missing samtools executable, please report        # ← binaire corrompu (stub), plus le vrai samtools
$ module load samtools/1.21 && samtools view aligned_mods.bam | head -2000 | grep -c "MM:Z"
2000                                                # ← 1,0 s avec le module
```

La corruption est une séquelle d'une **ancienne** version de `load_samtools_cluster` qui écrasait
`env/bin/samtools` par un wrapper ; les tâches d'array parallèles l'ont abîmé (la version actuelle
utilise un wrapper `mktemp` par-process et n'y touche plus).

### Cause secondaire — 07b relançait MicrobeMod sans garde-fou

Comme 07a ne produisait aucun TSV, l'ancien Bloc 2 de `07b_methyl_report.sh` **relançait**
`call_methylation` lui-même, sur le BAM plein, **sans `timeout`** → 4 h → TIMEOUT SLURM.

### Fix (2026-07-07)

| Fichier | Correctif |
|---|---|
| `steps/07a_methyl_call.sh` | `load_samtools_cluster` déplacé AVANT le pré-check (→ samtools/1.21). Pré-check borné `timeout 600` (timeout = « on suppose taggé », MicrobeMod validera). Sous-échantillonnage borné `timeout ${MICROBEMOD_SUBSAMPLE_TIMEOUT}`. Appel MicrobeMod en `timeout -k 120 ${MICROBEMOD_TIMEOUT}`. |
| `steps/07b_methyl_report.sh` | Bloc 2 ne lance PLUS JAMAIS MicrobeMod : réutilise seulement le TSV de 07a ; sinon skip bruyant + `_MM_MISSING`. Ne touche plus du tout à samtools. Le résumé final remonte les barcodes sans sortie (couvre le cas SIGKILL sans sentinelle). |
| `config.sh` | `MICROBEMOD_TIMEOUT` 3 h→8 h ; ajout `MICROBEMOD_SUBSAMPLE_TIMEOUT=3600` ; walltimes séparés `TIME_METHYL_CALL="10:00:00"` (07a array) et `TIME_METHYL="02:00:00"` (07b gather). |
| `submit_all.sh` | 07a soumis avec `TIME_METHYL_CALL` (07b garde `TIME_METHYL`). |

Résultat : chaque étape se **termine toujours** (sentinelle honnête si échec), plus aucun TIMEOUT SLURM muet.

### Vérification après relance

`logs/07a_methyl_call.0.*.out` doit montrer, dans l'ordre :
`samtools active: samtools 1.21` → `MM tag check OK (…)` → `subsampling to ~300x` →
`running MicrobeMod call_methylation (…, timeout 28800s)` → puis `microbemod_run.log` se remplit.

### Dette à traiter (non bloquant)

Réinstaller un vrai samtools dans l'env pour ne plus dépendre uniquement du module :
```bash
micromamba install -n nanopore_phage samtools=1.21
```

---

## [ENV] MicrobeMod échoue rc=1 : modkit trop récent (--only-tabs supprimé) — CAUSE RACINE CONFIRMÉE

### Symptômes

Après la correction du hang samtools (section précédente), 07a démarre enfin MicrobeMod, puis :

```
07a_methyl_call.*.out :
[FAIL] barcode01 (PAK_P1): MicrobeMod call_methylation FAILED (rc=1)
...
[DIAG] modkit self-check (modkit 0.6.2):
error: unexpected argument '--only-tabs' found
[DIAG] ↑ this is the real modkit error MicrobeMod hid (rc=2).

07b : sentinelle FAILED (rc=1) pour les 5 barcodes → rapport HTML avec stats globales
      seulement, puis die volontaire (ExitCode 1) car TOUS les barcodes ont échoué.
```

### Cause (confirmée 2026-07-07)

MicrobeMod 1.1.0 appelle en interne `modkit pileup … --only-tabs --filter-threshold 0.66`.
Le flag **`--only-tabs` a été supprimé dans modkit ≥ 0.3**. Le modkit trouvé dans le PATH
(celui de l'env micromamba `nanopore_phage`) est en **0.6.2** → `exit 2`, que MicrobeMod masque
derrière un `CalledProcessError` → rc=1. Le `config.sh` prépendait le bin de l'env en CROYANT y
trouver du modkit 0.2.x ; en réalité l'env est en 0.6.2.

### Fix (2026-07-07) — fournir un modkit 0.2.x épinglé et le forcer en tête de PATH

1. Installer une fois modkit 0.2.x dans son propre env (paquet bioconda = **`ont-modkit`**) :
```bash
micromamba create -n modkit_0.2 -c bioconda -c conda-forge ont-modkit=0.2.6
"${HOME}/.mamba/envs/modkit_0.2/bin/modkit" --version   # doit afficher 0.2.6
```
2. `config.sh` : `MODKIT_02X_BIN="${MAMBA_ROOT_PREFIX}/envs/modkit_0.2/bin/modkit"` + helper
   `use_modkit_02x()` (wrapper mktemp par-process, même principe que `load_samtools_cluster`).
3. `steps/07a_methyl_call.sh` : appelle `use_modkit_02x` après `load_samtools_cluster`, donc
   MicrobeMod voit `modkit` = 0.2.6 en priorité.
4. `submit_all.sh --check` : avertit si `MODKIT_02X_BIN` absent / mauvaise version.

### Vérification après relance

`07a_methyl_call.0.*.out` doit montrer `[INFO] modkit active: modkit 0.2.6 …`, puis MicrobeMod
qui produit un `*.tsv` dans `02_microbemod/barcode01_PAK_P1/` et `[INFO] Step 07a done`.

### Alternative (non retenue)

Patcher `run_modkit` (microbemod.py l.82) pour la syntaxe modkit 0.6.x : fragile, car le format
de sortie de `modkit pileup` a aussi changé entre 0.2 et 0.6 (le parser `read_modkit` casserait).
Épingler modkit 0.2.x est l'option robuste.

---

## [BIO][CODE] MicrobeMod réussit (rc=0) mais ne produit aucun TSV — résultat vide ≠ échec

### Symptôme

`07a_methyl_call.*.out` : `MicrobeMod done.` puis `[WARN] MicrobeMod produced no TSV.`
`microbemod_run.log` (fin) :
```
Potential sites (>33%) for 6mA: 36 bases | High quality (>66%): 0 | >90%: 0
0 sites for STREME ... Fewer than 10 methylated sites identified, no motif finding was run.
Note: No methylated sites identified, so there is no output.
```
Avant correctif, 07b faisait `die` (ExitCode 1) car il confondait « pas de TSV » avec un échec.

### Cause

Ce n'est PAS un bug d'outil : MicrobeMod a tourné entièrement et n'a trouvé **aucun site méthylé
au-dessus du seuil** (>66 % des reads). Résultat vide légitime — p.ex. un phage sans méthyltransférase
active. MicrobeMod n'écrit alors aucun fichier de sortie.

### Fix (2026-07-07)

`07b_methyl_report.sh` distingue désormais, quand il n'y a pas de TSV :
- **07a sentinelle = OK** → résultat vide valide : BED vide, 0 motif, rapport produit, **exit 0**
  (`_MM_EMPTY`, journalisé « No methylation found »).
- **07a sentinelle = FAILED/TIMEOUT, ou absente (killed)** → vrai échec (`_MM_MISSING`).
Le `die` final ne se déclenche que si **tous** les barcodes sont en vrai échec — jamais pour des vides.

### Vérification biologique à faire (résultat vide réel vs artefact)

Comparer le verdict MicrobeMod avec les stats globales du Bloc 1 (issues du `methylation.bed` de
l'étape 06), qui utilisent un seuil de couverture indépendant :
```bash
cat …/07_methyl/01_stats/barcode01_stats.txt        # colonnes MeanMeth% / %Hyper
```
- Si Bloc 1 montre aussi ~0 % hyper-méthylation → cohérent, le phage n'est pas (ou peu) méthylé.
- Si Bloc 1 montre de la méthylation mais MicrobeMod 0 → revoir le modèle de basecalling mod
  (MinKNOW/Dorado `--modified-bases`), la profondeur après sous-échantillonnage, ou les seuils.

### Résultat observé sur ce run (2026-07-07) — CONCORDANT

Les deux méthodes concordent sur les 5 barcodes : MicrobeMod 0 site >66 % ; Bloc 1 `pct_hyper` ≈ 0
(0,03–0,27 %). → **Aucune méthylation significative**, en particulier **aucun 6mA** (appelé
toutes-contextes → verdict solide). Le pipeline est fiable ; ce n'est pas un artefact.

⚠️ **Angle mort basecalling :** le modèle `5mCG_5hmCG` (section Méthodes) ne détecte le **5mC
qu'en contexte CpG**. La 5mC bactérienne/phagique non-CpG (Dcm `CCWGG`, motifs R-M) est donc
**invisible** ici — le « 0 5mC » n'est PAS concluant, contrairement au « 0 6mA ». Pour la 5mC
non-CpG, re-basecaller avec un modèle 5mC toutes-contextes (changement à l'étape 06, pas 07).
