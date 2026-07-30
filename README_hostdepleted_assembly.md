# Assemblage des phages par déplétion d'hôte

Correctif autonome pour l'assemblage des 5 phages du run `20260520_ONT_phages`.
Remplace l'étape 02 quand l'assemblage *de novo* échoue à cause de la
contamination par l'ADN génomique de l'hôte (*Pseudomonas aeruginosa*).

## Le problème

Les reads ONT contiennent beaucoup d'ADN génomique de l'hôte (résidu de prep).
L'assembleur assemble donc l'hôte **+** le phage, d'où les totaux aberrants vus
dans `06_summary` :

| Barcode | Échantillon | Contigs | Total | attendu |
|---|---|---|---|---|
| bc03 CHA_P1_rec | 1464 | 3,97 Mb | ~90 kb |
| bc04 PAK_P4_rec | 196 | 4,93 Mb | ~90 kb |
| bc05 PAK_P4_CHA | 257 | 6,6 Mb | ~90 kb |
| bc02 PAK_P4 | 16 | 0,93 Mb | ~90 kb |

(6,6 Mb ≈ un génome entier de *P. aeruginosa* → c'est bien de l'hôte.)

## La solution

Avant d'assembler, on **retire les reads d'hôte** par mapping compétitif, puis
on assemble les reads propres, et on garde Illumina pour le polissage des bases
(là où il existe). C'est l'ordre canonique « assemblage long-read → polish
short-read », appliqué à des reads décontaminés.

**Assembleur (objectif 1 contig de qualité, pour PhageTerm).** Par défaut
**Autocycler** (Wick 2025, successeur de Trycycler) : assemblage consensus
multi-assembleurs (flye, raven, miniasm, plassembler… selon ce qui est installé)
sur plusieurs sous-échantillons indépendants → meilleure exactitude structurale
et des extrémités. Autocycler gère lui-même le sous-échantillonnage, donc la
sur-couverture (5 000–13 000×) n'est plus un problème. **Repli automatique** sur
**Flye + Filtlong 100×** (downsampling pondéré qualité/longueur) si Autocycler ou
ses sous-assembleurs sont absents ou échouent. Choix via `--assembler autocycler|flye`.
Réorientation laissée à l'étape 03 (dnaapler), donc un assemblage **non
réorienté** reste disponible pour PhageTerm.

**Mapping compétitif (le point clé).** Chaque read est mappé à la fois sur
l'hôte et sur les références phages. Un read n'est retiré **que** s'il
ressemble *davantage* à l'hôte qu'à n'importe quelle référence phage (et qu'au
moins `HOST_COV_FRAC` de sa longueur s'aligne sur l'hôte). Conséquences :

- Les **séquences nouvelles/ingéniées** (aucun hit hôte) → **conservées**.
  C'est ce qui lève la tension « pas de génome hôte mais on veut garder le
  nouveau » : la déplétion d'hôte enlève l'hôte sans toucher au reste.
- Les **jonctions de recombinaison** (reads ONT longs, partiellement hôte) →
  **conservées** (fraction hôte sous le seuil).
- Protection contre un **prophage résident** de l'hôte similaire à ton phage :
  un read qui matche mieux ta référence phage que l'hôte est **gardé**.

Un garde-fou supplémentaire à l'assemblage sélectionne les contigs par
couverture / taille / circularité (le phage est à 800–13 000×, le débris hôte
résiduel bien plus bas), sans jamais vider l'assemblage.

## Compatibilité pipeline

Le script écrit `assembly.fasta` **à l'emplacement standard de l'étape 02**
(`results_nanopore/02_assembly/<barcode>_<sample>/assembly.fasta`) et sauvegarde
l'ancien assemblage contaminé en `.contaminated.bak.fasta`. Tu **reprends
ensuite le pipeline normal à l'étape 03** : Polypolish+Pypolca (bc03/04/05),
Medaka (bc01/02), réorientation dnaapler + renommage `<SAMPLE>_ctgNNN` →
`polished.fasta`. Tout l'aval (04 annotation, 07/08 méthylation, 09* comparaisons)
consomme `polished.fasta` **sans modification**. Bonus : ça répare aussi **bc02**,
ton point de comparaison pour bc04/bc05.

## Les 3 fichiers

- `fetch_host_refs.sh` — télécharge les génomes hôtes et fabrique
  `host_combined.fasta`. **À lancer sur une machine avec Internet** (les nœuds
  de calcul Maestro n'en ont pas).
- `assemble_phages_hostdepleted.sh` — le correctif : déplétion → Flye →
  sélection → finishing. **À lancer sur le cluster.**
- `README_hostdepleted_assembly.md` — ce fichier.

## Mode d'emploi

### 1. Récupérer le génome de l'hôte (machine avec Internet)

```bash
bash fetch_host_refs.sh
# -> host_refs/host_combined.fasta  (PAK CP020659.1 + PAO1 NC_002516.2)
```

Par défaut : **PAK** (`CP020659.1`, hôte de PAK_P1/PAK_P4 et de leurs
recombinants) + **PAO1** (`NC_002516.2`, référence universelle *P. aeruginosa*
qui capte la lignée CHA via le cœur génomique conservé à >99 %).

> Une assemblée publique « *P. aeruginosa* CHA » n'est pas identifiable sans
> ambiguïté sur NCBI. Si tu as le **génome CHA de ton labo** (le plus précis,
> puisque c'est l'hôte réel sur lequel tu as propagé les phages), ajoute-le :
> ```bash
> bash fetch_host_refs.sh /chemin/vers/CHA_host.fasta
> ```

Transfère ensuite le fichier sur le cluster :
```bash
scp host_refs/host_combined.fasta <user>@maestro:/pasteur/.../20260520_ONT_phages/refs/
```

### 2. Assembler (cluster)

```bash
cd /pasteur/.../phage_nanopore     # là où sont config.sh et les steps/

# les 3 prioritaires :
bash assemble_phages_hostdepleted.sh \
     --host-ref /pasteur/.../refs/host_combined.fasta \
     --barcodes barcode03,barcode04,barcode05 --force

# ou les 5 (recommandé, répare aussi bc02) :
bash assemble_phages_hostdepleted.sh \
     --host-ref /pasteur/.../refs/host_combined.fasta --force
```

Le script utilise les références phages de `config.sh` (`EXTERNAL_REF_1/2/3`,
les `.gbk`) pour le mapping compétitif. Tu peux aussi passer ta propre FASTA
phage via `--phage-ref`.

### 3. Nettoyer les résultats aval périmés, puis relancer

Chaque étape se **saute si sa sortie existe déjà**. Les résultats 03→10 du run
contaminé doivent donc être supprimés, sinon ils seraient conservés. Utilise le
script de nettoyage (dry-run par défaut) :

```bash
bash clean_for_reassembly.sh          # montre ce qui serait supprimé
bash clean_for_reassembly.sh --yes    # supprime réellement
```

Puis relance le pipeline (l'étape 02 se saute toute seule et garde ton nouvel
assemblage ; il n'y a **pas** de flag `--from`) :

```bash
bash submit_all.sh                    # chaîne complète 00→10
# ou une seule étape :  bash submit_all.sh --step 03
```

> `submit_all.sh` accepte `--check`, `--dry-run`, `--resume`, `--step N`. Pour
> tout régénérer sans rien supprimer, tu peux aussi forcer : `FORCE=true bash
> submit_all.sh` (dépend de l'export d'environnement SLURM ; la suppression est
> plus fiable).

## Sorties

- `02_assembly/<bc>_<sample>/assembly.fasta` — assemblage propre (entrée de 03).
- `02_assembly/<bc>_<sample>/assembly.contaminated.bak.fasta` — ancien, sauvegardé.
- `02b_hostdepleted/<bc>_<sample>/` — reads propres, PAF, logs intermédiaires.
- `02b_hostdepleted/decontamination_summary.tsv` — reads avant/après, % gardés,
  contigs Flye vs retenus, bp, circularité. **À vérifier en premier.**

## Paramètres réglables (env ou flags)

| Variable | Défaut | Rôle |
|---|---|---|
| `HOST_COV_FRAC` | 0.5 | un read est « hôte » si ≥ cette fraction s'aligne sur l'hôte **et** qu'il matche l'hôte mieux que toute réf phage |
| `COV_FRAC` | 0.10 | garde les contigs Flye de couverture ≥ `COV_FRAC` × couverture max |
| `MIN_PHAGE_LEN` / `MAX_PHAGE_LEN` | 20000 / 200000 | fenêtre de taille plausible d'un contig phage |
| `--force` | — | refait même si `assembly.fasta` existe |

## Points de vigilance

- **Vérifie `decontamination_summary.tsv`.** Un % de reads gardés très bas
  (<10 %) sur un parent (bc01/bc02) peut signaler un `HOST_COV_FRAC` trop
  agressif ou un mauvais génome hôte — ajuste et relance avec `--force`.
- **bc05 (PAK_P4 × CHA)** est une mosaïque de deux phages : le mapping
  compétitif utilise l'**union** des références (PAK_P1 + PAK_P4 + CHA_P1), donc
  les deux origines sont protégées. Garde quand même bc05 à l'œil.
- **Qualité de base (optionnel).** Pour les hybrides, l'étape 03 applique
  Polypolish+Pypolca (Illumina) sur l'assemblage Flye — très bon. Un Medaka
  (ONT) *avant* le polish Illumina gagnerait quelques bases ; non implémenté ici
  pour rester non intrusif, mais facile à ajouter si besoin.
- **Inserts nouveaux longs.** S'ils s'assemblent en contig séparé, le garde-fou
  les garde s'ils sont à haute couverture (≥ `COV_FRAC` × max) et dans la
  fenêtre de taille. Inspecte `for_bandage.gfa` en cas de doute.
```
