# WARNINGS — validité & limites du pipeline phage

*Dernière mise à jour : 2026-07-17*

Ce fichier liste les paramètres et hypothèses qui conditionnent la **validité
d'interprétation** des résultats. À lire avant de conclure sur un échantillon.

---

## Illumina ~50 bp (polissage hybride)

- Les lectures Illumina de ce projet font **~50 bp** (natives). Elles servent
  **uniquement** au polissage de bases dans les **régions uniques** (Polypolish
  `--careful` + Pypolca `--careful`, étape 03), **jamais** à définir la structure.
- Elles ne résolvent **pas** la structure à travers les répétitions / DTR : la
  structure reste définie par l'ONT.
- Le QV « hybride » dépend **fortement de la profondeur** Illumina locale ; une
  couverture faible/inégale = gain de qualité limité et hétérogène.

## Méthylation

- Le basecalling MinKNOW utilisé est `5mCG_5hmCG` = **5mC en contexte CpG uniquement**.
  → Le **5mC non-CpG** (Dcm `CCWGG`, motifs de systèmes R-M) est **invisible**.
- Un résultat « **0 5mC** » **n'est pas concluant** sans re-basecalling all-context
  (voir `rebasecall_dorado.sh` : modèle SUP + `6mA 4mC_5mC`).
- **MicrobeMod 1.1.0 exige `modkit` 0.2.x** (`--only-tabs`, supprimé en ≥ 0.3).
  Sans un modkit 0.2.x épinglé (`MODKIT_02X_BIN`), chaque barcode échoue en rc=1.
- Un TSV MicrobeMod **vide n'est pas un échec** : c'est « 07a OK mais 0 motif au seuil ».
  Croiser avec `pct_hyper` de l'étape 07b (Bloc 1) pour distinguer absence réelle de
  méthylation vs artefact.

## Assemblage

- **Autocycler = défaut** (consensus multi-assembleurs). Repli automatique
  **Flye + Filtlong (~100× d'entrée)** si Autocycler (ou tous ses sous-assembleurs)
  indisponibles.
- Les **DTR** peuvent **sur-fragmenter Flye** (5–20 petits contigs au lieu d'un
  génome unique) → utiliser `FLYE_META_BARCODES` pour ces barcodes.
- Un **contig unique est une métrique QC**, pas une cible : la sortie de l'assembleur
  est conservée telle quelle. Si le génome reste fragmenté, `for_bandage.gfa` est
  exposé pour finition manuelle (voir `MANUAL_FINISHING.md`).

## Complétude

- La complétude est jugée par la **détection de DTR** (`finish_assembly.py`) et par
  **CheckV** (étape 04), **pas** par une colonne « circular » (retirée : elle était
  alimentée par `assembly_info.txt` de Flye, absent sous Autocycler → toujours 0, et
  biologiquement trompeuse pour des phages à DTR).

## Termini / réorientation

- Ordre : **PhageTerm prioritaire**, **dnaapler (terL) en repli**. La réorientation
  **n'est pas garantie** (terminase non trouvée, terminus ambigu).
- L'affichage « Termini / reorientation » a été **retiré** du rapport final : la
  réorientation est toujours exécutée en amont (étape 03), simplement plus affichée.

## Environnement / cluster

- Le **samtools de l'env micromamba est cassé** (stub). Toutes les étapes qui
  l'utilisent (04, 06) passent par `load_samtools_cluster()` → module **`samtools/1.21`**.
- Scratch : toujours `TMPDIR=/local/scratch/tmp` (jamais `/tmp`).

## Filtre parasite

- `MIN_CONTIG_LEN=35000` (35 kb) supprime les petits contigs (fragments d'hôte /
  débris). **Sécurité** : le filtre est **ignoré si aucun contig n'atteint le seuil**
  (un génome fragmenté n'est jamais effacé). Mettre à 0 pour désactiver.

## Comparaison — RETIRÉE (2026-07-17)

- Plus de **variant calling** read-based (breseq / snippy / gdtools).
- Plus de **comparaison d'assemblage** (nucdiff / MUMmer, ex-étape 08).
- Plus de **comparaison de méthylation clone-vs-référence** (ex-Section 2 de l'étape 07b).
- Chaque barcode est désormais **assemblé et analysé indépendamment**.
- Le code retiré est archivé sous `scripts/legacy/` (`08_compare.sh`).
