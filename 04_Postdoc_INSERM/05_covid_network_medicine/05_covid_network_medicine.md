# CovMulNet19 – Workflow d’analyse multilayer et bootstrap

Ce dépôt contient un pipeline structuré pour analyser le réseau multilayer CovMulNet19 et comparer ses propriétés structurales à un ensemble de réseaux aléatoires (mock networks).  
Le réseau observé et les mocks ont été fournis par le CoMuNe Lab (M. De Domenico).

L’objectif du pipeline est d’estimer, pour différents types d’entités (GO, drugs, diseases, symptoms), leur degré structuré dans le réseau réel et leur position par rapport à un modèle nul obtenu via bootstrap multilayer.

---

## Structure du pipeline

Le workflow est organisé en six étapes :

1. **Import des réseaux**  
2. **Calcul des mesures observées**  
3. **Calcul des distributions mock (µ, σ, valeurs)**  
4. **Calcul des Z-scores**  
5. **Calcul des p-values (erf / Chebyshev, avec tests de normalité)**  
6. **Classement des entités et génération des figures**

---

## Étape 1 — Import des réseaux
Détection des fichiers nodes_XXXX.csv / edges_XXXX.csv et sélection des mock networks valides.

**Scripts :**
- [`5. covid_network_medicine/01_multilayer_pipeline/io_networks_clean.py`](5. covid_network_medicine/01_multilayer_pipeline/io_networks_clean.py)

```bash
python 01_multilayer_pipeline/io_networks_clean.py \
  --mock_basepath Mock_networks/
```
---

## Étape 2 — Mesures observées (réseau réel)
Calcul des degrés observés entre un type d’entités source et un type cible (directed ou undirected).

**Scripts :**
- [`01_multilayer_pipeline/compute_observed_degrees_clean.py`](01_multilayer_pipeline/compute_observed_degrees_clean.py)

```bash
python 01_multilayer_pipeline/compute_observed_degrees_clean.py \
  --nodes COVID19_GDDS_nodes.csv \
  --edges COVID19_GDDS_edges.csv \
  --source_type protein \
  --target_type GO \
  --directed True \
  --output results/observed/
```
---

## Étape 3 — Distributions mock (µ, σ, valeurs)
Agrégation des degrés à travers l’ensemble des mock networks pour obtenir les distributions structurales.

**Scripts :**
- [`02_bootstrap_pipeline/compute_mock_distributions_clean.py`](02_bootstrap_pipeline/compute_mock_distributions_clean.py)

```bash
python 02_bootstrap_pipeline/compute_mock_distributions_clean.py \
  --mock_basepath Mock_networks/ \
  --node_template_file COVID19_GDDS_nodes.csv \
  --focal_type GO \
  --neighbor_type protein \
  --outfile results/distributions/GO.tsv
```
---

## Étape 4 — Z-scores
Comparaison du réseau observé aux réseaux mock :  
Z = (observed − mean) / sd.

**Scripts :**
- [`02_bootstrap_pipeline/compute_zscores_clean.py`](02_bootstrap_pipeline/compute_zscores_clean.py)

```bash
python 02_bootstrap_pipeline/compute_zscores_clean.py \
  --focal_types GO \
  --observed_files results/observed/protein_to_GO_directed_withZeros.tsv \
  --distribution_files results/distributions/GO.tsv \
  --outdir results/zscores/
```
---

## Étape 5 — p-values
Tests de normalité (Shapiro, D’Agostino) → choix entre erf ou Chebyshev.

**Scripts :**
- [`02_bootstrap_pipeline/compute_pvalues_clean.py`](02_bootstrap_pipeline/compute_pvalues_clean.py)

```bash
python 02_bootstrap_pipeline/compute_pvalues_clean.py \
  --focal_types GO \
  --zscores_files results/zscores/GO_zscores.tsv \
  --distribution_files results/distributions/GO.tsv \
  --outdir results/pvalues/
```
---

## Étape 6 — Classement & figures
Classement des entités selon p-values ou Z-score et génération des figures associées.

**Scripts :**
- [`02_bootstrap_pipeline/ranking_and_plots_clean.py`](02_bootstrap_pipeline/ranking_and_plots_clean.py)

```bash
python 02_bootstrap_pipeline/ranking_and_plots_clean.py \
  --focal_types GO \
  --pvalue_files results/pvalues/GO_pvalues.tsv \
  --outdir results/final/ \
  --topN 20 \
  --ranking_metric p_shapiro \
  --plot_metric zscore
```
---

