# **Drug Repurposing for COVID-19 through Network Medicine*

**Affiliation:** INSERM U1037 – Centre de Recherches en Cancérologie de Toulouse (CRCT)  
**Period:** 2020–2021  
**Publication:** [Network and Systems Medicine, 2020](https://www.liebertpub.com/doi/10.1089/nsm.2020.0011)


## 🧭 Context
At the onset of the COVID-19 pandemic, identifying potential therapeutic candidates required integrative strategies beyond single-target screening. This project used **network medicine** approaches to explore interactions between SARS-CoV-2 proteins, host cellular pathways, and drug targets, with the goal of repositioning existing compounds.  


## 🎯 Objectives
- Integrate multi-omics and molecular interaction data to construct a **virus–host–drug network**.  
- Identify biologically plausible drug candidates through **topological proximity** and **pathway enrichment**.  
- Test robustness of network-based predictions using simulated perturbations.  


## 🧪 Methods
- **Data integration:** Host–virus interactome from public datasets (BioGRID, IntAct), drug–target relationships from DrugBank and ChEMBL.  
- **Network modeling:** Weighted graph representation of molecular associations.  
- **Simulation:** Random rewiring and node removal to assess prediction stability.  
- **Analysis:** Centrality and community detection to highlight key druggable modules.  
- **Validation:** Cross-checking candidate lists with published clinical data and ongoing trials.  

## 💡 Contributions
- Implemented random network simulations to evaluate robustness of predicted drug–disease associations.  
- Automated analysis of node connectivity and topological metrics for ranking candidate drugs.  
- Contributed to visualization and reporting of systemic network perturbations.  
- Participated in manuscript review and interpretation of results.  

## 🔗 Reference
*Verstraete N.*, et al. *CovMulNet19, Integrating Proteins, Diseases, Drugs, and Symptoms: A Network Medicine Approach to COVID-19.*  
*Network and Systems Medicine*, 2020. [DOI:10.1089/nsm.2020.0011](https://www.liebertpub.com/doi/10.1089/nsm.2020.0011)


# CovMulNet19 – Workflow d’analyse multilayer et bootstrap

Ce dépôt contient un pipeline structuré pour analyser le réseau multilayer CovMulNet19 et comparer ses propriétés structurales à un ensemble de réseaux aléatoires (mock networks). Le réseau observé et les mocks ont été fournis par le CoMuNe Lab (M. De Domenico).

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

