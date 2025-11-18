# **Drug Repurposing for COVID-19 through Network Medicine**

**Affiliation:** INSERM U1037 – Centre de Recherches en Cancérologie de Toulouse (CRCT)  
**Period:** 2020–2021  
**Publication:** [Network and Systems Medicine, 2020](https://www.liebertpub.com/doi/10.1089/nsm.2020.0011)


## Context
At the onset of the COVID-19 pandemic, identifying potential therapeutic candidates required integrative strategies beyond single-target screening. This project used **network medicine** approaches to explore interactions between SARS-CoV-2 proteins, host cellular pathways, and drug targets, with the goal of repositioning existing compounds.  


## Objectives
- Integrate multi-omics and molecular interaction data to construct a **virus–host–drug network**.  
- Identify biologically plausible drug candidates through **topological proximity** and **pathway enrichment**.  
- Test robustness of network-based predictions using simulated perturbations.  


## Methods
- **Data integration:** Host–virus interactome from public datasets (BioGRID, IntAct), drug–target relationships from DrugBank and ChEMBL.  
- **Network modeling:** Weighted graph representation of molecular associations.  
- **Simulation:** Random rewiring and node removal to assess prediction stability.  
- **Analysis:** Centrality and community detection to highlight key druggable modules.  
- **Validation:** Cross-checking candidate lists with published clinical data and ongoing trials.  

## Contributions
- Implemented random network simulations to evaluate robustness of predicted drug–disease associations.  
- Automated analysis of node connectivity and topological metrics for ranking candidate drugs.  
- Contributed to visualization and reporting of systemic network perturbations.  
- Participated in manuscript review and interpretation of results.  

## 🔗 Reference
*Verstraete N.*, et al. *CovMulNet19, Integrating Proteins, Diseases, Drugs, and Symptoms: A Network Medicine Approach to COVID-19.*  
*Network and Systems Medicine*, 2020. [DOI:10.1089/nsm.2020.0011](https://www.liebertpub.com/doi/10.1089/nsm.2020.0011)


# CovMulNet19 – Multilayer and Bootstrap Analysis Workflow

This repository provides a structured pipeline to analyze the multilayer CovMulNet19 network and compare its structural properties to a collection of random networks (mock networks). The observed network and the mock networks were provided by the CoMuNe Lab (M. De Domenico).

The goal of the pipeline is to estimate, for different types of entities (GO, drugs, diseases, symptoms), their structured degree in the real network and their position relative to a null model obtained through multilayer bootstrap.

---

## Pipeline structure

The workflow is organized into six steps:

1. **Network import**
2. **Computation of observed measures**
3. **Computation of mock distributions (µ, σ, values)**
4. **Computation of Z-scores**
5. **Computation of p-values (erf / Chebyshev, with normality tests)**
6. **Entity ranking and figure generation**

---

## Step 1 — Network import

Detection of nodes_XXXX.csv / edges_XXXX.csv files and selection of valid mock networks.

**Scripts:**

* [`5. covid_network_medicine/01_multilayer_pipeline/io_networks_clean.py`](5. covid_network_medicine/01_multilayer_pipeline/io_networks_clean.py)

```bash
python 01_multilayer_pipeline/io_networks_clean.py \
  --mock_basepath Mock_networks/
```

---

## Step 2 — Observed measures (real network)

Computation of observed degrees between a source entity type and a target entity type (directed or undirected).

**Scripts:**

* [`01_multilayer_pipeline/compute_observed_degrees_clean.py`](01_multilayer_pipeline/compute_observed_degrees_clean.py)

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

## Step 3 — Mock distributions (µ, σ, values)

Aggregation of degrees across all mock networks to obtain structural distributions.

**Scripts:**

* [`02_bootstrap_pipeline/compute_mock_distributions_clean.py`](02_bootstrap_pipeline/compute_mock_distributions_clean.py)

```bash
python 02_bootstrap_pipeline/compute_mock_distributions_clean.py \
  --mock_basepath Mock_networks/ \
  --node_template_file COVID19_GDDS_nodes.csv \
  --focal_type GO \
  --neighbor_type protein \
  --outfile results/distributions/GO.tsv
```

---

## Step 4 — Z-scores

Comparison of the observed network to mock networks:
Z = (observed − mean) / sd.

**Scripts:**

* [`02_bootstrap_pipeline/compute_zscores_clean.py`](02_bootstrap_pipeline/compute_zscores_clean.py)

```bash
python 02_bootstrap_pipeline/compute_zscores_clean.py \
  --focal_types GO \
  --observed_files results/observed/protein_to_GO_directed_withZeros.tsv \
  --distribution_files results/distributions/GO.tsv \
  --outdir results/zscores/
```

---

## Step 5 — p-values

Normality tests (Shapiro, D’Agostino) → selection between erf or Chebyshev.

**Scripts:**

* [`02_bootstrap_pipeline/compute_pvalues_clean.py`](02_bootstrap_pipeline/compute_pvalues_clean.py)

```bash
python 02_bootstrap_pipeline/compute_pvalues_clean.py \
  --focal_types GO \
  --zscores_files results/zscores/GO_zscores.tsv \
  --distribution_files results/distributions/GO.tsv \
  --outdir results/pvalues/
```

---

## Step 6 — Ranking & figures

Ranking of entities based on p-values or Z-scores and generation of associated figures.

**Scripts:**

* [`02_bootstrap_pipeline/ranking_and_plots_clean.py`](02_bootstrap_pipeline/ranking_and_plots_clean.py)

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

