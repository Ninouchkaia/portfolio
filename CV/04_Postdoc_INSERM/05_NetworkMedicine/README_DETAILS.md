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