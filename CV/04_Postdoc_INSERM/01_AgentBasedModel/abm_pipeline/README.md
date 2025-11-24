# ABM Pipeline – Command-Line Interface

This module provides the Python tools used to process, validate and analyse
the output of the NetLogo / OpenMOLE simulations for the Tumor–Immune ABM.

The pipeline consists of four components:

- `parameter_exploration/` – parsing BehaviorSpace outputs, Pareto extraction, best-set selection
- `model_validation/` – computing NRMSE and comparing model trajectories to experimental data
- `advanced_analysis/` – PCA, violin plots, parameter statistics
- `sensitivity/` – one-factor sensitivity analysis

All file paths are defined in `abm_pipeline/config.py`.

---

# 1. Running commands

All commands are executed via:

```bash
python -m abm_pipeline.cli <command> [options]
````

Typer converts internal names:

* `stats_params` → `stats-params`
* `make_behaviorspace` → `make-behaviorspace`
  etc.

---

# 2. Main commands

## 2.1. Pareto extraction

```bash
python -m abm_pipeline.cli pareto \
    --input results/behaviorspace/general_model/NSGAII_exploration_output.txt \
    --output data/pareto/general_pareto.tsv
```

## 2.2. Best-set extraction

```bash
python -m abm_pipeline.cli bestsets \
    --pareto data/pareto/general_pareto.tsv \
    --k 1 \
    --output data/pareto/general_bestset.tsv
```

For all patients:

```bash
python -m abm_pipeline.cli bestsets-all
```

---

## 2.3. BehaviorSpace generation

```bash
python -m abm_pipeline.cli make-behaviorspace \
    --template openmole/parameter_exploration.oms \
    --output openmole/generated_behaviorspace.xml
```

---

## 2.4. Validation

```bash
python -m abm_pipeline.cli validate-all
```

Outputs written to:

```
results/validation/
```

---

## 2.5. Advanced analysis

### PCA

```bash
python -m abm_pipeline.cli pca \
    --input data/pareto/general_pareto.tsv
```

### Violin plots

```bash
python -m abm_pipeline.cli violin-all
```

### Parameter statistics

```bash
python -m abm_pipeline.cli stats-params
```

Results saved to:

```
results/advanced_analysis/
```

---

## 2.6. Sensitivity analysis

Uses experiments defined in `abm_pipeline/sensitivity_config.py`.

```bash
python -m abm_pipeline.cli sensitivity-plot-all \
    --csv-dir data/sensitivity_csv \
    --save-dir results/sensitivity
```

---

# 3. Folder structure

```
data/
  experimental/        # experimental measurements
  pareto/              # extracted pareto fronts
results/
  validation/
  advanced_analysis/
  sensitivity/
openmole/
netlogo_model/
abm_pipeline/
```

Dossiers créés automatiquement par `ensure_dirs_exist()` dans `config.py`.

=======
# abm_pipeline — Internal Developer Documentation

## Overview

This document provides a complete technical description of the **abm_pipeline** Python package.  

The package implements the entire computational workflow used in the NLC–CLL agent-based modeling study:

1. Parameter exploration and optimization (OpenMOLE outputs processing)
2. Pareto front construction
3. Extraction of best parameter sets for each patient
4. Generation of BehaviorSpace XML files for NetLogo
5. Creation of shell scripts for headless simulation
6. Simulation validation (NRMSE, R², plots)
7. Sensitivity analyses
8. High-level parameter analyses (PCA, violin plots, stats)

---

# Package Structure

```
abm_pipeline/
|
├── cli.py
├── config.py
│
├── parameter_exploration/
│   ├── initial_ranges/
│   ├── nsga2_analysis/
│   ├── instantiate_models/
│   ├── shell_commands/
│   └── utils.py
│
├── model_validation/
│   ├── metrics.py
│   ├── plots.py
│   └── validator.py
│
├── sensitivity/
│   ├── plots.py
│   └── utils.py
│
└── advanced_analysis/
    ├── plots_advanced.py
    ├── stats_advanced.py
    └── utils.py
```

---

# config.py

Central configuration module containing:

- `PATIENTS_WITH_MONO`  
- `PATIENT_IDS`  
- Global constants and paths  

---

# CLI (cli.py)

Entry point for the entire pipeline.  
Run:

```
python -m abm_pipeline.cli --help
```

Provides commands for:

- Pareto construction  
- Best-set extraction  
- XML generation  
- Shell script generation  
- Validation plots  
- NRMSE calculations  
- Sensitivity plotting  
- PCA, violin plots, statistics  

---

# 1. Parameter Exploration

## 1.1 initial_ranges/

Tools to aggregate raw OpenMOLE outputs.

### aggregate_data.py  
Concatenates populationX.csv tables into a single aggregated file.

---

## 1.2 nsga2_analysis/

Handles the NSGA-II results from OpenMOLE.

### pareto_front.py  
Computes Pareto fronts and saves PNG plots.

### extract_best_sets.py  
Extracts best_via, knee_point, best_conc parameter sets  
and aggregates results across all patients.

### export_for_git.py  
Reformats Pareto tables for publication.

---

## 1.3 instantiate_models/

### xml_behavior_space.py  
Generates BehaviorSpace XML files for:

- best_via  
- best_conc  
- knee_point  
- class1/class2  
- sensitivity ranges  

### behavior_space_files.py  
Creates NetLogo-readable parameter files.

---

## 1.4 shell_commands/

Generates `.sh` scripts for NetLogo headless runs.

### patients.py  
Creates one script per patient.

### kneepoint.py  
Knee-point simulations (0, 1, 2).

### averaged.py  
Class-averaged simulations.

### sensitivity.py  
Shell scripts for sensitivity experiments.

---

# 2. Model Validation

## 2.1 metrics.py

Implements RMSE and NRMSE calculations exactly as in the original scripts.

Provides:

- rmse  
- nrmse_maxmin  
- nrmse_mean  
- nrmse_std  
- compute_viability_conc_nrmse  

---

## 2.2 plots.py

Produces simulation vs experimental comparison figures:

- Reads BehaviorSpace CSV  
- Aggregates runs  
- Computes NRMSE and R²  
- Two subplots (viability + concentration)  

---

## 2.3 validator.py

High-level module computing all metrics for all patients.  
Writes:

```
NRMSE_via_max_min.tsv
NRMSE_conc_mean.tsv
NRMSE_sum_std.tsv
```

---

# 3. Sensitivity Analysis

## sensitivity/utils.py  
Loads sensitivity CSV outputs.

## sensitivity/plots.py  
Plots viability and concentration curves for each perturbation.  
Batch mode supported.

---

# 4. Advanced Analysis

These analyses explore the global structure of optimized parameter space.

## Statistical tests between groups
plots_advanced.py  
Implements:
- Violin plots  
- PCA scatter plot  
- PCA scree plot  

### Violin plots
From a combined Pareto TSV file (e.g. data/pareto/pareto_ABM_2D_all.tsv):
```
python -m abm_pipeline.cli violin-all data/pareto/pareto_ABM_2D_all.tsv
```
This calls make_violinplots_all_parameters and generates a violin plot per parameter, saved under:
`results/advanced_analysis/violin/`

### 2. Principal Component Analysis (PCA)
```
python -m abm_pipeline.cli pca data/pareto/pareto_ABM_2D_all.tsv
```
This uses run_pca_analysis to standardize parameters, compute PCA and plot:
- 2D PCA scatter of parameter sets
- scree plot of explained variance

Outputs saved under:
`results/advanced_analysis/pca/`

## Statistical tests between groups
stats_advanced.py  
Mann–Whitney U tests on parameter distributions.
11.3 Statistical tests between groups

Given two Pareto TSV files (e.g. class1 vs class2):
```
python -m abm_pipeline.cli stats-params \
    data/pareto/pareto_class1.tsv \
    data/pareto/pareto_class2.tsv
```


This calls `run_parameter_stats_tests` to run Mann–Whitney U tests per parameter and write:
`results/advanced_analysis/stats/stats_comparison.tsv` 



---

# Example Workflow

### Pareto
```
python -m abm_pipeline.cli pareto PATIENT/filtered.txt PATIENT/pareto
```

### Best sets
```
python -m abm_pipeline.cli bestsets PATIENT/pareto.txt PATIENT/best.tsv
```

### XML generation
```
python -m abm_pipeline.cli make-behaviorspace PATIENT/best.tsv PATIENT/params.txt --patient-dict data/patient_dict.txt
```

### Shell scripts
```
python -m abm_pipeline.cli patient-shell .
```

### Simulations
```
bash patient_command_PATIENT.sh
```

### Validation
```
python -m abm_pipeline.cli plot-sim PATIENT stocha_knee_point
```

### NRMSE
```
python -m abm_pipeline.cli validate-all .
```

### Sensitivity
```
python -m abm_pipeline.cli sensitivity-plot-all sensitivity_csv results/sensitivity
```

### PCA
```
python -m abm_pipeline.cli pca pareto/all.tsv
```

---