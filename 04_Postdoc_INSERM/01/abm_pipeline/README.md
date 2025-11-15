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

## plots_advanced.py  
Implements:

- Violin plots  
- PCA scatter plot  
- PCA scree plot  

## stats_advanced.py  
Mann–Whitney U tests on parameter distributions.

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


