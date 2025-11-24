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

```
