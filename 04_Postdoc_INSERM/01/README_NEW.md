# ABM Parameter Exploration Pipeline

Ce dépôt contient les scripts utilisés pour l'exploration de l'espace de paramètres,
l'optimisation NSGA-II (OpenMOLE), la sélection des sets de paramètres et la
validation des modèles présentés dans l'article *[Titre du papier]*.

## Structure

- `abm_pipeline/config.py` : définition des patients et des chemins généraux
- `abm_pipeline/parameter_exploration/`
  - `initial_ranges/` : agrégation des sorties OpenMOLE (populationX.csv → outputs.txt)
  - `nsga2_analysis/` : calcul du Pareto front, extraction des best sets, export Git
  - `instantiate_models/` : fichiers de paramètres NetLogo (texte + XML BehaviorSpace)
  - `shell_commands/` : génération des scripts `.sh` pour NetLogo headless
- `abm_pipeline/model_validation/` *(à compléter)* : plots, scores, RMSE, PCA, violons
- `abm_pipeline/sensitivity_config.py` : définition des ranges et des expériences de sensibilité

Le CLI principal est accessible via :

```
python -m abm_pipeline.cli --help
```



## Pipeline principal (par patient)

Exemple pour le patient CAS1802 :

### 1. Agrégation des simulations OpenMOLE
```
python -m abm_pipeline.cli aggregate CAS1802/ABM_2D_CAS1802 CAS1802/outputs_ABM_2D_CAS1802.txt
```

### 2. Calcul Pareto + figure
```
python -m abm_pipeline.cli pareto \
    CAS1802/outputs_ABM_2D_CAS1802_duplicates_removed_filtered_only_samples_kept_50.0.txt \
    CAS1802/pareto_ABM_2D_CAS1802
```

### 3. Extraction des 3 sets (best_via, knee_point, best_conc)
```

python -m abm_pipeline.cli bestsets \
    CAS1802/pareto_ABM_2D_CAS1802.txt \
    CAS1802/best_param_sets_ABM_2D_CAS1802.tsv
```

### 4. Fichier de paramètres NetLogo pour le patient
```
python -m abm_pipeline.cli make-behaviorspace \
    CAS1802/best_param_sets_ABM_2D_CAS1802.tsv \
    CAS1802/netlogo_best_param_sets_ABM_2D_CAS1802.txt \
    --patient-dict patient_dict.txt
```

### 5. Génération du XML BehaviorSpace (3 expériences)

### 6. Génération du script .sh complet pour ce patient
```
python -m abm_pipeline.cli patient-shell .
```

### 7. Lancement du script sous Git Bash
```bash
bash patient_command_CAS1802.sh
```


## Analyses de sensibilité

### 1. Génération du XML de sensibilité à partir des paramètres du knee point
```
python -m abm_pipeline.cli sensi-xml \
    class1_processing/best_param_sets_pareto_ABM_2D_9patients_1_class1_50.tsv \
    sensitivity_analysis_experiment_file_class1.xml \
    1.28 4.55  # valeurs moyennes mono_init, apo_init
```

### 2. Génération du script .sh de sensibilité
```
python -m abm_pipeline.cli sensi-shell \
    /C/Users/Nina/.openmole/nina-windows/webui/projects/ABM_2D_9patients_1_class1.nlogo \
    sensitivity_analysis_experiment_file_class1.xml \
    sensitivity_experiments_class1.sh
```

### 3. Lancement (Git Bash)
```bash
bash sensitivity_experiments_class1.sh
```

