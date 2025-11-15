abm_pipeline/
│
├── __init__.py
├── cli.py
├── config.py
│
└── parameter_exploration/
    ├── __init__.py
    ├── utils.py
    │
    ├── initial_ranges/
    │   ├── __init__.py
    │   └── aggregate_data.py
    │
    ├── nsga2_analysis/
    │   ├── __init__.py
    │   ├── pareto_front.py
    │   ├── extract_best_sets.py
    │   └── export_for_git.py
    │
    └── instantiate_models/
        ├── __init__.py
        └── behavior_space_files.py


```python
python -m abm_pipeline.cli aggregate PATH/TO/ABM_2D_CAS1802 outputs_ABM_2D_CAS1802.txt
python -m abm_pipeline.cli pareto CAS1802/outputs_ABM_2D_CAS1802_duplicates_removed_filtered_only_samples_kept_50.0.txt CAS1802/pareto_ABM_2D_CAS1802
python -m abm_pipeline.cli bestsets CAS1802/pareto_ABM_2D_CAS1802.txt CAS1802/best_param_sets_ABM_2D_CAS1802.tsv
python -m abm_pipeline.cli make-behaviorspace CAS1802/best_param_sets_ABM_2D_CAS1802.tsv CAS1802/netlogo_best_param_sets_ABM_2D_CAS1802.txt --patient-dict patient_dict.txt

```


# Une figure pour un patient & un set
python -m abm_pipeline.cli plot-sim CAS1802 stocha_knee_point

# Tous les NRMSE pour tous les patients (3 param_sets par défaut)
python -m abm_pipeline.cli validate-all .


python -m abm_pipeline.cli sensitivity-plot perturb-gui-apo-mov \
    sensitivity_csv \
    figures_sensitivity

    

>>> faire un jupyternotebook !!!
Reproduire les figures de l'article 

Figure 3B–C :

Calcul du Pareto (commandes pareto, bestsets)

Génération de la figure via plot_pareto.py (à ajouter dans model_validation)

Figure 5A–C :

Lancer les scripts patient_command_{PATIENT}.sh

Utiliser plot_sim_vs_exp_with_scores.py dans abm_pipeline/model_validation/

Figure 4A–C, S4, S5… :

Utiliser les modules violinplots, pca, stats_tests dans abm_pipeline/model_validation/.