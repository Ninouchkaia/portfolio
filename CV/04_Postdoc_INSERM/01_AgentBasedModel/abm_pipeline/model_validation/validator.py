# abm_pipeline/model_validation/validator.py

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import pandas as pd

from abm_pipeline.model_validation.metrics import compute_viability_conc_nrmse
from abm_pipeline.parameter_exploration.utils import (
    get_patient_ids,
    logger,
)
from abm_pipeline.model_validation.plots import _load_patient_exp_data, _load_behaviorspace_csv


DEFAULT_PARAM_SETS = ["stocha_best_via", "stocha_best_conc", "stocha_knee_point"]


def compute_metrics_for_patient(
    patient: str,
    param_sets: List[str] | None = None,
    base_dir: str = ".",
) -> Dict[str, Dict[str, Dict[str, float]]]:
    """
    Calcule tous les NRMSE (via/conc/sum, maxmin/mean/std) pour un patient
    et plusieurs param_sets, comme dans RMSE.py.

    Retourne un dict:
    {
      'stocha_best_via': {
         'viability': {...},
         'concentration': {...},
         'sum': {...},
      },
      ...
    }
    """
    if param_sets is None:
        param_sets = DEFAULT_PARAM_SETS

    base_path = Path(base_dir)
    patient_data = _load_patient_exp_data(patient, base_dir=base_dir)

    all_results: Dict[str, Dict[str, Dict[str, float]]] = {}

    for param_set in param_sets:
        logger.info(f"[{patient}] metrics for {param_set}")
        simu_file = base_path / patient / "BehaviorSpace" / f"{param_set}.csv"
        df_viability_simu, df_remaining_simu = _load_behaviorspace_csv(str(simu_file))

        time_points = list(patient_data["Day"])
        filtered_viability_simu = df_viability_simu[
            df_viability_simu.index.isin(time_points)
        ]
        filtered_remaining_simu = df_remaining_simu[
            df_remaining_simu.index.isin(time_points)
        ]

        via_simu_mean = filtered_viability_simu.mean(axis=1)
        conc_simu_mean = filtered_remaining_simu.mean(axis=1)

        via_simu_mean.index = time_points
        conc_simu_mean.index = time_points

        metrics = compute_viability_conc_nrmse(
            patient_data=patient_data,
            viability_sim=via_simu_mean,
            conc_sim=conc_simu_mean,
            patient=patient,
        )
        all_results[param_set] = metrics

    return all_results


def compute_metrics_all_patients(
    param_sets: List[str] | None = None,
    base_dir: str = ".",
    output_prefix: str = "NRMSE",
) -> None:
    """
    Reproduit la structure de sortie de RMSE.py :

    - NRMSE_via_max_min.tsv
    - NRMSE_via_mean.tsv
    - NRMSE_via_stdev.tsv
    - ... idem pour concentration et somme

    Chaque fichier : index = patients, colonnes = param_sets
    """
    if param_sets is None:
        param_sets = DEFAULT_PARAM_SETS

    patients = get_patient_ids()

    # DataFrames vides
    df_via_maxmin = pd.DataFrame(columns=param_sets, index=patients)
    df_via_mean = pd.DataFrame(columns=param_sets, index=patients)
    df_via_std = pd.DataFrame(columns=param_sets, index=patients)

    df_conc_maxmin = pd.DataFrame(columns=param_sets, index=patients)
    df_conc_mean = pd.DataFrame(columns=param_sets, index=patients)
    df_conc_std = pd.DataFrame(columns=param_sets, index=patients)

    df_sum_maxmin = pd.DataFrame(columns=param_sets, index=patients)
    df_sum_mean = pd.DataFrame(columns=param_sets, index=patients)
    df_sum_std = pd.DataFrame(columns=param_sets, index=patients)

    for patient in patients:
        results = compute_metrics_for_patient(
            patient=patient,
            param_sets=param_sets,
            base_dir=base_dir,
        )
        for param_set in param_sets:
            m = results[param_set]
            via = m["viability"]
            conc = m["concentration"]
            sm = m["sum"]

            df_via_maxmin.loc[patient, param_set] = via["maxmin"]
            df_via_mean.loc[patient, param_set] = via["mean"]
            df_via_std.loc[patient, param_set] = via["std"]

            df_conc_maxmin.loc[patient, param_set] = conc["maxmin"]
            df_conc_mean.loc[patient, param_set] = conc["mean"]
            df_conc_std.loc[patient, param_set] = conc["std"]

            df_sum_maxmin.loc[patient, param_set] = sm["maxmin"]
            df_sum_mean.loc[patient, param_set] = sm["mean"]
            df_sum_std.loc[patient, param_set] = sm["std"]

    # écriture des TSV (comme dans RMSE.py)
    out_prefix = Path(output_prefix)
    df_via_maxmin.to_csv(f"{out_prefix}_via_max_min.tsv", sep="\t", index=True)
    df_via_mean.to_csv(f"{out_prefix}_via_mean.tsv", sep="\t", index=True)
    df_via_std.to_csv(f"{out_prefix}_via_stdev.tsv", sep="\t", index=True)

    df_conc_maxmin.to_csv(f"{out_prefix}_conc_max_min.tsv", sep="\t", index=True)
    df_conc_mean.to_csv(f"{out_prefix}_conc_mean.tsv", sep="\t", index=True)
    df_conc_std.to_csv(f"{out_prefix}_conc_stdev.tsv", sep="\t", index=True)

    df_sum_maxmin.to_csv(f"{out_prefix}_sum_max_min.tsv", sep="\t", index=True)
    df_sum_mean.to_csv(f"{out_prefix}_sum_mean.tsv", sep="\t", index=True)
    df_sum_std.to_csv(f"{out_prefix}_sum_stdev.tsv", sep="\t", index=True)

    logger.info("All NRMSE tables written.")
