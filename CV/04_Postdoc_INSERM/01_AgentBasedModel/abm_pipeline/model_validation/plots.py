# abm_pipeline/model_validation/plots.py

from __future__ import annotations

from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import mean_squared_error, r2_score

from abm_pipeline.model_validation.metrics import compute_viability_conc_nrmse
from abm_pipeline.parameter_exploration.utils import logger


def _load_patient_exp_data(patient: str, base_dir: str = ".") -> pd.DataFrame:
    """
    Lit les données expérimentales du patient.
    On suppose un fichier: {base_dir}/{patient}/{patient}.csv
    avec colonnes: 'Day', f'{patient}_viability', f'{patient}_concentration'
    """
    path = Path(base_dir) / patient / f"{patient}.csv"
    df = pd.read_csv(path, index_col=0)
    # dans le script RMSE.py, Day est multiplié par 24 → heures
    if "Day" in df.columns:
        df["Day"] = df["Day"].apply(lambda x: x * 24)
    return df


def _load_behaviorspace_csv(simu_file_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Lit le CSV BehaviorSpace et reconstruit deux dataframes :
    - df_viability_simu : viability par step × run
    - df_remaining_simu : concentration par step × run

    On colle à la logique du script historique :
      - step = colonne 21
      - viability = colonne 23
      - remainingCellRatio = colonne 24
      - on saute les 7 premières lignes (header BehaviorSpace).
    """
    viability_dict = {}
    remaining_dict = {}

    with open(simu_file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    for line in lines[7:]:
        line = line.replace('"', "").strip().split(",")
        if len(line) < 25:
            continue
        run_number = int(line[0])
        step = int(line[21])
        viability = float(line[23])
        remaining = float(line[24])

        viability_dict.setdefault(run_number, {})[step] = viability
        remaining_dict.setdefault(run_number, {})[step] = remaining

    df_viability = pd.DataFrame.from_dict(viability_dict)
    df_remaining = pd.DataFrame.from_dict(remaining_dict)
    return df_viability, df_remaining


def plot_sim_vs_exp_with_scores(
    patient: str,
    param_set: str,
    base_dir: str = ".",
    behaviorspace_dir: str | None = None,
    save_dir: str | None = None,
) -> None:
    """
    Reprend la logique de plot_sim_vs_exp_with_scores.py :

    - patient: ex 'CAS1802'
    - param_set: ex 'stocha_best_via', 'stocha_knee_point', ...
    - base_dir: racine du projet (dossiers {patient}/...).
    - behaviorspace_dir: sous-dossier des CSV BehaviorSpace (par défaut {patient}/BehaviorSpace).
    - save_dir: dossier où sauver la figure (par défaut {patient}/figures).

    La figure contient:
      - courbe simulée + points exp pour viabilité
      - courbe simulée + points exp pour concentration
      - NRMSE et r² affichés dans le titre de chaque subplot.
    """
    base_path = Path(base_dir)
    if behaviorspace_dir is None:
        behaviorspace_path = base_path / patient / "BehaviorSpace"
    else:
        behaviorspace_path = base_path / behaviorspace_dir

    simu_file = behaviorspace_path / f"{param_set}.csv"
    logger.info(f"Loading simulations from {simu_file}")

    patient_data = _load_patient_exp_data(patient, base_dir=base_dir)
    df_viability_simu, df_remaining_simu = _load_behaviorspace_csv(str(simu_file))

    # On garde uniquement les time points présents dans les données exp :
    # dans plot_sim_vs_exp_with_scores, la colonne 'Day' est directement utilisée.
    time_points = list(patient_data["Day"])
    filtered_viability_simu = df_viability_simu[
        df_viability_simu.index.isin(time_points)
    ]
    filtered_remaining_simu = df_remaining_simu[
        df_remaining_simu.index.isin(time_points)
    ]

    # Moyennes (sur les runs)
    via_simu_mean = filtered_viability_simu.mean(axis=1)
    conc_simu_mean = filtered_remaining_simu.mean(axis=1)

    # Recalage index des séries sur les time_points
    via_simu_mean.index = time_points
    conc_simu_mean.index = time_points

    # NRMSE (max-min, mean, std) + somme
    metrics = compute_viability_conc_nrmse(
        patient_data=patient_data,
        viability_sim=via_simu_mean,
        conc_sim=conc_simu_mean,
        patient=patient,
    )

    # NRMSE supplémentaire type sklearn (comme dans plot_sim_vs_exp_with_scores.py)
    y_via_exp = patient_data[f"{patient}_viability"]
    y_conc_exp = patient_data[f"{patient}_concentration"]
    nrms_via_sklearn = mean_squared_error(
        y_via_exp, via_simu_mean, squared=False
    ) / (y_via_exp.max() - y_via_exp.min())
    nrms_conc_sklearn = mean_squared_error(
        y_conc_exp, conc_simu_mean, squared=False
    ) / (y_conc_exp.max() - y_conc_exp.min())
    nrms_via_sklearn = round(float(nrms_via_sklearn), 2)
    nrms_conc_sklearn = round(float(nrms_conc_sklearn), 2)

    # r²
    r2_via = round(float(r2_score(y_via_exp, via_simu_mean)), 2)
    r2_conc = round(float(r2_score(y_conc_exp, conc_simu_mean)), 2)

    # --- Figure -------------------------------------------------------------
    fig, axes = plt.subplots(2, 1, figsize=(6, 8), sharex=True)

    # Viabilité
    axes[0].plot(
        via_simu_mean.index / 24.0, via_simu_mean, label="Simulation", linewidth=2
    )
    axes[0].plot(
        patient_data["Day"] / 24.0,
        patient_data[f"{patient}_viability"],
        label="Experimental",
        color="black",
        marker="o",
        linestyle="--",
    )
    axes[0].set_ylabel("Viability (%)")
    axes[0].set_xticks(range(0, 16, 2))
    axes[0].set_title(
        f"{patient} – {param_set}\n"
        f"NRMSE_via={nrms_via_sklearn}, R²_via={r2_via}"
    )
    axes[0].legend()

    # Concentration
    axes[1].plot(
        conc_simu_mean.index / 24.0,
        conc_simu_mean,
        label="Simulation",
        linewidth=2,
    )
    axes[1].plot(
        patient_data["Day"] / 24.0,
        patient_data[f"{patient}_concentration"],
        label="Experimental",
        color="black",
        marker="o",
        linestyle="--",
    )
    axes[1].set_xlabel("Day")
    axes[1].set_ylabel("Concentration Ratio (%)")
    axes[1].set_xticks(range(0, 16, 2))
    axes[1].set_title(
        f"NRMSE_conc={nrms_conc_sklearn}, R²_conc={r2_conc}"
    )
    axes[1].legend()

    fig.tight_layout()

    # savefig
    if save_dir is None:
        save_path = base_path / patient / "figures"
    else:
        save_path = Path(save_dir)
    save_path.mkdir(parents=True, exist_ok=True)
    out_file = save_path / f"{patient}_{param_set}_sim_vs_exp.png"
    fig.savefig(out_file, dpi=300, bbox_inches="tight")

    logger.info(f"Figure saved to {out_file}")
