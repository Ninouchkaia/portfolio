# abm_pipeline/model_validation/metrics.py

from __future__ import annotations

import math
from typing import Sequence, Dict, Tuple

import numpy as np
import pandas as pd


def _to_array(x: Sequence[float]) -> np.ndarray:
    return np.array(list(x), dtype=float)


def rmse(y_true: Sequence[float], y_pred: Sequence[float]) -> float:
    """Root Mean Square Error."""
    y_t = _to_array(y_true)
    y_p = _to_array(y_pred)
    return float(np.sqrt(np.mean((y_p - y_t) ** 2)))


def nrmse_maxmin(y_true: Sequence[float], y_pred: Sequence[float]) -> float:
    """
    Normalized RMSE = RMSE / (max - min)
    (c'est ce que tu fais dans RMSE.py).
    """
    y_t = _to_array(y_true)
    denom = float(y_t.max() - y_t.min())
    if denom == 0:
        return float("nan")
    return rmse(y_t, y_pred) / denom


def nrmse_mean(y_true: Sequence[float], y_pred: Sequence[float]) -> float:
    """
    Normalized RMSE = RMSE / mean(y_true).
    """
    y_t = _to_array(y_true)
    denom = float(y_t.mean())
    if denom == 0:
        return float("nan")
    return rmse(y_t, y_pred) / denom


def nrmse_std(y_true: Sequence[float], y_pred: Sequence[float]) -> float:
    """
    Normalized RMSE = RMSE / std(y_true).
    """
    y_t = _to_array(y_true)
    denom = float(y_t.std(ddof=1))
    if denom == 0:
        return float("nan")
    return rmse(y_t, y_pred) / denom


def compute_viability_conc_nrmse(
    patient_data: pd.DataFrame,
    viability_sim: pd.Series,
    conc_sim: pd.Series,
    patient: str,
) -> Dict[str, Dict[str, float]]:
    """
    Calcule les 3 NRMSE pour la viabilité et la concentration pour un patient.

    patient_data :
        dataframe avec colonnes 'Day', '{patient}_viability', '{patient}_concentration'.
    viability_sim / conc_sim :
        séries avec les mêmes index (en heures ou en jours) que patient_data['Day'].
    """
    via_exp = patient_data[f"{patient}_viability"]
    conc_exp = patient_data[f"{patient}_concentration"]

    # s'assurer que les index sont alignés
    via_exp = via_exp.loc[viability_sim.index]
    conc_exp = conc_exp.loc[conc_sim.index]

    via = {
        "maxmin": nrmse_maxmin(via_exp, viability_sim),
        "mean": nrmse_mean(via_exp, viability_sim),
        "std": nrmse_std(via_exp, viability_sim),
    }
    conc = {
        "maxmin": nrmse_maxmin(conc_exp, conc_sim),
        "mean": nrmse_mean(conc_exp, conc_sim),
        "std": nrmse_std(conc_exp, conc_sim),
    }

    sums = {
        "maxmin": via["maxmin"] + conc["maxmin"],
        "mean": via["mean"] + conc["mean"],
        "std": via["std"] + conc["std"],
    }

    return {"viability": via, "concentration": conc, "sum": sums}
