# abm_pipeline/sensitivity/utils.py

from __future__ import annotations
from pathlib import Path
import pandas as pd
from typing import Dict, Tuple


def load_sensitivity_csv(csv_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Lecture d’un fichier BehaviorSpace de sensibilité.
    Retourne :
      - df_viability : DataFrame index=step, columns=run_number
      - df_concentration : pareil
    """
    viability_dict: Dict[int, Dict[int, float]] = {}
    remaining_dict: Dict[int, Dict[int, float]] = {}

    with open(csv_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    for line in lines[7:]:  # header BehaviorSpace
        row = line.replace('"', "").strip().split(",")
        if len(row) < 25:
            continue

        run_number = int(row[0])
        step = int(row[21])
        viability = float(row[23])
        remaining = float(row[24])

        viability_dict.setdefault(run_number, {})[step] = viability
        remaining_dict.setdefault(run_number, {})[step] = remaining

    df_viability = pd.DataFrame.from_dict(viability_dict)
    df_remaining = pd.DataFrame.from_dict(remaining_dict)

    return df_viability, df_remaining

