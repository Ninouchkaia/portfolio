# abm_pipeline/parameter_exploration/nsga2_analysis/extract_best_sets.py

from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

from abm_pipeline.parameter_exploration.utils import (
    logger,
    get_patient_ids,
)


def extract_best_param_sets(pareto_file: str) -> pd.DataFrame:
    """
    Depuis un fichier Pareto (CSV), renvoie un DataFrame avec 3 lignes :
      - best_via_set
      - knee_point_set
      - best_conc_set
    """
    pf = pd.read_csv(pareto_file)
    pf["distances"] = np.sqrt(
        pf["delta_fitness_via"] ** 2 + pf["delta_fitness_conc"] ** 2
    )

    best_via = pf[pf.delta_fitness_via == pf.delta_fitness_via.min()]
    best_conc = pf[pf.delta_fitness_conc == pf.delta_fitness_conc.min()]
    knee_point = pf[pf.distances == pf.distances.min()]

    best_params = pd.concat([best_via, knee_point, best_conc])
    best_params = best_params.drop(columns=["distances"])
    best_params.insert(
        loc=0,
        column="set",
        value=["best_via_set", "knee_point_set", "best_conc_set"],
    )
    return best_params


def extract_best_param_sets_to_file(pareto_file: str, output_tsv: str | None = None):
    pf_path = Path(pareto_file)
    if output_tsv is None:
        output_tsv = f"best_param_sets_{pf_path.stem}.tsv"

    df = extract_best_param_sets(pareto_file)
    df.to_csv(output_tsv, sep="\t", index=False)
    logger.info(f"Best parameter sets written to {output_tsv}")


def make_best_sets_all_patients(
    set_name: Literal["best_via_set", "best_conc_set", "knee_point_set"],
    output_file: str,
) -> None:
    """
    Concatène, pour tous les patients, la ligne correspondante de
    best_param_sets_ABM_2D_{patient}.tsv dans un seul TSV.
    """
    patients = get_patient_ids()
    all_sets = pd.DataFrame(columns=patients)

    for patient in patients:
        path = Path(f"{patient}/best_param_sets_ABM_2D_{patient}.tsv")
        df = pd.read_csv(path, sep="\t", header=0, index_col=0)
        best_set = df.loc[set_name]
        all_sets[patient] = best_set.T

    all_sets.to_csv(output_file, index=True, sep="\t")
    logger.info(f"All patients {set_name} written to {output_file}")
