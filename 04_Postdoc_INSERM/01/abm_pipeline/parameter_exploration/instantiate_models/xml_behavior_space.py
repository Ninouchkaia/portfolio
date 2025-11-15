# abm_pipeline/parameter_exploration/instantiate_models/xml_behavior_space.py

from pathlib import Path
from typing import Dict

import pandas as pd

from abm_pipeline.parameter_exploration.utils import logger


def _load_knee_point_dict(best_sets_tsv: str) -> Dict[str, float]:
    """
    Lit le TSV best_param_sets_* et retourne un dict param_name -> valeur
    en prenant la ligne 'knee_point_set'.
    """
    df = pd.read_csv(best_sets_tsv, sep="\t")
    knee_row = df[df["set"] == "knee_point_set"].iloc[0]
    # les colonnes de param commencent après (set, delta_fitness_via, delta_fitness_conc)
    param_cols = knee_row.index[3:]
    knee_dict = {name: float(knee_row[name]) for name in param_cols}
    return knee_dict


def make_sensitivity_experiment_xml(
    best_sets_tsv: str,
    output_xml: str,
    gui_prop_mono_init: float,
    gui_prop_apo_init: float,
    ranges_dict: Dict[str, list],
    repetitions: int = 3,
    time_limit: int = 312,
) -> None:
    """
    Génère un fichier XML de type BehaviorSpace pour l'analyse de sensibilité.

    Équivalent conceptuel à tes scripts:
      - make_behavior_space_experiment_file_class1.py
      - make_behavior_space_experiment_file_class2.py

    Paramètres
    ----------
    best_sets_tsv : str
        Fichier TSV best_param_sets_... (contenant la ligne 'knee_point_set').
    output_xml : str
        Fichier XML de sortie.
    gui_prop_mono_init : float
        Proportion initiale de monocytes (valeur moyenne utilisée pour la classe).
    gui_prop_apo_init : float
        Proportion initiale apoptotique (valeur moyenne utilisée pour la classe).
    ranges_dict : dict
        Dictionnaire {param_name: [first, step, last]} pour la perturbation.
        C’est ici que tu recopies les valeurs de tes scripts d’origine.
    repetitions : int
        Nombre de répétitions BehaviorSpace.
    time_limit : int
        Nombre de pas de temps.
    """
    knee_point_dict = _load_knee_point_dict(best_sets_tsv)
    out_path = Path(output_xml)

    logger.info(f"Writing sensitivity BehaviorSpace XML → {out_path}")
    with out_path.open("w", encoding="utf-8") as fw:
        fw.write("<experiments>\n")

        for param_name in ranges_dict:
            remaining_params = list(ranges_dict.keys())
            remaining_params.remove(param_name)

            exp_name = f"perturb-{param_name}"
            fw.write(
                f'  <experiment name="{exp_name}" repetitions="{repetitions}" '
                'runMetricsEveryStep="true">\n'
            )
            fw.write("    <setup>setup</setup>\n")
            fw.write("    <go>go</go>\n")
            fw.write(f'    <timeLimit steps="{time_limit}"/>\n')
            fw.write("    <metric>getSeed</metric>\n")
            fw.write("    <metric>getViability</metric>\n")
            fw.write("    <metric>getRemainingCellRatio</metric>\n\n")

            # paramètres globaux mono/apo init
            fw.write(
                '    <enumeratedValueSet variable="gui-prop-mono-init">'
                f'<value value="{gui_prop_mono_init}"/></enumeratedValueSet>\n'
            )
            fw.write(
                '    <enumeratedValueSet variable="gui-prop-apo-init">'
                f'<value value="{gui_prop_apo_init}"/></enumeratedValueSet>\n'
            )

            # les autres paramètres sont fixés sur la valeur du knee_point
            for other_param in remaining_params:
                val = knee_point_dict[other_param]
                fw.write(
                    f'    <enumeratedValueSet variable="{other_param}">'
                    f'<value value="{val}"/></enumeratedValueSet>\n'
                )

            first, step, last = ranges_dict[param_name]
            fw.write(
                f'    <steppedValueSet variable="{param_name}" '
                f'first="{first}" step="{step}" last="{last}"/>\n'
            )
            fw.write("  </experiment>\n\n")

        fw.write("  </experiments>\n")

    logger.info("Sensitivity XML generated.")
