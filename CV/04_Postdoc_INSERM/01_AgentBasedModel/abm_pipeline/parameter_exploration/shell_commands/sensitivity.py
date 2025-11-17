# abm_pipeline/parameter_exploration/shell_commands/sensitivity.py

from pathlib import Path

from abm_pipeline.parameter_exploration.utils import logger

from .patients import NETLOGO_HEADLESS


def generate_sensitivity_shell_scripts(
    model_path: str,
    xml_file: str,
    out_sh: str,
    exp_list: list[str],
) -> None:
    """
    Reprend la logique de make_shell_commands*.py :

    - model_path : chemin complet vers le modèle NetLogo .nlogo
    - xml_file   : fichier XML pour l'analyse de sensibilité
    - out_sh     : nom du .sh (sensitivity_experiments_class1.sh, etc.)
    - exp_list   : ["perturb-gui-apo-mov", ...]
    """
    sh_path = Path(out_sh)
    logger.info(f"Writing {sh_path}")

    with sh_path.open("w", encoding="utf-8") as fw:
        fw.write("#!/bin/bash\n\n")
        for exp in exp_list:
            fw.write(f'echo "This is a shell script for exp {exp}"\n')
            fw.write(
                f'"{NETLOGO_HEADLESS}" --model {model_path} '
                f"--setup-file {xml_file} --experiment {exp} "
                f"--table ABM_2D_sensitivity_{exp}.csv --threads 4\n\n"
            )
