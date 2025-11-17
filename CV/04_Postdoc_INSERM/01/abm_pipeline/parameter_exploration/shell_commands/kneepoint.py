# abm_pipeline/parameter_exploration/shell_commands/kneepoint.py

from pathlib import Path
from typing import Literal

from abm_pipeline.parameter_exploration.utils import get_patient_ids, logger
from .patients import NETLOGO_HEADLESS


def generate_kneepoint_scripts(
    setup_xml: str,
    class_label: Literal["0", "1_class1", "1_class2", "2"],
    output_dir: str = ".",
) -> None:
    """
    Reprend la logique de:
      - kneepoint0_by_patient_simu_shell.py
      - kneepoint1_class1_by_patient_simu_shell.py
      - kneepoint1_class2_by_patient_simu_shell.py
      - kneepoint2_by_patient_simu_shell.py

    Un script .sh par patient, avec l'appel NetLogo correspondant.
    """
    patients = get_patient_ids()
    out_path = Path(output_dir)

    for patient in patients:
        sh_name = f"kneepoint{class_label}_by_patient_simu_{patient}.sh"
        sh_path = out_path / sh_name
        logger.info(f"Writing {sh_path}")

        exp_name = f"kneepoint{class_label}_simu_{patient}"
        behavior_csv = f"{patient}/BehaviorSpace/{exp_name}.csv"

        with sh_path.open("w", encoding="utf-8") as fw:
            fw.write("#!/bin/bash\n\n")
            fw.write(
                f'"{NETLOGO_HEADLESS}" --model {patient}/ABM_2D_{patient}.nlogo '
                f"--setup-file {setup_xml} "
                f"--experiment {exp_name} "
                f"--table {behavior_csv} --threads 4\n"
            )
