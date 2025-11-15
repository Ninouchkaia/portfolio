# abm_pipeline/parameter_exploration/shell_commands/averaged.py

from pathlib import Path
from abm_pipeline.parameter_exploration.utils import get_patient_ids, logger
from .patients import NETLOGO_HEADLESS


def generate_averaged_class_scripts(
    xml_file: str,
    class_label: str,
    output_dir: str = ".",
) -> None:
    """
    Équivalent à averaged_simu_shell.py :
    crée un patient_command_averaged_classX_simu.sh
    qui lance 'averaged_simu_{PATIENT}' pour tous les patients.
    """
    patients = get_patient_ids()
    out_path = Path(output_dir) / f"patient_command_averaged_{class_label}_simu.sh"
    logger.info(f"Writing {out_path}")

    with out_path.open("w", encoding="utf-8") as fw:
        fw.write("#!/bin/bash\n\n")
        for patient in patients:
            fw.write(
                f'"{NETLOGO_HEADLESS}" --model {patient}/ABM_2D_{patient}.nlogo '
                f"--setup-file {xml_file} "
                ff"--experiment averaged_simu_{patient} "
                f"--table {patient}/BehaviorSpace/averaged_{class_label}_simu_{patient}.csv "
                "--threads 4\n\n"
            )
