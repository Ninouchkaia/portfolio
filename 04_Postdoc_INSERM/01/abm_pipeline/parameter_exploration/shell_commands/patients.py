# abm_pipeline/parameter_exploration/shell_commands/patients.py

from pathlib import Path

from abm_pipeline.parameter_exploration.utils import (
    get_patient_ids,
    logger,
)


NETLOGO_HEADLESS = "/I/Program Files/NetLogo 6.1.0/netlogo-headless.bat"


def generate_patient_command_scripts(base_dir: str = ".") -> None:
    """
    Reprend la logique de list_of_commands.py :

    Pour chaque patient, crée un fichier patient_command_{PATIENT}.sh
    avec la chaîne d’appels :
      - aggregateData
      - remove_duplicates...
      - copy_for_git_paretoFront...
      - extract_param_sets_from_pareto...
      - parse_best_param_for_behavior_space_adapted...
      - make_behavior_space_experiment_file...
      - 3 simulations stocha_* via NetLogo
      - 3 plots plot_sim_vs_exp.py
    """
    patients = get_patient_ids()
    base_path = Path(base_dir)

    for patient in patients:
        sh_path = base_path / f"patient_command_{patient}.sh"
        logger.info(f"Writing {sh_path}")

        with sh_path.open("w", encoding="utf-8") as fw:
            fw.write("#!/bin/bash\n\n")
            fw.write(f'echo "This is a shell script for patient {patient}"\n\n')

            fw.write(
                f"python scripts/aggregateData.py {patient}/ABM_2D_{patient} "
                f"{patient}/outputs_ABM_2D_{patient}.txt\n\n"
            )
            fw.write(
                "python scripts/"
                f"remove_duplicates_generic_with_filtering_keeping_only_samples.py "
                f"{patient}/outputs_ABM_2D_{patient}.txt fitnessVia 50\n\n"
            )
            fw.write(
                "python scripts/copy_for_git_paretoFrontGenericStochastic.py "
                f"{patient}/outputs_ABM_2D_{patient}"
                "_duplicates_removed_filtered_only_samples_kept_50.0.txt "
                f"{patient}/pareto_ABM_2D_{patient}\n\n"
            )
            fw.write(
                "python scripts/extract_param_sets_from_pareto_adapted.py "
                f"{patient}/pareto_ABM_2D_{patient}.txt "
                f"{patient}/best_param_sets_ABM_2D_{patient}.tsv\n\n"
            )
            fw.write(
                "python scripts/parse_best_param_for_behavior_space_adapted.py "
                f"{patient}/best_param_sets_ABM_2D_{patient}.tsv "
                f"{patient}/netlogo_best_param_sets_ABM_2D_{patient}.txt\n\n"
            )
            fw.write(
                "python scripts/make_behavior_space_experiment_file.py "
                f"{patient}/best_param_sets_ABM_2D_{patient}.tsv "
                f"{patient}/experiment_file.xml\n\n"
            )

            for exp in ["stocha_best_via", "stocha_knee_point", "stocha_best_conc"]:
                fw.write(
                    f'"{NETLOGO_HEADLESS}" --model {patient}/ABM_2D_{patient}.nlogo '
                    f"--setup-file {patient}/experiment_file.xml "
                    f"--experiment {exp} "
                    f"--table {patient}/BehaviorSpace/{exp}.csv --threads 4\n\n"
                )

            for exp in ["stocha_best_via", "stocha_knee_point", "stocha_best_conc"]:
                fw.write(
                    f"python scripts/plot_sim_vs_exp.py {patient} {exp}\n"
                )

        # sur Git Bash : chmod +x patient_command_*.sh
