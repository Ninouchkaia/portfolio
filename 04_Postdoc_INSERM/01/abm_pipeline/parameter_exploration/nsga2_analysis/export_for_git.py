# abm_pipeline/parameter_exploration/nsga2_analysis/export_for_git.py

from pathlib import Path

from abm_pipeline.parameter_exploration.utils import (
    build_patient_mono_dict,
    get_patient_ids,
    logger,
)


def export_patient_data_for_git(output_root: str = "patient_data_for_git") -> None:
    """
    Reproduit le comportement de ton script make_files_for_git.py :

    - crée patient_data_for_git/patient_n
    - copie pareto_ABM_2D_{ID}.txt en pareto_front_patient_n.txt
    - met un header "propre" dans le fichier cible
    """
    output_root_path = Path(output_root)
    output_root_path.mkdir(parents=True, exist_ok=True)

    patients = get_patient_ids()
    patient_dict = build_patient_mono_dict()

    for idx, patient in enumerate(patients, start=1):
        patient_dir = output_root_path / f"patient_{idx}"
        patient_dir.mkdir(exist_ok=True)

        logger.info(f"Exporting pareto for {patient} → patient_{idx}")
        src_file = Path(f"{patient}/pareto_ABM_2D_{patient}.txt")
        dest_file = patient_dir / f"pareto_front_patient_{idx}.txt"

        with src_file.open("r", encoding="utf-8") as file_read:
            data = file_read.readlines()

        data[0] = (
            "delta_fitness_via,delta_fitness_conc, apoCellsMovementProba,"
            " needSigCellsMvtProba, layersAroundNLC, antiApoBoost, monoPhagoEff,"
            " nlcPhagoEff, macroPhagoEff, macroKillEff, cllSensingDistance,"
            " monocyteSensingDistance, nlcSensingDistance, macrophageSensingDistance,"
            " nlcThreshold, signalInitMean, signalInitStd, monoDiffThreshold,"
            " monoDiffTimeStd, gammaLifeInitRate, gammaLifeInitShape\n"
        )

        with dest_file.open("w", encoding="utf-8") as file_write:
            file_write.writelines(data)
