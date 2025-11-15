# abm_pipeline/parameter_exploration/instantiate_models/behavior_space_files.py

from pathlib import Path

from abm_pipeline.parameter_exploration.utils import logger


def _read_patient_dict(patient_dict_file: str, patient: str):
    mono_init = None
    apo_init = None
    with open(patient_dict_file, "r", encoding="utf-8") as f:
        data = f.readlines()
    for line in data[1:]:
        line = line.replace("\n", "").split(" ")
        patient_name = line[0]
        if patient_name == patient:
            mono_init = line[1]
            apo_init = line[2]
            break
    if mono_init is None or apo_init is None:
        raise ValueError(f"Patient {patient} non trouvé dans {patient_dict_file}")
    return mono_init, apo_init


def make_behavior_space_file_for_patient(
    best_sets_tsv: str,
    output_file: str,
    patient_dict_file: str = "patient_dict.txt",
) -> None:
    """
    Version "adapted" : ajoute les paramètres mono/apo init spécifiques au patient
    et écrit un fichier de paramètres lisible par NetLogo.
    """
    tsv_path = Path(best_sets_tsv)
    patient = tsv_path.parts[0]  # on suppose "PATIENT/..."
    mono_init, apo_init = _read_patient_dict(patient_dict_file, patient)

    out_path = Path(output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with tsv_path.open("r", encoding="utf-8") as file_read, out_path.open(
        "w", encoding="utf-8"
    ) as file_write:
        data = file_read.readlines()
        n = 1
        for line in data[1:]:
            cols = line.replace("\n", "").split("\t")
            (
                _set_name,
                _delta_via,
                _delta_conc,
                gui_apo_mov,
                gui_need_sig_mov,
                gui_layers,
                gui_alpha,
                gui_mono_phago_eff,
                gui_NLC_phago_eff,
                gui_M_phago_eff,
                gui_M_kill_eff,
                gui_cll_sens_dist,
                gui_mono_sens_dist,
                gui_nlc_sens_dist,
                gui_macro_sens_dist,
                gui_nlc_threshold,
                gui_sig_init,
                gui_sig_init_std,
                gui_diff_mean,
                gui_diff_std,
                gui_life_init_gamma,
                gui_alpha_distrib,
            ) = cols

            if n == 1:
                file_write.write("Best_via_set\n")
            elif n == 2:
                file_write.write("Knee_point_set\n")
            elif n == 3:
                file_write.write("Best_conc_set\n")

            file_write.write(f'["gui-prop-mono-init" {mono_init}]\n')
            file_write.write(f'["gui-prop-apo-init" {apo_init}]\n')
            file_write.write(f'["gui-apo-mov" {gui_apo_mov}]\n')
            file_write.write(f'["gui-need-sig-mov" {gui_need_sig_mov}]\n')
            file_write.write(f'["gui-layers" {gui_layers}]\n')
            file_write.write(f'["gui-alpha" {gui_alpha}]\n')
            file_write.write(f'["gui-mono-phago-eff" {gui_mono_phago_eff}]\n')
            file_write.write(f'["gui-NLC-phago-eff" {gui_NLC_phago_eff}]\n')
            file_write.write(f'["gui-M-phago-eff" {gui_M_phago_eff}]\n')
            file_write.write(f'["gui-M-kill-eff" {gui_M_kill_eff}]\n')
            file_write.write(f'["gui-cll-sens-dist" {gui_cll_sens_dist}]\n')
            file_write.write(f'["gui-mono-sens-dist" {gui_mono_sens_dist}]\n')
            file_write.write(f'["gui-nlc-sens-dist" {gui_nlc_sens_dist}]\n')
            file_write.write(f'["gui-macro-sens-dist" {gui_macro_sens_dist}]\n')
            file_write.write(f'["gui-nlc-threshold" {gui_nlc_threshold}]\n')
            file_write.write(f'["gui-sig-init" {gui_sig_init}]\n')
            file_write.write(f'["gui-sig-init-std" {gui_sig_init_std}]\n')
            file_write.write(f'["gui-diff-mean" {gui_diff_mean}]\n')
            file_write.write(f'["gui-diff-std" {gui_diff_std}]\n')
            file_write.write(f'["gui-life-init-gamma" {gui_life_init_gamma}]\n')
            file_write.write(f'["gui-alpha-distrib" {gui_alpha_distrib}]\n\n')
            n += 1

    logger.info(f"BehaviorSpace param file written to {out_path}")


def make_behavior_space_file_generic(best_sets_tsv: str, output_file: str) -> None:
    """
    Version sans mono/apo init spécifiques (ancienne version générique).
    """
    tsv_path = Path(best_sets_tsv)
    out_path = Path(output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with tsv_path.open("r", encoding="utf-8") as file_read, out_path.open(
        "w", encoding="utf-8"
    ) as file_write:
        data = file_read.readlines()
        n = 1
        for line in data[1:]:
            cols = line.replace("\n", "").split("\t")
            (
                _set_name,
                _delta_via,
                _delta_conc,
                gui_apo_mov,
                gui_need_sig_mov,
                gui_layers,
                gui_alpha,
                gui_mono_phago_eff,
                gui_NLC_phago_eff,
                gui_M_phago_eff,
                gui_M_kill_eff,
                gui_cll_sens_dist,
                gui_mono_sens_dist,
                gui_nlc_sens_dist,
                gui_macro_sens_dist,
                gui_nlc_threshold,
                gui_sig_init,
                gui_sig_init_std,
                gui_diff_mean,
                gui_diff_std,
                gui_life_init_gamma,
                gui_alpha_distrib,
            ) = cols

            if n == 1:
                file_write.write("Best_via_set\n")
            elif n == 2:
                file_write.write("Knee_point_set\n")
            elif n == 3:
                file_write.write("Best_conc_set\n")

            file_write.write(f'["gui-apo-mov" {gui_apo_mov}]\n')
            file_write.write(f'["gui-need-sig-mov" {gui_need_sig_mov}]\n')
            file_write.write(f'["gui-layers" {gui_layers}]\n')
            file_write.write(f'["gui-alpha" {gui_alpha}]\n')
            file_write.write(f'["gui-mono-phago-eff" {gui_mono_phago_eff}]\n')
            file_write.write(f'["gui-NLC-phago-eff" {gui_NLC_phago_eff}]\n')
            file_write.write(f'["gui-M-phago-eff" {gui_M_phago_eff}]\n')
            file_write.write(f'["gui-M-kill-eff" {gui_M_kill_eff}]\n')
            file_write.write(f'["gui-cll-sens-dist" {gui_cll_sens_dist}]\n')
            file_write.write(f'["gui-mono-sens-dist" {gui_mono_sens_dist}]\n')
            file_write.write(f'["gui-nlc-sens-dist" {gui_nlc_sens_dist}]\n')
            file_write.write(f'["gui-macro-sens-dist" {gui_macro_sens_dist}]\n')
            file_write.write(f'["gui-nlc-threshold" {gui_nlc_threshold}]\n')
            file_write.write(f'["gui-sig-init" {gui_sig_init}]\n')
            file_write.write(f'["gui-sig-init-std" {gui_sig_init_std}]\n')
            file_write.write(f'["gui-diff-mean" {gui_diff_mean}]\n')
            file_write.write(f'["gui-diff-std" {gui_diff_std}]\n')
            file_write.write(f'["gui-life-init-gamma" {gui_life_init_gamma}]\n')
            file_write.write(f'["gui-alpha-distrib" {gui_alpha_distrib}]\n\n')
            n += 1

    logger.info(f"Generic BehaviorSpace file written to {out_path}")
