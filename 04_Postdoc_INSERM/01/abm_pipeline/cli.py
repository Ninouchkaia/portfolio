# abm_pipeline/cli.py

import typer

from abm_pipeline.parameter_exploration.initial_ranges.aggregate_data import (
    aggregate_population_files,
)
from abm_pipeline.parameter_exploration.nsga2_analysis import (
    build_pareto_front,
    extract_best_param_sets_to_file,
    make_best_sets_all_patients,
    export_patient_data_for_git,
)
from abm_pipeline.parameter_exploration.instantiate_models import (
    make_behavior_space_file_for_patient,
)

app = typer.Typer(help="ABM parameter exploration pipeline")


@app.command()
def aggregate(
    rootdir: str = typer.Argument(..., help="Dossier contenant populationX.csv"),
    output_file: str = typer.Argument(..., help="Fichier de sortie agrégé"),
):
    """Concatène les fichiers populationX.csv (ancienne aggregateData.py)."""
    aggregate_population_files(rootdir, output_file)


@app.command()
def pareto(
    input_file: str = typer.Argument(..., help="CSV issu d'OpenMOLE"),
    output_prefix: str = typer.Argument(
        ..., help="Préfixe pour pareto_*.png et pareto_*.txt"
    ),
):
    """Construit le Pareto front + figure (copy_for_git_paretoFrontGenericStochastic)."""
    build_pareto_front(input_file, output_prefix)


@app.command()
def bestsets(
    pareto_file: str = typer.Argument(..., help="Fichier Pareto .csv ou .txt"),
    output_tsv: str = typer.Argument(
        ..., help="Fichier TSV best_param_sets_*.tsv"
    ),
):
    """Extrait best_via, knee_point, best_conc (extract_param_sets_from_pareto)."""
    extract_best_param_sets_to_file(pareto_file, output_tsv)


@app.command()
def bestsets_all(
    set_name: str = typer.Argument(
        ...,
        help="best_via_set | knee_point_set | best_conc_set",
    ),
    output_file: str = typer.Argument(..., help="TSV de sortie"),
):
    """Concatène les sets pour tous les patients."""
    make_best_sets_all_patients(set_name, output_file)


@app.command()
def export_git(
    output_root: str = typer.Argument(
        "patient_data_for_git", help="Dossier de sortie"
    ),
):
    """Prépare les fichiers pareto_* pour dépôt Git (make_files_for_git.py)."""
    export_patient_data_for_git(output_root)


@app.command()
def make_behaviorspace(
    best_sets_tsv: str = typer.Argument(
        ..., help="TSV best_param_sets_ABM_2D_patient.tsv"
    ),
    output_file: str = typer.Argument(
        ..., help="Fichier de paramètres pour NetLogo"
    ),
    patient_dict: str = typer.Option(
        "patient_dict.txt", help="Fichier patient_dict.txt"
    ),
):
    """Génère le fichier de paramètres NetLogo pour un patient."""
    make_behavior_space_file_for_patient(best_sets_tsv, output_file, patient_dict)


# dans abm_pipeline/cli.py (suite)

from abm_pipeline.parameter_exploration.shell_commands import (
    generate_patient_command_scripts,
    generate_kneepoint_scripts,
    generate_averaged_class_scripts,
    generate_sensitivity_shell_scripts,
)
from abm_pipeline.parameter_exploration.instantiate_models.xml_behavior_space import (
    make_sensitivity_experiment_xml,
)


@app.command()
def patient_shell(
    base_dir: str = typer.Argument(".", help="Dossier où écrire les scripts .sh"),
):
    """Génère patient_command_{PATIENT}.sh pour tous les patients."""
    generate_patient_command_scripts(base_dir)


@app.command()
def kneepoint_shell(
    setup_xml: str = typer.Argument(..., help="Fichier XML de setup (kneepoint*.xml)"),
    label: str = typer.Argument("0", help="0, 1_class1, 1_class2, 2"),
    out_dir: str = typer.Argument(".", help="Dossier de sortie"),
):
    """Génère un script .sh par patient pour les simulations kneepoint."""
    generate_kneepoint_scripts(setup_xml, label, out_dir)


@app.command()
def averaged_shell(
    xml_file: str = typer.Argument(..., help="XML class1_averaged.xml ou autre"),
    class_label: str = typer.Argument("class1", help="class1 ou class2"),
    out_dir: str = typer.Argument(".", help="Dossier de sortie"),
):
    """Génère un script .sh pour lancer les simulations averaged_* sur tous les patients."""
    generate_averaged_class_scripts(xml_file, class_label, out_dir)


@app.command()
def sensi_xml(
    best_sets_tsv: str = typer.Argument(...),
    output_xml: str = typer.Argument(...),
    mono_init: float = typer.Argument(..., help="gui-prop-mono-init"),
    apo_init: float = typer.Argument(..., help="gui-prop-apo-init"),
):
    """
    Génére un fichier XML BehaviorSpace pour l'analyse de sensibilité.
    (les ranges_dict sont à coder dans un module de config séparé).
    """
    from abm_pipeline.sensitivity_config import RANGES_DICT_CLASS1  # à créer

    make_sensitivity_experiment_xml(
        best_sets_tsv,
        output_xml,
        gui_prop_mono_init=mono_init,
        gui_prop_apo_init=apo_init,
        ranges_dict=RANGES_DICT_CLASS1,
    )


@app.command()
def sensi_shell(
    model_path: str = typer.Argument(..., help="Chemin complet vers ABM_2D_*.nlogo"),
    xml_file: str = typer.Argument(..., help="Fichier XML de sensibilité"),
    out_sh: str = typer.Argument("sensitivity_experiments.sh"),
):
    """Génère le script .sh de lancement des expériences de sensibilité."""
    from abm_pipeline.sensitivity_config import EXP_LIST  # ["perturb-gui-apo-mov", ...]

    generate_sensitivity_shell_scripts(
        model_path=model_path,
        xml_file=xml_file,
        out_sh=out_sh,
        exp_list=EXP_LIST,
    )


def main():
    app()


if __name__ == "__main__":
    main()
