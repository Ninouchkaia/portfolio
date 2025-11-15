# abm_pipeline/parameter_exploration/initial_ranges/aggregate_data.py

import os
from pathlib import Path

from abm_pipeline.parameter_exploration.utils import logger, debug


def aggregate_population_files(rootdir: str, output_file: str) -> None:
    """
    Concatène tous les fichiers CSV du dossier `rootdir` dans `output_file`.

    - 'population1.csv' est recopié avec son header.
    - tous les autres 'populationX.csv' sont recopiés en sautant la première ligne.

    Comportement équivalent à l'ancien aggregateData.py.
    """
    rootdir_path = Path(rootdir)
    output_path = Path(output_file)
    logger.info(f"Aggregating population files from {rootdir_path} → {output_path}")

    with output_path.open("w", encoding="utf-8") as file_write:
        for _, _, files in os.walk(rootdir_path):
            for fname in sorted(files):
                path_to_file = rootdir_path / fname
                debug(f"Reading {path_to_file}")
                with path_to_file.open("r", encoding="utf-8") as file_read:
                    lines = file_read.readlines()
                    if fname == "population1.csv":
                        file_write.writelines(lines)
                    else:
                        file_write.writelines(lines[1:])
