# abm_pipeline/advanced_analysis/utils.py

from pathlib import Path
import pandas as pd


def load_pareto_dataframe(file_path: str) -> pd.DataFrame:
    """Charge un fichier Pareto_ABM_2D_patient.txt ou TSV."""
    path = Path(file_path)
    if path.suffix == ".txt":
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def ensure_dir(path: str | Path):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path
