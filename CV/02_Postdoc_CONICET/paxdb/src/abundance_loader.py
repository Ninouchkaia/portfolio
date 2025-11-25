# paxdb/src/abundance_loader.py

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)


def load_abundance_table(
    path: Path,
    id_col: str = "protein_id",
    abundance_col: str = "abundance",
) -> pd.DataFrame:
    """
    Load a PaxDB-like protein abundance table.

    Parameters
    ----------
    path : Path
        TSV/CSV file with protein IDs and abundance values.
    id_col : str
        Name of the column containing protein IDs.
    abundance_col : str
        Name of the column containing abundance values.

    Returns
    -------
    DataFrame
        Columns: [id_col, abundance_col], filtered to non-null abundances.
    """
    logger.info("Loading abundance table from %s", path)
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    if id_col not in df.columns or abundance_col not in df.columns:
        raise ValueError(
            f"Required columns '{id_col}' and '{abundance_col}' not found in {path}"
        )

    df = df[[id_col, abundance_col]].copy()
    df[abundance_col] = pd.to_numeric(df[abundance_col], errors="coerce")
    df = df.dropna(subset=[abundance_col])

    logger.info("Loaded %d entries from %s", len(df), path)
    return df


def load_species_metadata(path: Path) -> pd.DataFrame:
    """
    Load species metadata TSV describing all species to analyze.

    Expected columns:
    - species_id
    - fasta_path
    - abundance_path
    - mapping_path (not used in the minimal example, but kept for compatibility)
    """
    logger.info("Loading species metadata from %s", path)
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str)

    required_cols = ["species_id", "fasta_path", "abundance_path", "mapping_path"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"Species metadata file {path} is missing required columns: {missing}"
        )

    return df
