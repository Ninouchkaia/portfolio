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
        Name of the column containing abundance scores.

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
    df = df[df[abundance_col] > 0]

    logger.info(
        "Loaded %d proteins with positive abundance from %s", len(df), path.name
    )
    return df


def load_mapping_table(
    path: Optional[Path],
    string_col: str = "string_id",
    uniprot_col: str = "uniprot_id",
) -> Optional[pd.DataFrame]:
    """
    Load an optional STRING→UniProt mapping table.

    Parameters
    ----------
    path : Path or None
        If None or empty, returns None.
    string_col : str
        Column name with STRING / PaxDB protein ids.
    uniprot_col : str
        Column name with UniProt ids.

    Returns
    -------
    DataFrame or None
        If path is provided: DataFrame with [string_col, uniprot_col].
        Otherwise: None.
    """
    if path is None:
        logger.info("No mapping table provided.")
        return None

    if not path.exists():
        raise FileNotFoundError(f"Mapping file not found: {path}")

    logger.info("Loading STRING→UniProt mapping from %s", path)
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    if string_col not in df.columns or uniprot_col not in df.columns:
        raise ValueError(
            f"Required columns '{string_col}' and '{uniprot_col}' not found in {path}"
        )

    df = df[[string_col, uniprot_col]].dropna()
    logger.info("Loaded %d mapping rows from %s", len(df), path.name)
    return df


def map_abundances_to_fasta_ids(
    abund_df: pd.DataFrame,
    fasta_ids: set[str],
    mapping_df: Optional[pd.DataFrame] = None,
    abundance_id_col: str = "protein_id",
    abundance_col: str = "abundance",
    string_col: str = "string_id",
    uniprot_col: str = "uniprot_id",
) -> pd.DataFrame:
    """
    Map abundance entries to FASTA sequence identifiers.

    Two cases:
    1. No mapping_df: we assume abundance_df[abundance_id_col] already matches
       FASTA IDs directly; we filter on intersection.
    2. mapping_df provided: we join on STRING/PaxDB IDs and map to UniProt IDs.

    Parameters
    ----------
    abund_df : DataFrame
        Abundance table.
    fasta_ids : set of str
        Identifiers present in the FASTA proteome (seq_id from fasta_parser).
    mapping_df : DataFrame, optional
        Mapping between STRING/PaxDB ids and UniProt ids.
    abundance_id_col : str
    abundance_col : str
    string_col : str
    uniprot_col : str

    Returns
    -------
    DataFrame
        Columns: ['seq_id', 'abundance']
    """
    if mapping_df is None:
        logger.info(
            "Mapping abundances directly to FASTA IDs (no external mapping table)."
        )
        df = abund_df.rename(columns={abundance_id_col: "seq_id"})[["seq_id", abundance_col]]
        df = df[df["seq_id"].isin(fasta_ids)].copy()
        df = df.rename(columns={abundance_col: "abundance"})
        logger.info(
            "Mapped %d abundance entries to %d FASTA ids (direct matching).",
            len(df),
            len(fasta_ids),
        )
        return df

    logger.info("Mapping abundances via STRING→UniProt mapping table.")
    merged = abund_df.merge(
        mapping_df,
        left_on=abundance_id_col,
        right_on=string_col,
        how="inner",
    )

    merged = merged.rename(columns={uniprot_col: "seq_id"})
    merged = merged[["seq_id", abundance_col]].copy()
    merged = merged[merged["seq_id"].isin(fasta_ids)]
    merged = merged.rename(columns={abundance_col: "abundance"})

    logger.info(
        "Mapped %d abundance entries to %d FASTA ids (via mapping table).",
        len(merged),
        len(fasta_ids),
    )
    return merged
