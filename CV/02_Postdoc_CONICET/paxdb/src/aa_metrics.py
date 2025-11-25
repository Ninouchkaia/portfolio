# paxdb/src/aa_metrics.py

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np
import pandas as pd

from paxdb.src.protein import Protein

logger = logging.getLogger(__name__)

STANDARD_AA = list("ACDEFGHIKLMNPQRSTVWY")


def compute_proteome_aa_usage(
    proteins: Iterable[Protein],
    amino_acids: Iterable[str] = STANDARD_AA,
) -> pd.DataFrame:
    """
    Compute amino-acid usage across a proteome.

    Parameters
    ----------
    proteins : Iterable[Protein]
        Collection of Protein objects with abundance.
    amino_acids : Iterable[str]
        List of amino acids to consider.

    Returns
    -------
    DataFrame
        Index: amino acids; columns: ['count', 'weighted_count', 'freq', 'weighted_freq'].
    """
    aa_list = list(amino_acids)
    counts = dict.fromkeys(aa_list, 0.0)
    weighted_counts = dict.fromkeys(aa_list, 0.0)

    total_residues = 0.0
    total_weighted_residues = 0.0

    for protein in proteins:
        seq = protein.sequence
        abundance = protein.abundance

        for aa in seq:
            if aa not in counts:
                continue
            counts[aa] += 1.0
            weighted_counts[aa] += abundance
            total_residues += 1.0
            total_weighted_residues += abundance

    data = []
    for aa in aa_list:
        c = counts[aa]
        wc = weighted_counts[aa]
        freq = c / total_residues if total_residues > 0 else 0.0
        wfreq = wc / total_weighted_residues if total_weighted_residues > 0 else 0.0
        data.append((aa, c, wc, freq, wfreq))

    df = pd.DataFrame(
        data,
        columns=["aa", "count", "weighted_count", "freq", "weighted_freq"],
    ).set_index("aa")

    return df


def load_amino_acid_costs(path: Path) -> Dict[str, float]:
    """
    Load amino-acid biosynthetic costs from a TSV file.

    Parameters
    ----------
    path : Path
        TSV file with columns: 'aa', 'cost_atp'.

    Returns
    -------
    dict
        Mapping single-letter AA -> cost.
    """
    df = pd.read_csv(path, sep="\t", comment="#", dtype={"aa": str})
    if "aa" not in df.columns or "cost_atp" not in df.columns:
        raise ValueError(
            f"Cost file {path} must contain columns 'aa' and 'cost_atp'."
        )

    costs = dict(zip(df["aa"], df["cost_atp"]))
    return costs


def aa_vector_from_costs(
    aa_costs: Dict[str, float],
    amino_acids: Iterable[str] = STANDARD_AA,
) -> np.ndarray:
    """
    Build a cost vector aligned with the STANDARD_AA list.
    """
    aa_list = list(amino_acids)
    return np.array([aa_costs.get(aa, np.nan) for aa in aa_list], dtype=float)


def aa_vector_from_freqs(
    aa_usage: pd.DataFrame,
    column: str = "freq",
    amino_acids: Iterable[str] = STANDARD_AA,
) -> np.ndarray:
    """
    Build a usage-frequency vector aligned with the STANDARD_AA list.
    """
    aa_list = list(amino_acids)
    return np.array(
        [aa_usage.loc[aa, column] if aa in aa_usage.index else np.nan for aa in aa_list],
        dtype=float,
    )
