# paxdb/src/aa_metrics.py

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np
import pandas as pd

from .protein import Protein

logger = logging.getLogger(__name__)

STANDARD_AA = list("ACDEFGHIKLMNPQRSTVWY")


def compute_proteome_aa_usage(
    proteins: Iterable[Protein],
    weighted: bool = True,
    allowed_aas: Iterable[str] = STANDARD_AA,
) -> Dict[str, float]:
    """
    Compute proteome-level amino acid relative abundances.

    Parameters
    ----------
    proteins : iterable of Protein
        Proteins with abundance and sequence information.
    weighted : bool
        If True, weight counts by abundance; otherwise each protein
        contributes equally (unweighted composition).
    allowed_aas : iterable of str
        Consider only this set of amino acids (typically 20 standard).

    Returns
    -------
    dict
        Mapping amino-acid -> relative frequency in proteome.
        Frequencies sum to 1.0 over allowed_aas (if any non-zero counts).
    """
    allowed = set(allowed_aas)
    total_counts = {aa: 0.0 for aa in allowed}

    for p in proteins:
        counts = p.weighted_amino_acid_counts() if weighted else p.amino_acid_counts()
        for aa, c in counts.items():
            if aa in allowed:
                total_counts[aa] += float(c)

    total = float(sum(total_counts.values())) or 1.0
    freqs = {aa: c / total for aa, c in total_counts.items()}
    logger.debug(
        "Computed %sweighted AA usage: %s",
        "" if weighted else "un",
        freqs,
    )
    return freqs


def load_amino_acid_costs(path: Path) -> Dict[str, float]:
    """
    Load amino acid synthesis costs (e.g. ATP/time units).

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
    df["cost_atp"] = pd.to_numeric(df["cost_atp"], errors="coerce")
    df = df.dropna(subset=["cost_atp"])
    costs = dict(zip(df["aa"].str.upper(), df["cost_atp"].astype(float)))
    logger.info("Loaded amino acid costs for %d residues from %s", len(costs), path)
    return costs


def aa_vector_from_freqs(
    freqs: Dict[str, float],
    ordered_aas: Iterable[str] = STANDARD_AA,
) -> np.ndarray:
    """
    Convert AA frequency dict into a numeric vector in fixed order.
    """
    aa_list = list(ordered_aas)
    return np.array([float(freqs.get(aa, 0.0)) for aa in aa_list], dtype=float)


def aa_vector_from_costs(
    costs: Dict[str, float],
    ordered_aas: Iterable[str] = STANDARD_AA,
) -> np.ndarray:
    """
    Convert AA cost dict into a numeric vector in fixed order.
    """
    aa_list = list(ordered_aas)
    return np.array([float(costs.get(aa, np.nan)) for aa in aa_list], dtype=float)
